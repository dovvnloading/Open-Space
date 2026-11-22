"""Lossless, per-cache packing of Chandra catalog response bodies in the research workflow."""
import argparse
from collections import defaultdict
from functools import lru_cache
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import tempfile
from urllib.parse import unquote, urlsplit
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

import campaign

BLOCK = 8 * 1024**2
PACK = 16 * 1024**2
MAX_ZIP = 24 * 1024**2
MAX_MEMBERS = 512
MANIFEST = "MEMBERS.json"
FORMAT = "chandra-catalog-zip-v1"
H = re.compile(r"[0-9a-f]{64}")


def need(test, msg):
    if not test:
        raise ValueError(msg)


def sha(p):
    with p.open("rb") as f:
        return hashlib.file_digest(f,"sha256").hexdigest()


def stamp(p):
    s=p.stat()
    return (s.st_dev,s.st_ino,s.st_size,s.st_mtime_ns,s.st_ctime_ns)


@lru_cache(maxsize=131072)
def checked_hash(path,identity):
    return sha(path)


def digest(b):
    return hashlib.sha256(b).hexdigest()


def encoded(x):
    return (json.dumps(x,sort_keys=True,separators=(",",":"),allow_nan=False)+"\n").encode()


def valid_name(name, namespace):
    p=PurePosixPath(name)
    need(not p.is_absolute() and ".." not in p.parts and str(p)==name and
         p.parts[:2]==("catalogs",namespace) and p.suffix==".bin","Invalid logical response path")
    return name


def safe(p, root):
    need(not p.is_symlink() and p.resolve().is_relative_to(root.resolve()),"Unsafe evidence path")
    return p


def source_map(rows):
    return digest(encoded([[n,*rows[n]] for n in sorted(rows)]))


@lru_cache(maxsize=32)
def read_index(path, namespace, identity):
    index = campaign.read(path)
    need(index["format"] == FORMAT and index["namespace"] == namespace, "Unexpected catalog index")
    need(source_map(index["memos"]) == index["memo_tree_sha256"], "Original request-memo index changed")
    return index


@lru_cache(maxsize=64)
def pack_manifest(path, namespace, checksum, identity):
    need(checked_hash(path, identity) == checksum, "ZIP size or content address changed")
    entries = defaultdict(list)
    with ZipFile(path) as z:
        infos = z.infolist()
        need(len(infos) <= MAX_MEMBERS + 1 and len({r.filename for r in infos}) == len(infos)
             and all(not r.is_dir() and r.file_size <= BLOCK for r in infos), "Invalid ZIP member size or inventory")
        rows = json.loads(z.read(MANIFEST))
        need({r["member"] for r in rows} | {MANIFEST} == set(z.namelist()) and len(rows) == len(infos) - 1,
             "ZIP member inventory changed")
        for r in rows:
            name = valid_name(r["path"], namespace)
            need(H.fullmatch(r["source_sha256"]) and H.fullmatch(r["sha256"]) and
                 r["bytes"] == z.getinfo(r["member"]).file_size and 0 <= r["bytes"] <= BLOCK,
                 "Invalid member receipt")
            entries[name].append((path.name, r))
    return entries


class Store:
    """The top-level index locates bounded ZIPs; each ZIP has its own member manifest."""
    def __init__(self, directory, namespace, only=None):
        self.dir=directory
        self.namespace=namespace
        self.sources={}
        self.entries=defaultdict(list)
        self.index=None
        p=directory/"index.json"
        if not p.is_file():
            return
        self.index=read_index(p,namespace,stamp(p))
        for loose in self.index["loose"]:
            if only is not None and loose["path"] != only:
                continue
            n=valid_name(loose["path"],namespace)
            need(n not in self.sources and isinstance(loose["bytes"],int) and loose["bytes"]>=0 and
                 H.fullmatch(loose["sha256"]),"Invalid or duplicate loose response")
            self.sources[n]=(loose["bytes"],loose["sha256"])
        for item in self.index["packs"]:
            if only is not None and not item["first_path"] <= only <= item["last_path"]:
                continue
            zpath=directory/item["file"]
            need(item["file"]==item["sha256"]+".zip" and H.fullmatch(item["sha256"])
                 and safe(zpath,directory).stat().st_size==item["bytes"]<MAX_ZIP,
                 "ZIP size or content address changed")
            members=pack_manifest(zpath,namespace,item["sha256"],stamp(zpath))
            need(members and (min(members),max(members))==(item["first_path"],item["last_path"]),
                 "Invalid ZIP logical-path range")
            if only is not None:
                if only in members:
                    self.entries[only].extend(members[only])
            else:
                for n, chunks in members.items():
                    self.entries[n].extend(chunks)
        for n, chunks in self.entries.items():
            chunks.sort(key=lambda x:x[1]["part"])
            raw=chunks[0][1]
            need(n not in self.sources and [r["part"] for _,r in chunks]==list(range(len(chunks))) and
                 all((r["source_sha256"],r["source_bytes"])==(raw["source_sha256"],raw["source_bytes"]) for _,r in chunks)
                 and sum(r["bytes"] for _,r in chunks)==raw["source_bytes"],"Incomplete or conflicting catalog parts")
            self.sources[n]=(raw["source_bytes"],raw["source_sha256"])
        if only is None:
            need(source_map(self.sources)==self.index["logical_tree_sha256"],"Catalog logical-path inventory changed")

    def blocks(self, name, evidence, opened):
        if name in self.entries:
            for pack,r in self.entries[name]:
                if pack not in opened:
                    if len(opened)>=8:
                        opened.pop(next(iter(opened))).close()
                    opened[pack]=ZipFile(self.dir/pack)
                data=opened[pack].read(r["member"])
                need(digest(data)==r["sha256"],"Catalog part checksum mismatch")
                yield data
        else:
            p=safe(evidence/name,evidence)
            with p.open("rb") as f:
                while b:=f.read(BLOCK):
                    yield b


def verify_sources(store,root,names=None,dest=None):
    opened={}
    try:
        for name in sorted(store.sources if names is None else names):
            size,h=store.sources[name]
            out=None
            if dest is not None:
                p=safe(dest/name,dest)
                p.parent.mkdir(parents=True,exist_ok=True)
                if p.is_file():
                    need(p.stat().st_size==size and sha(p)==h,"Conflicting restoration target")
                    continue
                out=tempfile.NamedTemporaryFile(dir=p.parent,delete=False)
            result=hashlib.sha256()
            seen=0
            try:
                for data in store.blocks(name,root,opened):
                    result.update(data)
                    seen+=len(data)
                    if out:
                        out.write(data)
                need((seen,result.hexdigest())==(size,h),f"Catalog body mismatch: {name}")
                if out:
                    out.close()
                    try:
                        os.link(out.name,p)
                    except FileExistsError:
                        need(p.stat().st_size==size and sha(p)==h,"Concurrent restore disagrees")
            finally:
                if out:
                    out.close()
                    Path(out.name).unlink(missing_ok=True)
    finally:
        for z in opened.values():
            z.close()


@lru_cache(maxsize=128)
def cached_store(directory,namespace,index_identity,name):
    return Store(directory,namespace,only=name)


@lru_cache(maxsize=131072)
def verified_packed(root,name,index_identity,pack_identities):
    namespace=PurePosixPath(name).parts[1]
    db=cached_store(root/"catalog_store"/namespace,namespace,index_identity,name)
    verify_sources(db,root,[name])


def memo_unchanged(root, path):
    """Check the original loose HTTP memo against a packed producer's receipt."""
    logical = safe(path, root).relative_to(root)
    need(len(logical.parts) >= 3 and logical.parts[0] == "catalogs", "Expected a catalog request memo")
    namespace = logical.parts[1]
    index = root / "catalog_store" / namespace / "index.json"
    if not index.is_file():
        return path.is_file()
    if not path.is_file():
        return False
    old = read_index(index, namespace, stamp(index))["memos"].get(logical.as_posix())
    current = stamp(path)
    return old == [current[2], checked_hash(path, current)]


def recorded_reply(root,path,receipt,destination=None):
    """Authenticate one logical HTTP body and optionally copy it to private staging."""
    path=safe(path,root)
    logical=path.relative_to(root).as_posix()
    parts=PurePosixPath(logical).parts
    need(len(parts)>=3 and parts[0]=="catalogs","Expected an original survey catalog response")
    name=valid_name(logical,parts[1])
    size,h=receipt["size_bytes"],receipt["sha256"]
    need(isinstance(size,int) and size>=0 and H.fullmatch(h),"Invalid recorded catalog byte identity")
    failure=f"Original catalog response is missing or altered: {name}"
    db=None
    if path.is_file():
        identity=stamp(path)
        need(identity[2]==size and checked_hash(path,identity)==h,failure)
    else:
        directory=root/"catalog_store"/parts[1]
        index=directory/"index.json"
        need(index.is_file(),failure)
        index_identity=stamp(index)
        db=cached_store(directory,parts[1],index_identity,name)
        need(db.sources.get(name)==(size,h) and name in db.entries,failure)
        packed=[]
        for pack in dict.fromkeys(pack for pack,_ in db.entries[name]):
            filename=directory/pack
            identity=stamp(filename)
            need(checked_hash(filename,identity)==pack.removesuffix(".zip"),failure)
            packed.append(identity)
        if destination is None or destination.is_file():
            verified_packed(root,name,index_identity,tuple(packed))
    if destination is None:
        return
    need(not destination.resolve().is_relative_to(root.resolve()),"Use independent reply staging")
    destination.parent.mkdir(parents=True,exist_ok=True)
    if destination.exists():
        need(destination.stat().st_size==size and sha(destination)==h,"Conflicting staged catalog reply")
        return
    with tempfile.NamedTemporaryFile(dir=destination.parent,prefix=".restore-",delete=False) as out:
        temp=Path(out.name)
    opened={}
    try:
        if db is None:
            shutil.copyfile(path,temp)
        else:
            with temp.open("wb") as out:
                for chunk in db.blocks(name,root,opened):
                    out.write(chunk)
        need(temp.stat().st_size==size and sha(temp)==h,failure)
        try:
            os.link(temp,destination)
        except FileExistsError:
            need(destination.stat().st_size==size and sha(destination)==h,"Concurrent reply restoration disagrees")
    finally:
        for z in opened.values():
            z.close()
        temp.unlink(missing_ok=True)


def request_memos(folder,root,sources):
    """Validate every successful or failed HTTP body identified by a stored memo."""
    memo_bytes={}
    for p in folder.rglob("*.json"):
        safe(p,root)
        memo=campaign.read(p)
        for item in [memo,*memo.get("attempts",[])]:
            if "body_file" not in item:
                continue
            name=item["body_file"]
            need(Path(name).name==name,"Nonlocal body in request memo")
            rel=(p.parent/name).relative_to(root).as_posix()
            need(sources.get(rel)==(item["size_bytes"],item["sha256"]),f"Original memo/body mismatch: {rel}")
        memo_bytes[p.relative_to(root).as_posix()]=(p.stat().st_size,sha(p))
    return memo_bytes


def protected(root):
    """Keep direct browser links and every raw real instrumental pass query accessible."""
    keep=set()
    for p in (root/"observations").glob("*.json.gz"):
        row=campaign.read_obs(p)
        for prod in row.get("products",[]):
            for t in (prod.get("search") or {}).get("triggers",[]):
                if not t.get("instrument_pass"):
                    continue
                cache=t.get("catalog_cache") or row.get("catalog_cache")
                for q in (t.get("vet") or {}).get("queries",{}).values():
                    for item in [q,*q.get("attempts",[])]:
                        if cache and item.get("body_file"):
                            keep.add(str(PurePosixPath(cache)/item["body_file"]))
    for site in {campaign.HERE,root}:
        for md in site.rglob("*.md"):
            for u in re.findall(r"\]\(([^)\n]+)\)",md.read_text(errors="replace")):
                ref=urlsplit(u.strip("<>"))
                if ref.scheme or ref.netloc or not ref.path.endswith(".bin"):
                    continue
                target=(md.parent/unquote(ref.path)).resolve()
                for base in (root,campaign.EVIDENCE):
                    if target.is_relative_to(base/"catalogs"):
                        keep.add(target.relative_to(base).as_posix())
    return keep


def info(name):
    obj=ZipInfo(name,date_time=(1980,1,1,0,0,0))
    obj.create_system=3
    obj.external_attr=0o100644<<16
    obj.compress_type=ZIP_DEFLATED
    return obj


class Writer:
    def __init__(self,path):
        self.path=path
        self.z=None
        self.list=[]
        self.size=0
        self.packs=[]

    def flush(self):
        if self.z is None:
            return
        contents=encoded(self.list)
        need(len(contents)<=BLOCK,"Oversized ZIP member manifest")
        self.z.writestr(info(MANIFEST),contents,compresslevel=6)
        self.z.close()
        path=self.path/"writing.zip"
        checksum=sha(path)
        need(path.stat().st_size<MAX_ZIP,"ZIP is too large")
        self.packs.append({"file":checksum+".zip","sha256":checksum,"bytes":path.stat().st_size,
                           "member_count":len(self.list),"native_chunk_bytes":self.size,
                           "first_path":self.list[0]["path"],"last_path":self.list[-1]["path"]})
        path.rename(self.path/(checksum+".zip"))
        self.z=None
        self.list=[]
        self.size=0

    def add(self,name,source_size,source_hash,index,b):
        if self.z is not None and (len(self.list)>=MAX_MEMBERS or self.size+len(b)>PACK):
            self.flush()
        if self.z is None:
            self.z=ZipFile(self.path/"writing.zip","x",compression=ZIP_DEFLATED,compresslevel=6)
        member=name if source_size<=BLOCK else f"pieces/{digest(name.encode())}/{index:08d}.part"
        self.z.writestr(info(member),b,compresslevel=6)
        self.list.append({"path":name,"source_bytes":source_size,"source_sha256":source_hash,"part":index,
                          "member":member,"sha256":digest(b),"bytes":len(b)})
        self.size+=len(b)


def install(original,target):
    # Transfer across runner volumes without exposing a partially copied pack.
    try:
        original.replace(target)
    except OSError as exc:
        import errno
        if exc.errno != errno.EXDEV:
            raise
        with tempfile.NamedTemporaryFile(dir=target.parent,prefix=".catalog-install-",delete=False) as out:
            scratch=Path(out.name)
        try:
            shutil.copyfile(original,scratch)
            need(sha(scratch)==sha(original),"Cross-volume pack transfer changed bytes")
            scratch.replace(target)
        finally:
            scratch.unlink(missing_ok=True)


def pack(root,namespace,keep,dry,prune):
    folder=root/"catalogs"/namespace
    directory=root/"catalog_store"/namespace
    prior=Store(directory,namespace)
    sources=dict(prior.sources)
    native={}
    for p in folder.rglob("*.bin"):
        safe(p,root)
        n=valid_name(p.relative_to(root).as_posix(),namespace)
        native[n]=(p.stat().st_size,sha(p))
        need(n not in sources or sources[n]==native[n],f"Unchanged logical name has different bytes: {n}")
        sources[n]=native[n]
    memos=request_memos(folder,root,sources)
    if prior.index:
        need(all(memos.get(n)==tuple(v) for n,v in prior.index["memos"].items()),
             "An original loose request memo was removed or altered")
    pinned=set(sources)&keep
    # Hydrate newly linked bodies before switching the canonical index.
    if not dry:
        previous_links=pinned-set(native)
        verify_sources(prior,root,previous_links,root)
    with tempfile.TemporaryDirectory(dir=os.environ.get("RUNNER_TEMP"),prefix="catalog-pack-") as tmp:
        stage=Path(tmp)
        writer=Writer(stage)
        opened={}
        try:
            for n in sorted(sources.keys()-pinned):
                size,h=sources[n]
                check=hashlib.sha256()
                seen=0
                if n in native:
                    def chunks(path):
                        with path.open("rb") as f:
                            while b:=f.read(BLOCK):
                                yield b
                    body=chunks(root/n)
                else:
                    body=prior.blocks(n,root,opened)
                for i,b in enumerate(body):
                    if not b:
                        continue
                    writer.add(n,size,h,i,b)
                    check.update(b)
                    seen+=len(b)
                if size==0:
                    writer.add(n,size,h,0,b"")
                need(seen==size and check.hexdigest()==h,"Catalog input changed during packing")
            writer.flush()
        finally:
            if writer.z:
                writer.z.close()
            for z in opened.values():
                z.close()
        index={"format":FORMAT,"namespace":namespace,"logical_tree_sha256":source_map(sources),
               "memos":memos,"memo_tree_sha256":source_map(memos),
               "packs":writer.packs,"loose":[{"path":n,"bytes":sources[n][0],"sha256":sources[n][1]} for n in sorted(pinned)]}
        (stage/"index.json").write_bytes(encoded(index))
        checked=Store(stage,namespace)
        verify_sources(checked,root,checked.entries)
        facts={"namespace":namespace,"bodies":len(sources),"native_bytes":sum(s for s,_ in sources.values()),
               "loose_bytes":sum(sources[n][0] for n in pinned),"loose_bodies":len(pinned),
               "zip_bytes":sum(p["bytes"] for p in writer.packs),"packs":len(writer.packs)}
        if not dry:
            verify_sources(checked,root,pinned)
            directory.mkdir(parents=True,exist_ok=True)
            for part in writer.packs:
                p=directory/part["file"]
                if p.exists():
                    need(sha(p)==part["sha256"],"Conflicting existing ZIP")
                else:
                    install(stage/part["file"],p)
            # Publish verified new packs before the index; only then remove
            # optional loose duplicates. Interrupted runs always keep all bytes.
            install(stage/"index.json",directory/"index.json")
            if prune:
                for n in sorted(set(checked.entries)&set(native)):
                    need(sha(root/n)==sources[n][1],"Catalog response changed before pruning")
                    (root/n).unlink()
            if prior.index:
                for old in {r["file"] for r in prior.index["packs"]}-{r["file"] for r in index["packs"]}:
                    (directory/old).unlink(missing_ok=True)
        return facts


def namespaces(root):
    return sorted({p.name for base in (root/"catalogs",root/"catalog_store") if base.is_dir()
                   for p in base.iterdir() if p.is_dir() and not p.name.startswith(".")})


def selected(root,args):
    if getattr(args,"obsid",None) is not None:
        row=campaign.read_obs(root/"observations"/f"{args.obsid}.json.gz")
        roots={row.get("catalog_cache")}
        for prod in row.get("products",[]):
            scan=prod.get("search") or {}
            for stream in [scan,*scan.get("controls",[])]:
                roots.update(t.get("catalog_cache") for t in stream.get("triggers",[]))
        return sorted({p.parts[1] for value in roots if value for p in [PurePosixPath(value)]
                       if p.parts[:1]==("catalogs",) and len(p.parts)>1})
    return [args.cache] if getattr(args,"cache",None) else namespaces(root)


def extract(root,ns,destination):
    safe(destination/"catalogs"/ns,destination).mkdir(parents=True,exist_ok=True)
    store=Store(root/"catalog_store"/ns,ns)
    # Copy the original HTTP metadata. A previously uncommitted
    # archive or overlapping working copy may contain additional loose bodies.
    for p in (root/"catalogs"/ns).rglob("*"):
        safe(p,root)
        if not p.is_file():
            continue
        target=safe(destination/p.relative_to(root),destination)
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():
            need(sha(p)==sha(target),"Conflicting temporary extraction")
        else:
            # Never let a consumer's ordinary write mutate the source checkout.
            shutil.copyfile(p,target)
    verify_sources(store,root,dest=destination)
    actual={p.relative_to(destination).as_posix():(p.stat().st_size,sha(p)) for p in
            (destination/"catalogs"/ns).rglob("*.bin")}
    memos=request_memos(destination/"catalogs"/ns,destination,actual)
    if store.index:
        need(source_map(memos)==store.index["memo_tree_sha256"],"Restored request metadata changed")


def lock(root,namespace,exclusive):
    import fcntl
    key=digest(str(root/namespace).encode())
    p=Path(os.getenv("RUNNER_TEMP",tempfile.gettempdir()))/("catalog-store-"+key+".lock")
    stream=p.open("a")
    fcntl.flock(stream,fcntl.LOCK_EX if exclusive else fcntl.LOCK_SH)
    return stream


README="""# Archived catalog responses

Every HTTP memo remains at its original path under `evidence/catalogs/`. The loose binary files include replies needed by instrument-passing real triggers and directly linked review evidence. Each remaining reply retains its original logical pathname, byte length and SHA-256 in its producer’s ZIP member manifest. Each producer index records its sorted logical file inventory checksum and the compressed hash of each ZIP. Open a producer directory and download a ZIP from GitHub to read its `MEMBERS.json`. Large logical replies are split into ordered members. The acquisition memos include failed HTTP responses and retry histories.

From the repository root in the Chandra workflow, locate one reply, verify all bundles, or reconstruct a cache in a separate directory.

```bash
python results/chandra_acis_fxt/catalog_pack.py where catalogs/CACHE/REPLY.bin
python results/chandra_acis_fxt/catalog_pack.py verify
python results/chandra_acis_fxt/catalog_pack.py extract --cache CACHE --destination "$RUNNER_TEMP/catalog-view"
python results/chandra_acis_fxt/catalog_pack.py replay OBSID
```

The replay command restores the caches used by that observation in a temporary directory and invokes the original checksummed observation replay. The logical filenames and response bytes supplied to the frozen classifier are identical to the original working files. [Operational details](../../CATALOG_STORAGE.md)
"""


def main():
    require_env=(os.getenv("GITHUB_ACTIONS")=="true" and os.getenv("GITHUB_REPOSITORY")=="dovvnloading/Open-Space"
                 and os.getenv("GITHUB_REF")=="refs/heads/research/chandra-acis-fast-transients")
    need(require_env,"Use the primary Chandra research workflow")
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--evidence",type=Path,default=campaign.EVIDENCE)
    commands=cli.add_subparsers(dest="mode",required=True)
    for mode in ("measure","pack","verify","extract","where","replay"):
        c=commands.add_parser(mode)
        if mode in {"measure","pack","verify"}:
            c.add_argument("--cache")
        if mode=="pack":
            c.add_argument("--prune",action="store_true",help="remove checked loose duplicates after verification")
        if mode=="extract":
            which=c.add_mutually_exclusive_group(required=True)
            which.add_argument("--all",action="store_true")
            which.add_argument("--obsid",type=int)
            which.add_argument("--cache")
            c.add_argument("--destination",type=Path,required=True)
        if mode=="where":
            c.add_argument("member")
        if mode=="replay":
            c.add_argument("obsid",type=int)
    args=cli.parse_args()
    root=args.evidence.resolve()
    if args.mode=="where":
        p=PurePosixPath(args.member)
        need(len(p.parts)>2,"Expected catalogs/CACHE/body.bin")
        n=valid_name(args.member,p.parts[1])
        db=Store(root/"catalog_store"/p.parts[1],p.parts[1])
        in_tree=(root/n).is_file()
        need(in_tree or n in db.entries,"Source path absent from the catalog store")
        print(json.dumps({"logical_path":n,"loose":in_tree,"native_size_sha256":db.sources.get(n),
                          "zip_members":db.entries.get(n,[])},indent=2))
        return
    choices=selected(root,args)
    need(all(re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*",n) for n in choices),"Invalid catalog namespace")
    if args.mode=="replay":
        import reproduce
        from unittest.mock import patch
        with tempfile.TemporaryDirectory() as tmp:
            home=Path(tmp)
            view=home/"view"
            view.mkdir()
            for f in ("inventory","observations","photometry"):
                if (root/f).exists():
                    (view/f).symlink_to(root/f,target_is_directory=True)
            for ns in choices:
                with lock(root,ns,False):
                    extract(root,ns,view)
            (home/"replay").mkdir()
            with patch.object(reproduce,"EVIDENCE",view):
                reproduce.reproduce(args.obsid,home/"replay")
        return
    keep=protected(root) if args.mode in {"pack","measure"} else set()
    if args.mode=="pack":
        (root/"catalog_store").mkdir(exist_ok=True)
    target=getattr(args,"destination",None)
    if target:
        target=target.resolve()
        need(not target.is_relative_to(root) and not root.is_relative_to(target),"Extract into an independent temporary view")
    totals=defaultdict(int)
    for ns in choices:
        with lock(root,ns,args.mode=="pack"):
            if args.mode in {"measure","pack"}:
                facts=pack(root,ns,keep,args.mode=="measure",getattr(args,"prune",False))
                for k,v in facts.items():
                    if isinstance(v,int):
                        totals[k]+=v
            elif args.mode=="extract":
                extract(root,ns,target)
            else:
                db=Store(root/"catalog_store"/ns,ns)
                need(db.index is not None,f"Missing pack index for {ns}")
                verify_sources(db,root)
                memos=request_memos(root/"catalogs"/ns,root,db.sources)
                need(memos=={n:tuple(v) for n,v in db.index["memos"].items()} and
                     source_map(memos)==db.index["memo_tree_sha256"],
                     "Loose request receipts differ from the verified deposit")
                for p in (root/"catalogs"/ns).rglob("*.bin"):
                    safe(p,root)
                    n=p.relative_to(root).as_posix()
                    need(db.sources.get(n)==(p.stat().st_size,sha(p)),"Unindexed or altered loose catalog response")
    if args.mode=="pack":
        (root/"catalog_store/README.md").write_text(README)
    if totals:
        print(json.dumps(dict(totals),sort_keys=True))


if __name__=="__main__":
    main()

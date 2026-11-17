"""Search a ranked series of independent DESI files, bounded by two input shards.

Linux can keep the compressed input products in anonymous memory while CFITSIO
decodes them. On other platforms use disk-backed temporary files (--disk).
Every source file is verified against the DESI SHA-256 list before detection.
"""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from contextlib import ExitStack
import csv
from dataclasses import asdict
import gzip
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time

from injection_window import pack_injection_window
from lithium_profile import profile
from desi_shard import BASE, digest, obtain, product_groups, targets
from laser_search import scan_target


def acquire_inputs(cache, stack, survey, program, hp, mem):
    records = []
    jobs = []
    for folder, names, check in product_groups(survey, program, hp):
        base = BASE + folder + '/'
        p = cache / check
        obtain(base + check, p)
        records.append({'url': base + check, 'sha256': digest(p), 'bytes': p.stat().st_size})
        expected = {l.split()[1].lstrip('*'): l.split()[0] for l in p.read_text().splitlines() if l.strip()}
        for name in names:
            path = cache / name
            fd = None
            if mem and (name.endswith('.fits.gz') or name.startswith('rvmod_')):
                fd = stack.enter_context(os.fdopen(os.memfd_create(name), 'w+b'))
                path.symlink_to(f"/proc/self/fd/{fd.fileno()}")
            jobs.append((path, fd, base + name, expected[name]))

    def fetch(job):
        path, fd, url, expected = job
        if fd is None:
            obtain(url, path, expected)
        else:
            for attempt in range(3):
                fd.seek(0)
                fd.truncate()
                result = subprocess.run(['curl','-fLsS','--connect-timeout','30','--speed-time','45','--speed-limit','1024','--max-time','900',
                                         '--max-filesize','950000000',url], stdout=fd)
                if result.returncode == 0 and digest(path) == expected:
                    break
                if result.returncode == 63:
                    raise ValueError(f"per-file RAM budget exceeded: {url}")
            else:
                raise ValueError(f"failed download or SHA-256 mismatch: {url}")
        return {'url':url, 'sha256':expected, 'bytes':path.stat().st_size}

    with ThreadPoolExecutor(max_workers=4) as pool:
        records.extend(pool.map(fetch,jobs))
    return sorted(records,key=lambda row:row['url'])


def process_shard(task):
    survey, program, hp, out, mem = task
    tag = f"{survey}-{program}-{hp}"
    index = out / f"scale-count-{tag}.json"
    if index.exists():
        return json.loads(index.read_text())
    started = time.monotonic()
    count = Counter()
    output = out / f"scale-scan-{tag}.jsonl.gz"
    window_dir = out / "windows"
    with tempfile.TemporaryDirectory(dir=out) as work, ExitStack() as stack:
        cache = Path(work)
        sources = acquire_inputs(cache,stack,survey,program,hp,mem)
        (out/f'scale-inputs-{tag}.json').write_text(json.dumps(sources,indent=2)+'\n')
        transferred = time.monotonic()
        print(f"{tag} input verified after {transferred-started:.1f}s",flush=True)
        with gzip.open(str(output)+'.partial','wt',compresslevel=3) as f:
            for meta, arms in targets(cache,survey,program,hp,count):
                ids = [r['expid'] for r in meta['exposures']]
                by_arm = {}
                for arm, recs in arms.items():
                    res = scan_target(recs,ids)
                    by_arm[arm] = asdict(res)
                    visits = {e.expid:e for e in recs}
                    for candidate in by_arm[arm]['candidates']:
                        for hit in candidate['hits']:
                            e = visits[hit['expid']]
                            raw = profile(e,hit['rest_A'],len(e.resolution)//2+4)
                            hit['raw_profile'] = raw
                            hit['rv_error_kms'] = e.velocity_error_kms

                    count['searched_arms'] += 1
                    count['triggers'] += sum(len(h) for h in res.exposure_hits.values())
                    count['persistent'] += len(res.candidates)
                f.write(json.dumps({'meta':meta,'arms':by_arm},allow_nan=False,separators=(',',':'))+'\n')
                pack_injection_window(arms['b'], meta, window_dir)
        os.replace(str(output)+'.partial',output)
    summary=dict(count, survey=survey,program=program,healpix=hp,
                 transfer_s=round(transferred-started,2),total_s=round(time.monotonic()-started,2),
                 scan_sha256=digest(output),window_files=count['eligible_stars'])
    tmp = index.with_suffix('.part')
    tmp.write_text(json.dumps(summary,indent=2)+'\n')
    os.replace(tmp,index)
    print(f"{tag} finished {count['eligible_stars']} stars, {count['persistent']} hits, {summary['total_s']}s",flush=True)
    return summary


def run_shard(task):
    """Serialize concurrent processes before the completion marker check."""
    import fcntl
    survey, program, hp, out, _ = task
    locks = out / "locks"
    locks.mkdir(exist_ok=True)
    with (locks / f"{survey}-{program}-{hp}.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        return process_shard(task)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('selection',type=Path)
    p.add_argument('output',type=Path)
    p.add_argument('--start',type=int,default=0)
    p.add_argument('--count',type=int,default=1)
    p.add_argument('--workers',type=int,default=2)
    p.add_argument('--disk',action='store_true')
    args=p.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    rows=[r for r in csv.DictReader(args.selection.open()) if (r['survey'],r['program'])==('main','backup')]
    tasks=[(r['survey'],r['program'],int(r['healpix']),args.output,not args.disk and hasattr(os,'memfd_create'))
           for r in rows[args.start:args.start+args.count]]
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for row in pool.map(run_shard,tasks):
            pass

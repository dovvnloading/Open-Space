"""Build flight destinations from the reports and pixels at one Git revision."""

from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path, PurePosixPath
import posixpath
import re
from urllib.parse import parse_qs, quote, unquote, urlsplit

from open_space.verify import markdown
from open_space.verify.tree import RepoTree
from observatory.check_media import fits_header
from observatory.flight_plates import VECTOR_PLOTS, render_plate


# Editorial names and observing units. Scientific values are read from the
# source tables, not retyped here. Keys identify individual public records.
STOPS = {
    "beta-pictoris-pulsations": ("β Pictoris", "star", "TESS · Stellar pulsations",
        "Pulsation-band amplitudes measured in one TESS observing sector."),
    "cassini-uvis-enceladus-clock": ("Enceladus", "moon", "Cassini UVIS · Stellar occultation",
        "Cassini FUV and high-speed photometer integrations during the 27 March 2017 occultation."),
    "borisov-704021": ("2I/Borisov", "comet", "CFHT · Interstellar comet encounter",
        "Four asteroid detections in CFHT images and two reconstructions of the 2019 encounter."),
    "epsilon-eridani-keck-jwst": ("ε Eridani · Source S2", "star", "Keck + JWST · Relative astrometry",
        "One image of the source and its reference. Offsets are measured relative to the field object."),
    "wr-117-2": ("WR 117-2", "star", "SPHEREx · Infrared variability",
        "Aperture flux at the catalogue position on matched detector-6 visits."),
    "at-2021aeuk": ("AT 2021aeuk", "transient", "Swift · Ultraviolet and X-ray monitoring",
        "Swift visits to the galaxy nucleus in the 2025 and 2026 observing seasons."),
    "m87-nova-astrometry": ("M87 · Nova 23", "transient", "Hubble ACS · Nova astrometry",
        "HST difference images recover three reported nova positions. The pointing marks ACS source 23."),
    "jwst-repeat-imaging-transients": ("GOODS-S 185", "transient", "JWST NIRCam · Nuclear brightening",
        "Repeat imaging at the location of JADES DR2 galaxy 103607. The flare remains unclassified."),
    "rgz-desi-giant-radio-galaxies": ("J000314.9+111544", "galaxy", "FIRST + VLASS + DESI · Radio galaxies",
        "Radio hotspot separations and DESI host spectra for two proposed giant systems. The optical field shows J000314.9+111544."),
    "archival-radio-supernovae": ("SN 2019njv", "transient", "VLASS · Late radio emission",
        "Beam-fitted 3 GHz amplitudes at two supernova positions. The pointing marks SN 2019njv."),
    "mars-crater-size-correlation": ("HiRISE · Martian impact field", "planet", "Mars · Crater mapping",
        "One impact field per population measurement. The displayed HiRISE field contains 76 mapped craters."),
}

SPOTLIGHTS = {"mars": ("planet", "Mars · Particle environment"), "sun": ("star", "Sun · Heliophysics"),
              "quasar": ("galaxy", "Gravitational lens · Time delay"),
              "dwarf": ("star", "Brown dwarf binary · Radio monitoring")}

GOODS = "results/jwst-repeat-imaging-transients/goods-s-oct2026/README.md"
HEAD = re.compile(r"^ {0,3}(#{1,6})\s+(.+?)\s*#*\s*$")
DEGREES = r"([+−-]?\d+\.\d+)°?,\s*([+−-]?\d+\.\d+)°?"
# A field in a paper with several positions always names the row and column.
POINTINGS = {
    "wr-117-2": ("J2000", "Published WR 117-2 position", rf"published J2000 position \({DEGREES}\)"),
    "at-2021aeuk": ("J2000", "Swift fixed optical aperture", rf"fixed optical position, {DEGREES} \(J2000\)"),
    "m87-nova-astrometry": ("ICRS", "ACS 23 · HST difference-image position",
                            rf"^\| 23 \| [^|]+ \| {DEGREES} \|"),
    "jwst-repeat-imaging-transients": ("J2000", "Site 185 · Nuclear aperture", rf"^\| 185 \| {DEGREES} \|"),
}

FITS_PLATES = {
    "borisov-704021": ("data/2453635p.fits", "CFHT / MegaCam", "(704021) 2007 XJ8 · CFHT exposure 2453635p",
        "CFHT exposure 2453635p at the predicted asteroid position in the Borisov field."),
    "archival-radio-supernovae": ("inputs/SN2019njv-3.fits", "NRAO / VLA Sky Survey", "3 GHz radio field · VLASS epoch 3",
        "VLASS 3 GHz field at SN 2019njv, in the third survey epoch. A variable host source remains possible."),
}


def plain(raw):
    for link in reversed(markdown.links(raw)):
        raw = raw[:link.start] + link.text + raw[link.end:]
    return " ".join(raw.replace("**", "").replace("`", "").strip().split())


def scope(document, anchor):
    """Return the exact section and its first source line, also for repeated headings."""
    if not anchor:
        return document, 1
    prose, _ = markdown.split_code(document)
    original = document.splitlines()
    start, level, counts = None, 0, {}
    for index, line in enumerate(prose.splitlines()):
        heading = HEAD.match(line)
        if not heading:
            continue
        depth, title = len(heading[1]), heading[2]
        if start is not None and depth <= level:
            return "\n".join(original[start:index]), start + 1
        if start is None:
            for slug in sorted({markdown.slug(title), markdown.slug(title.replace("_", "").replace("*", ""))}):
                seen = counts.get(slug, 0)
                counts[slug] = seen + 1
                if anchor == (f"{slug}-{seen}" if seen else slug):
                    start, level = index, depth
    return ("\n".join(original[start:]), start + 1) if start is not None else ("", 1)


def receipt(repository, tree, path, label, line=None):
    url = f"{repository}/blob/{tree.commit}/{quote(path)}"
    if line:
        url += f"?plain=1#L{line}"
    return {"path": path, "label": label, "url": url}


def report_table(document, line_offset, source):
    """Keep the published cell text, units, bounds and neighboring comparison rows."""
    original = document.splitlines()
    prose, _ = markdown.split_code(document)
    lines = prose.splitlines()
    heading = None
    for index, line in enumerate(lines[:-2]):
        match = HEAD.match(line)
        if match:
            heading = plain(match[2])
        if not line.strip().startswith("|") or not line.rstrip().endswith("|") or not re.fullmatch(
                r"\s*\|(?:\s*:?-{3,}:?\s*\|){2,12}\s*", lines[index + 1]):
            continue
        def cells(row):
            return [plain(value.replace("\\|", "|")) for value in re.split(r"(?<!\\)\|", row.strip()[1:-1])]
        columns = cells(original[index])
        rows = []
        for raw in original[index + 2:]:
            if not raw.strip().startswith("|") or not raw.strip().endswith("|"):
                break
            row = cells(raw)
            if len(row) != len(columns):
                rows = []
                break
            rows.append(row)
        if (not rows or not any(re.search(r"\d", value) for row in rows for value in row)
                or any(len(value) > 650 for row in rows[:12] for value in row)):
            continue
        return {"caption": heading or "Published measurements", "columns": columns, "rows": rows[:12],
                "totalRows": len(rows), "source": source("Published table", line_offset + index)}
    return None


def report_excerpt(document, first_line, summary):
    """Quote an opening scientific paragraph, preserving all of its qualifications."""
    for match in re.finditer(r"(?m)^(?!\s*$)\S[^\n]*(?:\n(?!\s*$)[^\n]+)*", document):
        paragraph = match[0].strip()
        if paragraph.startswith("#"):
            if re.search(r"\b(reproduc\w*|commands?|installation|running|how to|requirements|files)\b",
                         paragraph.splitlines()[0], re.I):
                return None
            continue
        if paragraph.startswith(("|", "!", "<", "- ", "* ", "> ", "```", "~~~")):
            continue
        # A trailing bibliography can have many links separated by middle dots.
        # Keep author prose intact when a link is itself the end of a sentence.
        for link in reversed(markdown.links(paragraph)):
            if link.end != len(paragraph) or not re.search(r"(?:[.!?]\**|·)\s*$", paragraph[:link.start]):
                break
            paragraph = paragraph[:link.start].rstrip(" ·\t")
        excerpt = plain(paragraph)
        excerpt = re.sub(r"^Noteworthy (?:result|correction) · [^.]+\.\s*", "", excerpt)
        if not re.search(r"\w", excerpt) or len(excerpt) < 50:
            continue
        if (len(excerpt) > 1_200 or excerpt == summary or
                re.match(r"(?i)^(run|use|download|install|execute|to (?:run|reproduce)|requirements)\b", excerpt)
                or not re.search(r"[.!?][’”)]*$", excerpt)):
            return None
        return excerpt, first_line + document.count("\n", 0, match.start())
    return None


def package_image(tree, out, repository, path, caption, media_type, credit, raw, suffix):
    digest = hashlib.sha256(raw).hexdigest()
    name = digest[:20] + suffix
    out.mkdir(parents=True, exist_ok=True)
    (out / name).write_bytes(raw)
    return {**receipt(repository, tree, path, "Original repository image"),
            "type": media_type, "src": f"data/flight/{name}", "alt": caption, "caption": caption,
            "sha256": digest, **({"credit": credit} if credit else {})}


def publish_image(tree, read_blob, out, repository, path, caption, media_type="research-figure", credit=None):
    # Publish only static image formats, with immutable names. File contents
    # come from the requested Git object, never a path from the working tree.
    suffix = PurePosixPath(path).suffix.lower()
    if path not in tree or (suffix not in {".png", ".jpg", ".jpeg", ".webp"} and path not in VECTOR_PLOTS):
        return None
    raw = read_blob(path)
    if not 16 < len(raw) <= 2_000_000:
        return None
    if path in VECTOR_PLOTS:
        valid = hashlib.sha256(raw).hexdigest() == VECTOR_PLOTS[path]
    elif suffix == ".png":
        valid = raw.startswith(b"\x89PNG\r\n\x1a\n\0\0\0\rIHDR")
    elif suffix == ".webp":
        valid = raw[:4] == b"RIFF" and raw[8:12] == b"WEBP"
    else:
        valid = raw.startswith(b"\xff\xd8\xff") and raw.endswith(b"\xff\xd9")
    if not valid:
        return None
    return package_image(tree, out, repository, path, caption, media_type, credit, raw, suffix)


def point(ra, dec, frame, label, source):
    ra, dec = (float(str(value).replace("−", "-")) for value in (ra, dec))
    if not (0 <= ra < 360 and -90 <= dec <= 90):
        raise ValueError(f"Invalid sky position at {source['path']}")
    return {"raDeg": ra, "decDeg": dec, "frame": frame, "label": label, "source": source}


def pointing(stop, report, path, tree, read_blob, repository):
    if stop in POINTINGS:
        frame, label, pattern = POINTINGS[stop]
        if frame not in report:
            return None
        match = re.search(pattern, report, re.M)
        if match:
            line = report.count("\n", 0, match.start()) + 1
            if match[0].startswith("|"):
                # The frame is defined in the column header of this table.
                # Other tables in the same paper can use a different frame.
                rows = report.splitlines()
                header = line - 1
                while header > 0 and rows[header - 1].lstrip().startswith("|"):
                    header -= 1
                if frame not in rows[header]:
                    return None
            return point(*match.groups(), frame, label, receipt(repository, tree, path, label, line))
    if stop == "archival-radio-supernovae" and "Coordinates are ICRS degrees" in report:
        table = "results/archival-radio-supernovae/followup-visits.csv"
        if table in tree:
            rows = csv.DictReader(io.StringIO(read_blob(table).decode("utf-8-sig")))
            positions = {(row["ra"], row["dec"]) for row in rows if row["name"] == "SN2019njv"}
            if len(positions) == 1:
                label = "SN 2019njv · Optical position in VLASS follow-up"
                return point(*positions.pop(), "ICRS", label, receipt(repository, tree, table, label))
    if stop == "quasar":
        plates = "observatory/assets/plates.json"
        wcs = "observatory/assets/j0826-r-wcs.fits"
        if plates in tree and wcs in tree:
            record = json.loads(read_blob(plates)).get("quasar", {})
            query = parse_qs(urlsplit(record.get("archiveUrl", "")).query)
            astrometry = record.get("astrometry", {})
            raw = read_blob(wcs)
            if hashlib.sha256(raw).hexdigest() == astrometry.get("sha256") and "ra" in query and "dec" in query:
                frame = fits_header(raw)[0].get("RADESYS")
                if frame in {"FK5", "ICRS"}:
                    label = "Pan-STARRS1 cutout center · between images A and B"
                    return point(query["ra"][0], query["dec"][0], frame, label,
                                 receipt(repository, tree, plates, "PS1 pointing and FITS WCS receipt"))
    return None


def build_destinations(tree: RepoTree, records, spotlights, repository, output: Path, read_blob):
    by_record = {item["recordId"]: item for item in spotlights}
    destinations = []
    out = output.parent / "flight"
    for record in records:
        spot = by_record.get(record["id"])
        path = record["path"]
        stop = record["id"].removeprefix("results-")
        anchor = unquote(urlsplit(record["url"]).fragment)
        entry = STOPS.get(stop)
        item = {"id": spot["id"] if spot else record["id"], "recordId": record["id"],
                "name": record["title"], "kind": "study", "setting": record["topic"], "featured": bool(spot or entry)}
        destinations.append(item)
        if spot:
            kind, setting = SPOTLIGHTS.get(spot["id"], ("study", record["topic"]))
            item.update(name=spot.get("target", record["title"]), kind=kind, setting=setting,
                        figure={**spot["figure"], "type": "research-figure"})
            for field in ("observingUnit", "limitation", "measurements"):
                if field in spot:
                    item[field] = spot[field]
            coords = pointing(spot["id"], "", "", tree, read_blob, repository)
            if coords:
                item["coordinates"] = coords
            continue  # Full spotlight tables stay in catalog.spotlights.
        if entry:
            item.update(zip(("name", "kind", "setting", "observingUnit"), entry))
        if stop == "jwst-repeat-imaging-transients":
            if GOODS in tree:
                path = GOODS
            else:
                item.update(name=record["title"], kind="study", setting=record["topic"], featured=False)
                item.pop("observingUnit", None)
        if path not in tree and path + "/README.md" in tree:
            path += "/README.md"
        if path not in tree or not path.endswith(".md"):
            continue
        document = tree.text(path)
        section, line = scope(document, anchor)
        source = lambda label, n: receipt(repository, tree, path, label, n)
        table = report_table(section, line, source)
        if table:
            item["table"] = table
        excerpt = report_excerpt(section, line, record["summary"])
        if excerpt:
            item["excerpt"] = excerpt[0]
            item["excerptSource"] = source("Investigation report", excerpt[1])
        position = pointing(stop, document, path, tree, read_blob, repository)
        if position:
            item["coordinates"] = position
        if path.startswith(("results/", "docs/")):
            origin = PurePosixPath(path)
            prefix = "/".join(origin.parts[:2]) + "/" if origin.parts[0] == "results" else "docs/"
            prose, _ = markdown.split_code(section)
            for link in markdown.links(prose):
                if not link.image or link.external or not link.text:
                    continue
                target = posixpath.normpath(posixpath.join(str(origin.parent), unquote(link.target)))
                if not target.startswith(prefix):
                    continue
                figure = publish_image(tree, read_blob, out, repository, target, plain(link.text))
                if figure:
                    item["figure"] = figure
                    break
        if stop == "rgz-desi-giant-radio-galaxies":
            image = publish_image(tree, read_blob, out, repository,
                "results/rgz-desi-giant-radio-galaxies/inputs/J000314_optical.jpg",
                "Legacy Surveys DR10 optical cutout centered on the DESI fiber in J000314.9+111544.",
                "archive-image", "DESI Legacy Imaging Surveys")
            if image:
                item["image"] = {**image, "fieldLabel": "Optical host field · Centered on DESI fiber"}
        if stop == "mars-crater-size-correlation":
            image = publish_image(tree, read_blob, out, repository,
                "results/mars-crater-size-correlation/pixel-check/field.jpg",
                "HiRISE red-channel crop of ESP_066902_1830 with 76 catalog crater positions circled.",
                "archive-image", "NASA/JPL-Caltech/University of Arizona")
            if image:
                item["image"] = {**image, "fieldLabel": "76-crater impact field · HiRISE red channel"}
        if stop in FITS_PLATES:
            filename, credit, field_label, caption = FITS_PLATES[stop]
            plate_path = f"results/{stop}/{filename}"
            if plate_path in tree:
                raw = read_blob(plate_path)
                png, epoch = render_plate(raw)
                caption += f" Observed {epoch[:10]} UTC. Native pixels, linear grayscale from the 1st to 99.5th percentile."
                image = package_image(tree, out, repository, plate_path, caption, "archive-image", credit, png, ".png")
                item["image"] = {**image, "fieldLabel": field_label, "inputSha256": hashlib.sha256(raw).hexdigest()}
    return destinations

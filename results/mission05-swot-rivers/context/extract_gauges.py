"""Extract current-period gauge values from saved ZRA HTML bulletins.

Pass the directory of .html files, then the output directory. To reacquire
individual source pages, see the observed URLs in zra_bulletin_sources.csv.
Requires lxml; Python 3.10+. Local report dates use UTC+02:00.
"""
import argparse
import calendar
import csv
import hashlib
import re
from collections import defaultdict
from datetime import date
from pathlib import Path

from lxml import html


def write_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as dest:
        writer = csv.writer(dest, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def main(raw, output):
    output.mkdir(parents=True, exist_ok=True)
    observations = defaultdict(dict)
    documents = []
    skipped = []
    for path in sorted(raw.glob("*.html")):
        data = path.read_bytes()
        document = html.fromstring(data)
        articles = document.xpath('//article[starts-with(@id,"node-riverflows-")]')
        if len(articles) != 1:
            continue
        article = articles[0]
        ident = article.get("id").rsplit("-", 1)[-1]
        url = document.xpath('//link[@rel="canonical"]/@href')[0]
        title = document.xpath('string(//title)').lower()
        station = "Chavuma" if "chavuma" in title else "Nanas Farm" if "nana" in title else ""
        if not station:
            continue
        dates = article.xpath('.//*[contains(@class,"date-display")]/@content')
        if len(dates) != 2:
            skipped.append((ident, "missing reporting period", ""))
            continue
        start, end = [date.fromisoformat(item[:10]) for item in dates]
        docs_key = "zra_" + ident
        if end < date(2023, 3, 29) or start > date(2025, 5, 2):
            continue
        documents.append((docs_key, station, start, end, url, hashlib.sha256(data).hexdigest()))
        rows = article.xpath('.//div[contains(@class,"field-name-field-riverflow-data")]//table//tr')
        for row in rows:
            values = [value.text_content().strip() for value in row.xpath('./td|./th')]
            if not values or not values[0] or values[0].lower() == 'date':
                continue
            match = re.fullmatch(r"(\d{1,2})\s*[-/ .]?\s*([A-Za-z]{3,9})(?:[-/ ](?:\d{2}|\d{4}))?", values[0])
            if len(values) < 2 or not match:
                skipped.append((ident, "unparsed row", "|".join(values)))
                continue
            try:
                day, month = int(match[1]), list(calendar.month_abbr).index(match[2][:3].title())
                candidates = []
                for year in sorted({start.year, end.year}):
                    try:
                        when = date(year, month, day)
                    except ValueError:
                        continue
                    if start <= when <= end:
                        candidates.append(when)
                if len(candidates) != 1:
                    raise ValueError("date outside reporting period")
                when = candidates[0]
                token = values[1].replace(",", "").replace(" ", "")
                if not re.fullmatch(r"\d+(?:\.\d+)?",token):
                    raise ValueError("unreported or unparsed flow")
                flow = float(token)
            except ValueError as err:
                skipped.append((ident, str(err), "|".join(values)))
                continue
            if date(2023, 3, 29) <= when <= date(2025, 5, 2):
                observations[(station,when)][docs_key] = flow
    records = []
    conflicts = []
    for (station, day), by_doc in sorted(observations.items()):
        if len(set(by_doc.values())) == 1:
            records.append((station, day, f"{next(iter(by_doc.values())):g}", "|".join(sorted(by_doc))))
        else:
            conflicts.extend((station,day,source,value) for source,value in sorted(by_doc.items()))
    write_csv(output/"zra_daily_gauges.csv", ["station","date","q_m3s","source_ids"], records)
    write_csv(output/"zra_bulletin_sources.csv", ["source_id","station","period_start","period_end","original_url","retrieved_html_sha256"], sorted(documents))
    write_csv(output/"zra_excluded_conflicting_values.csv", ["station","date","source_id","q_m3s"], conflicts)
    write_csv(output/"zra_excluded_source_rows.csv", ["node_id","reason","row"], skipped)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    main(args.raw, args.output)

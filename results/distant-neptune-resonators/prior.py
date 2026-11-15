import csv, re
from pathlib import Path

w = Path(__file__).resolve().parent / "work"
alphabet = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz"


def ids(s):
    m = re.search(r"(\d{4}) ([A-Z])([A-Z])(\d+)?", s)
    number = re.match(r"^(\d+) ", s)
    keys = []
    if number:
        keys.append(number[1])
    if m:
        yr, a, b, num = m.groups()
        num = int(num or 0)
        keys.append("K" + yr[2:] + a + alphabet[num // 10] + str(num % 10) + b)
    return keys


if __name__ == "__main__":
    prior = {
        r[0]: r[1:5]
        for r in list(csv.reader(open(w / "volk.csv", encoding="utf-8-sig")))[1:]
    }
    for r in sorted(
        csv.DictReader(open(w / "screen.csv")), key=lambda x: float(x["arc_deg"])
    ):
        if float(r["arc_deg"]) > 350:
            continue
        print(
            r["designation"],
            r["p"] + "/" + r["q"],
            [(i, prior.get(i)) for i in ids(r["designation"])],
        )

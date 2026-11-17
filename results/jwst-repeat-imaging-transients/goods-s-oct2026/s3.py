import urllib.parse as p, urllib.request as u, sys, xml.etree.ElementTree as ET

base = "https://stpubdata.s3.amazonaws.com/"


def keys(pre, delimiter=None, limit=1000):
    a = {"list-type": 2, "prefix": pre, "max-keys": limit}
    if delimiter:
        a["delimiter"] = delimiter
    for n in range(100):
        tree = ET.fromstring(u.urlopen(base + "?" + p.urlencode(a), timeout=45).read())
        ns = "{*}"
        for t in tree.findall(ns + "Contents"):
            yield t.findtext(ns + "Key"), int(t.findtext(ns + "Size"))
        for t in tree.findall(ns + "CommonPrefixes"):
            yield t.findtext(ns + "Prefix"), 0
        nt = tree.findtext(ns + "NextContinuationToken")
        if not nt:
            break
        a["continuation-token"] = nt


if __name__ == "__main__":
    for arg in sys.argv[1:]:
        print(list(keys(arg, "/"))[:45])

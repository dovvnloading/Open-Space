from query import *


def run(arg):
    name, a, d = arg
    q = f"SELECT ra,dec FROM ls_dr10.tractor WHERE ra BETWEEN {a-1.55} AND {a+1.55} AND dec BETWEEN {d-1.47} AND {d+1.47} AND flux_r>0.65 AND flux_ivar_r*flux_r*flux_r>25"
    url = "https://datalab.noirlab.edu/tap/sync?" + urlencode(
        dict(REQUEST="doQuery", LANG="ADQL", FORMAT="csv", MAXREC=600000, QUERY=q)
    )
    acquire("stationary_" + name, url)


if __name__ == "__main__":
    with ThreadPoolExecutor(2) as e:
        list(e.map(run, [("Chen1", 42.676, -15.021), ("Chen2", 45.297, -16.711)]))

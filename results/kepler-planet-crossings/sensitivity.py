from replay import *
import replay

sol = pd.read_csv(HERE / "isolated_fits.csv")
rows = []
for label in ["ShortCadence", "Block0", "Block1", "Block2", "Block3", "LDweak", "LDstrong"]:
    if label == "ShortCadence":
        geo = GEOMETRIES[label]
        u = (0.5256, 0.186)
    elif label.startswith("Block"):
        geo = isolated_geometry(sol, "c_" + label[-1], "d_" + label[-1])
        u = (0.5256, 0.186)
    else:
        u = (0.43, 0.23) if label == "LDweak" else (0.62, 0.12)
        iso = sol[(sol.u1 == u[0]) & (sol.u2 == u[1])].set_index("planet")
        geo = np.array(
            [
                *[
                    iso.loc[x, ["impact", "radius", "duration_hours"]].to_numpy(float)
                    for x in ["c", "d"]
                ],
                GEOMETRIES["DR25"][2],
            ]
        )
    replay.U1, replay.U2 = u
    replay.NORMALIZATION = 1 - u[0] / 3 - u[1] / 6
    ev = build_events(geo, "pdc", "best")
    for s, t, p, c in scan(geo, ev, np.arange(-180, 180, 2)):
        rows.append(
            dict(
                solution=label, sign_d=s, angle_deg=t, expected_delta_chi2=p, observed_delta_chi2=c
            )
        )
df = pd.DataFrame(rows)
df.to_csv(HERE / "sensitivity.csv", index=False, float_format="%.4f")
for sign in [1, -1]:
    p = df[df.sign_d == sign].pivot(
        index="angle_deg", columns="solution", values="observed_delta_chi2"
    )
    for cols in [list(p), ["ShortCadence", "LDstrong", "LDweak"], ["ShortCadence"]]:
        a = p[cols].min(axis=1)
        print(sign, ",".join(cols), a.max(), a.idxmax(), a[a > 25].index.tolist())

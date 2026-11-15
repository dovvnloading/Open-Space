"""Link Find_Orb's native radar projector using the already-built force adapter."""
import argparse
from pathlib import Path
import subprocess


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("build", type=Path, help="Force-enabled orbit packet build directory")
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    build = args.build.resolve()
    source = build / "find_orb"
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    text = (source / "elem_out.cpp").read_text()
    if "FO_ACCEL" not in text:
        raise RuntimeError("Rebuild Find_Orb using the force-enabled orbit packet")
    marker = "               set_locs( orbit, *orbit_epoch, obs, n_obs);"
    start = text.index('if (getenv("FO_ACCEL"))')
    block = text.index(marker, start)
    text = text[:block] + '''               for (int j = 0; j < n_obs; ++j)
                  obs[j].note2 = 'R';  /* Use radar's heliocentric light-time correction. */
''' + text[block:]
    start = text.index("                  const double dx = remainder(", block)
    stop = text.index("               }", start)
    text = text[:start] + '''                  char record[162];
                  memset(record, ' ', sizeof(record));
                  record[80] = record[161] = 0;
                  memcpy(record + 68, "253", 3);  /* Transmitter DSS-14. */
                  memcpy(record + 81 + 62, "085600", 6);  /* 8560 MHz. */
                  record[42] = '1';  /* Dummy positive sigma required by the projector. */
                  obs[j].second_line = record;
                  RADAR_INFO rinfo;
                  compute_radar_info(obs + j, &rinfo);
                  obs[j].second_line = NULL;
                  printf(" %.14g %.14g", rinfo.rtt_comp, rinfo.doppler_comp);
''' + text[stop:]
    patched = output / "radar_projection.cpp"
    patched.write_text(text)
    include = build / "prefix" / "include"
    lib = build / "prefix" / "lib"
    compiled = output / "radar_projection.o"
    subprocess.run(["g++", "-O3", "-std=c++17", "-I" + str(include), "-I" + str(source),
                    "-c", str(patched), "-o", str(compiled)], check=True)
    # Reuse the precise object list from the pinned upstream linker rule.
    make = subprocess.run(["make", "-Bn", "fo"], cwd=source, capture_output=True,
                          text=True, check=True)
    line = next(l for l in reversed(make.stdout.splitlines()) if " -o fo " in l)
    objects = [str(source / v) for v in line.split() if v.endswith(".o") and v != "elem_out.o"]
    subprocess.run(["g++", "-o", str(output / "fo"), *objects, str(compiled), "-L" + str(lib),
                    "-lm", "-llunar", "-ljpl", "-lsatell"], check=True)


if __name__ == "__main__":
    main()

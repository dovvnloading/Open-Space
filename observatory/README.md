# Open Space Observatory

The [public observatory](https://dovvnloading.github.io/Open-Space/) pairs an evidence archive with a procedural universe. On the observation deck, switch between a target's context image and its measured data. Use **Inspect image** to pan through the original plate, play eight hourly coronagraph frames, or locate the two quasar images. The featured plots have an observation cursor and a switch to the original published figure. The four investigation files show the observing unit, readout methods, source figures, browsable observation rows, limitations and next test. [Archive checks](https://dovvnloading.github.io/Open-Space/#audit) describe which checks cover the whole repository.

## Run locally

From a checkout, serve the bundled snapshot with the Python standard library.

```sh
python -m http.server 8000 --bind 127.0.0.1 --directory observatory
```

Open [the local observatory](http://127.0.0.1:8000/). To include current verification coverage in a fresh local snapshot, run both commands from the same committed revision. For a fork, pass `--repository owner/repo` to the generator. A full clone supplies first-parent publication dates for recent activity.

```sh
python -m open_space.verify --json runs/observatory-ledger.json --quiet
python observatory/build_catalog.py --verification-ledger runs/observatory-ledger.json
```

## Evidence and publication

The [Pages workflow](../.github/workflows/observatory.yml) builds a snapshot on relevant pushes to main. It verifies the committed bytes, checks the featured readouts against their source tables, validates the context images, and exercises every record’s displayed title, disposition and evidence links in Chromium, Firefox and WebKit. It also compares the rendered observation tables against all their published rows. Each classified record points to the exact README line that supplied its disposition. Additional research directories retain an unclassified research-record label. Folder checks cover the whole linked folder, which can contain several findings. The site keeps unresolved physical interpretations visible.

The curated target descriptions and figure sources live in [spotlights.json](data/spotlights.json). The MAVEN and STIX plots embed their source CSV digests. If the underlying tables change, inspect the new measurements, regenerate the plots using [render_figures.py](data/render_figures.py) with Matplotlib installed, and review the dossier text before publishing.

## Universe simulator

[![The Vesper Spiral in the running universe simulator](../docs/figures/observatory-universe.jpg)](https://dovvnloading.github.io/Open-Space/#universe)

The simulator generates 128 galaxies, each with 32,768 addressable systems. The destination atlas links directly to selected galaxies and worlds, and the sector map distinguishes five generated galaxy forms. Enter a galaxy, system or world using the scale controls, or travel to a system by address. Drag or use the arrow keys to pan and use the zoom controls, wheel or pinch to change scale. Schematic spacing is the default. **AU scale** shows relative orbital radii with a calibrated ruler. Advance a quarter of the selected orbit to inspect motion. The model uses circular orbits around one star. Body sizes, lighting and cloud patterns are illustrative.

## Checks

```sh
python -m unittest discover -s tests -p 'test_observatory_*.py' -v
node --test tests/observatory-simulation.test.mjs tests/observatory-viewer.test.mjs
python observatory/check_media.py
python -m pip install playwright==1.58.0
python -m playwright install chromium firefox webkit
python observatory/check_browser.py --directory observatory
python observatory/check_browser.py --directory observatory --engine firefox
python observatory/check_browser.py --directory observatory --engine webkit
```

The [featured-source replay](data/spotlight-replay.md) records which original data were recovered and where archived inputs were used. The site bundles its plates, figures, fonts and code. [Image sources and reuse terms](assets/CREDITS.md) · [Instrument Serif license](assets/fonts/instrument-serif-LICENSE.txt) · [DM Sans license](assets/fonts/dm-sans-LICENSE.txt) · [IBM Plex Mono license](assets/fonts/plex-mono-LICENSE.txt).

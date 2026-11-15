# Camera images for the 24 August 2026 quartet

All 18 distinct camera IDs in the six trajectory rows have a public night summary with a time-labelled detected-image thumbnail consistent with the associated trajectory time. The archive uses the night that began on 23 August. At least two cameras per geographic track show a visible streak in the corresponding time block. These are the original contributing instruments, so the images substantiate the presence of transients and allow visual checks for camera artifacts. They do not provide independent astrometry, subsecond arrival times or a second observing network.

| Track | Geographic region of meteor | Earliest trajectory start UTC across reductions | Cameras with especially legible features (block label UTC) |
| --- | --- | --- | --- |
| A (1) | Central Germany | 00:36:08.325 | CZ0006 (00:36:00), CZ000X (00:36:01) |
| B (2) | Northeastern Germany | 00:36:08.662 | DE000Y (00:36:07), DE0018 (00:36:01) |
| C (3) | Denmark | 00:36:17.035 | DK0007 (00:36:14), DK000A (00:36:14), DK000B (00:36:15) |
| D (4) | Northwestern Germany | 00:36:17.269 | DE000X (00:36:15), NL000C (00:36:12), NL0009 (00:36:07) |

The camera manifest gives each original SHA-256, source URL and zero-based row and column. Approximate feature positions are human-selected clues for reading a 320 × 180 image. The x coordinate increases rightward and y increases downward from the top-left, excluding the time label. Multiple transients can share a block. Target position and line-of-sight identification have not been checked against the astrometric calibration. The CZ000J thumbnail is included for completeness and has no unambiguous visible streak at this resolution.

RMS produces thumbnails from per-pixel maxima, writing the timestamp of the first four-frame-format (FF) file in the stack to whole seconds. A standard FF block comprises 256 frames. Nearby unrelated detections and fixed background structures can appear in these summaries. The source links retain the original public images. The inspection script displays the candidate, the median of the two temporally adjacent detected tiles on either side and their difference. Those comparison tiles are from the same field and may be minutes apart. See [RMS thumbnail generation](https://github.com/CroatianMeteorNetwork/RMS/blob/ae2f34f8e4a65bba8596dec2a9803a0cb4c2d40d/Utils/GenerateThumbnails.py) and [GMN observing methodology](https://doi.org/10.1093/mnras/stab2008).

Run from this directory with Python, Pillow and NumPy.

```sh
python inspect_tiles.py --cache downloaded --output inspected
```

Images are credited to the individual camera operators and Global Meteor Network. The repository distributes the manifest and the inspection method. The camera operators publish their images on the linked station pages.

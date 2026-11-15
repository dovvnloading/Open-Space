# Camera images for the 24 August 2026 quartet

All 18 distinct camera IDs in the six trajectory rows have a public night summary with a time-labelled detected-image thumbnail consistent with the associated trajectory time. The archive uses the night that began on 23 August. At least two cameras per geographic track show a visible streak in the corresponding time block. These are the original contributing instruments, so the images substantiate the presence of transients and allow visual checks for camera artifacts. They do not provide independent astrometry, subsecond arrival times or a second observing network.

| Track | Geographic region of meteor | Earliest trajectory start UTC across reductions | Cameras with especially legible features (block label UTC) |
| --- | --- | --- | --- |
| A (1) | Central Germany | 00:36:08.325 | CZ0006 (00:36:00), DK000J (00:36:05) |
| B (2) | Northeastern Germany | 00:36:08.662 | DE000Y (00:36:07), DE0018 (00:36:01) |
| C (3) | Denmark | 00:36:17.035 | DK0007 (00:36:14), DK000A (00:36:14), DK000F (00:36:11) |
| D (4) | Northwestern Germany | 00:36:17.269 | DE000X (00:36:15), NL000C (00:36:12), NL0009 (00:36:07) |

The camera manifest gives each original SHA-256, source URL and zero-based row and column. Approximate feature positions are human-selected clues for reading a 320 × 180 image. The x coordinate increases rightward and y increases downward from the top-left, excluding the time label. Multiple transients can share a block. The approximate pointing comparison below checks plausible feature positions. The CZ000J thumbnail is included for completeness and has no unambiguous visible streak at this resolution.

RMS produces thumbnails from per-pixel maxima, writing the timestamp of the first four-frame-format (FF) file in the stack to whole seconds. A standard FF block comprises 256 frames. Nearby unrelated detections and fixed background structures can appear in these summaries. The source links retain the original public images. The inspection script displays the candidate, the median of the two temporally adjacent detected tiles on either side and their difference. Those comparison tiles are from the same field and may be minutes apart. See [RMS thumbnail generation](https://github.com/CroatianMeteorNetwork/RMS/blob/ae2f34f8e4a65bba8596dec2a9803a0cb4c2d40d/Utils/GenerateThumbnails.py) and [GMN observing methodology](https://doi.org/10.1093/mnras/stab2008).

Run from this directory with Python, Pillow and NumPy.

```sh
python inspect_tiles.py --cache downloaded --output inspected
```

Images are credited to the individual camera operators and Global Meteor Network. The repository distributes the manifest and the inspection method. The camera operators publish their images on the linked station pages.

## Approximate image-plane check

For each camera, the nightly astrometry report prints central azimuth, altitude, rotation with respect to the horizon, and horizontal and vertical field of view. We transcribed those values in [pointing.csv](pointing.csv). The station position was estimated by intersecting sightlines from corresponding vertices in the publicly available 25 km and 100 km KML footprint outlines. Site longitudes and latitudes are rounded to 0.01° here. The KML is a later public footprint, so assumes the station remained at approximately the same site. The camera orientation is from the observing night itself.

Projecting the GMN begin and end coordinates into local east/north/up and applying a centered equidistant lens gives an approximate segment per camera. It omits lens-specific distortion, calibration residuals and the individual observation frames. The measured feature midpoint and orientation were identified from the images, so agreement cannot validate the catalog orbit independently.

For the two-camera northeastern Germany event, DE000Y shows a steep trail near (117, 62). DE0018 shows a nearly horizontal one near (273, 143). Both share the predicted orientation, with offsets of 3.0 and 11.4 pixels from the respective approximate segments. CZ000X and DK000B have conspicuous trails 69.9 and 64.5 pixels from their respective projected segments. They are unresolved contributors and do not support the visual count. CZ000J lacks a clear trail in the thumbnail. For the other 15 cameras, the indicated feature is within 11.4 pixels of its corresponding projected segment (14 within 6.1 pixels). Exact offsets from the rounded inputs can be obtained with the command below. The coordinates describe 320 × 180 pixel crops.

```sh
python project_tracks.py
```

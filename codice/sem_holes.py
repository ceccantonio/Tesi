"""Mean radius and period (with standard deviation) of the holes of a membrane, from its SEM photos.

    python sem_holes.py C4mem11            # the 4 photos C4mem11_00.tif ... C4mem11_03.tif of the folder FOLDER, pooled together

The scale (nm per pixel) is read from the 1 um bar of each photo. A hole is a dark spot; the holes cut by the border of the photo are ignored.
The radius is the one of the circle with the same area as the hole; the period is the distance between a hole and its nearest neighbour.
Holes that are not round, or whose size is far from the others (merged holes, dirt), are discarded.
It also saves <membrane>_check.png with the holes found, to verify the detection by eye.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from PIL import Image
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

FOLDER = Path(__file__).resolve().parent.parent / "INFRACHIP_v2" / "membranesceltesem"
BAR_UM = 1.0         # length of the scale bar written in the photo (um)
FOOTER = 713         # the data bar of the SEM starts at this row: it is cut away
THRESHOLD = 8        # a pixel belongs to a hole (dark inside, bright rim) if it differs from its surroundings by this many grey levels (0-255)
MIN_AREA = 300       # pixels: smaller dark spots are not holes
MAX_DEVIATION = 0.25 # a hole whose radius differs from the median of its photo by more than this fraction (0.25 = 25%) is discarded


def scale_nm_per_px(image):
    """Length in pixels of the longest horizontal white line of the data bar (the scale bar), right of the logo."""
    longest = 0
    for row in image[FOOTER:, 250:] >= 250:
        edges = np.diff(np.concatenate(([0], row.astype(int), [0])))
        runs = np.flatnonzero(edges == -1) - np.flatnonzero(edges == 1)
        longest = max(longest, runs.max(initial=0))
    return BAR_UM * 1000 / longest


def holes(path):
    """(radii in nm, centres in nm, nm per pixel, picture, mask of the holes found) of one photo."""
    image = np.array(Image.open(path).convert("L"), float)
    nm = scale_nm_per_px(image)
    picture = image[:FOOTER]
    contrast = ndi.gaussian_filter(picture, 1.5) - ndi.gaussian_filter(picture, 25)
    edge = ndi.binary_closing(abs(contrast) > THRESHOLD, iterations=4)             # dark inside + bright rim of the hole
    labels, n = ndi.label(ndi.binary_opening(ndi.binary_fill_holes(edge), iterations=3))
    ids, centres, radii = [], [], []
    for i, box in enumerate(ndi.find_objects(labels), 1):
        area = (labels[box] == i).sum()
        at_border = box[0].start == 0 or box[1].start == 0 or box[0].stop == FOOTER or box[1].stop == picture.shape[1]
        height, width = box[0].stop - box[0].start, box[1].stop - box[1].start
        round_enough = 0.7 < width / height < 1.4 and area / (width * height) > 0.65         # a circle fills 78% of its bounding box
        if area >= MIN_AREA and round_enough and not at_border:
            ids.append(i)
            radii.append(np.sqrt(area / np.pi) * nm)
            centres.append(np.array(ndi.center_of_mass(labels == i)) * nm)
    radii, centres, ids = np.array(radii), np.array(centres), np.array(ids)
    good = abs(radii / np.median(radii) - 1) < MAX_DEVIATION                       # merged holes and dirt are much bigger or smaller than the rest
    return radii[good], centres[good], nm, picture, np.isin(labels, ids[good])


def main(membrane):
    paths = sorted(FOLDER.glob(f"{membrane}_??.tif"))
    if not paths:
        sys.exit(f"no photos {membrane}_??.tif in {FOLDER}")
    all_radii, all_periods = [], []
    fig, axes = plt.subplots(1, len(paths), figsize=(5 * len(paths), 4), squeeze=False)
    for ax, path in zip(axes[0], paths):
        radii, centres, nm, picture, found = holes(path)
        periods = cKDTree(centres).query(centres, k=2)[0][:, 1] if len(centres) > 1 else np.array([])
        periods = periods[abs(periods / np.median(periods) - 1) < MAX_DEVIATION]       # a hole next to a missed or false one has a wrong neighbour distance
        all_radii += list(radii)
        all_periods += list(periods)
        print(f"{path.name}: {len(radii)} holes, radius {radii.mean():.0f} +- {radii.std():.0f} nm, period {periods.mean():.0f} +- {periods.std():.0f} nm ({nm:.2f} nm/px)")
        ax.imshow(picture, cmap="gray")
        ax.contour(found, [0.5], colors="r", linewidths=0.8)
        ax.set_title(path.name, fontsize=9)
        ax.axis("off")
    r, p = np.array(all_radii), np.array(all_periods)
    print(f"\n{membrane}: {len(r)} holes in {len(paths)} photos\n  radius = {r.mean():.0f} +- {r.std():.0f} nm\n  period = {p.mean():.0f} +- {p.std():.0f} nm")
    plt.tight_layout()
    plt.savefig(FOLDER / f"{membrane}_check.png", dpi=110)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "C4mem11")

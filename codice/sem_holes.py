"""Mean radius and period (with standard deviation) of the holes of a membrane, from its SEM photos.

    python sem_holes.py                    # all the membranes of the folder FOLDER: results in sem_holes_results.txt
    python sem_holes.py C4mem11            # only one: its 4 photos C4mem11_00.tif ... C4mem11_03.tif are pooled together

The scale (nm per pixel) is read from the 1 um bar of each photo. A hole is a dark spot; the holes cut by the border of the photo are ignored.
The radius is the one of the circle with the same area as the hole; the period is the distance between neighbouring holes (all the pairs); it is also measured with the 2D FFT of the photo (fft_lattice), which does not need the holes.
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
THRESHOLD = 6       # a pixel belongs to a hole (dark inside, bright rim) if it differs from its surroundings by this many grey levels (0-255)
MIN_AREA = 300       # pixels: smaller dark spots are not holes
MAX_DEVIATION = 0.25 # a hole whose radius differs from the median of its photo by more than this fraction (0.25 = 25%) is discarded
FFT_PAD = 4096       # the FFT of the photo is computed on a grid of this size (zero padding: finer sampling of the peaks)


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


def neighbour_distances(centres):
    """Distances between neighbouring holes (all the pairs of neighbours, each once): the period. Pairs that are not one period apart are discarded:
    diagonals (1.4 periods) and the pairs across a hole that was not detected (2 periods). The distance to the NEAREST neighbour only would
    underestimate the period and its standard deviation (the minimum of four noisy distances)."""
    if len(centres) < 2:
        return np.array([])
    tree = cKDTree(centres)
    pitch = np.median(tree.query(centres, k=2)[0][:, 1])                        # first guess of the period
    pairs = tree.query_pairs((1 + MAX_DEVIATION) * pitch, output_type="ndarray")
    distances = np.linalg.norm(centres[pairs[:, 0]] - centres[pairs[:, 1]], axis=1)
    return distances[distances > (1 - MAX_DEVIATION) * pitch]


def fft_lattice(picture, nm):
    """Lattice of the holes from the 2D FFT of the photo, without finding the holes: (|a1|, |a2|, angle between a1 and a2 in degrees), in nm.
    The photo (without the slow shading, with a Hann window) is transformed. The peaks of the spectrum are at the vectors g (cycles per pixel) of the
    reciprocal lattice: g1 is the strong peak of lowest frequency and g2 the next one that is not along the same direction (so that the diagonal peaks,
    at higher frequency, are not taken). Both are refined with the centroid of the peak. The real-space lattice vectors satisfy ai.gj = delta_ij.
    It measures the period of the whole photo at once, so it does not depend on the holes found or on the threshold."""
    image = picture - ndi.gaussian_filter(picture, 25)
    image = image * np.outer(np.hanning(image.shape[0]), np.hanning(image.shape[1]))
    power = np.abs(np.fft.fftshift(np.fft.fft2(image, s=(FFT_PAD, FFT_PAD)))) ** 2
    k = np.fft.fftshift(np.fft.fftfreq(FFT_PAD))                                 # cycles per pixel
    step = k[1] - k[0]
    kx, ky = np.meshgrid(k, k)                                                   # power[row, column] = power(ky, kx)
    power[np.hypot(kx, ky) < 0.004] = 0                                          # no DC and no slow structures (periods longer than 250 px)

    strong = (power == ndi.maximum_filter(power, size=15)) & (power > 0.02 * power.max())      # local maxima, not too weak
    rows, cols = np.nonzero(strong)
    order = np.argsort(np.hypot(kx[rows, cols], ky[rows, cols]))                 # lowest frequency first
    peaks = []
    for r, c in zip(rows[order], cols[order]):
        window = power[r - 2:r + 3, c - 2:c + 3]                                 # centroid of the peak
        w = window / window.sum()
        offsets = np.arange(-2, 3)
        g = np.array([kx[r, c] + (w.sum(axis=0) * offsets).sum() * step, ky[r, c] + (w.sum(axis=1) * offsets).sum() * step])
        if not peaks:
            peaks.append(g)
        elif abs(g[0] * peaks[0][1] - g[1] * peaks[0][0]) / (np.hypot(*g) * np.hypot(*peaks[0])) > np.sin(np.radians(25)):    # not along g1
            peaks.append(g)
            break
    g1, g2 = peaks
    for _ in range(10):                                                          # shortest basis of the lattice (Lagrange-Gauss): if g2 is a diagonal peak
        if np.linalg.norm(g1) > np.linalg.norm(g2):                              # (g1 + g2'), g2 - g1 is the fundamental one
            g1, g2 = g2, g1
        m = np.round(g1 @ g2 / (g1 @ g1))
        if m == 0:
            break
        g2 = g2 - m * g1
    a = np.linalg.inv(np.array([g1, g2])).T * nm                                 # rows: a1, a2 in nm
    cos = a[0] @ a[1] / (np.linalg.norm(a[0]) * np.linalg.norm(a[1]))
    return np.linalg.norm(a[0]), np.linalg.norm(a[1]), np.degrees(np.arccos(abs(cos)))


def main(membrane):
    paths = sorted(FOLDER.glob(f"{membrane}_??.tif"))
    if not paths:
        sys.exit(f"no photos {membrane}_??.tif in {FOLDER}")
    all_radii, all_periods, fft_periods = [], [], []
    fig, axes = plt.subplots(1, len(paths), figsize=(5 * len(paths), 4), squeeze=False)
    for ax, path in zip(axes[0], paths):
        radii, centres, nm, picture, found = holes(path)
        periods = neighbour_distances(centres)
        a1, a2, angle = fft_lattice(picture, nm)
        fft_periods.append((a1 + a2) / 2)
        all_radii += list(radii)
        all_periods += list(periods)
        print(f"{path.name}: {len(radii)} holes, radius {radii.mean():.0f} +- {radii.std():.0f} nm, period {periods.mean():.0f} +- {periods.std():.0f} nm ({nm:.2f} nm/px) | FFT: a1 {a1:.0f}, a2 {a2:.0f} nm, angle {angle:.0f} deg")
        ax.imshow(picture, cmap="gray")
        ax.contour(found, [0.5], colors="r", linewidths=0.8)
        ax.set_title(path.name, fontsize=9)
        ax.axis("off")
    r, p = np.array(all_radii), np.array(all_periods)
    print(f"{membrane}: {len(r)} holes in {len(paths)} photos, radius = {r.mean():.0f} +- {r.std():.0f} nm, period = {p.mean():.0f} +- {p.std():.0f} nm, period from the FFT = {np.mean(fft_periods):.0f} nm\n")
    plt.tight_layout()
    plt.savefig(FOLDER / f"{membrane}_check.png", dpi=110)
    plt.close()
    return len(paths), len(r), r.mean(), r.std(), p.mean(), p.std(), np.mean(fft_periods)


if __name__ == "__main__":
    membranes = sys.argv[1:] or sorted({path.name[:-7] for path in FOLDER.glob("*_??.tif")})       # C4mem11_00.tif -> C4mem11
    lines = [f"{'membrane':<12}{'photos':>7}{'holes':>7}{'radius_nm':>11}{'radius_std':>11}{'period_nm':>11}{'period_std':>11}{'period_fft':>11}"]
    for membrane in membranes:
        n_photos, n_holes, r_mean, r_std, p_mean, p_std, p_fft = main(membrane)
        lines.append(f"{membrane:<12}{n_photos:>7}{n_holes:>7}{r_mean:>11.0f}{r_std:>11.0f}{p_mean:>11.0f}{p_std:>11.0f}{p_fft:>11.0f}")
    (FOLDER / "sem_holes_results.txt").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines), f"\nSaved {FOLDER / 'sem_holes_results.txt'}")

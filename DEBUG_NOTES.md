# leaf_area.R — what was broken and what changed

Diagnosed against your own scans (`ADM_38_CABI`, `ADM_43_RHGR`, `ADM_37_VAUG`,
`ADM_38_DILA`, `ADM_38_ORTR`, `ADM_46_BECO`), both standards, and synthetic
scans with a known true area.

---

## The four bugs

### 1. The lid shadow was being measured instead of the leaf

`segment_leaf()` kept `ranked[1]` — the largest dark connected object. On your
400 dpi scans the darkest, largest object is not a leaf. It is the strip the
scanner lid casts down the right edge of the glass:

| scan | edge strip | actual leaf |
|---|---|---|
| ADM_38_CABI | 25,573 px (cols 3390–3399, all 4,676 rows) | 10,598 px |
| ADM_43_RHGR | 14,439 px | 10,005 px |

The strip runs from gray 240 down to 103 over the last ~12 columns, while the
leaves sit at 130–155 against paper at 245. So on any scan where the leaf is
small, the strip wins and the reported area is the strip.

**Fix:** trim `BORDER_CROP_IN` (0.10 in) off every side before anything looks at
the image, and prefer objects that do not run off the edge. `touches_edge` is
now a rejection rule, not just a note in the flags column.

### 2. Otsu's method is the wrong estimator for a leaf that is 0.1% of the scan

Otsu maximizes between-class variance, which assumes the two classes are of
broadly comparable size. Your leaves are 0.07–0.13% of the pixels, so the
optimum is set by the texture of the paper, not by the leaf. Concretely:
trimming 40 irrelevant pixels off the border of `ADM_38_CABI` moved Otsu's
cutoff from 190 to 205 and the measured area by **24%**. A method whose answer
moves 24% when you crop away something that is not the leaf is not measuring
the leaf.

Worse, `separability` stayed at 0.85–0.99 the whole time, so `MIN_SEPARABILITY`
never fired. Every row of the shipped `output/leaf_area.csv` has separability
> 0.98 and no flags.

**Fix:** find the leaf first, then place the cutoff halfway between the leaf's
own tone and the paper immediately around it, iterating until it stops moving.
For an edge blurred by the scanner's optics, the half-maximum crossing is the
unbiased position of the true boundary, and it does not care how much paper
surrounds the leaf. `separability` is replaced by `contrast` (paper − leaf, in
gray levels) and `area_sens_pct`, which is measured, not assumed — see below.

### 3. The size standard was scanned at the wrong resolution

`STANDARD_3inx3in.jpeg` is 2550 × 3508 = **300 dpi**. Your `ADM_*` scans are
3400 × 4676 = **400 dpi**. `CM2_PER_PX` came from the 300 dpi standard and was
applied to 400 dpi scans, inflating every `ADM_*` area by (400/300)² = **78%**,
on top of bug 1.

**Fix:** `STANDARDS` is now a table. Each scan is matched to the standard
scanned at its own resolution, so 300 and 400 dpi batches can sit in the same
folder. A scan with no matching standard is an error, not a silent 78%.

### 4. Four-connectivity shatters thin diagonal objects

`start[above] <= end[j] & start[j] <= end[above]` requires column overlap, which
is 4-connectivity. A needle lying diagonally on the glass is a staircase of
single pixels; under 4-connectivity it breaks into fragments and none of them
looks like a leaf. `ADM_38_CABI` and `ADM_38_ORTR` are exactly this case.

**Fix:** 8-connectivity for objects, 4-connectivity for background inside
`fill_holes()`, which is the standard pairing and keeps the two topologically
consistent.

### Smaller things
- `write.csv()` failed if `output/` did not exist. Now created.
- `holes_px` was reported but never explained; areas include filled holes, so
  `leaf_px - holes_px` is the outline minus insect damage.
- The masks are now written at the size of the original scan, so they overlay it
  directly.
- The multi-page TIFF warning was unreachable (`readTIFF(all = FALSE)` never
  returns a list).
- Typos: "How does it worl?", "the faction of total var". British spellings
  (grey, greyscale, metre, cut-off) changed to American throughout.

---

## Accuracy check

Synthetic 300 dpi scans with the same paper texture, shading gradient and lid
shadow as yours, so the true area is known exactly:

| test | true px | measured px | error |
|---|---|---|---|
| Low-contrast rectangle, 2.000 × 1.500 in (43 gray levels of contrast) | 270,000 | 269,999 | **−0.0004%** |
| Diagonal needle, 8 px wide | 17,003 | 17,001 | **−0.012%** |
| Green leaf, color scan | 300,000 | 299,997 | −0.001% |
| Leaf with an insect hole | 529,688 | 549,753 − 20,069 filled | −0.001% |
| Small high-contrast leaf | 18,000 | 17,996 | −0.022% |
| Leaf running off the edge | 560,000 | 538,998 | −3.75%, flagged `touches_edge` |
| Blank scan | — | error, run continues | — |
| Leaf plus a paper label | 490,000 | 489,998 | 0.000%, flagged `extra_object` |

The two independent standards cross-check each other: the post-it note at
400 dpi measures 402.7 dpi and the black square at 300 dpi measures 302.0 dpi —
the **same 0.67% scale offset**, from two different objects at two different
resolutions. That is a real difference between your scanner's nominal and actual
sampling pitch, and it is exactly what the physical standard is for. Using DPI
metadata alone would leave a 1.3% area bias in place.

---

## Before and after, on your scans

| scan | old cm² | new cm² | old / new |
|---|---|---|---|
| ADM_37_VAUG | 2.12 | 1.202 | 1.76× |
| ADM_38_CABI | 1.81 | 0.492 | 3.7× |
| ADM_38_DILA | 2.21 | 0.099 | 22× |
| ADM_38_ORTR | 1.29 | 0.196 | 6.6× |
| ADM_43_RHGR | 1.02 | 0.396 | 2.6× |
| ADM_46_BECO | 10.74 | 6.058 | 1.77× |

VAUG and BECO are off by 1.76× and 1.77×, essentially exactly (400/300)² = 1.778
— their segmentation was fine, only the calibration was wrong. The other four
picked up the lid shadow as well.

---

## What you need to check before trusting the numbers

1. **Measure `POST_IT_NOTE.tiff` with calipers.** It is now the ruler for every
   `ADM_*` area, and the script currently assumes 3 × 3 in = 58.0644 cm². If the
   note is really 76 mm (2.992 in), every ADM area is 0.5% high. Set
   `STANDARDS$area_cm2` to what you measure.
2. **Confirm the 300 dpi group** (`7a`, `8b`, `9a`, …) really belongs with
   `STANDARD_3inx3in.jpeg`. I could not stage those 26 MB files, so they are
   untested on real data, though the 300 dpi path is covered by the synthetics.
3. **Look at the masks.** `scans/processed_scans_jpeg/*.png` overlays the
   original scan pixel for pixel.

## Two new output columns worth reading

- `contrast` — gray levels between paper and leaf. Below `MIN_CONTRAST` (15) the
  boundary is guesswork.
- `area_sens_pct` — how much the area moves when the cutoff is nudged by ±10% of
  that contrast, with `area_cm2_lo` / `area_cm2_hi` as the bounds. This is a
  real error bar rather than a proxy. On your scans it runs 0.8% for the broad
  BECO leaf to 15% for the ORTR needle, which is honest: a needle's translucent
  margin genuinely does not have a well-defined edge at 400 dpi. Treat
  `soft_margin` rows as approximate.

## A note on dependencies

I kept this to base R plus `tiff`/`jpeg`/`png`. `imager` or `EBImage` would have
shortened the connected-components code, but the failures above were all
algorithmic, not computational — no library would have fixed any of them — and
the run is now ~2.5 s per scan (down from 4 s, despite doing considerably more
work), so speed is not a reason to add install friction for collaborators.

The one place a library would genuinely help is if you later want to separate
touching or overlapping leaves, which needs a watershed. Say the word if that
comes up.

## One workflow suggestion

Your scans are grayscale stored as RGB — R, G and B are byte-identical. Scanning
in **color** would help the low-contrast cases measurably: chlorophyll absorbs
red, so a green leaf on white paper separates far better in the red channel than
in any gray mix. The new `to_gray()` already picks the highest-contrast channel
automatically, so color scans would just work.

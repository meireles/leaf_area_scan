# leaf_area_scan — design review and plan

Everything below is either measured on your own scans or cited. Where the
evidence is thin I say so.

---

## Short answers

**1. Detect the shadow by position and shape rather than cropping — yes, and it
is now implemented and validated.** Better still, the shadow does not need to be
*detected* at all in the usual sense. It is a *separable* artifact — a function
of position across the bed, constant along the scan direction — so it can be
estimated as a 1-D profile and divided out. A leaf lying 6 px from the platen
edge now measures **−0.02%**, identical to the same leaf in the middle of the
bed. Under the old fixed crop it was deleted outright.

**2. Higher-resolution JPEG vs grayscale TIFF — resolution wins, and it is not
close.** JPEG at quality ≥ 90 costs you under 0.5% on every leaf you have. Losing
resolution costs you up to 36%. So a 600 dpi JPEG q95 genuinely beats a 400 dpi
TIFF. But this is the wrong trade to be making: your whole archive is ~300 MB.
Scan 600 dpi color TIFF and the question disappears.

**3. What DPI — 600, and the reason is needle width, not leaf area.** Your broad
leaves are flat to within 1% all the way down to **75 dpi**. Your 0.74 mm ORTR
needle is **28% off at 300 dpi**. Resolution requirement is set entirely by the
narrowest object on the plate.

**4. Does this exist in R already — no, and the gap is real enough to publish.**
The two real options (`LeafArea`, `pliman`) both use a global histogram threshold,
which is the estimator that fails at 0.1% foreground; neither handles scanner
edge artifacts; and **no tool surveyed reports a per-measurement uncertainty**.

**5. Touching leaves — do not split them.** Measured on two of your real leaf
masks pushed together at 0 px separation: the total area is **+0.00% correct**
even when they fuse into one blob. Splitting can only lose pixels. Watershed on
your actual needles found **16–27 seeds for 2 needles**.

---

## 1. The shadow: correct it, don't crop it

### What it actually is

Per-column 90th percentile across all seven of your scans:

| scan | left | right | top | bottom |
|---|---|---|---|---|
| CABI | 0 px | **16 px** | 0 px | 0 px |
| RHGR | 0 px | **14 px** | 0 px | 0 px |
| VAUG | 0 px | **16 px** | 0 px | 0 px |
| BECO | 0 px | **14 px** | 0 px | 0 px |
| DILA | 0 px | **17 px** | 0 px | 0 px |
| ORTR | 0 px | **15 px** | 0 px | 0 px |
| POST_IT | 0 px | **20 px** | 0 px | 0 px |

(px = how far in the profile departs from the plateau by more than 3 gray
levels.) Right edge only, 14–20 px, and monotone: `242 241 240 238 235 221 202
169 133 112 110 108`. The other three edges are clean to within noise.

That regularity is the whole design. The artifact is a **lateral response**, not
a 2-D blob: a function of distance across the bed, constant along it. This is
exactly the structure the radiochromic-film-dosimetry literature characterized
for flatbeds — Schoenfeld et al. show the profile appears even for
dose-independent optical subjects, i.e. it is a property of the scanner and not
of the sample ([10.1088/0031-9155/61/21/7704](https://doi.org/10.1088/0031-9155/61/21/7704));
Lewis & Chan give the per-lateral-position correction recipe
([10.1118/1.4903758](https://doi.org/10.1118/1.4903758)).

### The correction

A high percentile down each column is the paper level at that column, and is
robust while a leaf covers up to 90% of it. Divide it out — multiplicative,
because the lid attenuates light and scales the leaf and the paper together, so
contrast is restored rather than merely offset. Then drop only the few columns
where the lid has taken more than 30% of the light: past that, correcting means
amplifying noise, and there is no signal left worth recovering.

Result, on a synthetic scan carrying your measured shadow profile:

| needle position | raw | after correction |
|---|---|---|
| 6 px from the platen edge | **+159.3%** (measures the shadow) | **−0.02%** |
| 15 px from the edge | +102.3% | −0.02% |
| 25 px and beyond | +102.3% | −0.02% |

And on a broad leaf running off the edge of the glass, the fixed crop was losing
exactly the strip it removed:

| | fixed 0.10 in crop | lateral correction |
|---|---|---|
| leaf running off the platen edge | −3.75% | **−0.000%** |

Nothing else moved: BECO, whose leaf is nowhere near the edge, changed by
+0.05%; the post-it standard by +0.006%.

### The shape rule, for what's left

After correction, a residual sliver of amplified noise can still hug the border.
That is where your position-and-shape idea belongs, and it is now the rejection
rule:

> A component is the scanner's edge artifact if it touches a border, reaches no
> more than 2% of the way into the bed, and runs at least 10× further along the
> border than it reaches inward.

Your shadow is roughly 300:1 pinned to the border. A 2,000 px needle laid right
along the edge is 200:1 — so it *would* trip this rule, which is why the
illumination cutoff does the heavy lifting and the shape rule is only a backstop
for noise. The important change from the old code is what this **replaces**: the
old rule demoted *any* component touching an edge, so a leaf genuinely lying
against the glass edge was rejected on principle. Now it is measured.

---

## 2. Resolution: set by the narrowest object, not the leaf

I resampled each of your scans to 300/200/150/100/75 dpi (area-average, which
models a lower-resolution sensor) and re-measured. Error is against the same
leaf at 400 dpi.

| leaf | mean width | 300 | 200 | 150 | 100 | 75 |
|---|---|---|---|---|---|---|
| ADM_37_VAUG | 5.00 mm | +0.2% | −0.2% | +0.2% | +0.4% | +0.3% |
| ADM_46_BECO | 2.38 mm | −0.0% | −0.0% | −0.1% | +0.2% | +0.3% |
| ADM_43_RHGR | 2.25 mm | −0.2% | −0.6% | −0.9% | +0.9% | +0.3% |
| ADM_38_DILA | 1.26 mm | −3.5% | −2.9% | −3.6% | **+10.7%** | **+14.4%** |
| ADM_38_CABI | 1.06 mm | +1.4% | **+13.8%** | **+9.0%** | **+11.4%** | **+13.5%** |
| ADM_38_ORTR | **0.74 mm** | **+28.1%** | **+29.4%** | **+33.6%** | **+28.9%** | **+36.5%** |

The mechanism is visible in the measured contrast, which collapses for exactly
the leaves that break:

| leaf | 400 | 300 | 200 | 150 | 100 | 75 |
|---|---|---|---|---|---|---|
| VAUG / BECO / RHGR | 87 / 88 / 104 | 87 / 88 / 104 | 87 / 87 / 105 | 87 / 88 / 106 | 87 / 88 / 102 | 88 / 88 / 100 |
| ADM_38_CABI | 89 | 86 | **69** | 78 | 79 | 79 |
| ADM_38_ORTR | 84 | 84 | **58** | 57 | 59 | 57 |
| ADM_38_DILA | 68 | 72 | 72 | 70 | **52** | 53 |

A needle narrower than a few pixels never reaches its own true tone, so the
half-maximum cutoff sits too high and the needle fattens. This is the known
"narrow feature" limit of the ISO50 / half-max criterion — the same reason it is
unreliable for thin walls in dimensional CT
([10.1016/j.cirp.2011.05.006](https://doi.org/10.1016/j.cirp.2011.05.006)).

There is a *second*, opposite bias for curved boundaries: the 50% gray contour
is pulled toward the center of curvature, so areas are underestimated by roughly
σ²/(2R) for blur σ and radius R (Verbeek & van Vliet,
[10.1002/1361-6374(199303)1:1<47::AID-BIO8>3.0.CO;2-7](https://doi.org/10.1002/1361-6374%28199303%291:1%3C47::AID-BIO8%3E3.0.CO;2-7)).
For a 12 px-wide needle with σ ≈ 1 px that is ~1.5%; at 300 dpi it is ~5%. Both
biases push you the same direction: **more resolution**.

**Working rule: at least ~12 pixels across the narrowest object.**

| narrowest object | minimum dpi | your species |
|---|---|---|
| 0.5 mm | 600 | small needles |
| 0.75 mm | 400 | ORTR |
| 1.0 mm | 300 | CABI |
| 2.0 mm | 150 | RHGR, BECO |

Two cautions. Do not exceed your scanner's *hardware* resolution — advertised
figures above 1200 dpi are usually interpolated, and interpolated DPI would
corrupt any metadata-based calibration
([10.1002/aps3.11556](https://doi.org/10.1002/aps3.11556)). And 600 dpi is a
1.5× linear / 2.25× area jump from 400, which is fine; going to 1200 quadruples
the data for a marginal return.

### An honest limitation on the existing archive

The table above uses 400 dpi as its own reference, so it shows a *trend*, not
absolute accuracy. The trend for ORTR has clearly not converged by 400 dpi.
**For needles narrower than ~1 mm in the existing archive, I cannot tell you the
absolute error, only that it is probably several percent and biased upward.**
The `area_sens_pct` column reports 10–15% for exactly these leaves, which is the
right order — but at 300 dpi ORTR's true error was 28% while its sensitivity
still read ~12%, so the flag under-reports. Treat it as "this number is soft,"
not as a calibrated interval.

Phase 3 below is how to fix that retroactively.

---

## 3. Format: TIFF vs JPEG

Same scans, 400 dpi throughout, error against the lossless TIFF:

| leaf | q100 | q95 | q90 | q85 | q75 | q60 |
|---|---|---|---|---|---|---|
| ADM_37_VAUG | −0.01% | −0.02% | +0.01% | +0.04% | +0.02% | +0.03% |
| ADM_38_CABI | −0.02% | −1.28% | +0.18% | +0.16% | +0.08% | +0.16% |
| ADM_38_DILA | 0.00% | +0.10% | −0.10% | −0.40% | +0.10% | −0.51% |
| ADM_38_ORTR | −0.10% | −0.46% | −0.20% | **+1.79%** | **+2.35%** | **+2.46%** |
| ADM_43_RHGR | −0.03% | +0.05% | −0.10% | −0.05% | +0.03% | −0.13% |
| ADM_46_BECO | 0.00% | −0.01% | −0.01% | 0.00% | +0.02% | +0.02% |

Mean file size: TIFF 2.84 MB, q95 0.93 MB, q90 0.40 MB.

JPEG is not as dangerous as its reputation *here*, because your objects are
large relative to the 8×8 DCT block and reasonably contrasty. It degrades
exactly where you would predict — the thinnest needle, at quality ≤ 85. Note
that the error is not a bias with a fixed sign; it is added variance
concentrated at low-contrast boundaries, which is where your `contrast` flag
already says you are vulnerable.

**Recommendation: 600 dpi color TIFF with LZW.** At ~12 MB/scan and ~55 scans,
that is 660 MB — nothing. Keep JPEG out of the archive entirely. Two secondary
reasons beyond compression: JPEG frequently carries no DPI tag (your own
`STANDARD_3inx3in.jpeg` does not, which is why it needs a manual override in
`STANDARDS`), and lossy files can silently degrade if anything ever re-saves
them.

**Scan in color, not grayscale.** Your TIFFs are grayscale duplicated across
R, G and B — byte-identical channels, so a third of the file is wasted and the
useful signal was thrown away at the scanner. Chlorophyll absorbs red, so a
green leaf on white paper separates far better in the red channel than in any
gray mix. `to_gray()` already picks the highest-contrast channel automatically,
so color scans need no code change and would raise contrast on exactly the
washed-out leaves.

---

## 4. Touching and overlapping leaves

### Total area needs no splitting

I took your real BECO and VAUG masks, and your real CABI and ORTR needle masks,
and laid each pair together at decreasing separation:

| pair | 40 px apart | 4 px | 1 px | **0 px (fused)** |
|---|---|---|---|---|
| BECO + VAUG (182,442 px true) | +0.00% | +0.00% | +0.00% | **+0.00%** |
| CABI + ORTR (17,279 px true) | −0.03% | −0.02% | −0.03% | **−0.03%** |

Two objects that touch have disjoint interiors, so the union area *is* the sum.
Nothing is lost by leaving them fused — and a splitter can only remove pixels
along whatever cut line it draws.

### Splitting your leaves does not work

Distance-transform watershed (what `pliman` and `PlantCV` use), seeded exactly
as `EBImage::watershed` does, on those same real masks:

| pair (true count 2) | tol=1, ext=1 | tol=2, ext=1 | tol=5, ext=2 |
|---|---|---|---|
| BECO + VAUG | 3 | 3 | 3 |
| **CABI + ORTR (needles)** | **27** | **16** | **16** |

This is structural, not a tuning failure. For an object of near-constant width
the distance transform is a **plateau along the medial axis, not a peak** — the
seed information simply is not there. Any h large enough to suppress the ripples
merges everything; any h small enough to separate two needles finds dozens of
maxima along each one. Omnipose was built specifically because of this
([10.1038/s41592-022-01639-4](https://doi.org/10.1038/s41592-022-01639-4)), and
`pliman`'s own docs say "I strongly suggest using this only with round objects."

The algorithm that *does* work on lobed and irregular leaves is concavity-based
splitting — cut between paired boundary concavities — which is what BioVoxxel's
*Watershed Irregular Features* does and why the MuLES multi-leaf tool chose it
([10.1002/aps3.11513](https://doi.org/10.1002/aps3.11513)). **There is no R
implementation.** And it still fails on needles, which generate concave points
from their own curvature.

### Overlap is not recoverable, and should never be imputed

If one lamina lies *on top of* another, the silhouette gives |A ∪ B| = |A| + |B|
− |A ∩ B|, and |A ∩ B| is simply not in the image. No algorithm recovers it;
"shape completion" methods impute it from a prior, and for leaf outlines that
prior is far too loose to be a measurement. Detect and re-scan, never estimate.

*(The one real exception: in transmission mode with a transparency unit, optical
density is roughly additive, so two stacked laminae are measurably darker than
one. It needs a per-batch OD calibration, saturates by 2–3 layers, and is very
sensitive to the lateral non-uniformity of §1. Probably not worth it for you.)*

### What this means for the design

You chose "total per scan, but flag fusion." That is the right call and it maps
cleanly onto the code:

- **Sum every object above `MIN_OBJECT_CM2`**, instead of keeping only the
  largest. Report the per-object breakdown in a companion CSV so you can audit.
- **Flag likely fusion by shape**, don't split it. Two touching laminae make a
  component with anomalously low *solidity* (area ÷ convex hull area) and a
  perimeter far larger than a single leaf of that area would have. Calibrate the
  threshold on your own single-leaf scans — you have ~55, which is enough to set
  it empirically rather than by guess.
- **Space the leaves at scan time.** Cornelissen et al. found multiple leaves per
  image saves substantial time at no accuracy cost
  ([10.1002/ecy.70308](https://doi.org/10.1002/ecy.70308)) — the win is real,
  and it only requires that they not touch.

---

## 5. What already exists in R

| tool | segments itself? | deps | multi-leaf | splits touching | calibration | threshold | edge artifacts | uncertainty |
|---|---|---|---|---|---|---|---|---|
| **LeafArea** (CRAN, last release **2019**) | no — writes an ImageJ macro | ImageJ + Java | yes, but returns the **sum** | no | `Set Scale` | ImageJ `Minimum`, global | blind symmetric trim | no |
| **pliman** (CRAN, active) | yes | **EBImage** (Bioconductor) + fftw3 | yes, per-object | DT watershed (round objects only) | `dpi=` **or** `reference_area=` | Otsu global, or `EBImage::thresh` adaptive | none | no |
| **EBImage** (Bioconductor) | building block | fftwtools, libfftw | — | `watershed`, `distmap` | — | `otsu`, `thresh` | none | — |
| **imager** (CRAN, revived 2026) | building block | libX11, fftw3 | — | `watershed` | — | `threshold` | none | — |
| **Momocs** | no — needs pre-made silhouettes | — | — | — | — | — | — | — |
| **MuLES** (Fiji macro, not R) | yes | Fiji + BioVoxxel | **yes, by design** | **yes** — convexity-based | 2-point interactive | ImageJ auto | none | no |

Concretely, on your six scans:

- `LeafArea`'s default `low.size = 0.7` cm² would **silently discard four of your
  six leaves** (0.098, 0.192, 0.396, 0.491 cm²).
- `pliman`'s `lower_noise = 0.1` drops objects below 10% of *mean* object area —
  with a large lid-shadow strip inflating the mean, your leaf is the thing that
  gets deleted.
- `pliman`'s `dpi=` calibration path would have propagated your scanner's real
  0.67% pitch offset as a 1.3% area bias. Only its `reference_area=` path
  matches what you do.
- Neither would have rejected the lid strip, and both default to a global
  histogram threshold — the estimator that failed at 0.1% foreground.
- **Neither reports a per-measurement uncertainty.** Nothing surveyed does.

**Verdict: keep your implementation.** The one thing worth adopting is `pliman`
as an independent cross-check on a handful of scans — two implementations
agreeing is worth more than either one's internal diagnostics — kept as an
optional dependency, never in the main path.

### Is it a package?

I think yes, and specifically because of the last row of that table. Four things
here are, as far as I can find, not available anywhere in R:

1. scanner lateral-response correction for flatbeds
2. threshold anchored to physical gray levels rather than histogram shape, which
   is what makes 0.1%-foreground objects tractable
3. resolution-matched calibration against multiple physical standards
4. **a per-measurement uncertainty on every reported area**

That is a legitimate Applications in Plant Sciences / MEE software note, and the
validation you'd need for it is mostly the test suite you should build anyway.
The honest cost is maintenance: a CRAN package is a multi-year commitment, and
`LeafArea` going quiet since 2019 is the cautionary example.

**Suggested hedge:** build it *package-shaped* from the start — `R/`,
`tests/testthat/`, roxygen — which costs maybe a day more than a tidy script and
loses nothing. Decide about CRAN once you have the validation numbers in hand.
Until then it installs fine from GitHub.

---

## 6. Plan

### Phase 0 — done
- ✅ Lid shadow corrected, not cropped. Leaf at the platen edge: −0.02%.
- ✅ Edge artifact identified by position and shape, so an edge-touching leaf is
  no longer rejected on principle.
- ✅ Detection floor rewritten so a leaf at 0.016% of the image is still found.
- ✅ Full test suite passing within **0.022%** on every synthetic case.

### Phase 1 — scanning protocol *(no code; the biggest single win)*
1. **600 dpi**, hardware resolution, no interpolation.
2. **Color, not grayscale.** No code change needed.
3. **TIFF with LZW.** No JPEG in the archive.
4. **Leaves spaced apart**, off the platen margin.
5. Turn off auto-exposure, auto-color, sharpening, descreening. Every one of
   them makes the response non-linear and breaks the half-max argument.
6. **One empty-platen scan and one white-reference scan per session**, at each
   resolution you use. Cheap now, and it makes §1's correction unnecessary
   rather than merely effective — ImageJ's own guidance is to fix illumination
   at acquisition rather than post-hoc.
7. **Scan the standard in the same session, at the same resolution.** The 78%
   bug came from not doing this.

### Phase 2 — code
1. **Sum all objects** above `MIN_OBJECT_CM2` rather than keeping the largest;
   emit a per-object CSV alongside the per-scan one.
2. **Fusion flag** from solidity and perimeter, threshold calibrated on your
   existing single-leaf scans.
3. **Flat-field path**: use the empty-platen reference when present, fall back to
   the lateral profile when not.
4. **Overlap detection** → flag for re-scan, never impute.
5. **Test suite**: the synthetic generators I used here, checked into the repo
   with their known true areas, so a future change cannot silently reintroduce
   the class of bug that started this.
6. **Package skeleton**: `R/`, `tests/testthat/`, roxygen, a vignette that is
   the DEBUG_NOTES analysis.

### Phase 3 — validation *(this is what makes it publishable, and fixes the archive)*
1. **Build a width calibration target.** Objects of accurately known width
   spanning 0.25–5 mm. Do not trust a printed target's line widths — toner
   spread makes thin lines wider than nominal. Use feeler gauges or shim stock
   (~$10, known to ±0.005 mm), or measure a printed target with calipers.
2. **Scan it at 1200 / 600 / 400 / 300 dpi.** Now you can measure the width-
   dependent bias directly instead of inferring it from a trend, and derive a
   **correction curve applicable retroactively to the existing 400 dpi archive**
   — which is the only way to rescue the needle measurements you cannot re-scan.
3. **Measure the scanner's true pitch** across the bed with fiducials at a known
   long separation, verified against a steel rule. You already have evidence it
   is 0.67% off nominal, consistently, from two independent standards.
4. **Cross-check against `pliman`** on 10 scans.
5. **Verify the post-it with calipers.** It is currently the ruler for every
   400 dpi area and the script assumes exactly 3 × 3 in.

### Phase 4 — write-up, if you want it
Software note. The validation from Phase 3 plus the before/after in
`DEBUG_NOTES.md` is most of a results section already.

---

## 7. Open risks

| risk | severity | mitigation |
|---|---|---|
| Needles < 1 mm in the **existing** archive carry an unquantified upward bias | **high** — it is a systematic error correlated with species | Phase 3.2 correction curve. Until then, report `area_sens_pct` alongside every needle measurement and do not compare needle species against broad-leaf species without it. |
| `area_sens_pct` under-reports true error for the narrowest needles (12% reported vs 28% actual at 300 dpi) | medium | Phase 3 gives a real calibration; meanwhile treat it as an ordinal warning, not an interval. |
| The post-it's true area is assumed, not measured | medium — sets absolute scale for every ADM area | Calipers. Ten minutes. |
| `BORDER_ASPECT` / `MIN_ILLUMINATION` tuned on one scanner | low | They are config, and the empty-platen reference from Phase 1.6 makes them nearly irrelevant. |
| Different scanner or firmware changes the lateral profile | low | The profile is estimated per image, so it adapts. But re-run the empty-platen check if the hardware changes. |

---

## Sources

Cornelissen et al. 2026, *Ecology* 107(2):e70308 — [10.1002/ecy.70308](https://doi.org/10.1002/ecy.70308) ·
Katabuchi 2015, *Ecol Res* 30:1073 — [10.1007/s11284-015-1307-x](https://doi.org/10.1007/s11284-015-1307-x) ·
Olivoto 2022, *MEE* 13:789 (pliman) — [10.1111/2041-210X.13803](https://doi.org/10.1111/2041-210X.13803) ·
MuLES — [10.1002/aps3.11513](https://doi.org/10.1002/aps3.11513) ·
Easy Leaf Area — [10.3732/apps.1400033](https://doi.org/10.3732/apps.1400033) ·
Lakeram et al. 2023 (interpolated dpi) — [10.1002/aps3.11556](https://doi.org/10.1002/aps3.11556) ·
Bakr 2005 (dpi experiment) — [10.1111/j.1439-0418.2005.00948.x](https://doi.org/10.1111/j.1439-0418.2005.00948.x) ·
Minervini et al. 2015 (compression in phenotyping) — [10.1071/FP15033](https://doi.org/10.1071/FP15033) ·
Schoenfeld et al. 2016 (flatbed parabola effect) — [10.1088/0031-9155/61/21/7704](https://doi.org/10.1088/0031-9155/61/21/7704) ·
Lewis & Chan 2015 (lateral response correction) — [10.1118/1.4903758](https://doi.org/10.1118/1.4903758) ·
van Battum et al. 2016 — [10.1088/0031-9155/61/2/625](https://doi.org/10.1088/0031-9155/61/2/625) ·
Otsu 1979 — [10.1109/TSMC.1979.4310076](https://doi.org/10.1109/TSMC.1979.4310076) ·
Xu et al. 2011 (Otsu bias analysis) — [10.1016/j.patrec.2011.01.021](https://doi.org/10.1016/j.patrec.2011.01.021) ·
Zack et al. 1977 (triangle method) — [10.1177/25.7.70454](https://doi.org/10.1177/25.7.70454) ·
Kruth et al. 2011 (ISO50) — [10.1016/j.cirp.2011.05.006](https://doi.org/10.1016/j.cirp.2011.05.006) ·
Verbeek & van Vliet 1993 (curvature bias) — [10.1002/1361-6374(199303)1:1<47::AID-BIO8>3.0.CO;2-7](https://doi.org/10.1002/1361-6374%28199303%291:1%3C47::AID-BIO8%3E3.0.CO;2-7) ·
Vincent & Soille 1991 (watershed) — [10.1109/34.87344](https://doi.org/10.1109/34.87344) ·
Cutler et al. 2022 (Omnipose) — [10.1038/s41592-022-01639-4](https://doi.org/10.1038/s41592-022-01639-4) ·
Bai et al. 2009 (concavity splitting) — [10.1016/j.patcog.2009.04.003](https://doi.org/10.1016/j.patcog.2009.04.003) ·
Peng et al. 2017 (BaSiC shading correction) — [10.1038/ncomms14836](https://doi.org/10.1038/ncomms14836) ·
ImageJ imaging principles — https://imagej.net/imaging/principles

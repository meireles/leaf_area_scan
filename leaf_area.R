################################################################################
## Leaf area from flatbed scans
##
## Measures the one-sided area of a single leaf per scan and writes a tidy CSV
## plus a binary mask per scan, so every number can be audited by eye.
##
## Pipeline, in order:
##   1. read_scan()        read the file, reduce it to one 8-bit gray channel
##   2. flatten_scanner()  divide out the scanner's lateral response: the dark
##                         band the lid casts along the edge of the glass, plus
##                         any shading across the bed
##   3. segment_leaf()     find the leaf on a coarse copy, then measure it at
##                         full resolution inside its own neighborhood, with the
##                         cutoff placed halfway between the leaf's own tone and
##                         the paper immediately around it
##   4. calibrate_standards() / match_standard()
##                         convert pixels to cm2 using a physical standard
##                         scanned at the SAME resolution as the leaves
##
## Every function below is documented in roxygen2 form, so this file can be
## split into an R package's R/ directory without rewriting the comments.
## Nothing executes until the "Run" section at the very bottom; set
##
##   options(leaf_area.functions_only = TRUE)
##
## before source()ing to load the functions without processing anything. That is
## what the test suite does.
##
## Areas INCLUDE holes: insect damage and pale patches count as leaf, on the
## grounds that a hole is still inside the leaf outline. `holes_px` reports how
## much was filled, so `leaf_px - holes_px` is the outline minus the damage.
##
## Masks are written at the size of the original scan, so they overlay it
## directly, and losslessly, so re-counting a mask reproduces the area exactly.
##
## Flags in the output CSV:
##   low_contrast  the leaf is barely darker than the paper, so the boundary is
##                 guesswork. Rescan in color, or on a darker backing.
##   soft_margin   the area moves more than max_sensitivity percent when the
##                 cutoff is nudged. Typical of needles, hairy or translucent
##                 margins, or a leaf that was not lying flat on the glass.
##   extra_object  something else on the plate is a decent fraction of the size
##                 of the leaf. Check the mask; a second leaf will be ignored.
##   interior_recovered
##                 part of the leaf was too pale to threshold and was recovered
##                 by sealing the outline and accepting the enclosed region on
##                 the strength of its texture. `recovered_px` says how much.
##                 Check the mask.
##   touches_edge  the leaf reaches the edge of the measured region. Usually a
##                 lower bound, though after the lateral correction the loss is
##                 only the few columns the lid had effectively blinded.
##   dpi_mismatch  the standard and the file's own DPI tag disagree about scale.
##
## Requires the tiff, jpeg and png packages (only the ones matching your input
## formats are actually needed):
##
##   install.packages(c("tiff", "jpeg", "png"))
##
################################################################################


## =============================================================================
## Options
## =============================================================================

#' Tuning parameters for leaf area measurement
#'
#' Every threshold and magic number the pipeline uses, in one place, so that no
#' function reaches for a global. Call with no arguments for the defaults, which
#' are what the validation in `PROJECT_PLAN.md` was run against; override
#' individual entries by name.
#'
#' The defaults were tuned on 400 dpi grayscale scans of leaves 0.1 to 6 cm2 on
#' white paper. The ones most likely to need changing on a different scanner are
#' `min_illumination` and `border_aspect`, which describe the lid shadow, and
#' `min_object_cm2`, which sets what counts as debris.
#'
#' @param paper_percentile Percentile taken down each column (and along each
#'   row) as the paper level at that position. Robust while a leaf covers up to
#'   `100 - paper_percentile` percent of a column, so 90 tolerates a leaf
#'   spanning nine tenths of the bed.
#' @param profile_smooth_px Width of the running median applied to that profile.
#'   Damps sensor noise without smearing the shadow's ramp, which is only about
#'   15 px wide. Must be odd.
#' @param min_illumination Rows and columns receiving less than this fraction of
#'   the brightest light are dropped rather than corrected. Past this point
#'   correcting means amplifying noise and there is no signal worth recovering.
#'   Lowering it keeps more of the bed at the cost of noisier edges.
#' @param border_depth_frac,border_aspect Shape test that identifies the lid
#'   shadow. A dark component touching a border, reaching no more than
#'   `border_depth_frac` of the way into the bed, and running at least
#'   `border_aspect` times further along the border than inward, is the artifact
#'   rather than a specimen. Measured shadows here are roughly 300:1.
#' @param min_object_cm2 Smallest object reported as a leaf. Anything below is
#'   dust, a fiber or a speck. Also sets the sensitivity of the "is anything
#'   here at all?" test in [segment_leaf()].
#' @param detect_px Long side of the downsampled copy used to locate the leaf.
#'   Detection only has to answer "roughly where?", and doing it small
#'   suppresses paper noise and keeps the run-length labeling cheap.
#' @param detect_fraction How far from paper toward the darkest tone the
#'   detection cutoff sits. Deliberately permissive: it only has to catch the
#'   leaf, since [refine_cutoff()] then puts the real boundary in place.
#' @param core_erode_px Radius the leaf is shrunk by before its tone is read, so
#'   soft edge pixels do not drag the cutoff around.
#' @param recover_pale_interior Whether to run [recover_interior()], which
#'   looks inside the outline of a measured leaf for lamina too pale to
#'   threshold. Set `FALSE` to reproduce the behaviour of the pipeline without
#'   it.
#' @param texture_px Window radius for the local standard deviation that
#'   distinguishes a specimen from paper.
#' @param texture_ratio How many times the paper's own local SD a region must
#'   exceed to count as textured. Paper measures 0.4 to 0.6 and a tomentose
#'   lamina 5 to 7, so this is not a delicate threshold.
#' @param texture_frac Fraction of a candidate region that must read as textured
#'   before it is accepted as leaf.
#' @param interior_min_frac Smallest region [recover_interior()] will consider,
#'   as a fraction of the leaf already found. A lamina that thresholding missed
#'   is a substantial part of the leaf; anything smaller is a sliver of hull
#'   slack and is skipped, which is also what keeps the step affordable.
#' @param ring_inner_px,ring_outer_px The paper tone is read from a frame around
#'   the leaf standing off by `ring_inner_px` and extending `ring_outer_px`.
#'   Local rather than global, so shading across the glass cannot bias it.
#' @param max_iterations Cap on cutoff refinement passes. Convergence is
#'   normally reached in two or three.
#' @param sensitivity_frac Fraction of the leaf-to-paper contrast the cutoff is
#'   nudged by, in each direction, to produce the `area_cm2_lo` / `area_cm2_hi`
#'   error bar.
#' @param min_contrast Gray levels between paper and leaf below which the scan
#'   is too washed out to trust. Populates `low_contrast`.
#' @param max_sensitivity Percent change in area per nudge of the cutoff above
#'   which the margin is called soft. Populates `soft_margin`.
#' @param max_second_object Size of the runner-up object, as a fraction of the
#'   leaf, above which something else is clearly on the plate. Populates
#'   `extra_object`.
#' @param dpi_tolerance Relative disagreement allowed between the standard-based
#'   and DPI-based areas. Populates `dpi_mismatch`.
#' @param dpi_match_tolerance Relative difference in resolution within which a
#'   scan is considered to match a size standard.
#'
#' @return A named list of settings, passed as `opt` to the functions below.
#' @examples
#' opt <- leaf_area_options(min_object_cm2 = 0.05, max_sensitivity = 3)
leaf_area_options = function(
    ## Scanner lateral response
    paper_percentile   = 90,
    profile_smooth_px  = 5,
    min_illumination   = 0.7,
    border_depth_frac  = 0.02,
    border_aspect      = 10,
    ## Segmentation
    min_object_cm2     = 0.02,
    detect_px          = 800,
    detect_fraction    = 0.35,
    core_erode_px      = 3,
    recover_pale_interior = TRUE,
    texture_px         = 3,
    texture_ratio      = 4,
    texture_frac       = 0.5,
    interior_min_frac  = 0.05,
    ring_inner_px      = 12,
    ring_outer_px      = 60,
    max_iterations     = 12,
    sensitivity_frac   = 0.10,
    ## Quality control. These only populate `flags`; they never change an area.
    min_contrast       = 15,
    max_sensitivity    = 5.0,
    max_second_object  = 0.05,
    dpi_tolerance      = 0.05,
    dpi_match_tolerance = 0.02) {
  as.list(environment())
}


#' Tuning parameters in force for this run
#'
#' The default `opt` for every function below. Defined at the top level so that
#' a session which has sourced this file can inspect and change settings
#' interactively; pass `opt` explicitly to override it for one call.
#' @return A list; see [leaf_area_options()].
OPT = leaf_area_options()


## =============================================================================
## Image input
##
## Images come back from the tiff, jpeg and png packages as an array indexed
## [row, column, channel] with values in 0..1. Everything downstream works in
## 0..255 doubles instead, because the thresholds are easier to reason about and
## to compare against what you see in an image viewer.
## =============================================================================

#' Read a scan and reduce it to one gray channel
#'
#' @param path Path to a TIFF, JPEG or PNG.
#' @return A list with `gray`, a numeric matrix of 0..255 values indexed
#'   `[row, column]`, and `dpi`, the recorded resolution or `NA`.
#' @seealso [to_gray()] for the channel choice, [read_dpi()] for the metadata.
read_scan = function(path) {
  ext = tolower(tools::file_ext(path))

  img = switch(
    ext,
    tif  = ,
    tiff = read_with("tiff", tiff::readTIFF, path, info = TRUE),
    jpg  = ,
    jpeg = read_with("jpeg", jpeg::readJPEG, path),
    png  = read_with("png",  png::readPNG,   path),
    stop("unsupported file type: .", ext, call. = FALSE)
  )
  ## readTIFF(all = FALSE) returns a bare array, but a caller who changes that
  ## default gets a list of pages back; take the first and say so.
  if (is.list(img)) {
    warning(basename(path), " has several pages; using the first", call. = FALSE)
    img = img[[1L]]
  }

  list(gray = to_gray(img), dpi = read_dpi(path, img))
}

#' Call an image reader, with a useful error if its package is missing
#'
#' The three reader packages are only needed for the formats you actually have,
#' so they are loaded on demand rather than attached at the top of the file.
#'
#' @param pkg Package name, as a string.
#' @param reader The reader function itself.
#' @param path Path passed to `reader`.
#' @param ... Further arguments for `reader`.
#' @return Whatever `reader` returns.
read_with = function(pkg, reader, path, ...) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("package '", pkg, "' is needed to read ", basename(path), call. = FALSE)
  reader(path, ...)
}

#' Collapse an image to one 8-bit gray channel
#'
#' For a color scan this deliberately does NOT use luminance. Chlorophyll
#' absorbs red, so a green leaf on white paper separates far better in the red
#' channel than in a luminance mix that is 71 percent green. Rather than
#' hard-code that assumption, pick whichever channel puts the most distance
#' between the paper and the darkest few pixels; on a leaf that is not green,
#' or on a gray standard, the same rule still picks the most informative
#' channel. On a grayscale scan every channel is identical and this is a no-op.
#'
#' Channels are compared on a subsampled grid, since the comparison only has to
#' rank them, not measure them.
#'
#' @param img Numeric array `[row, column, channel]` with values in 0..1, or a
#'   bare `[row, column]` matrix for a grayscale file.
#' @return A numeric matrix of 0..255 values.
to_gray = function(img) {
  if (length(dim(img)) == 2L) return(round(img * 255))

  n = min(dim(img)[3L], 3L)           # ignore any alpha channel
  if (n == 1L) return(round(img[, , 1L] * 255))

  step = max(1L, as.integer(sqrt(prod(dim(img)[1:2]) / 2e5)))
  ri = seq(1L, dim(img)[1L], by = step)
  ci = seq(1L, dim(img)[2L], by = step)
  ## 0.2% guards against a dead pixel setting the dark end; 98% against a
  ## specular highlight setting the bright end.
  contrast = vapply(seq_len(n), function(k)
    as.numeric(diff(quantile(img[ri, ci, k], c(0.002, 0.98), names = FALSE))),
    numeric(1L))

  round(img[, , which.max(contrast)] * 255)
}

#' Horizontal resolution recorded in an image file
#'
#' Used two ways: to match a scan to a size standard taken at the same
#' resolution, and as an independent cross-check on that standard. It is
#' deliberately never the primary ruler, because a flatbed's true sampling pitch
#' can differ from its nominal resolution by a percent or two.
#'
#' @param path Path to the file, needed for the JPEG and PNG byte scans.
#' @param img The image as returned by the reader, whose attributes carry the
#'   TIFF resolution tags.
#' @return Dots per inch, or `NA_real_` if the file does not record it.
read_dpi = function(path, img) {
  dpi = switch(
    tolower(tools::file_ext(path)),
    tif = , tiff = {
      d = suppressWarnings(as.numeric(attr(img, "x.resolution")))
      if (identical(tolower(attr(img, "resolution.unit")), "centimeter"))
        d * 2.54 else d
    },
    jpg = , jpeg = jfif_dpi(path),
    png = png_dpi(path),
    NA_real_
  )
  if (!length(dpi) || !is.finite(dpi) || dpi <= 1) NA_real_ else as.numeric(dpi)
}

#' Find a byte marker near the start of a file
#'
#' Both the JFIF and pHYs headers sit within the first few hundred bytes, so
#' reading a 4 kB prefix and scanning it avoids pulling a multi-megabyte file
#' through memory just to read its resolution.
#'
#' @param path File to scan.
#' @param marker Raw vector to look for.
#' @param n Bytes to read from the head of the file.
#' @return A list with `bytes`, the prefix that was read, and `at`, the
#'   1-based index where `marker` starts, or `NA_integer_`.
find_marker = function(path, marker, n = 4096L) {
  bytes = readBin(path, "raw", n = n)
  for (i in which(bytes == marker[1L]))
    if (i + length(marker) - 1L <= length(bytes) &&
        identical(bytes[i + seq_along(marker) - 1L], marker))
      return(list(bytes = bytes, at = i))
  list(bytes = bytes, at = NA_integer_)
}

#' Resolution from a JPEG's JFIF header
#'
#' Layout after the marker `'J' 'F' 'I' 'F' 0x00`: version major, version minor,
#' density units, X density (2 bytes, big-endian), Y density (2 bytes). Units 1
#' means dots per inch, 2 means dots per centimeter.
#'
#' Many JPEGs carry no JFIF density at all, which is one reason to prefer TIFF
#' for archival scans.
#'
#' @param path Path to a JPEG.
#' @return Dots per inch, or `NA_real_`.
jfif_dpi = function(path) {
  f = find_marker(path, as.raw(c(0x4A, 0x46, 0x49, 0x46, 0x00)))
  if (is.na(f$at) || f$at + 9L > length(f$bytes)) return(NA_real_)
  units = as.integer(f$bytes[f$at + 7L])
  x     = as.integer(f$bytes[f$at + 8L]) * 256L + as.integer(f$bytes[f$at + 9L])
  if (x == 0L) NA_real_ else if (units == 1L) x else if (units == 2L) x * 2.54 else NA_real_
}

#' Resolution from a PNG's pHYs chunk
#'
#' Layout: the four bytes `pHYs`, then pixels per unit X (4 bytes, big-endian),
#' pixels per unit Y (4 bytes), then the unit specifier, where 1 means meters.
#'
#' @param path Path to a PNG.
#' @return Dots per inch, or `NA_real_`.
png_dpi = function(path) {
  f = find_marker(path, charToRaw("pHYs"))
  if (is.na(f$at) || f$at + 12L > length(f$bytes)) return(NA_real_)
  per_unit = sum(as.integer(f$bytes[f$at + 4:7]) * c(256^3, 256^2, 256, 1))
  if (as.integer(f$bytes[f$at + 12L]) == 1L && per_unit > 0)
    per_unit * 0.0254 else NA_real_
}


## =============================================================================
## Small matrix helpers
## =============================================================================

#' Shrink an image by taking the minimum of each block
#'
#' A mean or a median would average a one-pixel vein, or a needle narrower than
#' the block, away into the paper. The minimum keeps every dark structure
#' visible at the cost of pulling the background down slightly and uniformly,
#' which does not matter because the background level is re-estimated from the
#' shrunk image anyway.
#'
#' Implemented as a chain of `pmin()` calls rather than `apply()` over a
#' reshaped array, which is roughly fifty times faster on a full-page scan.
#'
#' @param m Numeric matrix.
#' @param f Integer shrink factor. `f <= 1` returns `m` unchanged.
#' @return A matrix with about `nrow(m)/f` rows and `ncol(m)/f` columns. Any
#'   partial block at the right or bottom edge is discarded.
downsample_min = function(m, f) {
  if (f <= 1L) return(m)
  h = (nrow(m) %/% f) * f
  w = (ncol(m) %/% f) * f
  if (h < f || w < f) return(m)
  m = m[seq_len(h), seq_len(w), drop = FALSE]

  top = seq(1L, h, by = f)                       # minimum down each block of f rows
  out = m[top, , drop = FALSE]
  for (k in seq_len(f - 1L)) out = pmin(out, m[top + k, , drop = FALSE])

  left = seq(1L, w, by = f)                      # then across each block of f columns
  m = out
  out = m[, left, drop = FALSE]
  for (k in seq_len(f - 1L)) out = pmin(out, m[, left + k, drop = FALSE])
  out
}

#' Erode a logical mask by a square structuring element
#'
#' `TRUE` only where the whole `(2r+1) x (2r+1)` window is `TRUE`. Separable, so
#' it runs as `2r` shifted ANDs rather than `(2r+1)^2` neighborhood lookups, and
#' stays vectorized in base R. Padding is `TRUE`, so the image border does not
#' eat into the mask.
#'
#' Used to find the leaf's interior, away from its soft edge, before reading its
#' tone.
#'
#' @param mask Logical matrix.
#' @param r Radius in pixels. Clamped so it cannot exceed the mask.
#' @return A logical matrix the same size as `mask`.
erode = function(mask, r) {
  h = nrow(mask); w = ncol(mask)
  r = min(as.integer(r), h - 1L, w - 1L)
  if (r <= 0L) return(mask)

  out = mask
  for (k in seq_len(r)) {                                   # horizontal pass
    left  = cbind(matrix(TRUE, h, k), mask[, seq_len(w - k), drop = FALSE])
    right = cbind(mask[, (k + 1L):w, drop = FALSE], matrix(TRUE, h, k))
    out = out & left & right
  }
  mask = out
  for (k in seq_len(r)) {                                   # vertical pass
    up   = rbind(matrix(TRUE, k, w), mask[seq_len(h - k), , drop = FALSE])
    down = rbind(mask[(k + 1L):h, , drop = FALSE], matrix(TRUE, k, w))
    out = out & up & down
  }
  out
}

#' Dilate a logical mask by a square structuring element
#'
#' The dual of [erode()], by De Morgan: growing the object is the same as
#' eroding its background.
#'
#' @param mask Logical matrix.
#' @param r Radius in pixels.
#' @return A logical matrix the same size as `mask`.
dilate = function(mask, r) !erode(!mask, r)

#' Close a logical mask: dilate, then erode by the same radius
#'
#' Bridges gaps up to about `2r` across without moving the outer boundary of
#' anything larger than the structuring element. Used to seal a leaf outline
#' that is interrupted, so that [fill_holes()] has something enclosed to fill.
#'
#' @param mask Logical matrix.
#' @param r Radius in pixels.
#' @return A logical matrix the same size as `mask`.
closing = function(mask, r) erode(dilate(mask, r), r)

#' Mean of each (2r+1) x (2r+1) window
#'
#' Separable, so it costs `2(2r+1)` shifted additions rather than `(2r+1)^2`
#' neighborhood lookups. Edge pixels are handled by replicating the border.
#'
#' @param m Numeric matrix.
#' @param r Window radius.
#' @return A numeric matrix the same size as `m`.
box_mean = function(m, r) {
  along = function(m, axis) {
    n = dim(m)[axis]
    acc = NULL
    for (k in -r:r) {
      j = pmin(pmax(seq_len(n) + k, 1L), n)          # replicate at the edges
      slab = if (axis == 2L) m[, j, drop = FALSE] else m[j, , drop = FALSE]
      acc = if (is.null(acc)) slab else acc + slab
    }
    acc / (2L * r + 1L)
  }
  along(along(m, 2L), 1L)
}

#' Local standard deviation over a square window
#'
#' The texture channel. Paper is optically smooth: on the scans this was built
#' for it measures a local SD of 0.4 to 0.6 gray levels, essentially the sensor
#' noise floor. Anything botanical is not: a leaf, and in particular a hairy or
#' tomentose one, measures 5 to 7, a separation of more than tenfold that holds
#' whether the surface is darker or brighter than the paper.
#'
#' That is why this exists. Intensity alone cannot find the white tomentum on
#' the underside of a *Rhododendron* leaf, because it is the same brightness as
#' the paper it is lying on. Texture can.
#'
#' @param m Numeric matrix.
#' @param r Window radius.
#' @return A numeric matrix of local standard deviations.
local_sd = function(m, r) {
  mu = box_mean(m, r)
  sqrt(pmax(box_mean(m * m, r) - mu * mu, 0))
}

#' Convex hull of a logical mask, as a filled region
#'
#' The outline of what was identified: the smallest convex region containing
#' every leaf pixel. Used by [recover_interior()] to define "inside the
#' outline", the region within which a pale lamina might be hiding.
#'
#' A convex hull is deliberately chosen over a morphological closing. A closing
#' needs a radius, and any radius large enough to bridge a gap in a leaf margin
#' is also large enough to bridge the notches of a serrate leaf, which squares
#' off the teeth of a *Betula* and quietly adds paper. The hull needs no radius
#' at all, and because it is only ever used as a place to LOOK -- never as
#' something to add -- its tendency to overshoot on a curved leaf costs nothing:
#' the overshoot is paper, and the texture test rejects paper.
#'
#' @param mask Logical matrix with at least one `TRUE`.
#' @return A logical matrix, `TRUE` inside the hull.
convex_hull_mask = function(mask) {
  edge = mask & !erode(mask, 1L)                    # boundary points suffice
  idx  = which(if (any(edge)) edge else mask, arr.ind = TRUE)
  h    = grDevices::chull(idx[, 2L], idx[, 1L])
  if (length(h) < 3L) return(mask)
  px = idx[h, 2L]; py = idx[h, 1L]; n = length(h)

  ## A convex polygon meets each scanline in a single interval, so the fill is
  ## just the running min and max of the edge crossings, row by row.
  ys  = seq.int(min(py), max(py))
  lo  = rep(Inf, length(ys)); hi = rep(-Inf, length(ys))
  for (i in seq_len(n)) {
    j = if (i == n) 1L else i + 1L
    if (py[i] == py[j]) {
      k = py[i] - ys[1L] + 1L
      lo[k] = min(lo[k], px[i], px[j]); hi[k] = max(hi[k], px[i], px[j])
    } else {
      yy = seq.int(min(py[i], py[j]), max(py[i], py[j]))
      xx = px[i] + (yy - py[i]) * (px[j] - px[i]) / (py[j] - py[i])
      k  = yy - ys[1L] + 1L
      lo[k] = pmin(lo[k], xx); hi[k] = pmax(hi[k], xx)
    }
  }
  out = matrix(FALSE, nrow(mask), ncol(mask))
  for (k in seq_along(ys))
    if (is.finite(lo[k]))
      out[ys[k], max(1L, floor(lo[k])):min(ncol(mask), ceiling(hi[k]))] = TRUE
  out
}

#' Bounding box of a logical mask
#'
#' @param mask Logical matrix with at least one `TRUE`.
#' @return An integer vector `c(first_row, last_row, first_column, last_column)`.
bounds = function(mask) {
  rows = which(rowSums(mask) > 0)
  cols = which(colSums(mask) > 0)
  c(rows[1L], rows[length(rows)], cols[1L], cols[length(cols)])
}

#' Estimate the paper level from the histogram mode
#'
#' The background is by far the most common tone in a leaf scan, so the mode is
#' a much steadier estimate of "paper" than a mean, which the leaf drags down,
#' or a median, which moves in whole gray levels. Refined to a counts-weighted
#' mean over the eleven levels around the peak, so the answer is continuous
#' rather than quantized.
#'
#' @param gray Numeric matrix of 0..255 values.
#' @return The paper level, as a double.
background_level = function(gray) {
  counts = tabulate(as.integer(gray) + 1L, nbins = 256L)
  peak   = which.max(counts) - 1L
  near   = max(0L, peak - 5L):min(255L, peak + 5L)
  sum(near * counts[near + 1L]) / sum(counts[near + 1L])
}


## =============================================================================
## Connected components
## =============================================================================

#' Label the connected components of a binary mask
#'
#' Run-length encoded with union-find: each row is reduced to its runs of
#' `TRUE`, runs are joined to overlapping runs in the row above, and the
#' resulting forest is flattened with path compression. Storing runs rather than
#' pixels is what keeps this tractable in base R on a 16 megapixel scan.
#'
#' Objects are labeled with 8-connectivity by default. This matters more than it
#' sounds: a needle lying diagonally across the glass is a staircase of single
#' pixels, and under 4-connectivity it shatters into hundreds of fragments, none
#' of which is recognizable as a leaf. Background, in [fill_holes()], is labeled
#' with 4-connectivity instead, which is the complement convention and keeps the
#' two topologically consistent.
#'
#' @param mask Logical matrix.
#' @param diagonal If `TRUE`, runs touching only at a corner are joined
#'   (8-connectivity). If `FALSE`, they must overlap in column (4-connectivity).
#' @return `NULL` if `mask` is entirely `FALSE`, otherwise a list describing the
#'   labeling:
#'   \describe{
#'     \item{row, start, end}{one entry per run: its row, and its first and last
#'       column}
#'     \item{dim}{dimensions of the original mask}
#'     \item{runs_of}{list of integer vectors, one per component, indexing into
#'       `row`/`start`/`end`}
#'     \item{size}{pixel count of each component, in the order of `runs_of`}
#'     \item{on_edge}{whether each component touches the image border}
#'   }
label_components = function(mask, diagonal = TRUE) {
  height = nrow(mask)
  width  = ncol(mask)
  slack  = if (diagonal) 1L else 0L

  ## Runs, from the transitions in each row. Padding with FALSE on both sides
  ## makes a run that starts at column 1 or ends at the last column ordinary.
  per_row = vector("list", height)
  for (y in seq_len(height)) {
    edges = diff(c(FALSE, mask[y, ], FALSE))
    opens = which(edges == 1L)
    if (length(opens))
      per_row[[y]] = cbind(y, opens, which(edges == -1L) - 1L)
  }
  runs = do.call(rbind, per_row)
  if (is.null(runs)) return(NULL)

  row = runs[, 1L]; start = runs[, 2L]; end = runs[, 3L]
  n = length(row)

  ## Union-find over runs, with path compression. `<<-` reaches the `parent`
  ## vector in this function's frame, so compression persists between calls.
  parent = seq_len(n)
  find = function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i = parent[i]
    }
    i
  }

  ## Join each run to any run in the row above it overlaps horizontally, or
  ## touches only at a corner when `diagonal` is TRUE. Rows with no runs break
  ## the chain, hence the consecutive-row check.
  by_row = split(seq_len(n), row)
  row_no = as.integer(names(by_row))
  for (k in seq_along(by_row)[-1L]) {
    if (row_no[k] != row_no[k - 1L] + 1L) next
    above = by_row[[k - 1L]]
    for (j in by_row[[k]]) {
      for (i in above[start[above] <= end[j] + slack & start[j] <= end[above] + slack]) {
        ri = find(i); rj = find(j)
        if (ri != rj) parent[ri] = rj
      }
    }
  }

  runs_of  = split(seq_len(n), vapply(seq_len(n), find, integer(1L)))
  on_edge  = row == 1L | row == height | start == 1L | end == width
  run_len  = end - start + 1L

  list(row = row, start = start, end = end, dim = c(height, width),
       runs_of = runs_of,
       size    = vapply(runs_of, function(k) sum(run_len[k]), numeric(1L)),
       on_edge = vapply(runs_of, function(k) any(on_edge[k]), logical(1L)))
}

#' Paint one component into a matrix just big enough to hold it
#'
#' [draw_runs()] allocates a matrix the size of the whole window, which is
#' ruinous when a hull leaves hundreds of small slivers to be examined one by
#' one. This allocates only the component's own bounding box.
#'
#' @param cc A labeling from [label_components()].
#' @param k Index of the component.
#' @return A list with the bbox-sized logical `mask` and the `rows` and `cols`
#'   it occupies in the original.
region_box = function(cc, k) {
  r  = cc$runs_of[[k]]
  r0 = min(cc$row[r]);   r1 = max(cc$row[r])
  c0 = min(cc$start[r]); c1 = max(cc$end[r])
  out = matrix(FALSE, r1 - r0 + 1L, c1 - c0 + 1L)
  for (i in r)
    out[cc$row[i] - r0 + 1L, (cc$start[i] - c0 + 1L):(cc$end[i] - c0 + 1L)] = TRUE
  list(mask = out, rows = r0:r1, cols = c0:c1)
}

#' Paint a set of runs into a fresh logical matrix
#'
#' @param cc A labeling from [label_components()].
#' @param which_runs Integer vector of run indices, typically
#'   `cc$runs_of[[k]]` for one component.
#' @return A logical matrix of `cc$dim`, `TRUE` on those runs.
draw_runs = function(cc, which_runs) {
  out = matrix(FALSE, cc$dim[1L], cc$dim[2L])
  for (k in which_runs) out[cc$row[k], cc$start[k]:cc$end[k]] = TRUE
  out
}

#' Fill enclosed holes in a mask
#'
#' Any background component that does not reach the image border is enclosed by
#' the object, so it is a hole. Insect damage and pale patches are therefore
#' counted as leaf, which is the right convention for one-sided leaf area: the
#' hole is still inside the outline. The caller can recover the unfilled area
#' from the reported `holes_px`.
#'
#' @param mask Logical matrix.
#' @return `mask` with enclosed background set to `TRUE`.
fill_holes = function(mask) {
  bg = label_components(!mask, diagonal = FALSE)
  if (is.null(bg)) return(mask)
  holes = unlist(bg$runs_of[!bg$on_edge], use.names = FALSE)
  if (length(holes)) mask | draw_runs(bg, holes) else mask
}

#' Is this component the scanner's edge artifact rather than a specimen?
#'
#' Identified by position and shape, deliberately not by "touches an edge", so
#' that a leaf lying against the edge of the glass is still measured. What the
#' lid shadow does, and a specimen does not, is hug one border while barely
#' penetrating: on the scans this was built for it runs the full 4,676 rows and
#' reaches 14 to 20 columns in, an aspect ratio of roughly 300:1 pinned to the
#' border.
#'
#' This is a backstop, not the main defense. [flatten_scanner()] has already
#' divided the shadow out and dropped the columns the lid effectively blinded;
#' what survives to be caught here is a residual sliver of amplified noise.
#'
#' @param cc A labeling from [label_components()].
#' @param k Index of the component to test.
#' @param opt Options from [leaf_area_options()].
#' @return `TRUE` if the component looks like the lid shadow.
edge_artifact = function(cc, k, opt) {
  r  = cc$runs_of[[k]]
  h  = cc$dim[1L]; w = cc$dim[2L]
  r0 = min(cc$row[r]);   r1 = max(cc$row[r])
  c0 = min(cc$start[r]); c1 = max(cc$end[r])
  tall = r1 - r0 + 1L
  wide = c1 - c0 + 1L

  (((c0 == 1L || c1 == w) &&                      # against a left/right border
    wide <= opt$border_depth_frac * w &&
    tall >= opt$border_aspect * wide) ||
   ((r0 == 1L || r1 == h) &&                      # against a top/bottom border
    tall <= opt$border_depth_frac * h &&
    wide >= opt$border_aspect * tall))
}

#' Largest component that is not the scanner's edge artifact
#'
#' @param cc A labeling from [label_components()], or `NULL`.
#' @param min_px Smallest component worth considering, in pixels.
#' @param opt Options from [leaf_area_options()].
#' @return The index of the chosen component, or `NA_integer_` if nothing
#'   reaches `min_px`. If every candidate looks like an artifact the largest is
#'   returned anyway, on the grounds that a wrong answer the caller can see in
#'   the mask beats no answer at all.
largest_specimen = function(cc, min_px, opt) {
  if (is.null(cc)) return(NA_integer_)
  ok = which(cc$size >= min_px)
  if (!length(ok)) return(NA_integer_)
  real = ok[!vapply(ok, function(k) edge_artifact(cc, k, opt), logical(1L))]
  if (!length(real)) real = ok
  real[which.max(cc$size[real])]
}

#' Component containing a given pixel
#'
#' Used to keep the refinement locked onto the same object between iterations,
#' rather than re-picking the largest each time and risking a jump to a
#' different one as the cutoff moves.
#'
#' @param cc A labeling from [label_components()].
#' @param seed Length-2 integer vector, `c(row, column)`.
#' @param min_px Smallest component worth considering, for the fallback.
#' @param opt Options from [leaf_area_options()].
#' @return The index of the component covering `seed`, or, if the cutoff moved
#'   far enough to drop that pixel, of the largest specimen.
component_at = function(cc, seed, min_px, opt) {
  hit = which(cc$row == seed[1L] & cc$start <= seed[2L] & cc$end >= seed[2L])
  if (length(hit)) {
    owner = which(vapply(cc$runs_of, function(k) hit[1L] %in% k, logical(1L)))
    if (length(owner)) return(owner[1L])
  }
  largest_specimen(cc, min_px, opt)
}

#' Count the objects on the plate and size up the runner-up
#'
#' Run on the coarse detection image rather than on the refinement window,
#' because refinement only ever sees a neighborhood around the leaf and would
#' miss a second leaf or a label sitting off in a corner.
#'
#' Note the approximation: the coarse image is a block minimum, so it is
#' systematically darker than the full-resolution scan, and applying a
#' full-resolution cutoff to it over-selects slightly. The census therefore errs
#' toward reporting too many objects rather than too few, which is the safe
#' direction for a warning flag.
#'
#' @param mask Logical matrix, the coarse image thresholded at the final cutoff.
#' @param min_px Smallest object to count, in coarse pixels.
#' @param opt Options from [leaf_area_options()].
#' @return A list with `n_objects` and `second_frac`, the size of the runner-up
#'   as a fraction of the largest.
object_census = function(mask, min_px, opt) {
  cc = label_components(mask)
  if (is.null(cc)) return(list(n_objects = 0L, second_frac = 0))
  big = which(cc$size >= max(1, min_px))
  big = big[!vapply(big, function(k) edge_artifact(cc, k, opt), logical(1L))]
  sizes = sort(cc$size[big], decreasing = TRUE)
  list(n_objects   = length(sizes),
       second_frac = if (length(sizes) > 1L) sizes[2L] / sizes[1L] else 0)
}


## =============================================================================
## Scanner lateral response
##
## A flatbed's lid casts a dark band along an edge of the platen, and the
## illumination is never perfectly even across the bed. On the scanner this was
## built for, the band is confined to the last 14 to 20 columns on the right and
## falls to gray 106 against paper at 245; the other three edges are clean to
## within 3 levels.
##
## The artifact is SEPARABLE: it is a function of position across the bed and is
## essentially constant along the scan direction. So it can be estimated as a
## 1-D profile instead of a 2-D surface, which is far better conditioned and
## needs no assumption about where the leaf is. This is the same structure the
## radiochromic film dosimetry literature characterizes for flatbeds, where the
## lateral response artifact appears even for dose-independent optical subjects
## and is therefore a property of the scanner rather than of the sample.
##
## Correcting rather than cropping keeps leaves lying near the edge of the
## glass. On a synthetic needle placed 6 px from the platen edge, the raw image
## measures the shadow instead of the needle (+159%); after this correction the
## needle comes back to within 0.02%. A broad leaf running off the edge of the
## bed goes from -3.75% under a fixed crop to -0.000% here.
## =============================================================================

#' Paper level along one axis of a scan
#'
#' A high percentile is the paper level at that position: it is untouched by a
#' leaf covering up to `100 - opt$paper_percentile` percent of the slice, which
#' is what makes this work without knowing where the leaf is. Subsampled along
#' the perpendicular axis, since this is a percentile and not a sum.
#'
#' @param gray Numeric matrix of 0..255 values.
#' @param axis 1 for a per-row profile, 2 for a per-column profile.
#' @param opt Options from [leaf_area_options()].
#' @return A numeric vector, one entry per row (`axis = 1`) or column
#'   (`axis = 2`).
paper_profile = function(gray, axis, opt = OPT) {
  along = dim(gray)[3L - axis]                    # the axis being subsampled
  step  = max(1L, as.integer(along %/% 1200L))
  keep  = seq(1L, along, by = step)
  slab  = if (axis == 2L) gray[keep, , drop = FALSE] else gray[, keep, drop = FALSE]

  prof = apply(slab, axis, quantile,
               probs = opt$paper_percentile / 100, names = FALSE)
  if (opt$profile_smooth_px >= 3L && length(prof) > opt$profile_smooth_px)
    prof = runmed(prof, opt$profile_smooth_px, endrule = "keep")
  prof
}

#' Divide out the lateral response along one axis
#'
#' The correction is multiplicative because the lid attenuates the light
#' reaching the sensor: it scales the leaf and the paper by the same factor, so
#' contrast is restored rather than merely offset. A subtractive correction
#' would flatten the background but leave a leaf sitting in the shadow with its
#' contrast still crushed.
#'
#' Slices receiving less than `opt$min_illumination` of the brightest light are
#' left uncorrected and marked dead. Past that point, correcting means
#' multiplying noise by more than 1.4 and there is no signal worth recovering;
#' dropping four columns is the data-driven replacement for a blind fixed crop
#' of forty.
#'
#' @param gray Numeric matrix of 0..255 values.
#' @param axis 1 for rows, 2 for columns.
#' @param opt Options from [leaf_area_options()].
#' @return A list with the corrected `gray` and a logical `live`, `TRUE` for
#'   slices bright enough to keep.
flatten_axis = function(gray, axis, opt = OPT) {
  prof  = paper_profile(gray, axis, opt)
  level = median(prof)
  live  = prof >= opt$min_illumination * level
  ## An all-dark image means the profile is meaningless, not that every slice
  ## should be discarded; leave it alone and let segmentation report the truth.
  if (!any(live)) return(list(gray = gray, live = rep(TRUE, length(prof))))

  gain = level / pmax(prof, 1)
  gain[!live] = 1
  gray = if (axis == 2L) gray * rep(gain, each = nrow(gray))
         else            gray * gain
  list(gray = pmin(gray, 255), live = live)
}

#' Correct the scanner's lateral response and drop the dead slices
#'
#' Columns first, then rows, so a scanner whose artifact runs the other way is
#' handled without configuration. Reports which rows and columns survived, so
#' the mask can be written back in the coordinates of the original scan.
#'
#' @param gray Numeric matrix of 0..255 values.
#' @param opt Options from [leaf_area_options()].
#' @return A list with the corrected and trimmed `gray`, the surviving `rows`
#'   and `cols` as indices into the original, and `dropped`, the total number of
#'   rows plus columns discarded.
flatten_scanner = function(gray, opt = OPT) {
  by_col = flatten_axis(gray, 2L, opt); gray = by_col$gray
  by_row = flatten_axis(gray, 1L, opt); gray = by_row$gray
  rows = which(by_row$live)
  cols = which(by_col$live)
  list(gray    = gray[rows, cols, drop = FALSE],
       rows    = rows,
       cols    = cols,
       dropped = (length(by_row$live) - length(rows)) +
                 (length(by_col$live) - length(cols)))
}

#' Put a mask back into the coordinates of the original scan
#'
#' @param mask Logical matrix measured on the trimmed image.
#' @param trim The list returned by [flatten_scanner()].
#' @param dims Dimensions of the original scan, `c(rows, columns)`.
#' @return A logical matrix of `dims`.
restore_mask = function(mask, trim, dims) {
  out = matrix(FALSE, dims[1L], dims[2L])
  out[trim$rows, trim$cols] = mask
  out
}


## =============================================================================
## Segmentation
##
## Otsu's method is the textbook choice here and it is the wrong tool for this
## picture. Otsu maximizes between-class variance, which assumes the two classes
## are of broadly comparable size; in these scans the leaf is 0.1% of the pixels,
## so the criterion is dominated by whatever split of the BACKGROUND maximizes
## the separation, and the cutoff slides into the paper. Measured: trimming 40
## irrelevant pixels off the border of one scan moved Otsu's cutoff from 190 to
## 205 and the measured area by 24%.
##
## Instead: find the leaf first, then put the cutoff halfway between the leaf's
## own tone and the paper immediately around it. For an edge blurred by a
## symmetric point spread function the half-maximum crossing sits at the true
## boundary, which is the same convention as ISO50 in dimensional metrology, and
## it is anchored to two physical levels rather than to the shape of a histogram
## -- so it does not care how much paper surrounds the leaf.
##
## The known limit is narrow features: when an object is not much wider than the
## blur, its interior never reaches its true tone, the cutoff sits too high and
## the object fattens. That is why resolution has to be set by the width of the
## narrowest object rather than by leaf area, and why `area_sens_pct` matters.
## =============================================================================

#' Paper tone in a frame around the leaf
#'
#' Read locally rather than from the whole scan, so that any shading left across
#' the glass cannot drag the cutoff around. Only pixels brighter than the
#' current cutoff count, so a second object sitting in the frame cannot darken
#' the estimate. Falls back to the global mode when the frame is too small to be
#' reliable, which happens when the leaf nearly fills its window.
#'
#' @param roi Numeric matrix, the refinement window.
#' @param mask Logical matrix, the leaf as currently segmented.
#' @param cut Current cutoff, used to exclude dark pixels from the frame.
#' @param opt Options from [leaf_area_options()].
#' @return The paper level, as a double.
paper_level = function(roi, mask, cut, opt = OPT) {
  b = bounds(mask)
  rows = max(1L, b[1L] - opt$ring_outer_px):min(nrow(roi), b[2L] + opt$ring_outer_px)
  cols = max(1L, b[3L] - opt$ring_outer_px):min(ncol(roi), b[4L] + opt$ring_outer_px)
  side = setdiff(cols, (b[3L] - opt$ring_inner_px):(b[4L] + opt$ring_inner_px))
  cap  = setdiff(rows, (b[1L] - opt$ring_inner_px):(b[2L] + opt$ring_inner_px))

  band = c(if (length(side)) roi[rows, side] else numeric(0),
           if (length(cap))  roi[cap,  cols] else numeric(0))
  band = band[band > cut]
  if (length(band) < 200) background_level(roi) else median(band)
}

#' Recover a leaf interior that thresholding could not see
#'
#' Some leaves are not uniformly darker than the paper. *Rhododendron
#' groenlandicum* scanned abaxial side up is the case this was written for: the
#' margins are revolute and roll under, so they read dark, while the lamina
#' between them is covered in dense white tomentum measuring 243 to 248 gray
#' against paper at 245. It is not darker than the paper at all, so no cutoff
#' can find it, and no amount of contrast adjustment will change that.
#'
#' Hole filling ought to rescue it, and sometimes does -- but only when the dark
#' margin happens to close into a loop. On the scan that prompted this the
#' margin is interrupted at the leaf base, so the pale interior drains to the
#' background and is not an enclosed hole. Worse, whether it closes turns on a
#' single gray level: one level either side of the converged cutoff is the
#' difference between 11,670 px and 18,355 px, a 57 percent swing. That
#' knife-edge is the defect; moving the cutoff would not remove it.
#'
#' What separates the tomentum from the paper is not brightness but **texture**.
#' Paper is optically smooth and measures a local standard deviation of 0.43 to
#' 0.49 gray levels, essentially the sensor noise floor. The tomentum measures
#' 5.2. That is more than a tenfold separation, and it holds whether the leaf
#' surface is darker or brighter than the paper.
#'
#' So: look inside the convex hull of what was already identified, take each
#' region there that is not yet leaf, and ask whether it is textured like a
#' specimen or smooth like paper.
#'
#' @section Why this cannot damage a correct measurement:
#' It only ever looks INSIDE the outline of what the validated pipeline already
#' found, and it only ever ADDS. Three independent guards decide what is added:
#' \itemize{
#'   \item Regions are judged well clear of their own rim, since a texture
#'     window straddling the leaf margin reads as textured whatever lies inside
#'     it. A region too thin to have an interior is refused rather than guessed
#'     at -- which is what keeps the notches between the teeth of a serrate
#'     *Betula* leaf from being swallowed.
#'   \item The verdict is texture, so a genuine hole in a leaf, a bay between
#'     two teeth, and the slack a convex hull leaves around a curved needle are
#'     all correctly read as paper and refused. Measured on the six scans that
#'     were already right, every such region reads 0.48 against a paper 0.45,
#'     and none is accepted.
#'   \item `recovered_px` is reported in the output and sets the
#'     `interior_recovered` flag, so an addition is never silent.
#' }
#'
#' @section Known limitation:
#' A leaf with a genuinely smooth pale interior -- glaucous, waxy, strongly
#' specular -- reads like paper to this feature and will be refused. Scanning in
#' color is the fix for that case: a pale lamina is rarely the same *hue* as
#' white paper even when it is the same brightness.
#'
#' @param roi Numeric matrix, the refinement window.
#' @param mask Logical matrix, the leaf as found by thresholding.
#' @param cut The converged cutoff, used to locate paper for the texture
#'   reference.
#' @param opt Options from [leaf_area_options()]. Set
#'   `recover_pale_interior = FALSE` to disable this step entirely.
#' @return A list with the possibly enlarged `mask` and `recovered_px`.
recover_interior = function(roi, mask, cut, opt = OPT) {
  none = list(mask = mask, recovered_px = 0)
  if (!isTRUE(opt$recover_pale_interior)) return(none)

  ## Work in the leaf's own bounding box. The refinement window can be the whole
  ## scan when a leaf runs off the bed, and there is nothing to find out there.
  b = bounds(mask)
  pad  = 2L * opt$texture_px + 2L
  rows = max(1L, b[1L] - pad):min(nrow(mask), b[2L] + pad)
  cols = max(1L, b[3L] - pad):min(ncol(mask), b[4L] + pad)
  win  = mask[rows, cols, drop = FALSE]
  roi  = roi[rows, cols, drop = FALSE]

  hull   = convex_hull_mask(win)
  inside = hull & !win

  ## A lamina worth recovering is a substantial part of the leaf, not a sliver.
  ## Requiring that up front is also what keeps this affordable: the hull of a
  ## curved needle leaves hundreds of small staircase gaps, and examining each
  ## one costs far more than the answer is worth.
  floor_px = max(50, opt$interior_min_frac * sum(win))
  if (sum(inside) < floor_px) return(none)

  sd = local_sd(roi, opt$texture_px)
  paper = !dilate(hull, 2L * opt$texture_px) & roi > cut
  if (sum(paper) < 500) return(none)
  hot = opt$texture_ratio * median(sd[paper])

  cc = label_components(inside)
  if (is.null(cc)) return(none)

  keep = matrix(FALSE, nrow(win), ncol(win))
  for (k in seq_along(cc$size)) {
    if (cc$size[[k]] < floor_px) next
    reg  = region_box(cc, k)
    core = erode(reg$mask, 2L * opt$texture_px)
    if (sum(core) < 50) next                        # too thin to judge: refuse
    patch = sd[reg$rows, reg$cols, drop = FALSE]
    if (mean(patch[core] > hot) > opt$texture_frac)
      keep[reg$rows, reg$cols] = keep[reg$rows, reg$cols] | reg$mask
  }
  if (!any(keep)) return(none)

  out = mask
  out[rows, cols] = win | keep
  list(mask = out, recovered_px = sum(keep))
}

#' Iterate the cutoff to the half-maximum between leaf and paper
#'
#' Read the leaf's tone and the paper around it, put the cutoff halfway between
#' them, re-cut, repeat. Converges in a handful of passes because neither tone
#' moves much once the mask is roughly right, and the fixed point does not
#' depend on where the iteration started.
#'
#' @param roi Numeric matrix, the window around the leaf.
#' @param cut Starting cutoff.
#' @param min_px Smallest object worth considering, in pixels.
#' @param opt Options from [leaf_area_options()].
#' @return `NULL` if the object was lost, otherwise a list with the converged
#'   `cutoff`, the filled `mask`, a `seed` pixel inside it, the `roi` itself,
#'   the `paper` and `leaf_tone` levels and their difference as `contrast`,
#'   `holes_px`, and whether the object touches the window border.
refine_cutoff = function(roi, cut, min_px, opt = OPT) {
  last = NULL
  for (iteration in seq_len(opt$max_iterations)) {
    cc = label_components(roi <= cut)
    if (is.null(cc)) return(NULL)

    ## First pass picks the biggest thing; later passes stay locked on it.
    pick = if (is.null(last)) largest_specimen(cc, min_px, opt) else
             component_at(cc, last$seed, min_px, opt)
    if (is.na(pick)) return(NULL)

    mask = fill_holes(draw_runs(cc, cc$runs_of[[pick]]))
    holes_px = sum(mask) - cc$size[[pick]]

    ## The leaf's own tone, read away from its soft edge. The middle element of
    ## the eroded interior is a point that stays inside between iterations.
    core = erode(mask, opt$core_erode_px)
    if (sum(core) < 25) core = mask
    seed_idx = which(core, arr.ind = TRUE)
    seed = seed_idx[ceiling(nrow(seed_idx) / 2), ]

    leaf_tone = median(roi[core])
    paper = paper_level(roi, mask, cut, opt)

    ## The leaf must be darker than the paper somewhere, so a cutoff at or above
    ## the paper level cannot be meaningful. Without this bound, a leaf whose
    ## interior is paler than its margins and whose outline happens to close can
    ## poison its own tone estimate: hole filling hands the pale interior to
    ## `leaf_tone`, the cutoff climbs toward the paper, and the next pass calls
    ## most of the plate leaf. The clamp costs nothing on a normal scan -- the
    ## measured cutoffs here sit 30 to 50 gray levels below it -- and turns a
    ## runaway into an honest `low_contrast` flag.
    new_cut = min(round((paper + leaf_tone) / 2), paper - opt$min_contrast)
    last = list(cutoff = new_cut, mask = mask, seed = seed, roi = roi,
                paper = paper, leaf_tone = leaf_tone,
                contrast = paper - leaf_tone, holes_px = holes_px,
                touches_edge = cc$on_edge[[pick]])
    if (new_cut == cut) break
    cut = new_cut
  }
  last
}

#' Re-measure the same object at a different cutoff
#'
#' The other half of the error bar. Locked to the same seed pixel so that a
#' nudged cutoff cannot silently switch to measuring a different object.
#'
#' @param roi Numeric matrix, the refinement window.
#' @param cut Cutoff to re-measure at.
#' @param seed Length-2 integer vector identifying the object.
#' @param min_px Smallest object worth considering, in pixels.
#' @param opt Options from [leaf_area_options()].
#' @return Pixel count, or `NA_real_` if the object could not be found.
recount = function(roi, cut, seed, min_px, opt = OPT) {
  cc = label_components(roi <= cut)
  if (is.null(cc)) return(NA_real_)
  pick = component_at(cc, seed, min_px, opt)
  if (is.na(pick)) return(NA_real_)
  mask = fill_holes(draw_runs(cc, cc$runs_of[[pick]]))
  sum(recover_interior(roi, mask, cut, opt)$mask)
}

#' Measure the leaf in one flattened scan
#'
#' Four stages: locate the leaf on a coarse copy, refine the cutoff at full
#' resolution inside its own neighborhood, take a census of everything else on
#' the plate, and measure how much the answer depends on the cutoff.
#'
#' Working coarse first is what makes this fast enough in base R: detection only
#' has to answer "roughly where?", the block-minimum shrink suppresses paper
#' noise without losing thin objects, and the expensive full-resolution labeling
#' then runs on a window a few hundred pixels across rather than on 16 megapixels.
#'
#' @param gray Numeric matrix of 0..255 values, already passed through
#'   [flatten_scanner()].
#' @param min_px Smallest object reported as a leaf, in pixels.
#' @param opt Options from [leaf_area_options()].
#' @return A list with the full-size `mask`, the converged `cutoff`, the `paper`
#'   and `leaf_tone` levels and their `contrast`, `leaf_px`, `holes_px`,
#'   `n_objects` and `second_frac` from the census, `touches_edge`, and
#'   `leaf_px_lo` / `leaf_px_hi` bracketing the area.
#' @section Errors:
#' Stops with "no object stands out from the background" on a blank scan, and
#' with "no object found above the minimum size" when everything present is
#' smaller than `min_px`. Both are caught per scan by [run_leaf_area()], so one
#' bad file does not end the run.
segment_leaf = function(gray, min_px = 50, opt = OPT) {

  ## --- 1. Locate the leaf on a coarse copy ----------------------------------
  shrink = max(1L, as.integer(ceiling(max(dim(gray)) / opt$detect_px)))
  small  = downsample_min(gray, shrink)

  ## How dark does this scan get? Not a fixed quantile: the smallest leaf in the
  ## validation set is 0.016% of the image, so any percentile coarse enough to be
  ## robust is also coarse enough to miss it entirely. Take instead the k-th
  ## darkest pixel, where k is a fraction of the smallest object worth
  ## reporting. That shrugs off a few dead pixels while staying sensitive to a
  ## real speck.
  bg_small   = background_level(small)
  k          = min(max(5L, as.integer(min_px / shrink^2 / 4)), length(small))
  dark_small = sort(small, partial = k)[k]
  contrast   = bg_small - dark_small
  if (!is.finite(contrast) || contrast < 8)
    stop("no object stands out from the background", call. = FALSE)

  seed_cut = bg_small - opt$detect_fraction * contrast
  cc_small = label_components(small <= seed_cut)
  pick     = largest_specimen(cc_small, min_px / shrink^2, opt)
  if (is.na(pick)) stop("no object found above the minimum size", call. = FALSE)

  ## Coarse bounding box, mapped back to full-resolution coordinates. The extra
  ## `shrink` on the top and left undoes the block-minimum's 1-based offset.
  seed_runs = cc_small$runs_of[[pick]]
  box = c(min(cc_small$row[seed_runs]),   max(cc_small$row[seed_runs]),
          min(cc_small$start[seed_runs]), max(cc_small$end[seed_runs])) * shrink
  box = c(box[1L] - shrink, box[2L], box[3L] - shrink, box[4L])

  ## --- 2. Refine at full resolution, inside the leaf's own neighborhood -----
  margin = opt$ring_outer_px + 10L
  cut    = bg_small - 0.5 * contrast
  repeat {
    rows = max(1L, box[1L] - margin):min(nrow(gray), box[2L] + margin)
    cols = max(1L, box[3L] - margin):min(ncol(gray), box[4L] + margin)
    roi  = gray[rows, cols, drop = FALSE]

    fit = refine_cutoff(roi, cut, min_px, opt)
    if (is.null(fit)) stop("lost the object while refining the cutoff", call. = FALSE)

    ## If the leaf reaches the edge of the window, the window was too small and
    ## part of the leaf is outside it. Grow and try again, unless the window is
    ## already the whole image, in which case the leaf really does run off the
    ## measured region and `touches_edge` will say so.
    grew = any(fit$mask[1L, ]) || any(fit$mask[nrow(fit$mask), ]) ||
           any(fit$mask[, 1L]) || any(fit$mask[, ncol(fit$mask)])
    at_image_edge = rows[1L] == 1L && cols[1L] == 1L &&
                    rows[length(rows)] == nrow(gray) &&
                    cols[length(cols)] == ncol(gray)
    if (!grew || at_image_edge) break

    b = bounds(fit$mask)
    box = c(rows[b[1L]], rows[b[2L]], cols[b[3L]], cols[b[4L]])
    margin = margin * 2L
    cut = fit$cutoff
  }

  ## --- 3. Recover an interior that thresholding could not see ---------------
  rec = recover_interior(fit$roi, fit$mask, fit$cutoff, opt)

  ## --- 4. What else is on the plate? ----------------------------------------
  census = object_census(small <= fit$cutoff, min_px / shrink^2, opt)

  ## --- 5. How much does the answer depend on the cutoff? --------------------
  ## Move the cutoff by a fraction of the leaf-to-paper contrast, each way, and
  ## re-measure. This is an honest error bar: it is large exactly when the leaf
  ## margin is soft, translucent or out of focus, which is when the area is
  ## least defensible. Treat it as an ordinal warning rather than a calibrated
  ## interval -- for objects near the resolution limit it under-reports.
  wiggle = max(1, round(opt$sensitivity_frac * fit$contrast))
  lo = recount(fit$roi, fit$cutoff - wiggle, fit$seed, min_px, opt)
  hi = recount(fit$roi, fit$cutoff + wiggle, fit$seed, min_px, opt)

  full = matrix(FALSE, nrow(gray), ncol(gray))
  full[rows, cols] = rec$mask

  list(mask         = full,
       cutoff       = fit$cutoff,
       paper        = fit$paper,
       leaf_tone    = fit$leaf_tone,
       contrast     = fit$contrast,
       leaf_px      = sum(rec$mask),
       holes_px     = fit$holes_px,
       recovered_px = rec$recovered_px,
       n_objects    = census$n_objects,
       second_frac  = census$second_frac,
       touches_edge = fit$touches_edge,
       leaf_px_lo   = lo,
       leaf_px_hi   = hi)
}

#' Write a mask as a lossless PNG, leaf black on white
#'
#' Black on white to match the scans it came from, and PNG rather than JPEG so
#' that re-counting the mask reproduces the reported area exactly.
#'
#' @param mask Logical matrix at the size of the original scan.
#' @param path Destination path.
#' @return `path`, invisibly.
write_mask = function(mask, path) {
  if (!requireNamespace("png", quietly = TRUE))
    stop("package 'png' is needed to write masks; set write_masks = FALSE ",
         "to skip them", call. = FALSE)
  png::writePNG(ifelse(mask, 0, 1), path)
  invisible(path)
}


## =============================================================================
## Calibration
##
## The physical standard, not the DPI tag, is the ruler. A flatbed's true
## sampling pitch can differ from its nominal resolution by a percent or two,
## and the standard absorbs that: on the scanner this was built for, two
## independent standards at two different resolutions both measured 0.67% off
## nominal, which is a real property of the hardware and not a measurement error.
##
## But a standard is only a ruler for scans taken at ITS OWN resolution. Using a
## 300 dpi standard on 400 dpi scans inflates every area by (400/300)^2 = 78%.
## That was a real bug in the first version of this script, so matching is now
## enforced rather than assumed.
## =============================================================================

#' Measure every size standard and derive its pixel scale
#'
#' Each standard goes through exactly the same pipeline as a leaf, so that any
#' systematic bias in boundary placement partly cancels between the two.
#'
#' @param standards Data frame with columns `path`, `area_cm2` and `dpi`. Use
#'   `dpi = NA` to read the resolution from the file, and an explicit value for
#'   files that do not record one, which is common for JPEGs.
#' @param opt Options from [leaf_area_options()].
#' @param quiet Suppress the per-standard progress messages.
#' @return A data frame with `standard`, `dpi` and `cm2_per_px`, one row per
#'   standard.
#' @section Errors:
#' Stops if a standard is missing, if its resolution is neither recorded nor
#' supplied, or if two standards share a resolution, since a scan could then not
#' be matched to one of them unambiguously.
calibrate_standards = function(standards, opt = OPT, quiet = FALSE) {
  if (!quiet) message("Calibrating")

  calibration = do.call(rbind, lapply(seq_len(nrow(standards)), function(i) {
    path = standards$path[[i]]
    if (!file.exists(path))
      stop("Size standard not found: ", path, call. = FALSE)

    std = read_scan(path)
    dpi = if (is.na(standards$dpi[[i]])) std$dpi else standards$dpi[[i]]
    if (is.na(dpi))
      stop(basename(path), " records no resolution; set `dpi` for it in the ",
           "standards table so scans can be matched to it", call. = FALSE)

    ## The standard is a large, solid object, so the minimum-size floor only has
    ## to exclude dust. Scale it with resolution so it means the same thing at
    ## 75 dpi as at 600.
    min_px = opt$min_object_cm2 / (2.54 / dpi)^2
    seg = segment_leaf(flatten_scanner(std$gray, opt)$gray, min_px = min_px, opt = opt)
    cm2_per_px = standards$area_cm2[[i]] / seg$leaf_px

    if (!quiet)
      message(sprintf(
        "  %-28s %10s px = %.4f cm2 -> %.4g cm2/px (%.1f dpi measured, %s tagged)",
        basename(path), format(seg$leaf_px, big.mark = ","),
        standards$area_cm2[[i]], cm2_per_px, 2.54 / sqrt(cm2_per_px), dpi))
    if (seg$second_frac > opt$max_second_object)
      warning("More than one large object on ", basename(path), call. = FALSE)

    data.frame(standard = basename(path), dpi = dpi, cm2_per_px = cm2_per_px,
               stringsAsFactors = FALSE)
  }))

  if (anyDuplicated(round(calibration$dpi)))
    stop("Two size standards share a resolution; a scan could not be matched ",
         "to one of them unambiguously", call. = FALSE)
  calibration
}

#' Pick the size standard scanned at the same resolution as this scan
#'
#' @param dpi Resolution of the scan, or `NA`.
#' @param calibration The table from [calibrate_standards()].
#' @param opt Options from [leaf_area_options()].
#' @return The matching one-row slice of `calibration`.
#' @section Errors:
#' Stops when no standard is within `opt$dpi_match_tolerance` of the scan, and
#' when a scan records no resolution but more than one standard is available.
#' Both are deliberate: a silent fallback here is how areas end up 78% wrong.
match_standard = function(dpi, calibration, opt = OPT) {
  if (is.na(dpi)) {
    if (nrow(calibration) == 1L) return(calibration[1L, ])
    stop("this scan records no resolution, so it cannot be matched to one of ",
         "the ", nrow(calibration), " size standards", call. = FALSE)
  }
  hit = which(abs(calibration$dpi / dpi - 1) <= opt$dpi_match_tolerance)
  if (!length(hit))
    stop(sprintf("no size standard was scanned at %.0f dpi (have: %s)", dpi,
                 paste(sprintf("%.0f", calibration$dpi), collapse = ", ")),
         call. = FALSE)
  calibration[hit[1L], ]
}


## =============================================================================
## Per-scan measurement and output
## =============================================================================

#' Columns of the output CSV, in order
#'
#' Kept as a constant so that a successful row and an error row have identical
#' shape and `rbind()` cannot silently reorder anything.
#' @return A character vector of column names.
output_columns = function() {
  c("scan_id", "file", "width_px", "height_px", "dpi", "trimmed_px", "standard",
    "cutoff", "paper", "leaf_tone", "contrast", "n_objects", "leaf_px",
    "holes_px", "recovered_px", "frac_of_image", "area_cm2", "area_m2", "area_cm2_lo",
    "area_cm2_hi", "area_sens_pct", "area_cm2_dpi", "flags")
}

#' An all-NA output row, for a scan that could not be measured
#'
#' A failed scan still gets a row, so the output has one line per input file and
#' a missing measurement is visible rather than absent.
#'
#' @param scan_id File name without its extension.
#' @param path Path to the file.
#' @param flags Text placed in the `flags` column, normally the error message.
#' @return A one-row data frame with the columns of [output_columns()].
blank_row = function(scan_id, path, flags) {
  cols = output_columns()
  row = as.data.frame(setNames(rep(list(NA), length(cols)), cols),
                      stringsAsFactors = FALSE)
  row$scan_id = scan_id
  row$file    = basename(path)
  row$flags   = flags
  row
}

#' Measure one scan end to end
#'
#' Read, flatten, segment, calibrate, and evaluate the quality-control flags.
#' The flags never change a measured area; they only tell you how much to trust
#' one.
#'
#' @param path Path to the scan.
#' @param calibration The table from [calibrate_standards()].
#' @param opt Options from [leaf_area_options()].
#' @param mask_dir Directory to write the audit mask into, or `NULL` for none.
#' @return A one-row data frame with the columns of [output_columns()].
measure_scan = function(path, calibration, opt = OPT, mask_dir = NULL) {
  scan_id = tools::file_path_sans_ext(basename(path))

  img  = read_scan(path)
  std  = match_standard(img$dpi, calibration, opt)
  trim = flatten_scanner(img$gray, opt)
  gray = trim$gray
  seg  = segment_leaf(gray, min_px = opt$min_object_cm2 / std$cm2_per_px, opt = opt)

  area_cm2 = seg$leaf_px * std$cm2_per_px
  ## An area computed from the file's own DPI tag, as a check on the standard.
  ## The two should agree to within the scanner's pitch error.
  area_dpi = if (is.na(img$dpi)) NA_real_ else seg$leaf_px * (2.54 / img$dpi)^2
  sens = 100 * (seg$leaf_px_hi - seg$leaf_px_lo) / (2 * seg$leaf_px)

  flags = c(
    if (seg$contrast    < opt$min_contrast)      "low_contrast",
    if (isTRUE(sens     > opt$max_sensitivity))  "soft_margin",
    if (seg$second_frac > opt$max_second_object) "extra_object",
    if (seg$touches_edge)                        "touches_edge",
    if (seg$recovered_px > 0)                    "interior_recovered",
    if (!is.na(area_dpi) &&
        abs(area_dpi / area_cm2 - 1) > opt$dpi_tolerance) "dpi_mismatch"
  )

  if (!is.null(mask_dir))
    write_mask(restore_mask(seg$mask, trim, dim(img$gray)),
               file.path(mask_dir, paste0(scan_id, ".png")))

  data.frame(
    scan_id       = scan_id,
    file          = basename(path),
    width_px      = ncol(img$gray),
    height_px     = nrow(img$gray),
    dpi           = img$dpi,
    trimmed_px    = trim$dropped,
    standard      = std$standard,
    cutoff        = seg$cutoff,
    paper         = round(seg$paper, 1),
    leaf_tone     = round(seg$leaf_tone, 1),
    contrast      = round(seg$contrast, 1),
    n_objects     = seg$n_objects,
    leaf_px       = seg$leaf_px,
    holes_px      = seg$holes_px,
    recovered_px  = seg$recovered_px,
    frac_of_image = round(seg$leaf_px / length(img$gray), 5),
    area_cm2      = round(area_cm2, 4),
    area_m2       = round(area_cm2 / 1e4, 8),
    area_cm2_lo   = round(seg$leaf_px_lo * std$cm2_per_px, 4),
    area_cm2_hi   = round(seg$leaf_px_hi * std$cm2_per_px, 4),
    area_sens_pct = round(sens, 2),
    area_cm2_dpi  = round(area_dpi, 4),
    flags         = paste(flags, collapse = "|"),
    stringsAsFactors = FALSE
  )
}

#' Measure a folder of scans and write the results
#'
#' The entry point. Calibrates against the size standards, measures every image
#' in `scan_dir`, and writes one CSV row per file plus one audit mask per
#' successfully measured scan.
#'
#' A scan that fails is reported as a warning and an all-NA row; the run
#' continues. A failure in calibration, by contrast, stops everything, because
#' every area downstream would be wrong.
#'
#' @param scan_dir Folder of leaf scans.
#' @param standards Data frame of size standards; see [calibrate_standards()].
#' @param output_csv Path for the results CSV. Its directory is created if
#'   needed.
#' @param mask_dir Folder for the audit masks, or `NULL` to skip writing them.
#' @param opt Options from [leaf_area_options()].
#' @param recursive Search subfolders of `scan_dir` as well.
#' @param quiet Suppress progress messages.
#' @return The results data frame, invisibly. Also written to `output_csv`.
#' @examples
#' \dontrun{
#' run_leaf_area(
#'   scan_dir  = "scans/raw_scans_tiff",
#'   standards = data.frame(path = "scans/POST_IT_NOTE.tiff",
#'                          area_cm2 = 58.0644, dpi = NA),
#'   output_csv = "output/leaf_area.csv")
#' }
run_leaf_area = function(scan_dir,
                         standards,
                         output_csv,
                         mask_dir  = NULL,
                         opt       = OPT,
                         recursive = FALSE,
                         quiet     = FALSE) {

  files = list.files(scan_dir, pattern = "\\.(tiff?|jpe?g|png)$",
                     ignore.case = TRUE, full.names = TRUE,
                     recursive = recursive)
  if (!length(files))
    stop("No images found in ", scan_dir, call. = FALSE)

  calibration = calibrate_standards(standards, opt, quiet = quiet)

  if (!is.null(mask_dir))
    dir.create(mask_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(dirname(output_csv), showWarnings = FALSE, recursive = TRUE)

  if (!quiet) message("Measuring ", length(files), " scans in ", scan_dir)
  rows = vector("list", length(files))

  for (i in seq_along(files)) {
    path    = files[[i]]
    scan_id = tools::file_path_sans_ext(basename(path))
    started = Sys.time()

    rows[[i]] = tryCatch(
      measure_scan(path, calibration, opt, mask_dir),
      error = function(e) {
        warning("Failed on ", basename(path), ": ", conditionMessage(e),
                call. = FALSE)
        blank_row(scan_id, path, paste0("error: ", conditionMessage(e)))
      })

    if (!quiet)
      message(sprintf("  [%*d/%d] %-16s %9.3f cm2 +-%4.1f%%  %4.1fs %s",
                      nchar(length(files)), i, length(files), scan_id,
                      rows[[i]]$area_cm2, rows[[i]]$area_sens_pct,
                      as.numeric(Sys.time() - started, units = "secs"),
                      rows[[i]]$flags))
  }

  leaf_area = do.call(rbind, rows)
  leaf_area = leaf_area[order(leaf_area$scan_id), output_columns()]
  write.csv(leaf_area, output_csv, row.names = FALSE)

  if (!quiet)
    message(sprintf("Wrote %s  (%d of %d scans flagged)", output_csv,
                    sum(leaf_area$flags != "", na.rm = TRUE), nrow(leaf_area)))
  invisible(leaf_area)
}


## =============================================================================
## Run
##
## Everything above is a definition; nothing has happened yet. Edit the settings
## here and source this file, or run it with Rscript.
##
## To load the functions WITHOUT processing anything -- for testing, or to call
## them from another script -- set options(leaf_area.functions_only = TRUE)
## before source()ing.
## =============================================================================

if (!isTRUE(getOption("leaf_area.functions_only"))) {

  ## Size standards. One row per physical object of known area that you scanned.
  ## A scan is calibrated with the standard whose resolution matches its own, so
  ## 300 and 400 dpi batches can sit in the same folder.
  ##
  ##   path      file holding the standard
  ##   area_cm2  its true area. 3 in x 3 in = 9 in2 = 58.0644 cm2.
  ##             MEASURE THIS. It sets the absolute scale of every area below.
  ##   dpi       resolution override, or NA to read it from the file
  STANDARDS = data.frame(
    path     = c("scans/POST_IT_NOTE.tiff", "scans/STANDARD_3inx3in.jpeg"),
    area_cm2 = c(58.0644,                    58.0644),
    dpi      = c(NA,                         300),
    stringsAsFactors = FALSE
  )

  leaf_area = run_leaf_area(
    scan_dir   = "scans/raw_scans_tiff",
    standards  = STANDARDS,
    output_csv = "output/leaf_area.csv",
    mask_dir   = "scans/processed_scans_jpeg",   # NULL to skip audit masks
    recursive  = FALSE
  )
}

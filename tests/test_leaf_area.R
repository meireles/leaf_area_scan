################################################################################
## Test suite for leaf_area.R
##
## Runs with plain Rscript, no test framework and no data files: every fixture
## is a synthetic scan generated here, carrying the same paper texture, shading
## gradient and lid shadow measured on the real scans, so the true area of every
## object is known exactly.
##
##   Rscript tests/test_leaf_area.R
##
## Exits non-zero if anything fails, so it drops straight into CI. The `expect`
## helper maps one-to-one onto testthat's `expect_*`, if this later becomes a
## package.
################################################################################

options(leaf_area.functions_only = TRUE)
source(if (file.exists("leaf_area.R")) "leaf_area.R" else "../leaf_area.R")

set.seed(20260728)

## --- tiny test harness -------------------------------------------------------

PASS = 0L; FAIL = 0L

#' Record one assertion
#' @param label What is being checked, shown in the output.
#' @param ok Logical; `TRUE` for a pass.
#' @param detail Extra text shown when the assertion fails.
expect = function(label, ok, detail = "") {
  ok = isTRUE(ok)
  if (ok) PASS <<- PASS + 1L else FAIL <<- FAIL + 1L
  cat(sprintf("  %s %-58s %s\n", if (ok) "PASS" else "FAIL", label, detail))
}

#' Assert a measured value is within a relative tolerance of the truth
#' @param label What is being checked.
#' @param got,want Measured and true values.
#' @param tol_pct Allowed relative difference, in percent.
expect_close = function(label, got, want, tol_pct) {
  err = 100 * (got / want - 1)
  expect(label, is.finite(err) && abs(err) <= tol_pct,
         sprintf("%.0f vs %.0f  (%+.3f%%, tol %.2f%%)", got, want, err, tol_pct))
}

#' Assert an expression raises an error
#' @param label What is being checked.
#' @param expr Expression expected to fail.
expect_error = function(label, expr) {
  got = tryCatch({ force(expr); NULL }, error = function(e) conditionMessage(e))
  expect(label, !is.null(got), if (is.null(got)) "no error raised" else got)
}


## --- synthetic scan generator ------------------------------------------------

#' Approximate a Gaussian blur with repeated 3x3 box passes
#' @param m Numeric matrix.
#' @param passes How many box passes; 2 gives sigma of about 0.9 px.
#' @return A blurred matrix the same size.
blur = function(m, passes = 2L) {
  for (i in seq_len(passes)) {
    m = (cbind(m[, 1L], m[, -ncol(m), drop = FALSE]) + m +
         cbind(m[, -1L, drop = FALSE], m[, ncol(m)])) / 3
    m = (rbind(m[1L, ], m[-nrow(m), , drop = FALSE]) + m +
         rbind(m[-1L, , drop = FALSE], m[nrow(m), ])) / 3
  }
  m
}

#' A blank scan carrying the artifacts measured on the real scanner
#'
#' Paper at gray 245 with correlated texture, a gentle shading gradient across
#' the bed, and a multiplicative lid shadow over the last 18 columns falling to
#' about 43 percent transmission, which is what the real scans show.
#'
#' @param h,w Image size in pixels.
#' @return A list with the `paper` matrix and the `shade` profile that produced
#'   its lid shadow, so objects can be attenuated consistently.
blank_scan = function(h = 1400L, w = 1000L) {
  paper = matrix(245, h, w) + blur(matrix(rnorm(h * w, 0, 6), h, w), 2L)
  paper = paper - 4 * rep(seq_len(w) / w, each = h)          # shading across the bed
  shade = 1 - 0.57 * pmin(pmax((seq_len(w) - (w - 18L)) / 17, 0), 1)^1.8
  list(paper = paper * rep(shade, each = h), shade = shade)
}

#' Stamp an object into a scan
#'
#' The object is attenuated by the same lid shadow as the paper it sits on,
#' which is what makes it physically consistent: a leaf lying in the shadow gets
#' darker along with its background, so the contrast between them is scaled, not
#' offset.
#'
#' @param bg A list from [blank_scan()].
#' @param shape Logical matrix, `TRUE` inside the object.
#' @param tone Gray level of the object where it is fully lit.
#' @return A numeric matrix of 0..255 values.
stamp = function(bg, shape, tone) {
  soft = blur(shape * 1.0, 2L)
  lit  = tone * rep(bg$shade, each = nrow(shape))
  pmin(pmax(round(bg$paper * (1 - soft) + lit * soft), 0), 255)
}

#' Rectangle as a logical matrix
rect_shape = function(h, w, r0, r1, c0, c1) {
  m = matrix(FALSE, h, w); m[r0:r1, c0:c1] = TRUE; m
}


## =============================================================================
cat("\nUnit tests: matrix helpers\n")
## =============================================================================

m = matrix(c(9, 1, 5, 3,
             7, 2, 8, 4), nrow = 2, byrow = TRUE)
expect("downsample_min takes the block minimum",
       identical(as.vector(downsample_min(m, 2L)), c(1, 3)))
expect("downsample_min is a no-op at factor 1", identical(downsample_min(m, 1L), m))
expect("downsample_min discards a partial block",
       identical(dim(downsample_min(matrix(1, 7, 7), 3L)), c(2L, 2L)))

sq = matrix(FALSE, 11, 11); sq[4:8, 4:8] = TRUE
expect("erode(r=1) shrinks a 5x5 square to 3x3", sum(erode(sq, 1L)) == 9)
expect("erode(r=2) shrinks it to a single pixel", sum(erode(sq, 2L)) == 1)
expect("erode(r=3) erases it", sum(erode(sq, 3L)) == 0)
expect("erode(r=0) is a no-op", identical(erode(sq, 0L), sq))
expect("bounds finds the box", identical(bounds(sq), c(4L, 8L, 4L, 8L)))

flat = matrix(200, 100, 100); flat[1:10, ] = 50           # 10% dark
expect("background_level ignores the dark minority",
       abs(background_level(flat) - 200) < 0.01)


## =============================================================================
cat("\nUnit tests: connected components\n")
## =============================================================================

## A one-pixel-wide diagonal: 8-connectivity sees one object, 4-connectivity
## sees twelve. This is the difference between measuring a pine needle and
## measuring nothing at all.
diag_mask = matrix(FALSE, 12, 12)
for (i in 1:12) diag_mask[i, i] = TRUE
expect("8-connectivity keeps a diagonal whole",
       length(label_components(diag_mask, diagonal = TRUE)$size) == 1L)
expect("4-connectivity shatters it",
       length(label_components(diag_mask, diagonal = FALSE)$size) == 12L)

two = matrix(FALSE, 20, 20); two[2:6, 2:6] = TRUE; two[12:18, 12:18] = TRUE
cc = label_components(two)
expect("two objects are found", length(cc$size) == 2L)
expect("their sizes are right", identical(unname(sort(cc$size)), c(25, 49)))
expect("neither touches the border", !any(cc$on_edge))
expect("draw_runs round-trips a component",
       sum(draw_runs(cc, cc$runs_of[[which.max(cc$size)]])) == 49)

ring = matrix(FALSE, 20, 20); ring[5:15, 5:15] = TRUE; ring[8:12, 8:12] = FALSE
expect("fill_holes fills an enclosed hole", sum(fill_holes(ring)) == 121)
expect("fill_holes leaves a hole that drains to the border open",
       { notch = ring; notch[8:12, 13:15] = FALSE   # cut through the wall
         sum(fill_holes(notch)) == sum(notch) })

## The lid shadow: pinned to a border, deep aspect ratio. A needle of the same
## length lying free of the border, or a broad leaf running off it, is not.
shadow = matrix(FALSE, 400, 300); shadow[, 297:300] = TRUE
free   = matrix(FALSE, 400, 300); free[10:390, 200:203] = TRUE
broad  = matrix(FALSE, 400, 300); broad[100:250, 1:120] = TRUE
expect("edge_artifact catches the lid shadow",
       edge_artifact(label_components(shadow), 1L, OPT))
expect("edge_artifact spares a needle away from the border",
       !edge_artifact(label_components(free), 1L, OPT))
expect("edge_artifact spares a broad leaf running off the border",
       !edge_artifact(label_components(broad), 1L, OPT))

both = shadow | broad
cc2 = label_components(both)
expect("largest_specimen picks the leaf over the bigger shadow",
       sum(draw_runs(cc2, cc2$runs_of[[largest_specimen(cc2, 10, OPT)]])) == sum(broad))


## =============================================================================
cat("\nUnit tests: standard matching\n")
## =============================================================================

cal = data.frame(standard = c("a", "b"), dpi = c(300, 400),
                 cm2_per_px = c(7.1e-5, 4.0e-5), stringsAsFactors = FALSE)
expect("a 400 dpi scan matches the 400 dpi standard",
       match_standard(400, cal)$standard == "b")
expect("a 300 dpi scan matches the 300 dpi standard",
       match_standard(300, cal)$standard == "a")
expect("a scan within tolerance still matches",
       match_standard(399.999, cal)$standard == "b")
expect_error("an unmatched resolution is refused", match_standard(600, cal))
expect_error("an untagged scan with two standards is refused",
             match_standard(NA_real_, cal))
expect("an untagged scan with one standard falls back to it",
       match_standard(NA_real_, cal[1L, ])$standard == "a")


## =============================================================================
cat("\nUnit tests: lateral response correction\n")
## =============================================================================

bg = blank_scan()
prof = paper_profile(bg$paper, 2L)
expect("the profile finds paper in the lit part of the bed",
       abs(median(prof[1:900]) - 245) < 6, sprintf("median %.1f", median(prof[1:900])))
expect("the profile follows the shadow down at the edge",
       prof[1000] < 0.6 * median(prof), sprintf("%.0f vs %.0f", prof[1000], median(prof)))

trim = flatten_scanner(bg$paper)
kept = trim$gray[, 1:(ncol(trim$gray) - 2)]
expect("correction flattens the bed to within a few gray levels",
       diff(range(apply(kept, 2L, median))) < 6,
       sprintf("range %.1f", diff(range(apply(kept, 2L, median)))))
expect("only a handful of blinded columns are dropped",
       trim$dropped > 0 && trim$dropped < 12, sprintf("dropped %d", trim$dropped))


## =============================================================================
cat("\nEnd-to-end: known areas on synthetic scans\n")
## =============================================================================

H = 1400L; W = 1000L
dir.create(fixtures <- file.path(tempdir(), "fixtures", "scans"),
           recursive = TRUE, showWarnings = FALSE)

#' Write one fixture scan and return its true pixel area
write_fixture = function(name, shape, tone) {
  invisible(tiff::writeTIFF(stamp(blank_scan(H, W), shape, tone) / 255,
                            file.path(fixtures, paste0(name, ".tif"))))
  sum(shape)
}

## The size standard: a large solid square, 200 x 200 px, declared as 4 cm2.
std_shape = rect_shape(H, W, 100L, 299L, 100L, 299L)
invisible(tiff::writeTIFF(stamp(blank_scan(H, W), std_shape, 200) / 255,
                          file.path(dirname(fixtures), "standard.tif")))
STD_TRUE_PX   = sum(std_shape)
STD_TRUE_CM2  = 4.0

truth = c(
  ## a big low-contrast block: the easy case, and the tightest tolerance
  block   = write_fixture("block",   rect_shape(H, W, 300L, 749L, 200L, 799L), 200),
  ## a thin needle in the middle of the bed
  needle  = write_fixture("needle",  rect_shape(H, W, 200L, 1199L, 500L, 509L), 150),
  ## the same needle deep in the lid shadow, down to 78% transmission, but
  ## still inside the region the correction can recover
  at_edge = write_fixture("at_edge", rect_shape(H, W, 200L, 1199L, W - 17L, W - 8L), 150),
  ## a leaf with an insect hole; areas include holes, so `holes_px` is subtracted
  holed   = write_fixture("holed",
              { s = rect_shape(H, W, 400L, 799L, 300L, 699L)
                s[550:649, 450:549] = FALSE; s }, 140)
)
## a needle running off the very edge of the bed, into the columns the lid has
## blinded. Not recoverable, and not meant to be: the test is that the loss is
## confined to those columns and that the row says so.
blinded = rect_shape(H, W, 200L, 1199L, W - 9L, W)
invisible(tiff::writeTIFF(stamp(blank_scan(H, W), blinded, 150) / 255,
                          file.path(fixtures, "blinded.tif")))
## a scan with nothing on it at all
invisible(tiff::writeTIFF(blank_scan(H, W)$paper / 255,
                          file.path(fixtures, "blank.tif")))
## a leaf plus a second, smaller, darker object elsewhere on the plate
two_obj = rect_shape(H, W, 300L, 699L, 200L, 599L)
invisible(tiff::writeTIFF(stamp(blank_scan(H, W),
                      two_obj | rect_shape(H, W, 1000L, 1199L, 700L, 899L), 160) / 255,
                          file.path(fixtures, "two_objects.tif")))

res = run_leaf_area(
  scan_dir   = fixtures,
  standards  = data.frame(path = file.path(dirname(fixtures), "standard.tif"),
                          area_cm2 = STD_TRUE_CM2, dpi = 400,
                          stringsAsFactors = FALSE),
  output_csv = file.path(tempdir(), "fixtures", "out.csv"),
  mask_dir   = NULL,
  quiet      = TRUE)
row = function(id) res[match(id, res$scan_id), ]

cat("\n")
expect_close("standard recovers its own pixel count",
             STD_TRUE_CM2 / row("block")$area_cm2 * row("block")$leaf_px, STD_TRUE_PX,
             1.0)
expect_close("block: large, low contrast",   row("block")$leaf_px,   truth[["block"]],   0.10)
expect_close("needle: thin, mid-bed",        row("needle")$leaf_px,  truth[["needle"]],  0.50)
expect_close("needle 6 px from the platen edge, in the lid shadow",
             row("at_edge")$leaf_px, truth[["at_edge"]], 0.50)
expect_close("holed leaf, once filled holes are subtracted",
             row("holed")$leaf_px - row("holed")$holes_px, truth[["holed"]], 0.10)
expect_close("two objects: the larger one is measured",
             row("two_objects")$leaf_px, sum(two_obj), 0.10)

cat("\n")
expect("a leaf in the shadow measures the same as one in the open",
       abs(row("at_edge")$leaf_px / row("needle")$leaf_px - 1) < 0.005,
       sprintf("%d vs %d px", row("at_edge")$leaf_px, row("needle")$leaf_px))
expect("holes are reported, not silently swallowed",
       abs(row("holed")$holes_px - 100 * 100) < 400,
       sprintf("holes_px = %d, true hole = 10000", row("holed")$holes_px))
expect("a second object on the plate is flagged",
       grepl("extra_object", row("two_objects")$flags), row("two_objects")$flags)
expect("a blank scan errors without killing the run",
       is.na(row("blank")$leaf_px) && grepl("^error:", row("blank")$flags),
       row("blank")$flags)
expect("every input file produced a row", nrow(res) == 7L, sprintf("%d rows", nrow(res)))
expect("a leaf running into the blinded columns is flagged, not silently short",
       grepl("touches_edge", row("blinded")$flags), row("blinded")$flags)
## Not a soft "most of it survived" check: it should recover exactly the part
## of the leaf that was still lit, and nothing more.
expect_close("and it keeps exactly the part of the leaf that was still lit",
             row("blinded")$leaf_px,
             sum(blinded[, flatten_scanner(blank_scan(H, W)$paper)$cols]), 1.0)
expect("the error bar brackets the measurement",
       with(row("needle"), area_cm2_lo <= area_cm2 && area_cm2 <= area_cm2_hi))
expect("a crisp block is less cutoff-sensitive than a thin needle",
       row("block")$area_sens_pct < row("needle")$area_sens_pct,
       sprintf("%.2f%% vs %.2f%%", row("block")$area_sens_pct, row("needle")$area_sens_pct))
expect("options are honored, not ignored",
       { o = leaf_area_options(min_object_cm2 = 1e6)
         inherits(try(segment_leaf(flatten_scanner(
           stamp(blank_scan(H, W), rect_shape(H, W, 300L, 320L, 300L, 320L), 150),
           o)$gray, min_px = 1e9, opt = o), silent = TRUE), "try-error") })


## =============================================================================
cat(sprintf("\n%d passed, %d failed\n", PASS, FAIL))
if (FAIL > 0L) quit(status = 1L)

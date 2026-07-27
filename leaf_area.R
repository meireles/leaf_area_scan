################################################################################
## Leaf area from scans
##
## Measures the one-sided area of a single leaf per scan and writes a tidy CSV
## plus a binary mask per scan so every number can be audited by eye.
##
## How does it worl?
##   1. Read the scan and reduce it to 8-bit greyscale
##   2. Pick the black/white cut-off automatically with Otsu's method.
##   3. Keep only the largest connected dark object and fill any holes in it.
##   4. Convert pixels to cm2 using a 3 in x 3 in reference square scanned on
##      the same machine, and cross-check that against the scan's DPI metadata.
##
## Requires the tiff, jpeg and png packages (only the ones matching your input
## formats are actually needed).
##
##   install.packages(c("tiff", "jpeg", "png"))
##
################################################################################


## =============================================================================
## Configuration
## =============================================================================

SCAN_DIR      = "scans/raw_scans_tiff"        # folder of leaf scans
STANDARD_PATH = "scans/STANDARD_3inx3in.jpeg" # size reference, same scanner
MASK_DIR      = "scans/processed_scans_jpeg"  # masks go here, written as PNG:
                                              # lossless, so re-counting a mask
                                              # reproduces the area exactly
OUTPUT_CSV    = "output/leaf_area.csv"

STANDARD_AREA_CM2 = 58.0644  # 3 in x 3 in = 9 in2. CHANGE THIS DEPENDING ON WHAT THE STANDARD IS
RECURSIVE         = FALSE    # also search subfolders of SCAN_DIR?
WRITE_MASKS       = TRUE     # write a mask per scan so you can audit it

## Quality-control thresholds. These only populate the `flags` column;
## they never change a measured area.
MIN_SEPARABILITY  = 0.90  # Otsu separability below this = poor leaf/background contrast
MAX_SECOND_OBJECT = 0.01  # 2nd largest object bigger than this x the leaf = something else on the scan
DPI_TOLERANCE     = 0.02  # relative disagreement allowed between standard- and DPI-based area


## =============================================================================
## Input
## Returns an array indexed [row, column, channel] with values in 0..1.
## Convention for tiff, jpeg and png packages.
## =============================================================================

read_grey = function(path) {
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
  if (is.list(img)) {
    warning(basename(path), " has several pages; using the first", call. = FALSE)
    img = img[[1L]]
  }

  list(grey = round(luma(img) * 255), dpi = read_dpi(path, img))
}

read_with = function(pkg, reader, path, ...) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("package '", pkg, "' is needed to read ", basename(path), call. = FALSE)
  reader(path, ...)
}

## Collapse an image to a single luminance channel, dropping alpha.
luma = function(img) {
  if (length(dim(img)) == 2L) return(img)
  if (dim(img)[3L] >= 3L)
    0.2126 * img[, , 1L] + 0.7152 * img[, , 2L] + 0.0722 * img[, , 3L]
  else
    img[, , 1L]
}

## Horizontal resolution in dots per inch, or NA if the file does not record it.
## Used only as an independent check on the size standard.
read_dpi = function(path, img) {
  dpi = switch(
    tolower(tools::file_ext(path)),
    tif = , tiff = {
      d = suppressWarnings(as.numeric(attr(img, "x.resolution")))
      if (identical(tolower(attr(img, "resolution.unit")), "centimeter")) d * 2.54 else d
    },
    jpg = , jpeg = jfif_dpi(path),
    png = png_dpi(path),
    NA_real_
  )
  if (!length(dpi) || !is.finite(dpi) || dpi <= 1) NA_real_ else as.numeric(dpi)
}

## Find the byte offset of a marker near the start of a file, or NA.
find_marker = function(path, marker, n = 4096L) {
  bytes = readBin(path, "raw", n = n)
  for (i in which(bytes == marker[1L]))
    if (i + length(marker) - 1L <= length(bytes) &&
        identical(bytes[i + seq_along(marker) - 1L], marker))
      return(list(bytes = bytes, at = i))
  list(bytes = bytes, at = NA_integer_)
}

## units 1 = dots per inch, 2 = dots per cm.
jfif_dpi = function(path) {
  f = find_marker(path, as.raw(c(0x4A, 0x46, 0x49, 0x46, 0x00)))
  if (is.na(f$at) || f$at + 9L > length(f$bytes)) return(NA_real_)
  units = as.integer(f$bytes[f$at + 7L])
  x     = as.integer(f$bytes[f$at + 8L]) * 256L + as.integer(f$bytes[f$at + 9L])
  if (x == 0L) NA_real_ else if (units == 1L) x else if (units == 2L) x * 2.54 else NA_real_
}

## PNG pHYs chunk: "pHYs", pixels per unit X(4), Y(4), unit(1; 1 = metre).
png_dpi = function(path) {
  f = find_marker(path, charToRaw("pHYs"))
  if (is.na(f$at) || f$at + 12L > length(f$bytes)) return(NA_real_)
  per_unit = sum(as.integer(f$bytes[f$at + 4:7]) * c(256^3, 256^2, 256, 1))
  if (as.integer(f$bytes[f$at + 12L]) == 1L && per_unit > 0) per_unit * 0.0254 else NA_real_
}


## =============================================================================
## Thresholding
## =============================================================================

## Otsu's method: pick the grey level that best splits the image in two.
## Returns the cut-off and a separability - scale 0..1 (the faction of total
## var the split explains)
otsu = function(grey) {
  counts = as.double(tabulate(as.integer(grey) + 1L, nbins = 256L))
  levels = as.double(0:255)
  total  = sum(counts)

  below = cumsum(counts)                 # pixels at or below each level
  above = total - below
  sums  = cumsum(counts * levels)        # their summed grey value
  grand = sums[256L] / total

  between = rep(-Inf, 256L)
  ok = below > 0 & above > 0
  between[ok] = (grand * below[ok] - sums[ok])^2 / (below[ok] * above[ok])

  ## `between` is already the between-class variance, so the separability is
  ## just its peak as a fraction of the total variance.
  variance = sum(counts * (levels - grand)^2) / total
  list(threshold    = levels[which.max(between)],
       separability = if (variance > 0) max(between) / variance else 0)
}


## =============================================================================
## Connected components (run-length encoded, union-find)
## Label objects in a binary mask.
## =============================================================================

label_components = function(mask) {
  height = nrow(mask)
  width  = ncol(mask)

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

  ## Union-find over runs, with path compression.
  parent = seq_len(n)
  find = function(i) {
    while (parent[i] != i) {
      parent[i] <<- parent[parent[i]]
      i = parent[i]
    }
    i
  }

  ## Join each run to any run in the row above that it overlaps horizontally.
  by_row  = split(seq_len(n), row)
  row_no  = as.integer(names(by_row))
  for (k in seq_along(by_row)[-1L]) {
    if (row_no[k] != row_no[k - 1L] + 1L) next
    above = by_row[[k - 1L]]
    for (j in by_row[[k]]) {
      for (i in above[start[above] <= end[j] & start[j] <= end[above]]) {
        ri = find(i); rj = find(j)
        if (ri != rj) parent[ri] = rj
      }
    }
  }

  runs_of  = split(seq_len(n), vapply(seq_len(n), find, integer(1L)))
  on_edge  = row == 1L | row == height | start == 1L | end == width
  lengths_ = end - start + 1L

  list(row = row, start = start, end = end, dim = c(height, width),
       runs_of = runs_of,
       size    = vapply(runs_of, function(k) sum(lengths_[k]), numeric(1L)),
       on_edge = vapply(runs_of, function(k) any(on_edge[k]), logical(1L)))
}

## Paint a set of runs into a fresh bool matrix.
draw_runs = function(cc, which_runs) {
  out = matrix(FALSE, cc$dim[1L], cc$dim[2L])
  for (k in which_runs) out[cc$row[k], cc$start[k]:cc$end[k]] = TRUE
  out
}

## Any background that does not reach the image border is a hole.
fill_holes = function(mask) {
  bg = label_components(!mask)
  if (is.null(bg)) return(mask)
  holes = unlist(bg$runs_of[!bg$on_edge], use.names = FALSE)
  if (length(holes)) mask | draw_runs(bg, holes) else mask
}


## =============================================================================
## Segmentation: one greyscale scan in, one leaf mask out
## =============================================================================

segment_leaf = function(grey) {
  cut  = otsu(grey)
  mask = grey <= cut$threshold            # the leaf is darker than the background

  cc = label_components(mask)
  if (is.null(cc)) stop("no dark object found", call. = FALSE)

  ranked = order(cc$size, decreasing = TRUE)
  leaf   = ranked[1L]
  filled = fill_holes(draw_runs(cc, cc$runs_of[[leaf]]))

  list(mask         = filled,
       threshold    = cut$threshold,
       separability = cut$separability,
       n_objects    = length(cc$size),
       leaf_px      = sum(filled),
       holes_px     = sum(filled) - cc$size[[leaf]],
       discarded_px = sum(cc$size[-leaf]),
       second_frac  = if (length(ranked) > 1L) cc$size[[ranked[2L]]] / cc$size[[leaf]] else 0,
       touches_edge = cc$on_edge[[leaf]])
}

## Write the mask leaf-black on white, matching the scans it came from.
write_mask = function(mask, path) {
  png::writePNG(ifelse(mask, 0, 1), path)
}


## =============================================================================
## Calibrate against the size standard
## =============================================================================

if (!file.exists(STANDARD_PATH))
  stop("Size standard not found: ", STANDARD_PATH, call. = FALSE)

message("Size standard: ", STANDARD_PATH)
standard = segment_leaf(read_grey(STANDARD_PATH)$grey)
CM2_PER_PX = STANDARD_AREA_CM2 / standard$leaf_px

message(sprintf("  %s px = %s cm2  ->  %.4g cm2/px, i.e. %.1f dpi",
                format(standard$leaf_px, big.mark = ","), STANDARD_AREA_CM2,
                CM2_PER_PX, 2.54 / sqrt(CM2_PER_PX)))
if (standard$second_frac > MAX_SECOND_OBJECT)
  warning("The size standard has more than one large object on it; check ",
          STANDARD_PATH, call. = FALSE)


## =============================================================================
## Measure every scan
## =============================================================================

files = list.files(SCAN_DIR, pattern = "\\.(tiff?|jpe?g|png)$",
                    ignore.case = TRUE, full.names = TRUE, recursive = RECURSIVE)
if (!length(files))
  stop("No images found in ", SCAN_DIR, call. = FALSE)
if (WRITE_MASKS)
  dir.create(MASK_DIR, showWarnings = FALSE, recursive = TRUE)

blank_row = function(scan_id, path, flags) data.frame(
  scan_id = scan_id, file = basename(path),
  width_px = NA_integer_, height_px = NA_integer_, dpi = NA_real_,
  threshold = NA_real_, separability = NA_real_, n_objects = NA_integer_,
  leaf_px = NA_real_, holes_px = NA_real_, discarded_px = NA_real_,
  frac_of_image = NA_real_, area_cm2 = NA_real_, area_m2 = NA_real_,
  area_cm2_dpi = NA_real_, flags = flags, stringsAsFactors = FALSE)

message("Measuring ", length(files), " scans in ", SCAN_DIR)
rows = vector("list", length(files))

for (i in seq_along(files)) {
  path    = files[[i]]
  scan_id = tools::file_path_sans_ext(basename(path))
  started = Sys.time()

  rows[[i]] = tryCatch({
    img = read_grey(path)
    seg = segment_leaf(img$grey)

    area_cm2 = seg$leaf_px * CM2_PER_PX
    area_dpi = if (is.na(img$dpi)) NA_real_ else seg$leaf_px * (2.54 / img$dpi)^2

    flags = c(
      if (seg$separability < MIN_SEPARABILITY)  "low_contrast",
      if (seg$second_frac  > MAX_SECOND_OBJECT) "extra_object",
      if (seg$touches_edge)                     "touches_edge",
      if (!is.na(area_dpi) &&
          abs(area_dpi / area_cm2 - 1) > DPI_TOLERANCE) "dpi_mismatch"
    )

    if (WRITE_MASKS)
      write_mask(seg$mask, file.path(MASK_DIR, paste0(scan_id, ".png")))

    data.frame(
      scan_id       = scan_id,
      file          = basename(path),
      width_px      = ncol(img$grey),
      height_px     = nrow(img$grey),
      dpi           = img$dpi,
      threshold     = seg$threshold,
      separability  = round(seg$separability, 4),
      n_objects     = seg$n_objects,
      leaf_px       = seg$leaf_px,
      holes_px      = seg$holes_px,
      discarded_px  = seg$discarded_px,
      frac_of_image = round(seg$leaf_px / length(img$grey), 4),
      area_cm2      = round(area_cm2, 3),
      area_m2       = round(area_cm2 / 1e4, 7),
      area_cm2_dpi  = round(area_dpi, 3),
      flags         = paste(flags, collapse = "|"),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    warning("Failed on ", basename(path), ": ", conditionMessage(e), call. = FALSE)
    blank_row(scan_id, path, paste0("error: ", conditionMessage(e)))
  })

  message(sprintf("  [%*d/%d] %-14s %9.2f cm2  %4.1fs %s",
                  nchar(length(files)), i, length(files), scan_id,
                  rows[[i]]$area_cm2,
                  as.numeric(Sys.time() - started, units = "secs"),
                  rows[[i]]$flags))
}

## =============================================================================
## Write out
## =============================================================================

leaf_area = do.call(rbind, rows)
leaf_area = leaf_area[order(leaf_area$scan_id), ]
write.csv(leaf_area, OUTPUT_CSV, row.names = FALSE)

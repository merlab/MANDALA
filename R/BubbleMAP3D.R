#' Interactive BubbleMAP viewer
#'
#' Pools cells into bubbles on a 2-D embedding exactly as \code{BubbleMAP()}
#' does, then lifts each bubble to the mean of a chosen gene or metadata
#' metric and opens an interactive viewer. Elevation gain, the contour shell and the
#' highlight mode are adjustable in the viewer.
#'
#' Because a bubble's height is a mean over the cells it pools, a marker that
#' is far too sparse for a per-cell Z axis still produces usable relief: the
#' zero-inflation that flattens single cells is averaged away at bubble level.
#'
#' @param object A Seurat object, or a data.frame with columns Dim1/Dim2 (or
#'   UMAP1/UMAP2), the label column, and one column per gene.
#' @param gene Character vector of genes. All are computed and selectable in
#'   the viewer; the first is shown on open unless \code{start_gene} is given.
#' @param metric Character vector of numeric \code{meta.data} columns to offer
#'   alongside the genes — percent.mt, a module score, an ssGSEA score, a
#'   distance. A bubble's height is the mean over its cells either way, so a
#'   metric behaves exactly as a gene does. At least one of \code{gene} or
#'   \code{metric} is required.
#' @param reduction Name of the reduction supplying the floor coordinates.
#' @param label Column in \code{meta.data} holding the cell type.
#' @param assay,layer Where to read expression from. \code{layer} is mapped to
#'   \code{slot} automatically on Seurat v3/v4. Use a normalized layer: raw
#'   counts make heights depend on sequencing depth.
#' @param radius,min_group,rare_relax Bubble grouping, as in \code{BubbleMAP()}.
#' @param min_class Cell types with fewer cells are collapsed into "Other".
#' @param target_cells Optional cap-based stratified downsample before grouping.
#' @param interactive_cells Number of cells shipped to the viewer so that
#'   \code{radius} and \code{min_group} can be re-applied live. Those two
#'   decide which cells share a bubble, so they cannot act on finished bubbles
#'   — the grouping has to be redone, and that needs the cells. A stratified
#'   subsample of this size is embedded; expression is stored sparsely, so the
#'   cost is modest (~1.5 MB at the default). Set NULL to omit it, which hides
#'   the two sliders and keeps the file small. The opening view always uses
#'   every cell; only the sliders work from the subsample.
#' @param palette Named vector of colours. The default reproduces
#'   \code{BubbleMAP()}'s automatic ramp, so the two figures agree as long as
#'   both calls see the same cells (same \code{min_class}, same
#'   \code{target_cells}). To guarantee a match regardless, capture the palette
#'   from the 2-D call and pass it in:
#'   \code{bu <- BubbleMAP(...); BubbleMAP3D(..., palette = bu$palette)}.
#' @param short Named vector of shorter display names for crowded labels.
#' @param point_mode Open with every cell drawn as its own mark instead of
#'   pooled into bubbles. Used by \code{\link{SpatialMAP3D}}; the viewer's
#'   "Original points" button toggles it either way.
#' @param axes Axis captions. Defaults to the reduction key.
#' @param elevation,shell,highlight_all,highlight_lifted,flat Initial state of the viewer
#'   controls. \code{flat = TRUE} opens overhead with height switched off,
#'   matching the \code{\link{BubbleMAP}} layout; the button toggles between
#'   them. Use \code{\link{BubbleMAP}} for a scripted, vector print figure.
#' @param start_gene Which marker to show on open — a gene or a metric.
#'   Defaults to the first one supplied.
#' @param file Where to write the HTML. Defaults to a temporary file.
#' @param launch Open the viewer. Set FALSE to only write the file.
#'
#' @param bw_frac,gridsize,levels_norm Bandwidth fraction, grid size and
#'   normalized levels for the floor density contours.
#' @param title,subtitle Text for the viewer masthead.
#' @param template Path to the viewer template. NULL uses the packaged copy.
#' @param seed Random seed for the greedy grouping pass.
#'
#' @return The path to the written HTML file, invisibly.
#' @export
BubbleMAP3D <- function(
    object,
    gene          = NULL,
    metric        = NULL,
    reduction     = "umap",
    label         = "cell_type",
    assay         = NULL,
    layer         = "data",

    # bubbles — same rules as the 2-D function
    radius        = 0.42,
    min_group     = 8,
    rare_relax    = 3,
    min_class     = 150,
    target_cells  = NULL,
    interactive_cells = 60000,

    # floor contours
    bw_frac       = 0.045,
    gridsize      = 180,
    levels_norm   = c(0.10, 0.25, 0.42, 0.60, 0.80),

    # presentation
    palette       = NULL,
    short         = NULL,
    title         = NULL,
    subtitle      = NULL,

    # initial state of the viewer controls
    elevation     = 1,
    shell         = TRUE,
    highlight_all = FALSE,
    highlight_lifted = TRUE,
    flat          = FALSE,
    point_mode    = FALSE,
    start_gene    = NULL,
    axes          = NULL,

    file          = NULL,
    launch        = TRUE,
    template      = NULL,
    seed          = 42
) {
  if (!requireNamespace("jsonlite", quietly = TRUE))
    stop("BubbleMAP3D() requires the jsonlite package")
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("BubbleMAP3D() requires the MASS package")
  if (is.null(gene) && is.null(metric))
    stop("Supply at least one of gene or metric", call. = FALSE)
  if (!is.null(gene))   stopifnot(is.character(gene))
  if (!is.null(metric)) stopifnot(is.character(metric))
  set.seed(seed)

  # ---- 1. coordinates, labels, expression --------------------------------
  if (inherits(object, "Seurat")) {
    if (!reduction %in% names(object@reductions))
      stop("Reduction '", reduction, "' not found in the Seurat object")
    co  <- as.data.frame(
      object@reductions[[reduction]]@cell.embeddings[, 1:2, drop = FALSE])
    if (!label %in% colnames(object@meta.data))
      stop("Label column '", label, "' not found in meta.data")
    lab <- as.character(object[[label]][, 1])

    if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
    is_v5  <- utils::packageVersion("Seurat") >= "5.0.0"
    gd <- if (is_v5)
      SeuratObject::GetAssayData(object, assay = assay, layer = layer)
    else
      SeuratObject::GetAssayData(object, assay = assay, slot = layer)

    # An empty layer yields no rownames, which otherwise surfaces as "gene not
    # found" for every gene and sends you looking in the wrong place
    if (!nrow(gd) || !ncol(gd)) {
      have <- tryCatch(SeuratObject::Layers(object[[assay]]),
                       error = function(e) character(0))
      stop("Layer '", layer, "' of assay '", assay, "' is empty.",
           if (length(have))
             paste0(" Available: ", paste(have, collapse = ", "),
                    ". Pass layer = \"", have[1], "\".")
           else "")
    }
    miss <- setdiff(gene, rownames(gd))
    if (length(miss))
      stop("Not found in assay '", assay, "': ", paste(miss, collapse = ", "))

    if (length(gene)) {
      expr <- gd[gene, , drop = FALSE]
      expr <- as.matrix(if (inherits(expr, "Matrix")) Matrix::t(expr) else t(expr))
    } else {
      expr <- matrix(numeric(0), nrow = ncol(object), ncol = 0)
    }
    if (length(metric)) {
      mmiss <- setdiff(metric, colnames(object@meta.data))
      if (length(mmiss))
        stop("Not found in meta.data: ", paste(mmiss, collapse = ", "))
      md <- object@meta.data[, metric, drop = FALSE]
      bad <- names(md)[!vapply(md, is.numeric, logical(1))]
      if (length(bad))
        stop("metric must be numeric; these are not: ", paste(bad, collapse = ", "))
      expr <- cbind(expr, as.matrix(md))
    }

    # Raw counts make a bubble's height track sequencing depth as much as
    # expression, so say so rather than silently plotting it
    chk <- expr[seq_len(min(nrow(expr), 2000)), , drop = FALSE]
    if (all(chk == floor(chk)) && max(chk) > 30)
      warning("Layer '", layer, "' looks like raw counts. Bubble heights will ",
              "partly reflect sequencing depth; a normalised layer is safer.",
              call. = FALSE)

    ax <- toupper(sub("_$", "", object@reductions[[reduction]]@key))
    ax <- paste(ax, 1:2)

  } else {
    xy_cols <- if (all(c("Dim1", "Dim2") %in% colnames(object)))
      c("Dim1", "Dim2") else c("UMAP1", "UMAP2")
    stopifnot(all(xy_cols %in% colnames(object)),
              label %in% colnames(object),
              all(c(gene, metric) %in% colnames(object)))
    co   <- object[, xy_cols]
    lab  <- as.character(object[[label]])
    expr <- as.matrix(object[, c(gene, metric), drop = FALSE])
    ax <- c("Dim 1", "Dim 2")
  }
  markers <- c(gene, metric)
  kinds   <- c(rep("gene", length(gene)), rep("metric", length(metric)))
  colnames(expr) <- markers
  names(co) <- c("Dim1", "Dim2")

  keep_row <- stats::complete.cases(co) & !is.na(lab)
  co <- co[keep_row, ]; lab <- lab[keep_row]; expr <- expr[keep_row, , drop = FALSE]

  # ---- 2. optional cap-based stratified downsample ------------------------
  if (!is.null(target_cells) && nrow(co) > target_cells) {
    cnt <- table(lab); lo <- 1L; hi <- max(cnt)
    while (lo < hi) {
      mid <- (lo + hi + 1L) %/% 2L
      if (sum(pmin(cnt, mid)) <= target_cells) lo <- mid else hi <- mid - 1L
    }
    idx <- unlist(lapply(split(seq_along(lab), lab), function(ix)
      if (length(ix) <= lo) ix else sample(ix, lo)))
    co <- co[idx, ]; lab <- lab[idx]; expr <- expr[idx, , drop = FALSE]
    message("Per-class cap: ", lo, " -> ", nrow(co), " cells")
  }

  # ---- 3. collapse rare classes ------------------------------------------
  cnt   <- table(lab)
  keep  <- setdiff(names(cnt)[cnt >= min_class], "Unknown")
  klass <- ifelse(lab %in% keep, lab, "Other")
  lv    <- names(sort(table(klass), decreasing = TRUE))

  # ---- 4. palette ---------------------------------------------------------
  # BubbleMAP() ramps over unique(df$class) in ROW order — order of first
  # appearance in the object — then moves "Other" to the end. Sorting by
  # frequency instead shifts every class one slot along the ramp and the two
  # figures no longer agree, so the ordering is replicated exactly here.
  # This only holds when both calls see the same cells: keep min_class and
  # target_cells identical between them, or pass the palette through
  # explicitly (see below).
  if (is.null(palette)) {
    base_cols <- if (requireNamespace("ggsci", quietly = TRUE))
      ggsci::pal_d3("category20")(20) else scales::hue_pal()(20)
    lv_pal <- unique(klass)
    lv_pal <- c(setdiff(lv_pal, "Other"), intersect("Other", lv_pal))
    palette <- stats::setNames(
      grDevices::colorRampPalette(base_cols)(length(lv_pal)), lv_pal)
    palette["Other"] <- "#D8DCE0"
  }
  miss <- setdiff(lv, names(palette))
  if (length(miss))
    palette <- c(palette, stats::setNames(rep("#B0B8BF", length(miss)), miss))

  # ---- 5. bubbles ---------------------------------------------------------
  # Greedy single pass per class: take an unused seed, absorb every unused
  # cell within `radius`, emit one bubble, repeat. Neighbours come from a
  # uniform spatial hash whose cell size equals `radius`, so the 3x3 block
  # around a seed provably contains every point within range — the same
  # answer a full radius search gives, without materialising the neighbour
  # graph, which is what runs a large object out of memory.
  group_thresh <- function(n) {
    if (is.null(rare_relax) || n >= rare_relax * min_group) return(min_group)
    max(3L, min(as.integer(min_group), as.integer(ceiling(n / rare_relax))))
  }

  bubbles_of <- function(xy, ex) {
    n <- nrow(xy); thr <- group_thresh(n)
    if (n < thr) return(NULL)

    gx <- as.integer(floor(xy[, 1] / radius)); gx <- gx - min(gx)
    gy <- as.integer(floor(xy[, 2] / radius)); gy <- gy - min(gy)
    ny <- as.numeric(max(gy)) + 3                  # numeric: avoids overflow
    key <- as.numeric(gx) * ny + as.numeric(gy)

    ord    <- order(key); skey <- key[ord]
    starts <- c(1L, which(diff(skey) != 0) + 1L)
    ends   <- c(starts[-1] - 1L, n)
    ukey   <- skey[starts]

    cell_members <- function(k) {
      p <- findInterval(k, ukey)
      if (p < 1L || p > length(ukey) || ukey[p] != k) return(integer(0))
      ord[starts[p]:ends[p]]
    }
    r2 <- radius^2
    nb_of <- function(i) {
      cand <- integer(0); kx <- gx[i]; ky <- gy[i]
      for (dx in -1:1) for (dy in -1:1)
        cand <- c(cand, cell_members((kx + dx) * ny + (ky + dy)))
      if (!length(cand)) return(integer(0))
      d2 <- (xy[cand, 1] - xy[i, 1])^2 + (xy[cand, 2] - xy[i, 2])^2
      cand[d2 <= r2]
    }

    used <- logical(n); cap <- n %/% thr + 1L
    b_x <- numeric(cap); b_y <- numeric(cap); b_n <- integer(cap)
    b_e <- matrix(0, nrow = cap, ncol = ncol(ex))
    m <- 0L
    for (i in sample.int(n)) {
      if (used[i]) next
      g <- unique(c(i, nb_of(i))); g <- g[!used[g]]
      if (length(g) < thr) next
      m <- m + 1L
      b_x[m] <- mean(xy[g, 1]); b_y[m] <- mean(xy[g, 2]); b_n[m] <- length(g)
      b_e[m, ] <- colMeans(ex[g, , drop = FALSE])
      used[g] <- TRUE
    }
    if (m == 0L) return(NULL)
    list(x = b_x[seq_len(m)], y = b_y[seq_len(m)],
         n = b_n[seq_len(m)], e = b_e[seq_len(m), , drop = FALSE])
  }

  bub <- list(); k <- 1L

  if (isTRUE(point_mode)) {
    # No aggregation: each cell becomes its own mark. Spot-based assays are
    # already at their measurement unit, and on a regular grid the greedy pass
    # finds nothing to pool anyway unless radius exceeds the spot spacing.
    if (nrow(co) > 200000)
      warning("point_mode with ", nrow(co), " cells will produce a very large ",
              "file and a slow viewer; consider target_cells.", call. = FALSE)
    for (i in seq_len(nrow(co)))
      bub[[i]] <- list(cl = klass[i],
                       x = round(co$Dim1[i], 3), y = round(co$Dim2[i], 3),
                       n = 1L, e = round(as.numeric(expr[i, ]), 4))
    message("Point mode: ", nrow(co), " marks (no pooling)")
  } else {
    for (cl in lv) {
      sel <- klass == cl
      b <- bubbles_of(as.matrix(co[sel, ]), expr[sel, , drop = FALSE])
      if (is.null(b)) {
        message("  ", cl, ": ", sum(sel), " cells -> no bubble (threshold ",
                group_thresh(sum(sel)), " not met within radius ", radius, ")")
        next
      }
      for (j in seq_along(b$n)) {
        bub[[k]] <- list(cl = cl,
                         x = round(b$x[j], 3), y = round(b$y[j], 3),
                         n = b$n[j],
                         e = round(as.numeric(b$e[j, ]), 4))
        k <- k + 1L
      }
      message("  ", cl, ": ", sum(sel), " cells -> ", length(b$n),
              " bubbles (", round(100 * sum(b$n) / sum(sel), 1), "% represented)")
    }
    message("TOTAL ", nrow(co), " cells -> ", length(bub), " bubbles (",
            round(nrow(co) / max(1, length(bub))), "x reduction)")
  }

  if (!length(bub))
    stop("No bubbles formed. Every class fell below min_group within radius = ",
         radius, ". Raise radius above the spacing between points, lower ",
         "min_group, or use point_mode = TRUE if the data is already ",
         "aggregated (Visium spots, segmented cells).", call. = FALSE)

  # ---- 5b. cells for live regrouping in the viewer ------------------------
  cells <- NULL
  if (!is.null(interactive_cells) && interactive_cells > 0) {
    ic <- min(interactive_cells, nrow(co))
    cnt2 <- table(klass); lo <- 1L; hi <- max(cnt2)
    while (lo < hi) {
      mid <- (lo + hi + 1L) %/% 2L
      if (sum(pmin(cnt2, mid)) <= ic) lo <- mid else hi <- mid - 1L
    }
    sel <- sort(unlist(lapply(split(seq_along(klass), klass), function(ix)
      if (length(ix) <= lo) ix else sample(ix, lo))))

    # Expression is mostly zero for most markers, so ship only the non-zeros
    # as interleaved (index, value) pairs. On this data that is ~70k numbers
    # instead of 480k, and the viewer expands them once into dense arrays.
    e_sparse <- lapply(seq_along(markers), function(g) {
      v  <- expr[sel, g]
      nz <- which(v != 0)
      if (!length(nz)) return(I(numeric(0)))
      I(as.numeric(rbind(nz - 1L, round(v[nz], 2))))
    })
    cells <- list(
      classes = I(lv),
      x = I(round(co$Dim1[sel], 2)),
      y = I(round(co$Dim2[sel], 2)),
      c = I(match(klass[sel], lv) - 1L),
      e = e_sparse)
    message("Viewer regrouping: ", length(sel), " cells shipped, ",
            sum(vapply(e_sparse, length, integer(1))) / 2,
            " non-zero expression entries")
  }

  # ---- 6. floor contours and label anchors --------------------------------
  pad <- 0.6
  xr <- range(co$Dim1) + c(-pad, pad)
  yr <- range(co$Dim2) + c(-pad, pad)
  bw <- c(diff(xr), diff(yr)) * bw_frac

  contours <- list(); peaks <- list()
  for (cl in lv) {
    sel <- which(klass == cl)
    if (length(sel) < 200) next
    if (length(sel) > 25000) sel <- sample(sel, 25000)
    kd <- MASS::kde2d(co$Dim1[sel], co$Dim2[sel], h = bw,
                      n = gridsize, lims = c(xr, yr))
    kd$z <- kd$z / max(kd$z)
    ij <- which(kd$z == max(kd$z), arr.ind = TRUE)[1, ]
    peaks[[cl]] <- c(round(kd$x[ij[1]], 2), round(kd$y[ij[2]], 2))

    polys <- list(); p <- 1L
    for (li in seq_along(levels_norm)) {
      cls <- grDevices::contourLines(kd$x, kd$y, kd$z, levels = levels_norm[li])
      for (q in cls) {
        if (length(q$x) < 8) next
        step <- max(1L, length(q$x) %/% 90)
        ix <- seq(1, length(q$x), by = step)
        polys[[p]] <- list(
          l = li - 1L,
          p = unname(round(cbind(q$x[ix], q$y[ix]), 2)))
        p <- p + 1L
      }
    }
    contours[[cl]] <- polys
  }

  # ---- 7. assemble and write ---------------------------------------------
  if (is.null(template)) {
    template <- system.file("extdata", "bubblemap3d_template.html",
                            package = "MANDALA")
    if (!nzchar(template) || !file.exists(template))
      stop("Viewer template not found in the installed MANDALA package. ",
           "Pass template = <path to bubblemap3d_template.html>.",
           call. = FALSE)
  } else if (!file.exists(template)) {
    stop("Viewer template not found at:\n  ", template,
         "\nCheck the spelling and use forward slashes. ",
         "file.exists() on that path returns FALSE.", call. = FALSE)
  }

  payload <- list(
    genes     = as.list(markers),
    kinds     = as.list(kinds),
    n_cells   = nrow(co),
    radius    = radius,
    min_group = min_group,
    cells     = cells,
    xr        = round(xr, 2),
    yr        = round(yr, 2),
    axes      = if (is.null(axes)) ax else axes,
    levels    = levels_norm,
    class_n   = as.list(table(klass)[lv]),
    peaks     = peaks,
    contours  = contours,
    bubbles   = bub,
    palette   = as.list(palette[lv]),
    short     = if (is.null(short)) stats::setNames(list(), character(0))
                else as.list(short),
    title     = if (is.null(title)) "BubbleMAP3D" else title,
    subtitle  = if (is.null(subtitle))
      paste0(format(nrow(co), big.mark = ","),
             " cells pooled into bubbles, lifted by mean expression.")
      else subtitle,
    elevation     = elevation,
    shell         = shell,
    highlight_all = highlight_all,
    highlight_lifted = highlight_lifted,
    flat          = flat,
    point_mode    = point_mode,
    start_gene    = if (is.null(start_gene)) markers[1] else start_gene
  )

  json <- jsonlite::toJSON(payload, auto_unbox = TRUE, digits = 6,
                           null = "null", na = "null")
  html <- paste(readLines(template, warn = FALSE), collapse = "\n")
  html <- sub("__DATA__", json, html, fixed = TRUE)

  if (is.null(file)) file <- tempfile("BubbleMAP3D_", fileext = ".html")
  writeLines(html, file, useBytes = TRUE)

  if (launch) {
    # RStudio's pane cannot reach tempdir() unless the file sits under the
    # session temp directory, which tempfile() guarantees
    viewer <- getOption("viewer")
    if (!is.null(viewer) && startsWith(normalizePath(file, mustWork = FALSE),
                                       normalizePath(tempdir()))) {
      viewer(file)
    } else {
      utils::browseURL(paste0("file://", normalizePath(file)))
    }
  }
  invisible(file)
}

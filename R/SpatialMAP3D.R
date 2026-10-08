#' Interactive BubbleMAP viewer over spatial coordinates
#'
#' The same viewer as \code{BubbleMAP3D()}, but the floor is the tissue rather
#' than an embedding. Coordinates come from a Seurat image / FOV, from named
#' \code{meta.data} columns, or from a plain data.frame.
#'
#' Whether to pool the data into bubbles depends on the assay. Visium spots and
#' segmented cells from imaging platforms are already the unit of measurement,
#' so aggregating them discards the resolution the experiment was run for
#' \code{mode = "point"} draws one mark per spot. For dense single-cell spatial
#' with hundreds of thousands of cells, pooling is what makes the view legible.
#' Both are available at once in the viewer via the "Original points" button,
#' as long as \code{interactive_cells} covers the data.
#'
#' Coordinates are rescaled internally so the larger of the two spans is 36,
#' matching the range a UMAP occupies. Every length in the viewer: bubble
#' size, kernel bandwidths, the radius slider — is expressed relative to that,
#' so micron and pixel coordinates behave the same. \code{radius} is given in
#' the ORIGINAL units and converted for you; the scale factor is reported.
#'
#' @param object A Seurat object, or a data.frame carrying the coordinate
#'   columns, the label column and one column per gene.
#' @param gene Character vector of genes, selectable in the viewer.
#' @param metric Numeric \code{meta.data} columns to offer as height markers
#'   alongside the genes such as an ssGSEA score, a deconvolution proportion, a
#'   distance to a landmark. Often the more interesting axis for spatial data.
#' @param image Name of the image / FOV to take coordinates from. Defaults to
#'   the first one present. Ignored when \code{coords} is given.
#' @param coords Length-2 character vector naming coordinate columns in
#'   \code{meta.data} (or in the data.frame). Use this when the coordinates
#'   are not stored as a Seurat image.
#' @param mode "bubble" pools neighboring spots; "point" draws each one.
#' @param flip_y Image coordinate systems usually run top-down while plots run
#'   bottom-up. TRUE (the default) flips y so the tissue is not upside down.
#' @param radius Grouping radius in the ORIGINAL coordinate units.
#' @param axes Axis captions. Defaults to the coordinate column names.
#' @param ... Passed to \code{BubbleMAP3D()}: label, assay, layer, min_group,
#'   min_class, palette, elevation, shell, highlight_all, flat, template, file,
#'   launch, interactive_cells, and the rest.
#'
#' @param label Column in \code{meta.data} holding the region or cell type.
#' @param assay,layer Where to read expression from. Use a normalized layer.
#' @param title Title shown in the viewer masthead.
#' @param interactive_cells Number of cells shipped to the viewer for live
#'   regrouping and point mode. Defaults to every spot in point mode. NULL
#'   omits the cell table, which hides the two grouping sliders.
#'
#' @return Path to the written HTML file, invisibly.
#' @export
SpatialMAP3D <- function(
    object,
    gene    = NULL,
    metric  = NULL,
    image   = NULL,
    coords  = NULL,
    label   = "cell_type",
    mode    = c("bubble", "point"),
    flip_y  = TRUE,
    radius  = NULL,
    assay   = NULL,
    layer   = "data",
    axes    = NULL,
    title   = NULL,
    interactive_cells = NULL,
    ...
) {
  mode <- match.arg(mode)

  # ---- 1. coordinates ------------------------------------------------------
  if (inherits(object, "Seurat")) {
    if (!is.null(coords)) {
      stopifnot(length(coords) == 2, all(coords %in% colnames(object@meta.data)))
      xy <- object@meta.data[, coords, drop = FALSE]
      nm <- coords
    } else {
      imgs <- names(object@images)
      if (!length(imgs))
        stop("No images in the object. Pass coords = c(\"<x column>\", ",
             "\"<y column>\") to read coordinates from meta.data instead.")
      if (is.null(image)) image <- imgs[1]
      if (!image %in% imgs)
        stop("Image '", image, "' not found. Available: ",
             paste(imgs, collapse = ", "))
      tc <- Seurat::GetTissueCoordinates(object, image = image)

      # v4 Visium returns imagerow/imagecol; v5 FOV objects return x/y plus a
      # cell column, and the row order is not guaranteed to match the object
      if (all(c("imagerow", "imagecol") %in% colnames(tc))) {
        xy <- tc[, c("imagecol", "imagerow")]; nm <- c("x", "y")
      } else {
        cc <- intersect(c("x", "y"), colnames(tc))
        if (length(cc) != 2) cc <- colnames(tc)[1:2]
        if ("cell" %in% colnames(tc)) rownames(tc) <- tc$cell
        xy <- tc[, cc, drop = FALSE]; nm <- cc
      }
      keep <- intersect(colnames(object), rownames(xy))
      if (!length(keep))
        stop("Tissue coordinates and cell names do not overlap for image '",
             image, "'.")
      if (length(keep) < ncol(object))
        message("Coordinates cover ", length(keep), " of ", ncol(object),
                " cells; the rest are dropped")
      object <- subset(object, cells = keep)
      xy <- xy[colnames(object), , drop = FALSE]
    }

    if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)
    is_v5 <- utils::packageVersion("Seurat") >= "5.0.0"
    gd <- if (is_v5)
      SeuratObject::GetAssayData(object, assay = assay, layer = layer)
    else
      SeuratObject::GetAssayData(object, assay = assay, slot = layer)
    if (!nrow(gd) || !ncol(gd))
      stop("Layer '", layer, "' of assay '", assay, "' is empty.")
    miss <- setdiff(gene, rownames(gd))
    if (length(miss))
      stop("Not found in assay '", assay, "': ", paste(miss, collapse = ", "))
    if (length(gene)) {
      ex <- gd[gene, , drop = FALSE]
      ex <- as.matrix(if (inherits(ex, "Matrix")) Matrix::t(ex) else t(ex))
    } else {
      # metric-only: an empty matrix still has to carry one row per cell, or
      # the cbind below has nothing to align against
      ex <- matrix(numeric(0), nrow = ncol(object), ncol = 0)
    }
    if (length(metric)) {
      mmiss <- setdiff(metric, colnames(object@meta.data))
      if (length(mmiss))
        stop("Not found in meta.data: ", paste(mmiss, collapse = ", "))
      ex <- cbind(ex, as.matrix(object@meta.data[, metric, drop = FALSE]))
    }
    lab <- as.character(object[[label]][, 1])

  } else {
    if (is.null(coords)) {
      cand <- list(c("x", "y"), c("imagecol", "imagerow"),
                   c("col", "row"), c("X", "Y"))
      hit <- Filter(function(p) all(p %in% colnames(object)), cand)
      if (!length(hit))
        stop("Could not find coordinate columns. Pass coords = c(\"x\", \"y\").")
      coords <- hit[[1]]
    }
    stopifnot(all(coords %in% colnames(object)),
              label %in% colnames(object),
              all(c(gene, metric) %in% colnames(object)))
    xy  <- object[, coords, drop = FALSE]
    nm  <- coords
    ex  <- as.matrix(object[, c(gene, metric), drop = FALSE])
    lab <- as.character(object[[label]])
  }

  colnames(ex) <- c(gene, metric)
  x <- as.numeric(xy[[1]]); y <- as.numeric(xy[[2]])
  if (flip_y) y <- -y                       # image origin is top-left

  # ---- 2. rescale to the range the viewer's constants assume ---------------
  sc <- 36 / max(diff(range(x)), diff(range(y)))
  message("Coordinate scale: ", signif(sc, 4), " (spans ",
          signif(diff(range(x)), 4), " x ", signif(diff(range(y)), 4),
          " -> ", signif(diff(range(x))*sc, 3), " x ",
          signif(diff(range(y))*sc, 3), ")")

  df <- data.frame(Dim1 = x * sc, Dim2 = y * sc,
                   stringsAsFactors = FALSE)
  df[[label]] <- lab
  for (g in c(gene, metric)) df[[g]] <- ex[, g]

  # radius is quoted in the caller's units; convert
  if (is.null(radius)) {
    # 0.42 is a UMAP number and means nothing on a spot grid, where points sit
    # at a fixed pitch. Derive the default from the actual spacing instead:
    # anything below it groups each spot with itself and no bubble forms.
    ns  <- min(1500L, nrow(df))
    si  <- sample(seq_len(nrow(df)), ns)
    dm  <- as.matrix(stats::dist(df[si, c("Dim1", "Dim2")]))
    diag(dm) <- Inf
    pitch <- stats::median(apply(dm, 1, min))
    rad <- 1.6 * pitch
    message("Point spacing ~", signif(pitch/sc, 4), " (original units); ",
            "radius defaulting to ", signif(rad/sc, 4))
  } else {
    rad <- radius * sc
    message("Radius ", radius, " -> ", signif(rad, 3), " in scaled units")
  }

  # point mode needs every spot in the viewer, not a subsample
  if (is.null(interactive_cells))
    interactive_cells <- if (mode == "point") nrow(df) else
      min(60000L, nrow(df))
  if (mode == "point" && interactive_cells < nrow(df))
    warning("mode = 'point' with interactive_cells < the number of spots: ",
            "only the subsample will be drawn as points.", call. = FALSE)

  BubbleMAP3D(
    df, gene = gene, metric = metric, label = label,
    radius = rad,
    axes = if (is.null(axes)) nm else axes,
    title = if (is.null(title)) "SpatialMAP3D" else title,
    point_mode = (mode == "point"),
    interactive_cells = interactive_cells,
    ...
  )
}

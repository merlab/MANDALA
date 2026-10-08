#' Density-contour UMAP with local-neighbourhood bubbles
#'
#' Pools cells into bubbles by a greedy spatial-hash pass over a 2-D embedding
#' and draws them over per-class kernel density contours. Bubble area always
#' encodes cell count; gene expression and metadata metrics are carried by
#' shape, alpha or color so that size stays comparable across modes.
#'
#' @param object A Seurat object, or a data.frame with Dim1/Dim2 (or
#'   UMAP1/UMAP2) columns plus the label column.
#' @param reduction Reduction supplying the coordinates.
#' @param label Column in `meta.data` holding the cell type.
#' @param highlight Cell types to emphasise.
#' @param palette Named vector of colours; NULL builds one automatically.
#' @param target_cells Cap-based stratified downsample before grouping.
#' @param min_class Types with fewer cells are collapsed into "Other".
#' @param radius,min_group,rare_relax Bubble grouping parameters.
#' @param size_by "count", "gene" or "metric".
#' @param gene,metric What to encode when `size_by` is "gene" or "metric".
#' @param assay,layer Where to read expression from. `layer` follows the
#'   Seurat v5 naming and is mapped to `slot` automatically for v3/v4 objects.
#' @param slot Deprecated alias for `layer`, retained for backward
#'   compatibility. If supplied it overrides `layer` and emits a warning.
#' @param bubble_engine,bubble_max_size,bubble_cap_q Neighbour engine for the
#'   grouping pass, the largest drawn bubble size, and the quantile at which
#'   cell counts are capped so a few dense cores do not flatten the size scale.
#' @param bubble_alpha,bubble_alpha_hi,bubble_alpha_oth Opacity for normal,
#'   highlighted and "Other" bubbles.
#' @param bubble_min_dens Drop bubbles sitting below this normalized density in
#'   their own class's envelope. NULL keeps all of them.
#' @param bw_frac,gridsize,levels_norm Bandwidth fraction, grid size and
#'   normalized levels for the per-class density contours.
#' @param alpha_top,alpha_top_hi,alpha_other Contour fill opacity for normal,
#'   highlighted and "Other" classes.
#' @param edge_all,edge_width,edge_alpha,edge_hi Thin outline drawn on each
#'   class's contour, which is what makes the blobs read as distinct shapes.
#' @param label_classes,label_exclude Which classes carry a label.
#' @param label_size,label_size_hi,label_halo,label_halo_col Label type size
#'   and the halo drawn behind it.
#' @param label_method,label_declutter,declutter_iter,declutter_step,declutter_pull
#'   Placement: "fixed" pins each label to its density peak and breaks overlaps
#'   with a bounded tethered nudge; "repel" hands over to ggrepel.
#' @param label_force,label_pull,label_lines,label_seg_min ggrepel physics and
#'   leader lines, used only when \code{label_method = "repel"}.
#' @param label_rename,label_nudge Shorter display names, and manual offsets
#'   for labels whose density peak sits inside another class.
#' @param axis_arrows Draw small UMAP axis arrows at the bottom left.
#' @param expr_cutoff,expr_shapes Threshold and the shapes used above and below
#'   it when \code{size_by = "gene"}.
#' @param metric_alpha,metric_scale,metric_trans How a metric is encoded:
#'   bubble alpha or a continuous viridis colour, on the raw values or on their
#'   empirical percentile.
#' @param metric_palette,metric_direction viridisLite option name and ramp
#'   direction for the metric colour scale.
#' @param metric_edge,metric_edge_width In metric colour mode, ring each bubble
#'   in its cell-type colour so class identity stays readable, and how thick
#'   that ring is drawn.
#' @param bg,title,subtitle Panel background colour and plot text.
#' @param celltype_legend,celltype_legend_title,legend_ncol,legend_key_size,legend_title
#'   Whether to place a full cell-type key outside the panel, and how it is
#'   titled and laid out.
#' @param keep_data Return the full per-cell data frame in the result. FALSE by
#'   default, since it keeps one row per cell alive after plotting.
#' @param seed Random seed for the greedy grouping pass.
#'
#' @return Invisibly, a list with the plot, bubbles, contours, peaks and
#'   palette.
#' @export

BubbleMAP <- function(
    object,
    reduction        = "umap",
    label            = "RefinedLabel",
    highlight        = NULL,
    palette          = NULL,
    target_cells     = 120000,
    min_class        = 200,

    # bubbles
    radius           = 0.30,    #radius of each bubble
    min_group        = 8,       #bubble being valid with this minimum of cells
    rare_relax       = 3,       # classes with n < rare_relax * min_group get
                                # a reduced threshold (NULL = uniform rule)
    bubble_engine    = c("hash", "legacy"),
    bubble_max_size  = 5.0,
    bubble_cap_q     = 0.98,
    bubble_alpha     = 0.55,
    bubble_alpha_hi  = 0.95,
    bubble_alpha_oth = 0.15,
    bubble_min_dens  = 0.10,    # drop bubbles below this normalized density
                                # (NULL = keep all)

    # contours
    bw_frac          = 0.022,
    gridsize         = 400,
    levels_norm      = c(0.10, 0.25, 0.42, 0.60, 0.80),
    alpha_top        = 0.30,
    alpha_top_hi     = 0.52,
    alpha_other      = 0.08,
    edge_all         = TRUE,    # thin outline on every class (s
    edge_width       = 0.35,
    edge_alpha       = 0.35,
    edge_hi          = TRUE,

    # labels
    label_classes    = NULL,
    label_exclude    = NULL,
    label_size       = 2.9,
    label_size_hi    = 4.6,
    label_halo       = 0.14,
    label_halo_col   = NULL,    # NULL = bg. Use RGBA hex for a translucent
                                # halo, e.g. "#FBFBFCB0" (B0 ~ 69% opaque).
                                # Note: the halo is drawn as offset copies of
                                # the text, so alpha compounds where they
                                # overlap — the result reads more opaque than
                                # the value suggests.
    label_method     = c("fixed", "repel"),   # "fixed" = no ggrepel physics
    label_declutter  = TRUE,    # tethered nudge to break overlaps (fixed mode)
    declutter_iter   = 1200,
    declutter_step   = 0.060,   # push per iteration, data units
    declutter_pull   = 0.006,   # tether back toward the density peak
    axis_arrows      = TRUE,    # small UMAP1/UMAP2 arrows, bottom-left
    label_force      = 0,       # (repel only)
    label_pull       = 8,       # high = strong tether to the peak
    label_lines      = FALSE,   # TRUE = allow leader lines when pushed far
    label_seg_min    = 0.8,     # (only used when label_lines = TRUE)
    label_rename     = NULL,    # named chr vector of shorter display names,
                                # e.g. c("EoBasoMast Precursor" = "EoBasoMast")
    label_nudge      = NULL,    # named list of c(dx, dy), e.g.
                                # list(MPP = c(-1.2, -0.8))

    # gene / metric encoding (from BubbleMAP)
    size_by          = c("count", "gene", "metric"),
    gene             = NULL,
    assay            = "SCT",
    layer            = "data",     # Seurat v5 name; mapped to slot for v3/v4
    slot             = NULL,       # deprecated alias for `layer`
    metric           = NULL,
    expr_cutoff      = 0,
    expr_shapes      = c(19, 21),  # above / below cutoff
    metric_alpha     = c(0.2, 1),
    metric_scale     = c("alpha", "colour"),
    metric_trans     = c("identity", "rank"),
    metric_palette   = "viridis",   # any viridisLite option name
    metric_direction = 1,           # -1 reverses the colour ramp
    metric_edge      = TRUE,        # colour mode: ring bubbles by cell type
    metric_edge_width = 0.6,

    bg               = "#FBFBFC",
    title            = NULL,
    subtitle         = NULL,
    celltype_legend  = FALSE,   # TRUE = full cell-type key, right of the panel
    legend_ncol      = 1,
    legend_key_size  = 4,
    legend_title     = "HSC subpopulations",
    keep_data        = FALSE,   # TRUE = return the full per-cell data frame
    celltype_legend_title = "Cell type",
    seed             = 42
) {

  requireNamespace("ggplot2"); requireNamespace("dplyr")
  requireNamespace("ggrepel"); requireNamespace("KernSmooth")
  label_method  <- match.arg(label_method)
  size_by       <- match.arg(size_by)
  bubble_engine <- match.arg(bubble_engine)
  metric_scale  <- match.arg(metric_scale)
  metric_trans  <- match.arg(metric_trans)

  # `slot` was this function's original argument name; `layer` matches Seurat
  # v5 and the other functions in this package. Both are accepted.
  if (!is.null(slot)) {
    warning("`slot` is deprecated; use `layer` instead.", call. = FALSE)
    layer <- slot
  }
  metric_colour <- (size_by == "metric" && metric_scale == "colour")
  metric_ring   <- (metric_colour && isTRUE(metric_edge))
  if (metric_colour && !requireNamespace("viridisLite", quietly = TRUE))
    stop("metric_scale = 'colour' requires the viridisLite package")
  if (is.null(label_halo_col)) label_halo_col <- bg

  # legend title and tick formatter for the metric
  metric_label <- if (size_by != "metric") NULL else if (metric_trans == "rank")
    paste0(if (identical(metric, "percent.mt")) "Percent MT" else metric,
           " (percentile)")
  else if (identical(metric, "percent.mt")) "Percent MT" else metric

  metric_fmt <- if (size_by == "metric" && metric_trans == "identity" &&
                    identical(metric, "percent.mt"))
    function(x) paste0(round(x * 100, 1), "%")
  else if (size_by == "metric" && metric_trans == "rank")
    function(x) paste0(round(x * 100), "th")
  else waiver()

  # opacity of metric-coloured bubbles (kept high so the ramp reads cleanly)
  bubble_alpha_metric <- 0.9
  set.seed(seed)

  # ---- 1. Coordinates + labels ---------------------------------------------
  expr_vec <- NULL; metric_vec <- NULL

  if (inherits(object, "Seurat")) {
    if (!reduction %in% names(object@reductions))
      stop("Reduction '", reduction, "' not found in the Seurat object")
    co  <- as.data.frame(object@reductions[[reduction]]@cell.embeddings[, 1:2])
    lab <- object[[label]][, 1]

    # Seurat v3/v4/v5 compatible accessor
    is_v5 <- utils::packageVersion("Seurat") >= "5.0.0"
    get_ad <- function(obj, assay, layer) {
      if (is_v5) {
        SeuratObject::GetAssayData(obj, assay = assay, layer = layer)
      } else {
        SeuratObject::GetAssayData(obj, assay = assay, slot = layer)
      }
    }

    if (size_by == "gene") {
      if (is.null(gene)) stop("size_by = 'gene' requires a gene name")
      gd <- get_ad(object, assay, layer)
      if (!gene %in% rownames(gd))
        stop("Gene '", gene, "' not found in assay '", assay, "'")
      expr_vec <- as.numeric(gd[gene, ])
    }
    if (size_by == "metric") {
      if (is.null(metric)) stop("size_by = 'metric' requires a metric name")
      if (!metric %in% colnames(object@meta.data))
        stop("Metric '", metric, "' not found in meta.data")
      metric_vec <- as.numeric(object@meta.data[[metric]])
      if (metric == "percent.mt") metric_vec <- metric_vec / 100
    }

  } else {
    xy_cols <- if (all(c("Dim1", "Dim2") %in% colnames(object)))
      c("Dim1", "Dim2") else c("UMAP1", "UMAP2")
    stopifnot(all(xy_cols %in% colnames(object)))
    co  <- object[, xy_cols]
    lab <- object[[label]]

    # data.frame input: gene / metric must already be columns
    if (size_by == "gene") {
      if (is.null(gene) || !gene %in% colnames(object))
        stop("For a data.frame, 'gene' must name an existing column")
      expr_vec <- as.numeric(object[[gene]])
    }
    if (size_by == "metric") {
      if (is.null(metric) || !metric %in% colnames(object))
        stop("For a data.frame, 'metric' must name an existing column")
      metric_vec <- as.numeric(object[[metric]])
      if (metric == "percent.mt") metric_vec <- metric_vec / 100
    }
  }

  names(co) <- c("Dim1", "Dim2")
  df <- data.frame(co, class = as.character(lab), stringsAsFactors = FALSE)
  if (!is.null(expr_vec))   df$expr   <- expr_vec
  if (!is.null(metric_vec)) df$metric <- metric_vec
  df <- df[stats::complete.cases(df[, c("Dim1", "Dim2", "class")]), ]

  # ---- 2. Cap-based stratified downsample ----------------------------------
  if (!is.null(target_cells) && nrow(df) > target_cells) {
    cnt <- table(df$class)
    lo <- 1L; hi <- max(cnt)
    while (lo < hi) {
      mid <- (lo + hi + 1L) %/% 2L
      if (sum(pmin(cnt, mid)) <= target_cells) lo <- mid else hi <- mid - 1L
    }
    cap <- lo
    idx <- unlist(lapply(split(seq_len(nrow(df)), df$class), function(ix)
      if (length(ix) <= cap) ix else sample(ix, cap)))
    df <- df[idx, ]
    message("Per-class cap: ", cap, " -> ", nrow(df), " cells")
  }

  # ---- 3. Collapse rare classes --------------------------------------------
  cnt  <- sort(table(df$class), decreasing = TRUE)
  keep <- setdiff(names(cnt)[cnt >= min_class], "Unknown")
  df$class <- ifelse(df$class %in% keep, df$class, "Other")

  lv <- unique(df$class)
  lv <- c(setdiff(lv, "Other"), intersect("Other", lv))

  # ---- 4. Palette -----------------------------------------------------------
  if (is.null(palette)) {
    base_cols <- if (requireNamespace("ggsci", quietly = TRUE))
      ggsci::pal_d3("category20")(20) else scales::hue_pal()(20)
    palette <- setNames(colorRampPalette(base_cols)(length(lv)), lv)
    palette["Other"] <- "#D8DCE0"
  }
  miss <- setdiff(lv, names(palette))
  if (length(miss)) {
    warning("Not in palette, filling grey: ", paste(miss, collapse = ", "))
    palette <- c(palette, setNames(rep("grey60", length(miss)), miss))
  }
  if (is.null(highlight)) highlight <- character(0)
  highlight <- intersect(highlight, lv)

  lv <- c(setdiff(lv, highlight), highlight)   # highlight last => drawn on top
  lv_legend <- c(highlight, setdiff(lv, c(highlight, "Other")))  # key order
  df$class <- factor(df$class, levels = lv)

  # ---- 5. Binned KDE + contour polygons ------------------------------------
  pad <- 0.6
  xr <- range(df$Dim1) + c(-pad, pad)
  yr <- range(df$Dim2) + c(-pad, pad)
  bw <- c(diff(xr), diff(yr)) * bw_frac

  poly_list <- list(); edge_list <- list(); peaks <- list(); kde_store <- list()
  pk <- 1L; ek <- 1L

  for (cl in lv) {
    sub <- df[df$class == cl, ]
    if (nrow(sub) < 15) next

    k <- KernSmooth::bkde2D(as.matrix(sub[, c("Dim1", "Dim2")]),
                            bandwidth = bw,
                            gridsize  = c(gridsize, gridsize),
                            range.x   = list(xr, yr))
    k$fhat <- k$fhat / max(k$fhat)

    kde_store[[cl]] <- k

    ij <- which(k$fhat == max(k$fhat), arr.ind = TRUE)[1, ]
    peaks[[cl]] <- data.frame(class = cl,
                              Dim1 = k$x1[ij[1]], Dim2 = k$x2[ij[2]])

    is_hi <- cl %in% highlight
    top <- if (is_hi) alpha_top_hi else if (cl == "Other") alpha_other else alpha_top
    alphas <- seq(top * 0.16, top, length.out = length(levels_norm))

    for (li in seq_along(levels_norm)) {
      cls <- grDevices::contourLines(k$x1, k$x2, k$fhat, levels = levels_norm[li])
      if (!length(cls)) next
      fc <- grDevices::adjustcolor(palette[[cl]], alpha.f = alphas[li])
      for (p in cls) {
        poly_list[[pk]] <- data.frame(Dim1 = p$x, Dim2 = p$y,
                                      grp = sprintf("%s_%d_%d", cl, li, pk),
                                      fillcol = fc, stringsAsFactors = FALSE)
        pk <- pk + 1L
      }
    }

    draw_edge <- (edge_all && cl != "Other") || (edge_hi && is_hi)
    if (draw_edge) {
      cls <- grDevices::contourLines(k$x1, k$x2, k$fhat, levels = levels_norm[2])
      for (p in cls) {
        edge_list[[ek]] <- data.frame(
          Dim1 = p$x, Dim2 = p$y,
          grp  = sprintf("e_%s_%d", cl, ek),
          col  = palette[[cl]],
          lw   = if (is_hi) edge_width * 3 else edge_width,
          al   = if (is_hi) 0.85 else edge_alpha,
          stringsAsFactors = FALSE)
        ek <- ek + 1L
      }
    }
  }
  poly_df <- dplyr::bind_rows(poly_list)
  edge_df <- if (length(edge_list)) dplyr::bind_rows(edge_list) else NULL
  peak_df <- dplyr::bind_rows(peaks)

  # ---- 6. Greedy bubble grouping -------------------------------------------
  # Cells are grouped into bubbles by a single greedy pass: pick an unused
  # seed, absorb every unused cell within `radius`, emit one bubble, repeat.
  #
  # The neighbour lookup uses a uniform spatial hash whose cell size is exactly
  # `radius`. Any point within `radius` of a seed must therefore lie in the
  # seed's own grid cell or one of the 8 adjacent ones, so scanning that 3x3
  # block and applying the true distance test is EXACT -- it returns the same
  # neighbours a full radius search would, without ever materialising the
  # complete neighbour list. Neighbours are computed only for cells that
  # actually become seeds (roughly n / mean_group_size of them), which is the
  # other half of the saving.

  # threshold per class: relaxed for small populations
  group_thresh <- function(n) {
    if (is.null(rare_relax) || n >= rare_relax * min_group) return(min_group)
    max(3L, min(as.integer(min_group), as.integer(ceiling(n / rare_relax))))
  }

  bubbles_of <- function(sub) {
    xy <- as.matrix(sub[, c("Dim1", "Dim2")])
    n  <- nrow(xy)
    thr <- group_thresh(n)
    if (n < thr) return(NULL)

    has_expr   <- !is.null(sub$expr)
    has_metric <- !is.null(sub$metric)

    # ---- neighbour lookup ---------------------------------------------------
    if (bubble_engine == "legacy") {
      kk <- min(80L, n)
      z  <- RANN::nn2(xy, xy, k = kk, searchtype = "radius", radius = radius)
      nb_of <- function(i) { v <- z$nn.idx[i, ]; v[v > 0] }
    } else {
      gx <- as.integer(floor(xy[, 1] / radius))
      gy <- as.integer(floor(xy[, 2] / radius))
      gx <- gx - min(gx)
      gy <- gy - min(gy)
      ny <- as.numeric(max(gy)) + 3            # numeric: avoids int overflow
      key <- as.numeric(gx) * ny + as.numeric(gy)

      ord    <- order(key)
      skey   <- key[ord]
      starts <- c(1L, which(diff(skey) != 0) + 1L)
      ends   <- c(starts[-1] - 1L, n)
      ukey   <- skey[starts]

      cell_members <- function(k) {
        p <- findInterval(k, ukey)
        if (p < 1L || p > length(ukey)) return(integer(0))
        if (ukey[p] != k) return(integer(0))
        ord[starts[p]:ends[p]]
      }
      r2 <- radius^2
      nb_of <- function(i) {
        cand <- integer(0)
        kx <- gx[i]; ky <- gy[i]
        for (dx in -1:1) for (dy in -1:1)
          cand <- c(cand, cell_members((kx + dx) * ny + (ky + dy)))
        if (!length(cand)) return(integer(0))
        d2 <- (xy[cand, 1] - xy[i, 1])^2 + (xy[cand, 2] - xy[i, 2])^2
        cand[d2 <= r2]
      }
    }

    # ---- greedy pass --------------------------------------------------------
    used <- logical(n)
    cap  <- n %/% thr + 1L
    b_x <- numeric(cap); b_y <- numeric(cap); b_n <- integer(cap)
    b_e <- if (has_expr)   numeric(cap) else NULL
    b_m <- if (has_metric) numeric(cap) else NULL
    m <- 0L

    for (i in sample.int(n)) {
      if (used[i]) next
      g <- nb_of(i)
      g <- unique(c(i, g))
      g <- g[!used[g]]
      if (length(g) < thr) next
      m <- m + 1L
      b_x[m] <- mean(xy[g, 1])
      b_y[m] <- mean(xy[g, 2])
      b_n[m] <- length(g)
      if (has_expr)   b_e[m] <- mean(sub$expr[g],   na.rm = TRUE)
      if (has_metric) b_m[m] <- mean(sub$metric[g], na.rm = TRUE)
      used[g] <- TRUE
    }

    if (m == 0L) return(NULL)
    ix  <- seq_len(m)
    res <- data.frame(Dim1 = b_x[ix], Dim2 = b_y[ix], count = b_n[ix])
    if (has_expr)   res$expr_value   <- b_e[ix]
    if (has_metric) res$metric_value <- b_m[ix]
    res
  }

  bub <- dplyr::bind_rows(lapply(lv, function(cl) {
    sub <- df[df$class == cl, ]
    b   <- bubbles_of(sub)
    if (is.null(b)) {
      message("  ", cl, ": ", nrow(sub),
              " cells -> no bubble (threshold ", group_thresh(nrow(sub)),
              " not met within radius ", radius, ")")
      return(NULL)
    }
    cbind(b, class = cl)
  }))
  if (is.null(bub) || !nrow(bub))
    stop("No bubbles formed. Try a larger radius or a smaller min_group.")
  bub$class <- factor(bub$class, levels = lv)
  bub$hi    <- bub$class %in% highlight

  # drop bubbles sitting outside their own class's density envelope
  if (!is.null(bubble_min_dens)) {
    n0 <- nrow(bub)
    dens_at <- vapply(seq_len(nrow(bub)), function(i) {
      k <- kde_store[[as.character(bub$class[i])]]
      if (is.null(k)) return(1)
      ix <- max(1L, min(length(k$x1), findInterval(bub$Dim1[i], k$x1)))
      iy <- max(1L, min(length(k$x2), findInterval(bub$Dim2[i], k$x2)))
      k$fhat[ix, iy]
    }, numeric(1))
    bub <- bub[dens_at >= bubble_min_dens, ]
    message("Bubbles outside density envelope removed: ", n0 - nrow(bub))
  }

  # optional percentile transform of the metric (bubble means)
  if (size_by == "metric" && metric_trans == "rank" &&
      !is.null(bub$metric_value)) {
    bub$metric_raw   <- bub$metric_value
    bub$metric_value <- rank(bub$metric_value,
                             na.last = "keep") / sum(!is.na(bub$metric_value))
  }

  # precompute the viridis fill as hex: the fill aesthetic is already taken by
  # the contours (scale_fill_identity), so the ramp cannot go through a scale
  if (metric_ring && !is.null(bub$metric_value)) {
    ramp <- viridisLite::viridis(256, option = metric_palette,
                                 direction = metric_direction)
    rng  <- range(bub$metric_value, na.rm = TRUE)
    idx  <- if (diff(rng) == 0) rep(128L, nrow(bub)) else
      as.integer(round(1 + 255 * (bub$metric_value - rng[1]) / diff(rng)))
    idx[is.na(idx)] <- 128L
    bub$metric_fill <- ramp[idx]
  }

  cap_n <- stats::quantile(bub$count, bubble_cap_q, names = FALSE)
  bub$size_val <- pmin(bub$count, cap_n)
  message("Bubbles: ", nrow(bub), " | counts ", min(bub$count), "-",
          max(bub$count), " (size capped at ", round(cap_n), ")")

  # ---- 7. Plot --------------------------------------------------------------
  if (is.null(label_classes)) label_classes <- setdiff(lv, "Other")
  label_classes <- setdiff(label_classes, label_exclude)
  lab_df <- peak_df[peak_df$class %in% label_classes, ]
  lab_df$hi <- lab_df$class %in% highlight
  lab_df <- lab_df[order(!lab_df$hi), ]   # highlight first = placed first

  if (!is.null(label_rename)) {
    i <- match(lab_df$class, names(label_rename))
    lab_df$class <- ifelse(is.na(i), as.character(lab_df$class),
                           label_rename[i])
  }

  # manual offsets for labels whose density peak sits inside another class
  if (!is.null(label_nudge)) {
    for (nm in names(label_nudge)) {
      i <- which(lab_df$class == nm)
      if (length(i)) {
        lab_df$Dim1[i] <- lab_df$Dim1[i] + label_nudge[[nm]][1]
        lab_df$Dim2[i] <- lab_df$Dim2[i] + label_nudge[[nm]][2]
      }
    }
  }

  # tethered declutter: push overlapping labels apart, pull back to the peak.
  # Bounded and deterministic — unlike ggrepel, labels stay on their blob.
  if (label_method == "fixed" && label_declutter && nrow(lab_df) > 1) {
    anchor <- as.matrix(lab_df[, c("Dim1", "Dim2")])
    pos    <- anchor
    # width scales with the actual font size, not just character count
    sz     <- ifelse(lab_df$hi, label_size_hi, label_size)
    chw    <- diff(xr) / 95 * (sz / 2.9)
    w      <- nchar(as.character(lab_df$class)) * chw * 0.62 + chw * 2
    h      <- diff(yr) / 55 * (sz / 2.9) * 1.35

    for (it in seq_len(declutter_iter)) {
      dx <- outer(pos[, 1], pos[, 1], "-")
      dy <- outer(pos[, 2], pos[, 2], "-")
      needx <- outer(w, w, "+") / 2
      needy <- outer(h, h, "+") / 2
      ov <- (abs(dx) < needx) & (abs(dy) < needy)
      diag(ov) <- FALSE
      if (!any(ov)) break
      d  <- sqrt(dx^2 + dy^2) + 1e-9
      px <- rowSums(ifelse(ov, dx / d, 0)) * declutter_step
      py <- rowSums(ifelse(ov, dy / d, 0)) * declutter_step
      pos[, 1] <- pos[, 1] + px + (anchor[, 1] - pos[, 1]) * declutter_pull
      pos[, 2] <- pos[, 2] + py + (anchor[, 2] - pos[, 2]) * declutter_pull
    }
    lab_df$Dim1 <- pos[, 1]
    lab_df$Dim2 <- pos[, 2]
  }

  b_oth <- bub[bub$class == "Other", ]
  b_bg  <- bub[!bub$hi & bub$class != "Other", ]
  b_hi  <- bub[bub$hi, ]

  p <- ggplot() +
    geom_polygon(data = poly_df, aes(Dim1, Dim2, group = grp, fill = fillcol)) +
    scale_fill_identity()

  if (!is.null(edge_df))
    p <- p + geom_path(data = edge_df, aes(Dim1, Dim2, group = grp),
                       colour = edge_df$col, linewidth = edge_df$lw,
                       alpha = edge_df$al, lineend = "round")

  if (nrow(b_oth))
    p <- p + geom_point(data = b_oth, aes(Dim1, Dim2, size = size_val),
                        colour = palette[["Other"]],
                        alpha = bubble_alpha_oth, stroke = 0,
                        show.legend = FALSE)
  # gene / metric encoding follows BubbleMAP: bubble AREA always encodes the
  # number of cells; gene expression is shown as shape (above / below cutoff)
  # and a metric as alpha.
  if (size_by == "gene") {
    lab_hi <- paste0("\u2265 ", expr_cutoff)
    lab_lo <- paste0("< ", expr_cutoff)
    bub$expr_grp <- factor(ifelse(!is.na(bub$expr_value) &
                                    bub$expr_value >= expr_cutoff,
                                  lab_hi, lab_lo),
                           levels = c(lab_hi, lab_lo))
    b_bg <- bub[!bub$hi & bub$class != "Other", ]
    b_hi <- bub[bub$hi, ]
  }

  if (nrow(b_bg)) {
    p <- p + if (size_by == "gene") {
      geom_point(data = b_bg,
                 aes(Dim1, Dim2, size = size_val, colour = class,
                     shape = expr_grp),
                 alpha = bubble_alpha, stroke = 0.7, fill = "white")
    } else if (metric_ring) {
      geom_point(data = b_bg,
                 aes(Dim1, Dim2, size = size_val, colour = class),
                 shape = 21, fill = b_bg$metric_fill,
                 alpha = bubble_alpha_metric, stroke = metric_edge_width)
    } else if (metric_colour) {
      geom_point(data = b_bg,
                 aes(Dim1, Dim2, size = size_val, colour = metric_value),
                 alpha = bubble_alpha_metric, stroke = 0)
    } else if (size_by == "metric") {
      geom_point(data = b_bg,
                 aes(Dim1, Dim2, size = size_val, colour = class,
                     alpha = metric_value), stroke = 0)
    } else {
      geom_point(data = b_bg,
                 aes(Dim1, Dim2, size = size_val, colour = class),
                 alpha = bubble_alpha, stroke = 0)
    }
  }

  if (nrow(b_hi)) {
    p <- p +
      geom_point(data = b_hi, aes(Dim1, Dim2, size = size_val),
                 colour = "white", alpha = 0.9, stroke = 0,
                 show.legend = FALSE)
    p <- p + if (size_by == "gene") {
      geom_point(data = b_hi,
                 aes(Dim1, Dim2, size = size_val * 0.72, colour = class,
                     shape = expr_grp),
                 alpha = bubble_alpha_hi, stroke = 0.8, fill = "white")
    } else if (metric_ring) {
      geom_point(data = b_hi,
                 aes(Dim1, Dim2, size = size_val * 0.72, colour = class),
                 shape = 21, fill = b_hi$metric_fill,
                 alpha = 1, stroke = metric_edge_width * 1.6)
    } else if (metric_colour) {
      geom_point(data = b_hi,
                 aes(Dim1, Dim2, size = size_val * 0.72,
                     colour = metric_value),
                 alpha = 1, stroke = 0)
    } else if (size_by == "metric") {
      geom_point(data = b_hi,
                 aes(Dim1, Dim2, size = size_val * 0.72, colour = class,
                     alpha = metric_value), stroke = 0)
    } else {
      geom_point(data = b_hi,
                 aes(Dim1, Dim2, size = size_val * 0.72, colour = class),
                 alpha = bubble_alpha_hi, stroke = 0)
    }
  }

  p <- p +
    (if (label_method == "fixed" &&
         requireNamespace("shadowtext", quietly = TRUE)) {
       shadowtext::geom_shadowtext(
         data = lab_df, aes(Dim1, Dim2, label = class),
         size     = ifelse(lab_df$hi, label_size_hi, label_size),
         fontface = ifelse(lab_df$hi, "bold", "plain"),
         colour   = ifelse(lab_df$hi, "#0E1A24", "#3A4650"),
         bg.colour = label_halo_col,
         bg.r      = ifelse(lab_df$hi, label_halo * 1.3, label_halo),
         check_overlap = FALSE
       )
     } else if (label_method == "fixed") {
       # no shadowtext installed: plain text, no halo
       geom_text(
         data = lab_df, aes(Dim1, Dim2, label = class),
         size     = ifelse(lab_df$hi, label_size_hi, label_size),
         fontface = ifelse(lab_df$hi, "bold", "plain"),
         colour   = ifelse(lab_df$hi, "#0E1A24", "#3A4650")
       )
     } else {
       ggrepel::geom_text_repel(
         data = lab_df, aes(Dim1, Dim2, label = class),
         size     = ifelse(lab_df$hi, label_size_hi, label_size),
         fontface = ifelse(lab_df$hi, "bold", "plain"),
         colour   = ifelse(lab_df$hi, "#0E1A24", "#3A4650"),
         bg.color = label_halo_col,
         bg.r     = ifelse(lab_df$hi, label_halo * 1.3, label_halo),
         segment.color = "#98A2AB", segment.size = 0.22, segment.alpha = 0.7,
         min.segment.length = if (label_lines) label_seg_min else Inf,
         point.padding = 0.05,
         box.padding = ifelse(lab_df$hi, 0.45, 0.14),
         force = label_force, force_pull = label_pull,
         max.overlaps = Inf, max.iter = 8000, seed = seed
       )
     }) +
    (if (metric_ring) {
       # cell type stays on colour (the ring); the metric ramp is on fill,
       # precomputed, with its colourbar added separately below
       scale_colour_manual(
         values = palette,
         breaks = if (celltype_legend) lv_legend
                  else if (length(highlight)) highlight else waiver(),
         name   = if (celltype_legend) celltype_legend_title else legend_title,
         guide  = if (celltype_legend || length(highlight))
           guide_legend(order = 1, ncol = legend_ncol,
                        override.aes = list(size = legend_key_size, alpha = 1,
                                            shape = 21, fill = "white",
                                            stroke = 1))
         else "none"
       )
     } else if (metric_colour) {
       # metric on a continuous colour ramp; cell identity stays in the
       # contours and labels
       ggplot2::scale_colour_viridis_c(
         name   = metric_label,
         option = metric_palette,
         direction = metric_direction,
         labels = metric_fmt,
         guide  = guide_colourbar(order = 1, barheight = unit(45, "pt"))
       )
     } else {
       scale_colour_manual(
         values = palette,
         breaks = if (celltype_legend) lv_legend
                  else if (length(highlight)) highlight else waiver(),
         name   = if (celltype_legend) celltype_legend_title else legend_title,
         guide  = if (celltype_legend || length(highlight))
           guide_legend(order = 1, ncol = legend_ncol,
                        override.aes = list(size = legend_key_size, alpha = 1))
         else "none"
       )
     }) +
    scale_size_area(
      max_size = bubble_max_size,
      breaks   = unique(round(stats::quantile(bub$size_val,
                                             c(0.30, 0.70, 0.97),
                                             names = FALSE))),
      name     = "Cells per bubble",
      guide    = guide_legend(order = 2,
                              override.aes = list(colour = "#9AA5B1", alpha = 0.65))
    ) +
    coord_fixed(xlim = xr, ylim = yr, expand = FALSE) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_void(base_size = 12) +
    theme(
      plot.background   = element_rect(fill = bg, colour = NA),
      panel.background  = element_rect(fill = bg, colour = NA),
      plot.title        = element_text(face = "bold", size = 17,
                                       colour = "#0E1A24", hjust = 0),
      plot.subtitle     = element_text(size = 10, colour = "#67727C", hjust = 0,
                                       margin = margin(b = 10)),
      legend.background     = element_blank(),
      legend.box.background = element_rect(fill = "white", colour = "#DDE1E5"),
      legend.box.margin     = margin(4, 4, 4, 4),
      legend.box            = "vertical",
      legend.spacing.y      = unit(2, "pt"),
      legend.margin         = margin(6, 10, 6, 10),
      legend.title      = element_text(size = 10, colour = "#3A4650"),
      legend.text       = element_text(size = 9,  colour = "#3A4650"),
      plot.margin       = margin(14, 14, 14, 14)
    )

  if (axis_arrows) {
    ax_len <- diff(xr) * 0.11
    x0 <- xr[1] + diff(xr) * 0.015
    y0 <- yr[1] + diff(yr) * 0.015
    p <- p +
      annotate("segment", x = x0, y = y0, xend = x0 + ax_len, yend = y0,
               colour = "#8A939B", linewidth = 0.5,
               arrow = grid::arrow(length = unit(6, "pt"), type = "closed")) +
      annotate("segment", x = x0, y = y0, xend = x0, yend = y0 + ax_len,
               colour = "#8A939B", linewidth = 0.5,
               arrow = grid::arrow(length = unit(6, "pt"), type = "closed")) +
      annotate("text", x = x0 + ax_len / 2, y = y0 + diff(yr) * 0.012,
               label = "UMAP 1", colour = "#8A939B", size = 2.7, vjust = 0) +
      annotate("text", x = x0 + diff(xr) * 0.012, y = y0 + ax_len / 2,
               label = "UMAP 2", colour = "#8A939B", size = 2.7,
               angle = 90, vjust = 0)
  }

  if (size_by == "gene") {
    shp <- expr_shapes
    names(shp) <- levels(bub$expr_grp)
    p <- p + scale_shape_manual(
      name   = paste0(gene, " expression"),
      values = shp,
      guide  = guide_legend(order = 3,
                            override.aes = list(size = legend_key_size,
                                                colour = "#3A4650")))
  }
  if (metric_ring) {
    # zero-size, fully transparent layer purely to register a continuous fill
    # scale, so the colourbar renders without competing for the fill aesthetic
    # used by the contours
    p <- p +
      ggnewscale::new_scale_fill() +
      geom_point(data = bub,
                 aes(Dim1, Dim2, fill = metric_value),
                 size = 0, alpha = 0, shape = 21, stroke = 0,
                 inherit.aes = FALSE) +
      ggplot2::scale_fill_viridis_c(
        name      = metric_label,
        option    = metric_palette,
        direction = metric_direction,
        labels    = metric_fmt,
        guide     = guide_colourbar(order = 3,
                                    barheight = grid::unit(45, "pt")))
  }

  if (size_by == "metric" && !metric_colour) {
    p <- p + scale_alpha_continuous(
      name   = metric_label,
      range  = metric_alpha,
      labels = metric_fmt,
      guide = guide_legend(order = 3,
                           override.aes = list(size = legend_key_size,
                                               colour = "#3A4650")))
  }

  p <- p + if (celltype_legend) {
    # outside the panel, right-hand side white space
    theme(legend.position      = "right",
          legend.justification = "center",
          legend.box           = "vertical",
          legend.box.background = element_rect(fill = NA, colour = NA),
          legend.key           = element_blank(),
          legend.spacing.y     = unit(6, "pt"))
  } else if (utils::packageVersion("ggplot2") >= "3.5.0") {
    theme(legend.position = "inside",
          legend.position.inside = c(0.99, 0.01),
          legend.justification.inside = c(1, 0))
  } else {
    theme(legend.position = c(0.99, 0.01), legend.justification = c(1, 0))
  }

  print(p)
  invisible(list(plot = p, bubbles = bub, contours = poly_df,
                 peaks = peak_df, palette = palette,
                 data = if (keep_data) df else NULL,
                 class_n = table(df$class),
                 size_by = size_by, gene = gene, metric = metric,
                 seurat_version = if (inherits(object, "Seurat"))
                   as.character(utils::packageVersion("Seurat")) else NA))
}

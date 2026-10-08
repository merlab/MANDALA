#' Point-and-click front end for the MANDALA methods
#'
#' Opens a local Shiny app over a Seurat object already in your session. All
#' five methods are available, the grouping column can be switched without
#' re-running anything by hand, and figures are saved to a folder you pick.
#'
#' Runs locally, so there is no upload step and no size limit beyond your own
#' memory — the object is passed in directly rather than serialised.
#'
#' @param object A Seurat object.
#' @param launch.browser Open in the system browser rather than the RStudio
#'   pane. TRUE by default: the interactive viewers want the room.
#' @param ... Passed to \code{shiny::runApp()}.
#'
#' @return Invisibly, the last figure rendered.
#' @export
MANDALA <- function(object, launch.browser = TRUE, ...) {

  for (p in c("shiny", "Seurat")) if (!requireNamespace(p, quietly = TRUE))
    stop("MANDALA() requires the ", p, " package", call. = FALSE)
  if (!inherits(object, "Seurat"))
    stop("object must be a Seurat object", call. = FALSE)

  METHODS <- c("BubbleMAP", "BubbleMAP3D", "SpatialMAP3D",
               "HexMAP", "MagnicMAP")
  LIVE    <- c("BubbleMAP3D", "SpatialMAP3D")

  md          <- object@meta.data
  reductions  <- names(object@reductions)
  images      <- names(object@images)
  assays      <- names(object@assays)
  n_cells     <- ncol(object)

  # Grouping candidates. Cluster assignments are routinely stored as integers
  # — kmeans()$cluster and most leiden wrappers return them that way — so a
  # type check alone would hide them. Any column with few enough distinct
  # values is offered instead, which is what the level cap is for.
  n_uniq   <- vapply(md, function(x) length(unique(x)), integer(1))
  is_cat   <- vapply(md, function(x)
    is.factor(x) || is.character(x) || is.logical(x), logical(1))
  is_disc  <- vapply(md, function(x)
    is.numeric(x) && all(x == round(x), na.rm = TRUE), logical(1))
  # QC counts are whole numbers with few enough distinct values to look like a
  # grouping, but nobody wants to color by nFeature_RNA. Exclude them.
  is_qc    <- grepl("^nCount_|^nFeature_|^percent\\.|_prob$|score$|Score$",
                    names(md))
  cat_cols <- names(md)[(is_cat | is_disc) & !is_qc]
  n_lev    <- n_uniq[cat_cols]

  # An annotation column is almost always what you want first.
  default_label <- {
    pat <- c("^cell_?type$", "cell_?type", "^ident", "annotation",
             "seurat_clusters", "cluster", "compartment", "identity")
    hit <- NULL
    for (p in pat) {
      h <- grep(p, cat_cols, ignore.case = TRUE, value = TRUE)
      if (length(h)) { hit <- h[1]; break }
    }
    if (is.null(hit)) {
      chr <- cat_cols[is_cat[cat_cols]]
      hit <- if (length(chr)) chr[1] else cat_cols[1]
    }
    hit
  }
  num_cols <- names(md)[vapply(md, is.numeric, logical(1))]

  # Hand a function only arguments it can actually receive. Most are declared
  # formals, but SpatialMAP3D() takes several through `...` and forwards them,
  # so its allowed set is its own formals plus those of the function it calls.
  # Filtering on declared formals alone silently drops those — which is how
  # min_class stopped reaching the grouping step.
  FORWARDS <- list(SpatialMAP3D = "BubbleMAP3D")

  # The five methods name the same ideas differently — seurat_obj vs object,
  # label_col vs label vs group.by, gene_name vs gene. Rather than hard-code a
  # table per function, each canonical name carries a list of candidates and
  # the first one the function actually declares is used.
  # Order matters: the specific names come first. MagnicMAP declares both
  # `group.by` (the column) and `label` (a logical, whether to draw text), so
  # matching the generic name first would push a column name into a flag.
  SYNONYMS <- list(
    label_size = c("label_size", "label.size"),
    label   = c("label_col", "group.by", "group_by", "labels", "label"),
    gene    = c("gene_name", "features", "feature", "gene"),
    palette = c("palette", "custom_colors", "custom_colours", "cols"),
    layer   = c("layer", "slot"),
    metric  = c("metric", "metric_col")
  )
  rename_args <- function(fn, args) {
    fm <- names(formals(fn))
    for (canon in names(SYNONYMS)) {
      if (is.null(args[[canon]])) next
      hit <- intersect(SYNONYMS[[canon]], fm)
      if (!length(hit)) next
      if (hit[1] != canon) { args[[hit[1]]] <- args[[canon]]; args[[canon]] <- NULL }
    }
    args
  }
  allowed_args <- function(fname) {
    fn <- get(fname, envir = asNamespace("MANDALA"))
    a <- names(formals(fn))
    if ("..." %in% a && !is.null(FORWARDS[[fname]]))
      a <- union(a, names(formals(get(FORWARDS[[fname]],
                                      envir = asNamespace("MANDALA")))))
    setdiff(a, "...")
  }
  # a single gene for the functions that take one, a vector for those that do not
  one_gene <- function(fname)
    !"metric" %in% names(formals(get(fname, envir = asNamespace("MANDALA"))))
  call_method <- function(fname, args) {
    fn <- get(fname, envir = asNamespace("MANDALA"))
    args <- rename_args(fn, args)
    # the object goes in under whatever the first formal is called
    obj <- args$object; args$object <- NULL
    args <- args[intersect(names(args), allowed_args(fname))]
    args[[names(formals(fn))[1]]] <- obj
    do.call(fn, args)
  }

  # a native folder dialog, tried in the order most likely to work from a
  # running Shiny process
  choose_dir <- function() {
    d <- NULL
    if (requireNamespace("tcltk", quietly = TRUE))
      d <- tryCatch(tcltk::tk_choose.dir(caption = "Choose an output folder"),
                    error = function(e) NULL)
    if ((is.null(d) || is.na(d)) &&
        requireNamespace("rstudioapi", quietly = TRUE) &&
        rstudioapi::isAvailable())
      d <- tryCatch(rstudioapi::selectDirectory("Choose an output folder"),
                    error = function(e) NULL)
    if (is.null(d) || (length(d) == 1 && is.na(d))) NULL else d
  }

  # the object never changes within a session, so the key is just the settings
  digest_args <- function(a) {
    a$object <- NULL
    paste(vapply(a[order(names(a))],
                 function(x) paste(as.character(x), collapse = "|"),
                 character(1)), collapse = ";")
  }

  # A light scatter of the embedding, used as the pick target for MagnicMAP.
  # Downsampled: this is for aiming, not for reading.
  picker_data <- function(red, lab, n = 25000) {
    em <- object@reductions[[red]]@cell.embeddings[, 1:2, drop = FALSE]
    i  <- if (nrow(em) > n) sample(nrow(em), n) else seq_len(nrow(em))
    data.frame(x = em[i, 1], y = em[i, 2],
               g = as.character(md[[lab]])[i], stringsAsFactors = FALSE)
  }

  # The levels that survive min_class, in the order BubbleMAP() builds its
  # own ramp over — row order of first appearance, "Other" last. Matching that
  # here means the app's palette is the one BubbleMAP would have chosen.
  level_set <- function(lab, min_class) {
    v  <- as.character(md[[lab]])
    tb <- table(v)
    keep <- setdiff(names(tb)[tb >= min_class], "Unknown")
    k <- ifelse(v %in% keep, v, "Other")
    lv <- unique(k)
    c(setdiff(lv, "Other"), intersect("Other", lv))
  }
  default_palette <- function(lv) {
    base <- if (requireNamespace("ggsci", quietly = TRUE))
      ggsci::pal_d3("category20")(20) else scales::hue_pal()(20)
    p <- stats::setNames(grDevices::colorRampPalette(base)(length(lv)), lv)
    if ("Other" %in% lv) p[["Other"]] <- "#D8DCE0"
    p
  }

  serve_dir <- tempfile("scdataviz_")
  dir.create(serve_dir)
  shiny::addResourcePath("scdvout", serve_dir)

  ui <- shiny::fluidPage(
    shiny::tags$style(shiny::HTML("
      body { font-family: system-ui, sans-serif; }
      .container-fluid { max-width: none; padding: 0 16px; }
      /* the controls are a reference column, not the subject: keep them
         narrow so the figure gets the room */
      .well { padding: 10px; }
      .blk { border:1px solid #E2E7EB; border-radius:4px; padding:10px 12px;
             margin-bottom:12px; background:#FBFBFC; }
      .blk h5 { margin:0 0 8px; font-size:12px; font-weight:600; }
      .note { font-size:11.5px; color:#5A6670; margin-top:6px; }
      .shiny-input-container { margin-bottom:8px; }
    ")),
    shiny::titlePanel("MANDALA"),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        width = 3,
        shiny::div(class = "blk",
          shiny::h5("Object"),
          shiny::div(class = "note",
                     paste0(format(n_cells, big.mark = ","), " cells \u00b7 ",
                            length(assays), " assay(s) \u00b7 ",
                            length(reductions), " reduction(s)",
                            if (length(images))
                              paste0(" \u00b7 ", length(images), " image(s)")
                            else "")),
          shiny::fluidRow(
            shiny::column(6, shiny::selectInput("assay", "Assay", assays,
                            selected = Seurat::DefaultAssay(object))),
            shiny::column(6, shiny::uiOutput("layerui"))
          )
        ),

        # grouping, level cap and min_class sit together: they interact, and
        # the consequence of the combination is reported below them
        shiny::div(class = "blk",
          shiny::h5("Grouping"),
          shiny::fluidRow(
            shiny::column(12, shiny::numericInput("maxlev", "Max levels",
                                                 50, min = 2, step = 5)),
            shiny::column(12, shiny::numericInput("minclass", "Min cells/class",
                                                 150, min = 1, step = 10))
          ),
          shiny::uiOutput("labelui"),
          shiny::div(class = "note", shiny::textOutput("groupnote"))
        ),

        shiny::div(class = "blk",
          shiny::h5("Method"),
          shiny::selectInput("method", NULL, METHODS),
          shiny::uiOutput("methodui")
        ),

        shiny::div(class = "blk",
          shiny::h5("Bubbles"),
          shiny::fluidRow(
            shiny::column(6, shiny::numericInput("radius", "Radius",
                                                 0.42, min = 0.05, step = 0.05)),
            shiny::column(6, shiny::numericInput("mingroup", "Min cells/bubble",
                                                 8, min = 3, step = 1))
          ),
          shiny::numericInput("target", "Cap cells (blank = all)",
                              if (n_cells > 150000) 100000 else NA,
                              min = 5000, step = 10000),
          shiny::div(class = "note",
            "Stratified downsample before grouping. The bubble layout is a
             spatial summary, so a cap well above the bubble count changes the
             figure very little and cuts the wait a long way.")
        ),

        shiny::div(class = "blk",
          shiny::h5("Colours"),
          shiny::div(class = "note",
                     "One palette, used by every method. Click a swatch to
                      change a type's colour."),
          shiny::uiOutput("colourui"),
          shiny::actionButton("palreset", "Reset colours", width = "100%",
                              style = "margin-top:6px;")
        ),

        shiny::div(class = "blk",
          shiny::h5("Text"),
          shiny::sliderInput("textscale", "Size", min = 0.5, max = 2.5,
                             value = 1, step = 0.1),
          shiny::checkboxInput("textbold", "Bold", FALSE),
          shiny::div(class = "note",
                     "Live \u2014 no re-render needed. Applies to the static
                      methods; the interactive viewers carry their own.")
        ),

        shiny::actionButton("go", "Render", class = "btn-primary",
                            width = "100%"),

        shiny::div(class = "blk", style = "margin-top:12px",
          shiny::h5("Save"),
          shiny::actionButton("pickdir", "Choose folder\u2026", width = "100%"),
          shiny::div(class = "note", shiny::textOutput("dirnote")),
          shiny::fluidRow(
            shiny::column(6, shiny::numericInput("w", "Width (in)", 8, step = 1)),
            shiny::column(6, shiny::numericInput("h", "Height (in)", 6, step = 1))
          ),
          shiny::fluidRow(
            shiny::column(6, shiny::selectInput("fmt", "Format",
                            c("PDF (vector)" = "pdf", "PNG" = "png",
                              "Both" = "both"))),
            shiny::column(6, shiny::numericInput("dpi", "PNG dpi", 300,
                                                 min = 72, step = 50))
          ),
          shiny::actionButton("save", "Save figure", width = "100%"),
          shiny::div(class = "note",
                     "Interactive viewers save as HTML here; use their own
                      PDF and PNG buttons for the figure at the current angle.")
        )
      ),

      shiny::mainPanel(
        width = 9,
        shiny::uiOutput("out")
      )
    )
  )

  server <- function(input, output, session) {

    pal    <- shiny::reactiveVal(NULL)       # one palette, shared by every method
    zoom   <- shiny::reactiveValues()        # display zoom, per panel
    cache  <- new.env(parent = emptyenv())   # keyed on the full argument set
    # Layers differ per assay, and an object may carry only counts — the
    # default of "data" then reads an empty layer and fails on every gene.
    output$layerui <- shiny::renderUI({
      shiny::req(input$assay)
      ly <- tryCatch(SeuratObject::Layers(object[[input$assay]]),
                     error = function(e) c("counts", "data"))
      ly <- ly[!grepl("^scale", ly)]
      shiny::selectInput("layer", "Layer", ly,
                         selected = if ("data" %in% ly) "data" else ly[1])
    })

    saved  <- shiny::reactiveVal(NULL)   # chosen output folder
    result <- shiny::reactiveVal(NULL)   # last figure: list(kind, value, name)

    # ---- grouping column choices, filtered by the level cap ---------------
    labchoices <- shiny::reactive({
      ok <- cat_cols[n_lev[cat_cols] <= (input$maxlev %||% 50)]
      if (!length(ok)) return(character(0))
      stats::setNames(ok, paste0(ok, "  (", n_lev[ok], ")"))
    })

    output$labelui <- shiny::renderUI({
      ch <- labchoices()
      if (!length(ch))
        return(shiny::div(class = "note",
                          "No column has few enough levels. Raise Max levels."))
      sel <- if (!is.null(input$label) && input$label %in% ch) input$label
             else if (default_label %in% ch) default_label else ch[1]
      shiny::selectInput("label", "Column", ch, selected = sel)
    })

    # A new grouping column means new levels, so the palette is rebuilt. Colors
    # already chosen for a level of the same name are carried over, so switching
    # away and back does not lose your edits.
    shiny::observe({
      shiny::req(input$label, input$minclass)
      lv  <- level_set(input$label, input$minclass)
      old <- shiny::isolate(pal())
      new <- default_palette(lv)
      keep <- intersect(names(new), names(old))
      if (length(keep)) new[keep] <- old[keep]
      pal(new)
    })

    output$colourui <- shiny::renderUI({
      p <- pal()
      shiny::req(p)
      shiny::tagList(lapply(seq_along(p), function(i) {
        shiny::div(style = "display:flex;align-items:center;gap:8px;margin:2px 0;",
          shiny::tags$input(type = "color", id = paste0("col_", i),
                            value = tolower(p[[i]]),
                            class = "shiny-bound-input",
                            style = "width:22px;height:22px;padding:0;border:0;
                                     background:none;cursor:pointer;",
                            onchange = sprintf(
                              "Shiny.setInputValue('col_%d', this.value);", i)),
          shiny::span(style = "font-size:11.5px;", names(p)[i]))
      }))
    })

    # collect edits back into the shared palette
    lapply(1:40, function(i) {
      shiny::observeEvent(input[[paste0("col_", i)]], {
        p <- shiny::isolate(pal())
        if (!is.null(p) && i <= length(p)) {
          p[[i]] <- input[[paste0("col_", i)]]
          pal(p)
        }
      }, ignoreInit = TRUE)
    })

    output$groupnote <- shiny::renderText({
      shiny::req(input$label)
      # empty factor levels survive subsetting and would otherwise be counted
      tb   <- table(droplevels(as.factor(md[[input$label]])))
      mc   <- input$minclass %||% 150
      low  <- tb[tb < mc]
      paste0(length(tb), " levels present \u00b7 ",
             if (!length(low)) "none collapsed"
             else paste0(length(low), " below ", mc, " cells \u2192 \"Other\" (",
                         format(sum(low), big.mark = ","), " cells: ",
                         paste(names(low), collapse = ", "), ")"))
    })

    # ---- per-method controls ----------------------------------------------
    output$methodui <- shiny::renderUI({
      m <- input$method
      met_choices <- setdiff(num_cols,
                             names(md)[is_disc & n_uniq <= (input$maxlev %||% 50)])
      gene_input <- function(multi)
        shiny::selectizeInput("gene", "Gene(s)", choices = NULL, multiple = multi,
                              options = list(placeholder = "type to search",
                                             maxOptions = 200))
      switch(m,
        "BubbleMAP" = shiny::tagList(
          shiny::selectInput("reduction", "Reduction", reductions),
          shiny::selectInput("size_by", "Encode", c("count", "gene", "metric")),
          gene_input(FALSE),
          shiny::selectInput("metric", "Metric", c("", met_choices)),
          shiny::checkboxInput("celltype_legend", "Cell-type legend", TRUE)
        ),
        "BubbleMAP3D" = shiny::tagList(
          shiny::selectInput("reduction", "Reduction", reductions),
          gene_input(TRUE),
          shiny::selectInput("metric", "Metric(s)", met_choices, multiple = TRUE),
          shiny::checkboxInput("flat", "Open flat (2-D)", FALSE),
          shiny::checkboxInput("shell", "Contour shell", TRUE)
        ),
        "SpatialMAP3D" = shiny::tagList(
          if (length(images))
            shiny::selectInput("image", "Image", images)
          else shiny::div(class = "note", "No image in this object."),
          gene_input(TRUE),
          shiny::selectInput("metric", "Metric(s)", met_choices, multiple = TRUE),
          shiny::selectInput("mode", "Marks", c("point", "bubble"))
        ),
        "HexMAP" = shiny::tagList(
          shiny::selectInput("reduction", "Reduction", reductions),
          shiny::selectInput("plot_type", "Bin summary",
                             c("majority", "entropy", "gene")),
          shiny::conditionalPanel("input.plot_type == 'gene'",
                                  gene_input(FALSE)),
          shiny::fluidRow(
            shiny::column(6, shiny::numericInput("hex_density", "Hex density",
                                                 70, min = 10, step = 10)),
            shiny::column(6, shiny::numericInput("min_cells", "Min cells/bin",
                                                 3, min = 1, step = 1))
          )
        ),
        shiny::tagList(                       # MagnicMAP
          shiny::selectInput("reduction", "Reduction", reductions),
          shiny::div(class = "note",
                     "Click to place the centre, or drag a box to set centre
                      and zoom together."),
          shiny::plotOutput("picker", height = "220px",
                            click = "pick_click", brush = shiny::brushOpts(
                              "pick_brush", resetOnNew = TRUE)),
          shiny::fluidRow(
            shiny::column(6, shiny::numericInput("poi_x", "Centre x", 0, step = 1)),
            shiny::column(6, shiny::numericInput("poi_y", "Centre y", 0, step = 1))
          ),
          shiny::sliderInput("zoom_factor", "Zoom", min = 1.5, max = 20,
                             value = 4, step = 0.5),
          shiny::checkboxInput("autorender", "Re-render on each pick", TRUE)
        )
      )
    })

    shiny::observe({
      shiny::req(input$method, input$assay)
      g <- rownames(object[[input$assay]])
      shiny::updateSelectizeInput(session, "gene", choices = g, server = TRUE)
    })

    output$picker <- shiny::renderPlot({
      shiny::req(input$method == "MagnicMAP", input$reduction, input$label)
      d <- picker_data(input$reduction, input$label)
      ggplot2::ggplot(d, ggplot2::aes(x, y, colour = g)) +
        ggplot2::geom_point(size = 0.2, alpha = 0.5, show.legend = FALSE) +
        ggplot2::annotate("point", x = input$poi_x %||% 0, y = input$poi_y %||% 0,
                          shape = 3, size = 4, colour = "#0E1A24") +
        ggplot2::coord_fixed() +
        ggplot2::theme_void() +
        ggplot2::theme(plot.background =
                         ggplot2::element_rect(fill = "#FBFBFC", colour = NA))
    })

    shiny::observeEvent(input$pick_click, {
      shiny::updateNumericInput(session, "poi_x",
                                value = round(input$pick_click$x, 2))
      shiny::updateNumericInput(session, "poi_y",
                                value = round(input$pick_click$y, 2))
      if (isTRUE(input$autorender)) render_now()
    })

    # a dragged box gives the centre and the zoom at once: zoom is the ratio
    # of the whole embedding to the box, so the box becomes the magnified view
    shiny::observeEvent(input$pick_brush, {
      b <- input$pick_brush
      d <- picker_data(input$reduction, input$label, n = 5000)
      z <- max(1.5, min(20, round(diff(range(d$x)) / max(b$xmax - b$xmin, 1e-6), 1)))
      shiny::updateNumericInput(session, "poi_x", value = round((b$xmin+b$xmax)/2, 2))
      shiny::updateNumericInput(session, "poi_y", value = round((b$ymin+b$ymax)/2, 2))
      shiny::updateSliderInput(session, "zoom_factor", value = z)
      if (isTRUE(input$autorender)) render_now()
    })

    pal_settled <- shiny::debounce(shiny::reactive(pal()), 600)
    shiny::observeEvent(pal_settled(), {
      if (!is.null(result())) render_now()
    }, ignoreInit = TRUE)

    zoom_settled <- shiny::debounce(shiny::reactive(input$zoom_factor), 400)
    shiny::observeEvent(zoom_settled(), {
      if (isTRUE(input$autorender) && identical(input$method, "MagnicMAP"))
        render_now()
    }, ignoreInit = TRUE)

    # ---- render ------------------------------------------------------------
    render_now <- function() {
      shiny::req(input$label, input$method)
      m <- input$method

      args <- list(
        object      = object,
        label       = input$label,
        assay       = input$assay,
        layer       = input$layer,
        slot        = input$layer,
        palette     = pal(),
        min_class   = input$minclass,
        radius      = input$radius,
        min_group   = input$mingroup,
        target_cells = if (is.null(input$target) || is.na(input$target))
                         NULL else input$target,
        reduction   = input$reduction,
        gene        = if (length(input$gene) && any(nzchar(input$gene)))
                        input$gene else NULL,
        metric      = if (length(input$metric) && any(nzchar(input$metric)))
                        input$metric else NULL,
        size_by     = input$size_by,
        plot_type   = input$plot_type,
        hex_density = input$hex_density,
        min_cells   = input$min_cells,
        zoom_factor = input$zoom_factor,
        # always request the panels individually: a combined patchwork has no
        # layers to scale and takes a theme on its last panel only, so text
        # controls would silently do nothing. The app combines them instead.
        combine_plots = FALSE,
        point_of_interest = if (!is.null(input$poi_x) && !is.null(input$poi_y))
                              c(input$poi_x, input$poi_y) else NULL,
        celltype_legend = input$celltype_legend,
        image       = input$image,
        mode        = input$mode,
        flat        = input$flat,
        shell       = input$shell,
        launch      = FALSE
      )
      keep_null <- "target_cells"
      args <- args[!vapply(args, is.null, logical(1)) |
                     names(args) %in% keep_null]

      key <- paste(m, digest_args(args))
      if (!is.null(cache[[key]])) { result(cache[[key]]); return(NULL) }

      shiny::withProgress(message = paste("Running", m), value = 0.3, {
        res <- tryCatch(call_method(m, args),
                        error = function(e) e)
        if (inherits(res, "error")) {
          # one id, so a repeated failure replaces its predecessor instead of
          # stacking a toast per attempt
          shiny::showNotification(conditionMessage(res), type = "error",
                                  duration = 10, id = "renderfail")
          return(NULL)
        }
        stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        if (m %in% LIVE) {
          f <- file.path(serve_dir, paste0(m, "_", stamp, ".html"))
          file.copy(res, f, overwrite = TRUE)
          out <- list(kind = "html", value = f,
                      name = paste0(m, "_", input$label, "_", stamp, ".html"))
          cache[[key]] <- out; result(out)
        } else {
          # a method may hand back one plot, or several under its own names
          p <- if (inherits(res, "ggplot")) list(Figure = res)
               else if (!is.null(res$plot)) list(Figure = res$plot)
               else Filter(function(x) inherits(x, "ggplot"), res)
          if (!length(p)) stop("no plot returned by ", m)
          out <- list(kind = "plot", value = p,
                      name = paste0(m, "_", input$label, "_", stamp, ".pdf"))
          cache[[key]] <- out; result(out)
        }
      })
    }

    shiny::observeEvent(input$go, render_now())

    # Rendered separately, and reading its own value only inside isolate():
    # a renderUI that both creates an input and depends on it re-runs every
    # time the input changes, which loops forever.
    output$panelui <- shiny::renderUI({
      r <- result()
      if (is.null(r) || r$kind != "plot" || length(r$value) < 2) return(NULL)
      shiny::radioButtons("panel", "Panels", inline = TRUE,
                          choices = c("Both", names(r$value)),
                          selected = shiny::isolate(input$panel) %||% "Both")
    })

    output$out <- shiny::renderUI({
      r <- result()
      if (is.null(r))
        return(shiny::div(class = "note",
                          "Choose a grouping column and a method, then Render."))
      if (r$kind == "html")
        return(shiny::tags$iframe(src = file.path("scdvout", basename(r$value)),
                 style = "width:100%;height:90vh;border:1px solid #E2E7EB;"))
      keep <- seq_along(r$value)   # all panels exist; visibility is below
      n <- length(keep)
      shiny::tagList(
        shiny::uiOutput("panelui"),
        shiny::div(class = "note",
                   "Drag a box and double-click inside it to zoom;
                    double-click again to reset."),
        lapply(keep, function(i)
          shiny::conditionalPanel(
            condition = if (n == 1) "true" else
              sprintf("input.panel == 'Both' || input.panel == '%s'",
                      names(r$value)[i]),
            shiny::plotOutput(paste0("gg", i),
                              height = paste0(round(88/n), "vh"),
                              brush = shiny::brushOpts(paste0("gg", i, "_brush"),
                                                       resetOnNew = TRUE),
                              dblclick = paste0("gg", i, "_dbl"))))
      )
    })

    # which panels the radio selection leaves visible
    panel_keep <- function(r) {
      if (length(r$value) <= 1) return(1L)
      sel <- input$panel %||% "Both"
      if (identical(sel, "Both")) seq_along(r$value)
      else which(names(r$value) == sel)
    }

    # In-plot labels are a geom size, not a theme element, so the theme cannot
    # reach them. Walk the finished plot's text layers instead — the cached
    # original is untouched, since modifying a copy is what R does here, so the
    # slider stays live and reversible rather than triggering a re-render.
    scale_layer_text <- function(p, m, bold) {
      if (inherits(p, "patchwork")) return(p)   # no layers of its own
      for (i in seq_along(p$layers)) {
        g <- class(p$layers[[i]]$geom)[1]
        if (!grepl("text|label", g, ignore.case = TRUE)) next
        sz <- p$layers[[i]]$aes_params$size
        if (!is.null(sz)) p$layers[[i]]$aes_params$size <- sz * m
        if (bold) p$layers[[i]]$aes_params$fontface <- "bold"
      }
      p
    }

    # Scaling the base size alone does nothing where a theme sets an element
    # explicitly, which these do, so each of those elements is set as well.
    text_theme <- shiny::reactive({
      m <- input$textscale %||% 1
      f <- if (isTRUE(input$textbold)) "bold" else "plain"
      ggplot2::theme(
        text         = ggplot2::element_text(size = 11 * m, face = f),
        plot.title   = ggplot2::element_text(size = 17 * m, face = f),
        plot.subtitle= ggplot2::element_text(size = 10 * m, face = f),
        axis.text    = ggplot2::element_text(size = 9  * m, face = f),
        axis.title   = ggplot2::element_text(size = 11 * m, face = f),
        legend.text  = ggplot2::element_text(size = 9  * m, face = f),
        legend.title = ggplot2::element_text(size = 10 * m, face = f),
        strip.text   = ggplot2::element_text(size = 10 * m, face = f))
    })

    # one renderer and one zoom handler per possible panel
    lapply(1:3, function(i) {
      output[[paste0("gg", i)]] <- shiny::renderPlot({
        r <- result()
        shiny::req(r, r$kind == "plot", length(r$value) >= i)
        p <- r$value[[i]]
        p <- scale_layer_text(p, input$textscale %||% 1, isTRUE(input$textbold))
        # `&` reaches every panel of a patchwork; `+` would reach only the last
        p <- if (inherits(p, "patchwork")) p & text_theme() else p + text_theme()
        z <- zoom[[paste0("p", i)]]
        if (!is.null(z))
          # these are all embeddings, so hold the aspect ratio while zooming
          p <- p + ggplot2::coord_fixed(xlim = z$x, ylim = z$y, expand = FALSE)
        print(p)
      })

      shiny::observeEvent(input[[paste0("gg", i, "_dbl")]], {
        b <- input[[paste0("gg", i, "_brush")]]
        zoom[[paste0("p", i)]] <-
          if (is.null(b)) NULL else list(x = c(b$xmin, b$xmax),
                                         y = c(b$ymin, b$ymax))
      })
    })

    # ---- save --------------------------------------------------------------
    shiny::observeEvent(input$palreset, {
      shiny::req(input$label, input$minclass)
      pal(default_palette(level_set(input$label, input$minclass)))
    })

    shiny::observeEvent(input$pickdir, {
      d <- choose_dir()
      if (is.null(d)) {
        shiny::showNotification(
          "No folder chosen. Install the tcltk package if no dialog appeared.",
          type = "warning")
      } else saved(d)
    })

    output$dirnote <- shiny::renderText({
      if (is.null(saved())) "No folder chosen yet." else saved()
    })

    shiny::observeEvent(input$save, {
      r <- result()
      if (is.null(r)) {
        shiny::showNotification("Render something first.", type = "warning")
        return(NULL)
      }
      d <- saved()
      if (is.null(d)) {
        shiny::showNotification("Choose an output folder first.",
                                type = "warning")
        return(NULL)
      }
      dest <- file.path(d, r$name)
      if (r$kind == "html") {
        file.copy(r$value, dest, overwrite = TRUE)
        msg <- dest
      } else {
        nm  <- names(r$value)
        vis <- panel_keep(r)
        prep <- function(i) {
          p <- scale_layer_text(r$value[[i]], input$textscale %||% 1,
                                isTRUE(input$textbold))
          p <- p + text_theme()
          z <- zoom[[paste0("p", i)]]
          if (!is.null(z))
            p <- p + ggplot2::coord_fixed(xlim = z$x, ylim = z$y, expand = FALSE)
          p
        }
        # one helper, so every path honours the chosen format
        write_fig <- function(path_noext, p, hmul = 1) {
          out <- character(0)
          if (input$fmt %in% c("pdf", "both")) {
            f <- paste0(path_noext, ".pdf")
            ggplot2::ggsave(f, p, width = input$w, height = input$h * hmul,
                            device = grDevices::cairo_pdf)
            out <- c(out, f)
          }
          if (input$fmt %in% c("png", "both")) {
            f <- paste0(path_noext, ".png")
            ggplot2::ggsave(f, p, width = input$w, height = input$h * hmul,
                            dpi = input$dpi %||% 300, bg = "white")
            out <- c(out, f)
          }
          out
        }
        stem <- sub("\\.pdf$", "", dest)

        if (length(vis) > 1 && requireNamespace("patchwork", quietly = TRUE)) {
          # both panels visible: save the pair as one figure, stacked as shown
          f <- write_fig(stem, patchwork::wrap_plots(lapply(vis, prep), ncol = 1),
                         hmul = length(vis))
          msg <- paste(basename(f), collapse = ", ")
          shiny::showNotification(paste("Saved:", msg), type = "message",
                                  duration = 8)
          return(NULL)
        }
        files <- unlist(lapply(vis, function(i)
          write_fig(if (length(r$value) == 1) stem else paste0(stem, "_", nm[i]),
                    prep(i))))
        msg <- paste(basename(files), collapse = ", ")
      }
      shiny::showNotification(paste("Saved:", msg), type = "message",
                              duration = 8)
    })

    session$onSessionEnded(function() unlink(serve_dir, recursive = TRUE))
  }

  shiny::runApp(shiny::shinyApp(ui, server),
                launch.browser = launch.browser, ...)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

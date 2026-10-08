#' Create Magnified View Plots for Single-Cell Data Visualization (Seurat Compatible)
#'
#' @description
#' `MagnicMAP` creates side-by-side plots showing both an overview and a magnified view of 
#' a specific region of interest in single-cell dimensionality reduction plots.
#' Compatible with Seurat v3, v4, and v5. Automatically detects version and uses appropriate data access methods.
#'
#' @param seurat_obj A Seurat object containing single-cell data with dimensionality reduction computed.
#' @param point_of_interest Numeric vector of length 2 specifying the x,y coordinates of the 
#'   center point for magnification. These should be coordinates in the reduced dimensional space.
#' @param reduction Character string specifying the dimensionality reduction to use for plotting.
#'   Default is "umap". Common options include "umap", "tsne", "pca".
#' @param group.by Character string specifying the metadata column to use for cell grouping/coloring.
#'   If NULL (default), uses the default identity from the Seurat object.
#' @param zoom_factor Numeric value specifying the magnification level. Default is 4.
#'   Higher values create more zoomed-in views (smaller regions), lower values show larger regions.
#' @param colors Character vector of custom colors to use for cell groups. If NULL (default),
#'   uses the specified color_palette. This parameter is deprecated in favor of custom_colors
#'   for consistency with other functions.
#' @param color_palette Character string specifying the color palette to use when custom_colors
#'   is not provided. Default is "npg". Available options include "npg", "aaas", "nejm", 
#'   "jama", "jco", "lancet", "d3".
#' @param custom_colors Character vector of custom colors to use instead of predefined palettes.
#'   If provided, overrides both the `color_palette` and `colors` parameters.
#' @param pt.size Numeric value specifying the size of points in the plots. Default is 1.
#' @param label Logical value indicating whether to add text labels to clusters. Default is FALSE.
#' @param label.size Numeric value specifying the size of cluster labels when label = TRUE. Default is 4.
#' @param marker_color Character string specifying the color of the point of interest marker.
#'   Default is "black".
#' @param marker_size Numeric value specifying the size of the point of interest marker. Default is 3.
#' @param marker_shape Numeric value specifying the shape of the point of interest marker.
#'   Default is 3 (plus sign). See ggplot2 shape codes for options.
#' @param box_color Character string specifying the color of the zoom region box in the main plot.
#'   Default is "black".
#' @param box_linetype Character string specifying the line type of the zoom region box.
#'   Default is "dashed". Options include "solid", "dashed", "dotted", "dotdash", "longdash", "twodash".
#' @param combine_plots Logical value indicating whether to combine the main and zoom plots
#'   into a single figure. Default is TRUE. If FALSE, returns a list with separate plots.
#' @param plot_ratio Numeric vector of length 2 specifying the height ratio between main and zoom plots
#'   when combine_plots = TRUE. Default is c(2,2) for equal heights.
#'
#' @return
#' If combine_plots = TRUE (default), returns a combined patchwork plot object with version info.
#' If combine_plots = FALSE, returns a list containing the separate plots and version information.
#'
#' @export
MagnicMAP <- function(
    seurat_obj,
    point_of_interest,
    reduction = "umap",            
    group.by = NULL,               
    zoom_factor = 4,               
    colors = NULL,                 # Deprecated, use custom_colors
    color_palette = "npg",         
    custom_colors = NULL,
    pt.size = 1,                   
    label = FALSE,                 
    label.size = 4,               
    marker_color = "black",        
    marker_size = 3,              
    marker_shape = 3,             
    box_color = "black",          
    box_linetype = "dashed",      
    combine_plots = TRUE,         
    plot_ratio = c(2,2)          
) {
  # ========== SEURAT VERSION COMPATIBILITY ==========
  seurat_version <- packageVersion("Seurat")
  is_v5 <- seurat_version >= "5.0.0"
  
  # Compatible embeddings extraction
  get_embeddings_compat <- function(seurat_obj, reduction) {
    if (reduction %in% names(seurat_obj@reductions)) {
      return(seurat_obj@reductions[[reduction]]@cell.embeddings)
    } else {
      stop(paste("Reduction", reduction, "not found in Seurat object"))
    }
  }
  
  # Inform user about version compatibility
  message(paste("MagnicMAP: Using Seurat version", as.character(seurat_version)))
  if (!is_v5) {
    message("Running in v4/v3 compatibility mode")
  }
  # ===============================================
  
  # Handle deprecated colors parameter
  if (!is.null(colors) && is.null(custom_colors)) {
    warning("The 'colors' parameter is deprecated. Use 'custom_colors' instead.")
    custom_colors <- colors
  }
  
  # Input validation
  if (!reduction %in% names(seurat_obj@reductions)) {
    stop(sprintf("Reduction '%s' not found in Seurat object", reduction))
  }
  
  if (!is.null(group.by) && !group.by %in% names(seurat_obj@meta.data)) {
    stop(sprintf("group.by column '%s' not found in meta.data", group.by))
  }
  
  # Check required packages
  required_packages <- c("ggplot2", "Seurat")
  if (combine_plots) {
    required_packages <- c(required_packages, "patchwork")
  }
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is required but not installed. Please install it with install.packages('", pkg, "')"))
    }
  }
  
  # Function definitions inside MagnicMAP
  get_color_scale <- function(palette_name, custom_colors = NULL) {
    # Helper function for adaptive palette
    adaptive_pal_inner <- function(values) {
      force(values)
      n_colors <- length(values)
      function(n) {
        if (n <= n_colors) {
          values[seq_len(n)]
        } else {
          colorRampPalette(values, alpha = TRUE)(n)
        }
      }
    }
    
    # Use custom colors if provided
    if (!is.null(custom_colors)) {
      raw_cols <- custom_colors
      names(raw_cols) <- paste0("color", seq_along(raw_cols))
    } else {
      # Use built-in palettes from ggsci if available
      if (requireNamespace("ggsci", quietly = TRUE) && 
          palette_name %in% c("npg", "aaas", "nejm", "jama", "jco", "lancet", "d3")) {
        
        raw_cols <- switch(palette_name,
                           "npg"    = ggsci::pal_npg("nrc")(10),
                           "aaas"   = ggsci::pal_aaas("default")(10),
                           "nejm"   = ggsci::pal_nejm("default")(8),
                           "jama"   = ggsci::pal_jama("default")(7),
                           "jco"    = ggsci::pal_jco("default")(10),
                           "lancet" = ggsci::pal_lancet("lanonc")(9),
                           "d3"     = ggsci::pal_d3("category20")(20)
        )
      } else {
        # Default colors from ggplot2
        warning("Using default ggplot2 colors. Either ggsci is not available, or the specified palette name is invalid.")
        raw_cols <- scales::hue_pal()(8)
        names(raw_cols) <- paste0("color", 1:8)
      }
    }
    
    # Process colors
    raw_cols_rgb <- col2rgb(raw_cols)
    alpha_cols <- rgb(
      raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
      alpha = 255L, names = names(raw_cols),
      maxColorValue = 255L
    )
    
    # Return discrete scale
    discrete_scale("colour", "custom", adaptive_pal_inner(unname(alpha_cols)))
  }
  
  # Rest of the MagnicMAP function using compatible embeddings
  umap_coords <- get_embeddings_compat(seurat_obj, reduction)
  xlim <- range(umap_coords[,1])
  ylim <- range(umap_coords[,2])
  
  x_range <- diff(xlim)
  y_range <- diff(ylim)
  aspect_ratio <- y_range / x_range
  
  zoom_width <- x_range / zoom_factor
  zoom_height <- y_range / zoom_factor
  
  zoom_coords <- list(
    xmin = point_of_interest[1] - (zoom_width / 2),
    xmax = point_of_interest[1] + (zoom_width / 2),
    ymin = point_of_interest[2] - (zoom_height / 2),
    ymax = point_of_interest[2] + (zoom_height / 2)
  )
  
  plot_theme <- theme_minimal() +
    theme(aspect.ratio = aspect_ratio,
          panel.border = element_rect(fill = NA, color = "black"))
  
  plot_params <- list(
    object = seurat_obj,
    reduction = reduction,
    label = label,
    label.size = label.size,
    pt.size = pt.size
  )
  
  if (!is.null(group.by)) {
    plot_params$group.by <- group.by
  }
  
  # Create main plot
  main_plot <- tryCatch({
    do.call(Seurat::DimPlot, plot_params) +
      plot_theme +
      coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE)
  }, error = function(e) {
    stop(paste("Error creating main plot:", e$message))
  })
  
  # Apply color scheme
  if (!is.null(custom_colors)) {
    main_plot <- main_plot + scale_color_manual(values = custom_colors)
  } else {
    main_plot <- main_plot + get_color_scale(color_palette, custom_colors)
  }
  
  main_plot <- main_plot +
    annotate("rect",
             xmin = zoom_coords$xmin, xmax = zoom_coords$xmax,
             ymin = zoom_coords$ymin, ymax = zoom_coords$ymax,
             fill = NA, color = box_color, linetype = box_linetype) +
    annotate("point",
             x = point_of_interest[1],
             y = point_of_interest[2],
             color = marker_color,
             size = marker_size,
             shape = marker_shape) +
    labs(title = paste0("Overview (Seurat ", as.character(seurat_version), ")"))
  
  # Create zoom plot
  zoom_plot <- tryCatch({
    do.call(Seurat::DimPlot, plot_params) +
      plot_theme +
      coord_cartesian(
        xlim = c(zoom_coords$xmin, zoom_coords$xmax),
        ylim = c(zoom_coords$ymin, zoom_coords$ymax),
        expand = FALSE
      )
  }, error = function(e) {
    stop(paste("Error creating zoom plot:", e$message))
  })
  
  # Apply color scheme to zoom plot
  if (!is.null(custom_colors)) {
    zoom_plot <- zoom_plot + scale_color_manual(values = custom_colors)
  } else {
    zoom_plot <- zoom_plot + get_color_scale(color_palette, custom_colors)
  }
  
  zoom_plot <- zoom_plot +
    annotate("point",
             x = point_of_interest[1],
             y = point_of_interest[2],
             color = marker_color,
             size = marker_size,
             shape = marker_shape) +
    labs(title = paste0("Magnified View (", zoom_factor, "x)"))
  
  # Return results based on combine_plots parameter
  if (combine_plots) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      result <- main_plot / zoom_plot + patchwork::plot_layout(heights = plot_ratio)
      attr(result, "seurat_version") <- as.character(seurat_version)
      return(result)
    } else {
      warning("patchwork package not available. Returning separate plots.")
      combine_plots <- FALSE
    }
  }
  
  if (!combine_plots) {
    result <- list(
      main_plot = main_plot, 
      zoom_plot = zoom_plot,
      seurat_version = as.character(seurat_version),
      zoom_coordinates = zoom_coords
    )
    return(result)
  }
}

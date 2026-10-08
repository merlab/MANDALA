#' Create Hexagonal Bin Plots for Single-Cell Data Visualization (Seurat Compatible)
#'
#' @description
#' `HexMAP` creates hexagonal bin plots for single-cell RNA-seq data visualization.
#' Compatible with Seurat v3, v4, and v5. Automatically detects version and uses appropriate data access methods.
#'
#' @param seurat_obj A Seurat object containing single-cell data with dimensionality reduction computed.
#' @param reduction Character string specifying the dimensionality reduction to use for plotting.
#'   Default is "umap". Common options include "umap", "tsne", "pca".
#' @param label_col Character string specifying the metadata column containing cell type labels
#'   or classifications. This column must exist in the Seurat object's metadata.
#' @param hex_density Numeric value controlling the density of hexagonal bins. Default is 50.
#'   Higher values create more (smaller) bins, lower values create fewer (larger) bins.
#' @param plot_type Character string specifying the type of visualization. Must be one of:
#'   \itemize{
#'     \item "majority" - Color bins by majority cell type with transparency showing proportion
#'     \item "entropy" - Color bins by Shannon entropy (cell type diversity)
#'     \item "gene" - Color bins by majority cell type with transparency showing gene expression
#'   }
#' @param gene_name Character string specifying the gene name for expression visualization.
#'   Required when `plot_type = "gene"`. Must match a gene name in the specified assay.
#' @param assay Character string specifying the assay to use for gene expression data.
#'   Default is "SCT". Common options include "RNA", "SCT".
#' @param layer Character string specifying the data layer/slot to use from the assay.
#'   Default is "scale.data". Options include "counts", "data", "scale.data".
#' @param min_cells Integer specifying the minimum number of cells required per hexagonal bin
#'   for inclusion in the plot. Default is 3. Bins with fewer cells are filtered out.
#' @param color_palette Character string specifying the color palette for cell type coloring.
#'   Default is "npg". Available options include "npg", "aaas", "nejm", "jama", "jco", 
#'   "lancet", "d3".
#' @param custom_colors Character vector of custom colors to use instead of predefined palettes.
#'   If provided, overrides the `color_palette` parameter.
#' @param alpha Numeric value specifying the minimum alpha (transparency) level for bins.
#'   Default is 0.3. Alpha ranges from this value to 1.0 based on the mapped variable.
#' @param viridis Character string specifying the viridis color palette for entropy plots.
#'   Default is "magma". Options include "viridis", "plasma", "inferno", "magma", "cividis".
#'
#' @return
#' A ggplot2 object containing the hexagonal bin plot with version information in the title.
#'
#' @export
HexMAP <- function(
    seurat_obj,
    reduction = "umap",
    label_col,
    hex_density = 50,
    plot_type = c("majority", "entropy", "gene"),
    gene_name = NULL,
    assay = "SCT",
    layer = "data",
    min_cells = 3,
    color_palette = "npg",
    custom_colors = NULL,
    alpha = 0.3,
    viridis = "magma"
) {
  # ========== SEURAT VERSION COMPATIBILITY ==========
  seurat_version <- packageVersion("Seurat")
  is_v5 <- seurat_version >= "5.0.0"
  
  # Create compatible GetAssayData wrapper
  get_assay_data_func <- if (is_v5) {
    function(object, assay = NULL, layer = "data", ...) {
      SeuratObject::GetAssayData(object = object, assay = assay, layer = layer, ...)
    }
  } else {
    function(object, assay = NULL, layer = "data", ...) {
      # Map v5 layer names to v4/v3 slot names
      slot_mapping <- c(
        "data" = "data",
        "counts" = "counts", 
        "scale.data" = "scale.data",
        "scaled.data" = "scale.data"  # Alternative naming
      )
      
      slot_name <- slot_mapping[layer]
      if (is.na(slot_name)) {
        warning(paste("Unknown layer:", layer, "- defaulting to 'data'"))
        slot_name <- "data"
      }
      
      SeuratObject::GetAssayData(object = object, assay = assay, slot = slot_name, ...)
    }
  }
  
  # Compatible embeddings extraction
  get_embeddings_compat <- function(seurat_obj, reduction) {
    if (reduction %in% names(seurat_obj@reductions)) {
      return(seurat_obj@reductions[[reduction]]@cell.embeddings)
    } else {
      stop(paste("Reduction", reduction, "not found in Seurat object"))
    }
  }
  
  # Inform user about version compatibility
  message(paste("HexMAP: Using Seurat version", as.character(seurat_version)))
  if (!is_v5) {
    message("Running in v4/v3 compatibility mode")
  }
  # ===============================================
  
  plot_type <- match.arg(plot_type)
  if (plot_type == "gene" && is.null(gene_name)) {
    stop("gene_name must be provided when plot_type is 'gene'")
  }
  
  # Check required packages
  required_packages <- c("ggplot2", "dplyr")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(paste0("Package '", pkg, "' is required but not installed. Please install it with install.packages('", pkg, "')"))
    }
  }
  
  get_color_scale <- function(palette_name, custom_colors = NULL) {
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
    
    raw_cols_rgb <- col2rgb(raw_cols)
    alpha_cols <- rgb(
      raw_cols_rgb[1L, ], raw_cols_rgb[2L, ], raw_cols_rgb[3L, ],
      alpha = 255L, names = names(raw_cols),
      maxColorValue = 255L
    )
    
    discrete_scale("fill", "custom", adaptive_pal_inner(unname(alpha_cols)))
  }
  
  # Use compatible embeddings extraction
  umap_coords <- get_embeddings_compat(seurat_obj, reduction)
  
  hb_umap_df <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    class = seurat_obj@meta.data[[label_col]]
  )
  
  x_range <- diff(range(hb_umap_df$UMAP1))
  y_range <- diff(range(hb_umap_df$UMAP2))
  hex_size <- min(x_range, y_range) / hex_density
  
  hex_coords <- function(x, y, size) {
    sqrt3 <- sqrt(3)
    col <- round(x / (size * 1.5))
    row <- round((y - (col %% 2) * size * sqrt3/2) / (size * sqrt3))
    list(
      x = col * size * 1.5,
      y = row * size * sqrt3 + (col %% 2) * size * sqrt3/2
    )
  }
  
  coords <- hex_coords(hb_umap_df$UMAP1, hb_umap_df$UMAP2, hex_size)
  hb_umap_df$hex_x <- coords$x
  hb_umap_df$hex_y <- coords$y
  
  if (plot_type == "gene") {
    # Use compatible GetAssayData function
    tryCatch({
      gene_expr <- get_assay_data_func(seurat_obj, assay = assay, layer = layer)[gene_name,]
      hb_umap_df$gene_expr <- gene_expr[rownames(hb_umap_df)]
    }, error = function(e) {
      stop(paste("Error getting gene expression data:", e$message, 
                 "\nMake sure the gene exists and the assay/layer is correct for your Seurat version"))
    })
  }
  
  # Helper function to safely get majority class
  safe_majority <- function(class_vec) {
    if (length(class_vec) == 0) {
      return(NA_character_)
    }
    tab <- table(class_vec)
    if (length(tab) == 0) {
      return(NA_character_)
    }
    names(sort(tab, decreasing = TRUE))[1]
  }
  
  # Helper function to safely calculate entropy
  safe_entropy <- function(class_vec) {
    if (length(class_vec) == 0) {
      return(NA_real_)
    }
    props <- table(class_vec) / length(class_vec)
    if (length(props) == 0) {
      return(NA_real_)
    }
    -sum(props * log2(props + .Machine$double.xmin))
  }
  
  hex_summary <- switch(
    plot_type,
    "majority" = {
      hb_umap_df %>%
        dplyr::group_by(hex_x, hex_y) %>%
        dplyr::summarize(
          cell_count = dplyr::n(),
          majority = safe_majority(class),
          prop_majority = if (length(class) > 0) max(table(class)) / dplyr::n() else NA_real_,
          .groups = "drop"
        )
    },
    "entropy" = {
      hb_umap_df %>%
        dplyr::group_by(hex_x, hex_y) %>%
        dplyr::summarize(
          cell_count = dplyr::n(),
          majority = safe_majority(class),
          entropy = safe_entropy(class),
          .groups = "drop"
        )
    },
    "gene" = {
      hb_umap_df %>%
        dplyr::group_by(hex_x, hex_y) %>%
        dplyr::summarize(
          cell_count = dplyr::n(),
          majority = safe_majority(class),
          mean_expr = if (length(gene_expr) > 0) mean(gene_expr) else NA_real_,
          .groups = "drop"
        )
    }
  ) %>%
    dplyr::filter(cell_count >= min_cells) %>%
    dplyr::filter(!is.na(majority))  # Remove hexagons with no valid majority class
  
  p <- switch(
    plot_type,
    "majority" = {
      ggplot(hex_summary, aes(x = hex_x, y = hex_y)) +
        geom_hex(aes(fill = majority, alpha = prop_majority), stat = "identity") +
        get_color_scale(color_palette, custom_colors) +
        scale_alpha_continuous(range = c(alpha, 1)) +
        labs(alpha = "Proportion\nMajority", fill = "Cell Type")
    },
    "entropy" = {
      ggplot(hex_summary, aes(x = hex_x, y = hex_y)) +
        geom_hex(aes(fill = entropy), stat = "identity") +
        scale_fill_viridis_c(option = viridis) +
        labs(fill = "Shannon\nEntropy")
    },
    "gene" = {
      ggplot(hex_summary, aes(x = hex_x, y = hex_y)) +
        geom_hex(aes(fill = majority, alpha = mean_expr), stat = "identity") +
        get_color_scale(color_palette, custom_colors) +
        scale_alpha_continuous(range = c(alpha, 1)) +
        labs(alpha = paste0(gene_name, "\nExpression"), fill = "Cell Type")
    }
  ) +
    coord_equal() +
    theme_minimal() +
    theme(panel.grid = element_blank()) +
    labs(title = paste0("HexMAP - ", plot_type, " (Seurat ", as.character(seurat_version), ")"),
         x = paste0(reduction, "_1"),
         y = paste0(reduction, "_2"))
  
  return(p)
}

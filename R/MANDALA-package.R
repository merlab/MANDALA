#' @keywords internal
"_PACKAGE"

#' @import ggplot2
#' @importFrom grDevices col2rgb colorRampPalette rgb adjustcolor contourLines
#' @importFrom stats setNames complete.cases quantile median dist
#' @importFrom utils packageVersion
#' @importFrom dplyr %>%
#' @importFrom hexbin hexbin
NULL

utils::globalVariables(c(
  "Dim1","Dim2","grp","fillcol","size_val","expr_grp","metric_value",
  "hex_x","hex_y","cell_count","majority","prop_majority","entropy",
  "mean_expr","bin_x","bin_y","DimPlot","x","y"
))

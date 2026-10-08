#' Pan-cancer tumour microenvironment atlas (subset)
#'
#' Marker genes across cells from a pan-cancer atlas, with a UMAP reduction and
#' cell type labels in \code{cell_type}. Used by \code{\link{BubbleMAP}} and
#' \code{\link{BubbleMAP3D}}.
#'
#' The counts layer holds log-normalized values despite its name, so
#' \code{layer = "counts"} is the correct choice for plotting.
#'
#' @format A \link[SeuratObject]{Seurat} object with one \code{umap} reduction.
#' @source Rapozo Guimaraes et al. (2024) Nat Commun 15:5694, via CZ CELLxGENE.
"sample_scRNA"

#' Bone marrow Visium slide (subset)
#'
#' COL1A1, HBB, MPO and NKG7 across 1,714 spots, with RCTD deconvolution labels,
#' ssGSEA scores and distance-to-landmark columns in \code{meta.data}, plus the
#' \code{slice1} image. Used by \code{\link{SpatialMAP3D}}.
#'
#' The ssGSEA scores and \code{RCTD_*} proportions can be passed to
#' \code{SpatialMAP3D(metric = ...)} to put them on the height axis.
#'
#' @format A \link[SeuratObject]{Seurat} object with an SCT assay and one image.
"sample_spatial"

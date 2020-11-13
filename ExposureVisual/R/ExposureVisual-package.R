#' ExposureVisual - Heatmap visualization of mutational signatures exposure
#' of tumors genomes from MPF file
#'
#' The ExposureVisual package uses \code{decomposeTumorGenomes()} function of
#' \code{decompTumor2Sig} to determine the exposures reflecting the
#' contributions of a set of given mutational signatures (either of the
#' "Alexandrov model" by Alexandrov et al, Nature 500(7463):415-421,2013, or the
#' "Shiraishi model" by Shiraishi et al, PLoS Genet 11(12):e1005657, 2015) to
#' the mutation load of the individual tumor genomes, obtained using
#' \code{readGenomesFromMPF()} function of \code{\link{decompTumor2Sig}},
#' from a MPF file.
#' The resulting exposures are then hierarchically clustered based on
#' correlation distance between values and complete linkage between clusters,
#' scaled and visualized as a heatmap.
#'
#'
#' \tabular{ll}{
#' Package: \tab ExposureVisual\cr
#' Type: \tab Package\cr
#' Version: \tab 0.99.0\cr
#' Date: \tab 2020-08-23\cr
#' License: \tab GPL-2 \cr
#' }
#'
#' The package provides the following function:
#'
#' \tabular{ll}{
#'
#' ExposureHeatVisual():\tab read a set of genomes from a \cr
#'                      \tab Mutation Position Format (MPF) file,\cr
#'                      \tab determine the exposures of a set of\cr
#'                      \tab signatures to each individual\cr
#'                      \tab tumor genome, cluster and visualize\cr
#'                      \tab them with an heatmap.\cr
#' }
#'
#' @name ExposureVisual-package
#' @aliases ExposureVisual-package ExposureVisual
#' @docType package
#' @import decompTumor2Sig
#' @import pheatmap
#' @author Veronica Erconi, Politecnico di Milano\cr
#' Maintainer: Veronica Erconi\cr
#' E-Mail: <erconi.veronica@@gmail.com>
NULL

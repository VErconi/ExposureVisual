% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExposureHeatVisual.R
\name{ExposureHeatVisual}
\alias{ExposureHeatVisual}
\title{Heatmap visualization of the exposures of signatures in tumor genomes.}
\usage{
ExposureHeatVisual (MPFfile,  type="Alexandrov", numBases=3,
refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
transcriptAnno=
TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
trDir=FALSE, enforceUniqueTrDir=TRUE, signatures,
heat_main = "Exposure Heatmap",
cluster_rows=TRUE, cluster_cols=FALSE, cellwidth= 15, cellheight=15)
}
\arguments{
\item{MPFfile}{(Mandatory) The name of the MPF file which can also
be compressed with \code{gzip}.}

\item{type}{(Mandatory) Signature model or type (\code{"Alexandrov"} or
\code{"Shiraishi"}). Default: \code{"Alexandrov"}}

\item{numBases}{(Mandatory) Total number of bases (mutated base and
flanking bases) to be used for sequence patterns. Must be odd. Default: 3}

\item{refGenome}{(Mandatory) The reference genome (\code{BSgenome}) needed
to extract sequence patterns. Default: \code{BSgenome} object for hg19.}

\item{transcriptAnno}{(Optional) Transcript annotation (\code{TxDb} object)
used to determine the transcription direction. This is required only if
\code{trDir} is \code{TRUE}. Default: \code{TxDb} object for hg19.}

\item{trDir}{(Mandatory) Specifies whether the transcription direction is
taken into account in the signature model. If so, only mutations within
genomic regions with a defined transcription direction can be considered.
Default: \code{FALSE}}

\item{enforceUniqueTrDir}{(Optional) Used only if \code{trDir} is
\code{TRUE}. If \code{enforceUniqueTrDir} is TRUE (default), then mutations
which map to a region with multiple overlapping genes with opposing
transcription directions will be excluded from the analysis. If \code{FALSE},
the transcript direction encountered first in the transcript database (see
\code{transcriptAnno}) is assigned to the mutation. The latter was the
behavior until version 1.3.5 of \code{decompTumor2Sig} and is also the
behavior of \code{pmsignature}. However, it is preferable to exclude
these mutations from the count (default) because from mutation data alone
it cannot be inferred which of the two genes has the higher transcriptional
activity which might potentially be linked to the occurrence of the mutation.
(If you are unsure, use the default setting; this option exists mostly for
backward compatibility with older versions.)}

\item{signatures}{(Mandatory) A list of vectors, data frames or matrices.
Each of the objects represents one mutational signature. Vectors are
used for Alexandrov signatures, data frames or matrices for Shiraishi
signatures.}

\item{heat_main}{(Optional) A string representing the title of the heatmap,
if \code{heat_main} is not specified "Exposure Heatmap" (default)
will be shown as title.}

\item{cluster_rows}{(Optional) boolean value determining if rows should be
clustered or \code{hclust} object. If not specified, \code{TRUE} (default),
rows will be clustered.
Check \url{https://cran.r-project.org/web/packages/pheatmap/pheatmap.pdf}
and
\url{https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust}
for further information.}

\item{cluster_cols}{(Optional) similar to \code{cluster_rows} but for columns.
If not specified, columns will not be clustered, \code{FALSE} (default)}

\item{cellwidth}{(Optional) individual cell width in points.
If NA, then the values depend on the size of plotting window, if not specified,
the width will be of 15 (default).}

\item{cellheight}{(Optional) similar to \code{cellwidth} but for cell height.}
}
\value{
Exposure Heatmap
}
\description{
\code{ExposureHeatVisual()}allows to visualize as a heatmap the exposures,
determined with \code{decomposeTumorGenomes()} function of
\code{\link{decompTumor2Sig}}, of a set of signatures
(Alexandrov or Shiraishi-type) to each of a set of individual tumor genomes,
obtained using \code{readGenomesFromMPF()} function of
\code{\link{decompTumor2Sig}},from a MPF (mutation position format) file.
}
\seealso{
\url{https://www.bioconductor.org/packages/devel/bioc/html/decompTumor2Sig.html}\cr
\url{https://CRAN.R-project.org/package=pheatmap}\cr
}
\author{
Veronica Erconi\cr Bioinformatics for Computational Genomics -
Politecnico di Milano\cr Maintainer: Veronica Erconi\cr
E.mail: \href{mailto:erconi.veronica@gmail.com}{erconi.veronica@gmail.com}
}

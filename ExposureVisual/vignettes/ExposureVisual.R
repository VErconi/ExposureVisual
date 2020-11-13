## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----style, echo = FALSE, results = 'asis'------------------------------------
library(BiocStyle)

## ---- echo = FALSE------------------------------------------------------------
library(knitr)

## ----setup--------------------------------------------------------------------
library(ExposureVisual)

## ---- eval=FALSE--------------------------------------------------------------
#  ExposureHeatVisual (MPFfile,  type="Alexandrov", numBases=3,
#                  refGenome=BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19,
#                  transcriptAnno=
#                  TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
#                  trDir=FALSE, enforceUniqueTrDir=TRUE, signatures,
#                  heat_main = "Exposure Heatmap",
#                  cluster_rows=TRUE, cluster_cols=FALSE, cellwidth= 15, cellheight=15)
#  


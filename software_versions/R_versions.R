# moule load conda_R/4.0.x

library(VariantAnnotation)
library(tidyverse)
library(GenomicRanges)
library(httr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(data.table)
library(GenomicAlignments)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(rtracklayer)
library(foreach)
library(doParallel)
library(SummarizedExperiment)

writeLines(capture.output(sessionInfo()), 
           "./R_versions.txt")


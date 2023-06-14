# Getting the Duke blacklisted regions and storing as a GRanges

library(tidyverse)
library(GenomicRanges)
library(httr)

outDir <- "../outDir/duke"

dir.create(outDir, showWarnings = F)

blacklisted.file <- httr::content(GET("http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDukeMapabilityRegionsExcludable.bed.gz"))
blacklisted.tib <- read_tsv(gzcon(rawConnection(blacklisted.file)),
                            col_names=c("seqnames", "start", "end", "name", "score", "z")) %>%
                     select(-z)
blacklisted.tib <- blacklisted.tib %>% mutate(start=start+1)
blacklist.gr <- makeGRangesFromDataFrame(blacklisted.tib,
                                         keep.extra.columns=TRUE)
blacklist.gr <- keepSeqlevels(blacklist.gr, paste0("chr", c(1:22, "X", "Y")),
                              pruning.mode="coarse")

saveRDS(blacklist.gr, file = file.path(outDir, "duke-blacklist-hg19.rds"))

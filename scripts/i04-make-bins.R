# Makes 100kb bins and groups them for merging into larger 2.5Mb bins

library(tidyverse)
library(BSgenome.Hsapiens.UCSC.hg19)

dir.create("../outDir/bins", showWarnings = F)

bin.size <- 1e5
bins_100kb <- tileGenome(seqlengths = seqlengths(BSgenome.Hsapiens.UCSC.hg19)[paste0("chr", 1:22)], tilewidth = bin.size, cut.last.tile.in.chrom = TRUE)
names(bins_100kb) <- paste0("bin", 1:length(bins_100kb))

saveRDS(bins_100kb, file = "../outDir/bins/bins_100kb.rds")

# Identifying the 100kb bins that will be merged to create 2.5Mb bins
bins_100kb.by.chr <- split(bins_100kb, seqnames(bins_100kb))

n.merge <- 25
merge.grps <- NULL
for (i in 1:length(bins_100kb.by.chr)) {
  tmp.bins <- bins_100kb.by.chr[[i]]
  if (is.null(merge.grps)) {
    bin.grp <- rep(1:ceiling(length(tmp.bins)/n.merge), each = n.merge) %>% head(length(tmp.bins))
  } else {
    st <- max(merge.grps, na.rm = T)
    bin.grp <- rep((1 + st):(ceiling(length(tmp.bins)/n.merge) + st), each = n.merge) %>% head(length(tmp.bins))
  }
  bins.in.grp <- table(bin.grp)
  full.length.bins <- names(bins.in.grp)[bins.in.grp == n.merge] %>% as.numeric
  not.full.length.bins <- setdiff(bin.grp, full.length.bins)
  if (length(not.full.length.bins) > 0) {
    bin.grp[bin.grp == not.full.length.bins] <- NA
  }
  merge.grps <- c(merge.grps, bin.grp)
}

saveRDS(merge.grps, file = "../outDir/bins/merge_grps_100kb_to_2.5mb_bins.rds")

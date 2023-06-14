# Makes a GRanges objects with SBS variants from gnomad-v3

library(VariantAnnotation)

args <- commandArgs(trailingOnly = T)
tmpDir <- args[1]
outDir <- args[2]

#tf <- TabixFile(file.path(tmpDir, "gnomad.genomes.r3.0.sites.vcf.bgz"), yieldSize = 50000000)
tf <- TabixFile(file.path(tmpDir, "gnomad.genomes.r3.0.sites.vcf.gz"), yieldSize = 50000000)

i <- 1
open(tf)
while(length(gnomad <- readVcf(tf, param = ScanVcfParam(info = c("AF", "variant_type")), row.names = FALSE))) {
  rowRanges(gnomad)$AF <- unlist(info(gnomad)$AF)
  rowRanges(gnomad)$variant_type <- info(gnomad)$variant_type
  gnomad <- rowRanges(gnomad)
  gnomad$REF <- as.character(gnomad$REF)
  gnomad$ALT <- as.character(unlist(gnomad$ALT))
  gnomad <- gnomad[gnomad$variant_type == "snv"]
  gnomad <- gnomad[,-c(1, 3, 6)]
  gnomad$FILTER <- ifelse(gnomad$FILTER == "PASS", "P", "F")
  mcols(gnomad) <- mcols(gnomad)[, c(2, 3, 1, 4)]
  
  saveRDS(gnomad, file = file.path(tmpDir, paste0("gnomad_", i, ".sbs.rds")))
  i <- i + 1 
}
close(tf)


# Combining the individual GRanges objects into a single GRanges
gnomad <- GRanges()
for (i in 1:15) {
  tmp <- readRDS(file.path(tmpDir, paste0("gnomad_", i, ".sbs.rds")))
  gnomad <- c(gnomad, tmp)
}

# Saving each chromosome separately
chrs <- unique(as.character(seqnames(gnomad)))
for (i in 1:length(chrs)) {
  saveRDS(gnomad[seqnames(gnomad) == chrs[i]],
          file = file.path(outDir, paste0("gnomad.genomes.r3.0.sites.sbs.", chrs[i], ".rds")))
}

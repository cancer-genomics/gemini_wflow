# Note: bam file is assumed to be sorted and indexed with duplicate fragments marked or removed

args <- commandArgs(trailingOnly = TRUE)
bamFile <- args[1]
outDir <- args[2] 
tmpDir <- args[3]
nProcesses <- as.numeric(args[4])

library(dplyr)
library(tidyr)
library(readr)
library(tibble)
library(data.table)
library(GenomicAlignments)
library(BSgenome.Hsapiens.UCSC.hg19)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(VariantAnnotation)
library(rtracklayer)
library(foreach)
library(doParallel)

# Setting up the parallel computing environment
cl <- makeCluster(nProcesses)
registerDoParallel(cl)

if (!dir.exists(file.path(outDir, "variants"))) dir.create(file.path(outDir, "variants"))
if (!dir.exists(file.path(outDir, "cts"))) dir.create(file.path(outDir, "cts"))

# Sourcing helper functions
source("../functions/sbs_counter.R")

# Defining the path to the python file for computing the pileups
pyFile <- "../functions/pileup_runner.py"

# Defining the sample name
sampleName <- gsub(".bam$", "", basename(bamFile))

# Identifying the reference genome used for alignment
bamGenome <- ifelse ("1" %in% seqnames(seqinfo(BamFile(bamFile))), "GRCh37", "hg19")

# Reading in 100kb bins to generate position/variant counts in
bins <- readRDS("../outDir/bins/bins_100kb.rds")

# Initializing counts
nms <- c("pyR1", "pyR1_known_somatic", "pyR1_known_germline", "puR1", "puR1_known_somatic", "puR1_known_germline")
ct.names <- c("cg.pyR1", "cg.puR1", paste0("cg2at.", nms), paste0("cg2gc.", nms), paste0("cg2ta.", nms),
              "ta.pyR1", "ta.puR1", paste0("ta2at.", nms), paste0("ta2cg.", nms), paste0("ta2gc.", nms))

mcols(bins) <- matrix(ncol = length(ct.names), nrow = length(bins), data = 0, dimnames = list(NULL, ct.names))

# Defining chunks of genome to process together, each chunk is adjaent bins from the same chromosome.
# Reducing nParallelBins will save memory but increase runtime which may be necessary for >1-2x coverage WGS
# data (depending on the amount of memory available).
nParallelBins <- 100
bins.by.chr <- bins %>% split(., seqnames(.))
st <- 1
bin.grp <- NULL
for (i in 1:length(bins.by.chr)) {
 n.bins <- length(bins.by.chr[[i]])
 chr.grp <- rep(st:(st + ceiling(n.bins / nParallelBins)), each = nParallelBins) %>% head(., n.bins)
 bin.grp <- c(bin.grp, chr.grp)
 st <- max(bin.grp) + 1
}
bins <- split(bins, bin.grp)

# Defining the paths to gnomAD variants for each chromosome separately
gnomadDir <- "../outDir/gnomad"
chrs <- paste0("chr", c(1:22, "X", "Y"))
gnomadVariants <- file.path(gnomadDir, paste0("gnomad.genomes.r3.0.sites.sbs.", chrs, ".rds"))
names(gnomadVariants) <- chrs

# Defining the path to the hg19 to hg38 chain for lifting over variants
hg19ToHg38ChainPath <- "../outDir/chain/hg19ToHg38.over.chain"
hg19ToHg38Chain <- rtracklayer::import.chain(hg19ToHg38ChainPath)

# Defining a GRanges object that contains blacklisted regions (e.g. Duke blacklist)
blacklistedRegions <- readRDS("../outDir/duke/duke-blacklist-hg19.rds")

# Defining a GRanges object that contains known germline variants from matched normal sequencing in the sample being analyzed
knownGermline <- GRanges()

# Defining a GRanges object that contains known somatic variants from tumor sequencing in the sample being analyzed
knownSomatic <- GRanges()

# Defining the minimum Phred score for a base to be considered
minPhred <- 30

# Defining the minimum MAPQ score for a base to be considered
minMAPQ <- 40

# Defining the set of bins for each process to process
bins <- split(bins, cut(seq_along(bins), nProcesses, labels = FALSE))

ct.gr <- foreach(i = 1:nProcesses, .combine = "c", .packages = c("dplyr", "tidyr", "readr", "tibble", "data.table", "GenomicRanges", "BSgenome.Hsapiens.UCSC.hg19", "BSgenome.Hsapiens.UCSC.hg38", "Biostrings", "rtracklayer")) %dopar% {

  setDTthreads(threads = 1)

  # Defining a file path to write variants. 
  # If the file already exists (e.g. from a failed run) then remove it
  variantFile <- file.path(outDir, "variants", paste0(sampleName, "_variants_process", i, ".txt"))
  if (file.exists(variantFile)) file.remove(variantFile)

  # Defining a file path to write the pysam pileup to
  # # If the file already exists (e.g. from a failed run) then remove it
  pileupFile <- file.path(tmpDir, paste0(sampleName, "_pileup_process", i, ".txt"))
  if (file.exists(pileupFile)) file.remove(pileupFile)

  # Initializing a GRanges to store gnomad variants
  gnomad <- GRanges()

  # Getting the regions that will be processed by process N
  binsN <- bins[[i]]

  for (j in 1:length(binsN)) {

    # Reducing binsN[[j]] to a single chunk
    chunk <- reduce(binsN[[j]])
    stopifnot(length(chunk) == 1)

    region <- as.character(chunk)
    if (bamGenome == "GRCh37") {
      region <- gsub("^chr", "", region)
    } 
  
    # Pileup using pysam
    pyCmd <- paste("python3", pyFile, bamFile, region, minPhred, minMAPQ, pileupFile)
    system(pyCmd)

    cts <- count_sbs(pileupFile = pileupFile,
                     bamGenome = bamGenome,
                     region.gr = binsN[[j]],
                     bothMates = TRUE,
                     gnomadFilter = TRUE,
                     gnomadVariants = gnomadVariants,
                     hg19ToHg38Chain = hg19ToHg38Chain,
                     blacklistedRegions = blacklistedRegions,
                     knownGermline = knownGermline,
                     knownSomatic = knownSomatic,
                     writeVariants = TRUE,
                     variantFile = variantFile)

    mcols(binsN[[j]]) <- cts

  } # End of j loop

  # Removing temporary files
  if (file.exists(pileupFile)) invisible(file.remove(pileupFile))

  return(unlist(binsN, use.names = F))

} # End of i loop

# Saving the number of evaluable/variant positions to file
saveRDS(ct.gr, file = file.path(outDir, "cts", paste0(sampleName, "_cts.rds")))




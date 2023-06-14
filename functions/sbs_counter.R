#-------------------------------------------------------------
# Defining a function to annotate variants with gnomAD-v3
#--------------------------------------------------------------------------------------------------------------------------------------------
gnomADannotate <- function(chrom, pos, ref, alt, gnomadVariants, hg19ToHg38Chain) {

  # Creating a data.frame to store results
  out.df <- data.table(liftover.success = rep("Y", length(chrom)),
                       in.gnomad = rep("N", length(chrom)),
                       gnomad.pass = as.character(rep(NA, length(chrom))),
                       gnomad.af = as.numeric(rep(NA, length(chrom))), stringsAsFactors = FALSE)

  # Storing the variants as a GRanges 'mut.gr'
  mut.gr <- GRanges(seqnames = chrom, ranges = IRanges(start = pos, end = pos))

  # Lifting over variants from hg19 to hg38 to match gnomAD
  hg38.mut.grl <- rtracklayer::liftOver(mut.gr, hg19ToHg38Chain) # Note that some positions will not lift over (will mark as N)
  names(hg38.mut.grl) <- 1:nrow(out.df) # The name will be the index of the variant passed to this function
  n.liftover <- elementNROWS(hg38.mut.grl)

  no.liftover <- which(n.liftover == 0)
  if (length(no.liftover) > 0) {
    out.df$liftover.success[no.liftover] <- "N" # Marking variants that don't lift over to hg38
  }

  many.liftover <- which(n.liftover > 1)
  if (length(many.liftover) > 0) {
    out.df$liftover.success[many.liftover] <- "M" # Marking variants that lift over multiple times to hg38
    hg38.mut.grl <- hg38.mut.grl[-many.liftover]  # and excluding them from the gnomAD annotation stage
  }

  hg38.mut.gr <- unlist(hg38.mut.grl)
  if (length(hg38.mut.gr) == 0) return(out.df)
  rm(list = c("hg38.mut.grl", "n.liftover", "no.liftover", "many.liftover"))

  # Marking variants if they lifted over to one position but the reference base is different between the two builds
  hg38.mut.gr$REF <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, hg38.mut.gr))
  different.ref <- hg38.mut.gr$REF != ref[as.numeric(names(hg38.mut.gr))]
  if (any(different.ref)) {
    different.ref.ind <- which(different.ref)
    out.df$liftover.success[as.numeric(names(hg38.mut.gr))[different.ref.ind]] <- "R"
    hg38.mut.gr <- hg38.mut.gr[-different.ref.ind]
  }
  if (length(hg38.mut.gr) == 0) return(out.df)
  rm(list = c("different.ref"))

  hg38.mut.gr$ALT <- alt[as.numeric(names(hg38.mut.gr))]

  # Getting gnomAD variant positions for any chromosomes in chrom
  new_chrom <- unique(as.character(seqnames(hg38.mut.gr)))
  gnomad <<- gnomad[(seqnames(gnomad) %in% new_chrom)] # Only keeping chrommosomes that are still needed
  for (i in 1:length(new_chrom)) { # If a new chromosome is encountered, add it to gnomad
    if (!(new_chrom[i] %in% as.character(seqnames(gnomad)))) {
      gnomad <<- c(gnomad, readRDS(gnomadVariants[new_chrom[i]]))
    }
  }

  olaps <- findOverlaps(hg38.mut.gr, gnomad, type = "equal")
  hg38.mut.gr <- hg38.mut.gr[queryHits(olaps)]
  tmp.gnomad <- gnomad[subjectHits(olaps)]

  is.same.alt <- hg38.mut.gr$ALT == tmp.gnomad$ALT

  if (any(is.same.alt)) {
    hits <- as.numeric(names(hg38.mut.gr))[is.same.alt]
    out.df$in.gnomad[hits] <- "Y"
    out.df$gnomad.pass[hits] <- ifelse(tmp.gnomad$FILTER[is.same.alt] == "P", "Y", "N")
    out.df$gnomad.af[hits] <- tmp.gnomad$AF[is.same.alt]
  }

  return(out.df)

} # End of gnomADannotate
#--------------------------------------------------------------------------------------------------------------------------------------------


#----------------------------------------------------------------------------------
# Defining a function to call C:G and T:A evaluable positions and all SBS types.
#-----------------------------------------------------------------------------------------------------------------------
count_sbs <- function(pileupFile, bamGenome, region.gr, bothMates, gnomadFilter, gnomadVariants, hg19ToHg38Chain,
                      blacklistedRegions, knownGermline, knownSomatic, writeVariants, variantFile) {

  # Reading in the pileup file
  # Change to use fread if memory leak issue becomes fixed
  eval.df <- read_tsv(pileupFile, col_types = cols(), progress = FALSE)
  eval.df <- data.table(eval.df)
  if (nrow(eval.df) == 0) return(mcols(region.gr))

  # Converting the fragment identifier to an integer
  eval.df[, fragment := as.integer(factor(fragment))]

  # If the bamGenome is GRCh37 adding the chr prefix to chromosome names
  if (bamGenome == "GRCh37") {
    eval.df$chrom <- paste0("chr", eval.df$chrom)
  }

  # Getting the reference base for each position in eval.df
  pos.df <- eval.df[, .(chrom, pos)]
  pos.df <- pos.df[!duplicated(pos.df)]
  pos.df$ref <- getSeq(BSgenome.Hsapiens.UCSC.hg19, names = pos.df$chrom, start = pos.df$pos, end = pos.df$pos, as.character = TRUE)
  eval.df <- eval.df[pos.df, on = .(chrom, pos)]
  rm(list = c("pos.df"))

  # Removing any potential positions where the reference base is an 'N'
  eval.df <- eval.df[ref %in% c("A", "T", "C", "G"),]
  if (nrow(eval.df) == 0) return(mcols(region.gr))

  # Only keeping evaluable positions when the adjacent 2 bases in the read are also evaluable and are the reference allele
  tmp.eval.df <- eval.df[ref == obs]
  eval.df[, pos := pos - 1]
  eval.df[, is.left.eval := 0][tmp.eval.df, is.left.eval := 1, on = .(fragment, pos, read)]
  eval.df[, pos := pos + 2]
  eval.df[, is.right.eval := 0][tmp.eval.df, is.right.eval := 1, on = .(fragment, pos, read)]
  eval.df[, pos := pos - 1]
  eval.df <- eval.df[is.left.eval == 1 & is.right.eval == 1]

  eval.df <- eval.df[, .(fragment, read, strand, chrom, pos, ref, obs)]
  if (nrow(eval.df) == 0) return(mcols(region.gr))
  rm(list = c("tmp.eval.df"))

  # Filtering for positions in eval.df that overlap region.gr (this should be all or almost all of them)
  eval.gr <- GRanges(eval.df$chrom, IRanges(eval.df$pos, eval.df$pos))
  eval.df <- eval.df[overlapsAny(eval.gr, region.gr)]
  if (nrow(eval.df) == 0) return(mcols(region.gr))
  rm(list = c("eval.gr"))

  if (bothMates == TRUE) {
    eval.df <- eval.df[duplicated(eval.df[, .(fragment, chrom, pos, ref, obs)]), ] # Keeps the 2nd entry when there are two entries that have the same observed base and no entries if there is one entry
  } else if (bothMates == FALSE) {
    # Keeping mutations and evaluable positions that are only in one of the two reads or if in both of the reads selecting one of the reads to keep at random
    # Randomizing the rows in eval.df so that read 1 and read 2 are intermixed
    set.seed(0)
    eval.df <- eval.df[sample(1:nrow(eval.df), replace = FALSE),]
    eval.df <- eval.df[!duplicated(eval.df[,.(fragment, pos)]), ] # Keeps the first entry when there are two entries
  }

  if (nrow(eval.df) == 0) return(mcols(region.gr))
  #--------------------------------------------------------------------------------------------------------------------------------

  # Marking evaluable positions as variants or not
  eval.df[, is.variant := ifelse(obs == ref, 0, 1)]

  # Removing variants in the gnomAD database with an AF > 1/100,000
  if (gnomadFilter == TRUE) {
    eval.df$liftover.success <- as.character(NA)
    eval.df$in.gnomad <- as.character(NA)
    eval.df$gnomad.pass <- as.character(NA)
    eval.df$gnomad.af <- as.numeric(NA)
    eval.df.variant.ind <- which(eval.df$is.variant == 1)
    if (length(eval.df.variant.ind) > 0) {
      m <- gnomADannotate(chrom = eval.df$chrom[eval.df.variant.ind],
                          pos = eval.df$pos[eval.df.variant.ind],
                          ref = eval.df$ref[eval.df.variant.ind],
                          alt = eval.df$obs[eval.df.variant.ind],
                          gnomadVariants = gnomadVariants,
                          hg19ToHg38Chain = hg19ToHg38Chain)
      eval.df[eval.df.variant.ind, c("liftover.success", "in.gnomad", "gnomad.pass", "gnomad.af")] <- m
      rm(list = c("m"))
    }

    # Removing variant positions that: did not lift over to hg38, or that lifted over to hg38 but the
    # reference base was different in hg38 compared to hg19, or the variant didn't pass gnomAD filters, 
    # or the gnomAD allele frequency is > 1/100,000  
    eval.df <- eval.df[(is.na(liftover.success) | liftover.success == "Y") &
                             (is.na(gnomad.pass) | gnomad.pass == "Y") &
                             (is.na(gnomad.af) | gnomad.af <= 1/1e5)]

    eval.df <- eval.df[,1:8]
  }
  if (nrow(eval.df) == 0) return(mcols(region.gr))


  # Removing any evaluable positions that are in a blacklisted region
  if (length(blacklistedRegions) > 0) {
    eval.gr <- GRanges(eval.df$chrom, IRanges(eval.df$pos, eval.df$pos))
    blacklistedRegions.olaps <- findOverlaps(eval.gr, blacklistedRegions)
    if (length(blacklistedRegions.olaps) > 0) {
      eval.df <- eval.df[-unique(queryHits(blacklistedRegions.olaps)),]
    }
  }
  if (nrow(eval.df) == 0) return(mcols(region.gr))


  # Adding a column for whether a pyrimidine (y) or purine (u) is in R1
  eval.df[, r1.base := "y"]
  pur1.ind <- with(eval.df, which((read == 1 & ((ref %in% c("A", "G") & strand == "+") | (ref %in%  c("C", "T") & strand == "-"))) |
                                  (read == 2 & ((ref %in% c("A", "G") & strand == "-") | (ref %in% c("C", "T") & strand == "+")))))
  eval.df[pur1.ind, r1.base := "u"]

  # Assigning the counting group to each evaluable position / variant
  eval.df[, eval_grp := ifelse(ref %in% c("C", "G"), "cg", "ta")]
  eval.df[(ref == "C" & obs == "A") | (ref == "G" & obs == "T"), mut_grp := "cg2at"]
  eval.df[(ref == "C" & obs == "G") | (ref == "G" & obs == "C"), mut_grp := "cg2gc"]
  eval.df[(ref == "C" & obs == "T") | (ref == "G" & obs == "A"), mut_grp := "cg2ta"]
  eval.df[(ref == "T" & obs == "A") | (ref == "A" & obs == "T"), mut_grp := "ta2at"]
  eval.df[(ref == "T" & obs == "C") | (ref == "A" & obs == "G"), mut_grp := "ta2cg"]
  eval.df[(ref == "T" & obs == "G") | (ref == "A" & obs == "C"), mut_grp := "ta2gc"] 

  # Marking variants that are known germline or known somatic
  eval.df[, c("known_germline","known_somatic") := 0]
  variant.ind <- which(eval.df$is.variant == 1)
  if (length(variant.ind) > 0) {
    variant.gr <- with(eval.df[variant.ind,], GRanges(chrom, IRanges(pos, pos), ref = ref, mut = obs))
    variant.df <- as.data.table(variant.gr)

    # Marking known germline variants
    if (length(knownGermline) > 0) {
      variant.df[, known_germline := 0][as.data.table(knownGermline), known_germline := 1, on = .(seqnames, start, end, ref, mut)]
      eval.df[variant.ind, "known_germline"] <- variant.df[, "known_germline"]
    }

    # Marking known somatic variants
    if (length(knownSomatic) > 0) {
      variant.df[, known_somatic := 0][as.data.table(knownSomatic), known_somatic := 1, on = .(seqnames, start, end, ref, mut)]
      eval.df[variant.ind, "known_somatic"] <- variant.df[, "known_somatic"]
    }
  }

  # Getting the bins that each evaluable position overlaps
  eval.gr <- GRanges(eval.df$chrom, IRanges(eval.df$pos, eval.df$pos))
  eval.df$bin_id <- findOverlaps(eval.gr, region.gr) %>% subjectHits %>% names(region.gr)[.]

  # Recording the number of evaluable positions in each bin
  eval.bin.cts <- eval.df %>% 
                     group_by(bin_id, eval_grp, r1.base) %>% 
                     summarise(n = n(), .groups = "drop") %>%
                     mutate(tally_grp = paste0(eval_grp, ifelse(r1.base == "y", ".pyR1", ".puR1"))) %>%
                     dplyr::select(bin_id, tally_grp, n) %>%
                     pivot_wider(names_from = tally_grp, values_from = n) %>%
                     replace(is.na(.), 0) %>%
                     column_to_rownames(var = "bin_id") 
  if (nrow(eval.bin.cts) > 0) mcols(region.gr)[rownames(eval.bin.cts), colnames(eval.bin.cts)] <- eval.bin.cts

  # Recording the number of mutated positions in each bin
  mut.bin.cts <- eval.df %>%
                    filter(is.variant == 1) %>% 
                    group_by(bin_id, mut_grp, r1.base) %>%
                    summarise(n = n(), .groups = "drop") %>%
                    mutate(tally_grp = paste0(mut_grp, ifelse(r1.base == "y", ".pyR1", ".puR1"))) %>%
                    dplyr::select(bin_id, tally_grp, n) %>%
                    pivot_wider(names_from = tally_grp, values_from = n) %>%
                    replace(is.na(.), 0) %>%
                    column_to_rownames(var = "bin_id")
  if (nrow(mut.bin.cts) > 0) mcols(region.gr)[rownames(mut.bin.cts), colnames(mut.bin.cts)] <- mut.bin.cts

  # Recording the number of germline variant positions in each bin
  germline.bin.cts <- eval.df %>%
                        filter(known_germline == 1) %>%
                        group_by(bin_id, mut_grp, r1.base) %>%
                        summarise(n = n(), .groups = "drop") %>%
                        mutate(tally_grp = paste0(mut_grp, ifelse(r1.base == "y", ".pyR1", ".puR1"), "_known_germline")) %>%
                        dplyr::select(bin_id, tally_grp, n) %>%
                        pivot_wider(names_from = tally_grp, values_from = n) %>%
                        replace(is.na(.), 0) %>%
                        column_to_rownames(var = "bin_id")
  if (nrow(germline.bin.cts) > 0) mcols(region.gr)[rownames(germline.bin.cts), colnames(germline.bin.cts)] <- germline.bin.cts

  somatic.bin.cts <- eval.df %>%
                        filter(known_somatic == 1) %>%
                        group_by(bin_id, mut_grp, r1.base) %>%
                        summarise(n = n(), .groups = "drop") %>%
                        mutate(tally_grp = paste0(mut_grp, ifelse(r1.base == "y", ".pyR1", ".puR1"), "_known_somatic")) %>%
                        dplyr::select(bin_id, tally_grp, n) %>%
                        pivot_wider(names_from = tally_grp, values_from = n) %>%
                        replace(is.na(.), 0) %>%
                        column_to_rownames(var = "bin_id")
  if (nrow(somatic.bin.cts) > 0) mcols(region.gr)[rownames(somatic.bin.cts), colnames(somatic.bin.cts)] <- somatic.bin.cts

  # Writing the coordinates of identified variants to file
  if (writeVariants == TRUE) {
    mut.df <- eval.df[is.variant == 1]
    if (nrow(mut.df) > 0) {
     mut.df <- mut.df[, .(chrom, pos, ref, obs)]
     if (!file.exists(variantFile)) {
        write.table(mut.df, file = variantFile,
                    row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
      } else {
        write.table(mut.df, file = variantFile, append = TRUE,
                    row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
      }
    }
    rm(list = c("mut.df"))
  }

  rm(list = c("eval.df", "pur1.ind", "variant.ind", "eval.gr", "eval.bin.cts", "mut.bin.cts", "germline.bin.cts", "somatic.bin.cts"))

  return(mcols(region.gr))

} # End of count_sbs
#-----------------------------------------------------------------------------------------------------------------------



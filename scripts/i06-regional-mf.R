library(tidyverse)
library(SummarizedExperiment)

# Reading in the 100kb bins and creating 2.5Mb bins
bins <- readRDS("../outDir/bins/bins_100kb.rds")
merge.grps <- readRDS("../outDir/bins/merge_grps_100kb_to_2.5mb_bins.rds")
bins <- split(bins, merge.grps) %>% reduce %>% unlist
names(bins) <- paste0("bin", 1:length(bins))

#-------------------------------------------------------------------------------
# Creating a RangedSumamrizedExperiment to hold bins, bin counts, and metadata
#-------------------------------------------------------------------------------------------
ctDir <- "../outDir/cts"
ctFiles <- list.files(ctDir)
ct.list <- lapply(ctFiles, function(t) readRDS(file.path(ctDir, t)))

for (i in 1:length(ct.list)) {
  split.meta <- split(data.frame(mcols(ct.list[[i]])), merge.grps)
  meta <- lapply(split.meta, function(t) colSums(t)) %>% do.call("rbind", .)
  tmp.bins <- bins
  mcols(tmp.bins) <- meta
  ct.list[[i]] <- tmp.bins
}

ct_grps <- colnames(mcols(ct.list[[1]])) %>% .[-grep("known", .)]
assay.list <- vector("list", length(ct_grps)) %>% setNames(., ct_grps)
for (i in 1:length(ct_grps)) {
  x <- lapply(ct.list, function(t) mcols(t)[, ct_grps[i]]) %>% do.call("cbind", .)
  dimnames(x) <- list(names(bins),  gsub(".sorted_cts.rds", "", ctFiles))
  assay.list[[ct_grps[i]]] <- x
}

# Reading in sample metadata
meta.df <- readRDS("../data/metadata.rds")

# Creating a RangedSummarizedExperiment to hold the data
se <- SummarizedExperiment(assays = assay.list,
                           rowRanges = bins,
                           colData = meta.df)
#-------------------------------------------------------------------------------------------


#--------------------------------------------------------------------
# Generating regional differences in mutation frequencies using LOO
#----------------------------------------------------------------------------------------------------------------------
# Defining the types of mutations to analyze
types <- list(cg2at = c("cg", "cg2at"), cg2gc = c("cg", "cg2gc"), cg2ta = c("cg", "cg2ta"),
              ta2at = c("ta", "ta2at"), ta2cg = c("ta", "ta2cg"), ta2gc = c("ta", "ta2gc"))

# Creating empty columns in 'se' for each mutation type
feature.names <- names(types)
m <- matrix(nrow = ncol(se), ncol = length(feature.names), dimnames = list(NULL, paste0("multimf_", feature.names)))
colData(se) <- cbind(colData(se), m)

# Creating an empty list to save bins sets
metadata(se) <- vector("list", ncol(se)) %>% setNames(., colnames(se))
bin.list <- list(pos.bins = NA, neg.bins = NA)
for (i in 1:length(metadata(se))) {
  metadata(se)[[i]] <- list(cg2at = bin.list, cg2gc = bin.list, cg2ta = bin.list, ta2at = bin.list, ta2cg = bin.list, ta2gc = bin.list)
}

for (i in 1:ncol(se)) {
  print(noquote(paste0(i, "/", ncol(se))))
  se.loo <- se[,-i]
  se.lo <- se[,i]

  for (j in 1:length(types)) {
    e.pyR1 <- paste0(types[[j]][1], ".pyR1")
    e.puR1 <- paste0(types[[j]][1], ".puR1")
    v.pyR1 <- paste0(types[[j]][2], ".pyR1")
    v.puR1 <- paste0(types[[j]][2], ".puR1")

    # For each bin, identifying the difference in mutation frequency between the grp1 and grp2
    grp1.ind <- which(se.loo$patient_type == "Case")
    grp2.ind <- which(se.loo$patient_type == "Control")

    if (names(types[j]) == "cg2at") {
      grp1.eval <- rowSums(assays(se.loo)[[e.pyR1]][,grp1.ind])
      grp2.eval <- rowSums(assays(se.loo)[[e.pyR1]][,grp2.ind])
      grp1.variant <- rowSums(assays(se.loo)[[v.pyR1]][,grp1.ind])
      grp2.variant <- rowSums(assays(se.loo)[[v.pyR1]][,grp2.ind])
    } else {
      grp1.eval <- rowSums(assays(se.loo)[[e.pyR1]][,grp1.ind]) + rowSums(assays(se.loo)[[e.puR1]][,grp1.ind])
      grp2.eval <- rowSums(assays(se.loo)[[e.pyR1]][,grp2.ind]) + rowSums(assays(se.loo)[[e.puR1]][,grp2.ind])
      grp1.variant <- rowSums(assays(se.loo)[[v.pyR1]][,grp1.ind]) + rowSums(assays(se.loo)[[v.puR1]][,grp1.ind])
      grp2.variant <- rowSums(assays(se.loo)[[v.pyR1]][,grp2.ind]) + rowSums(assays(se.loo)[[v.puR1]][,grp2.ind])
    }

    grp1.mpm <- grp1.variant / grp1.eval * 1e6
    grp2.mpm <- grp2.variant / grp2.eval * 1e6

    # Computing the difference in mutation frequency between grp1 and grp2
    rd <- grp1.mpm - grp2.mpm

    # Ordering the bins by difference in mutation frequency between grp1 and grp2
    ord.bins <- names(sort(rd))
    pos.bins <- tail(ord.bins, 114)
    neg.bins <- head(ord.bins, 114)

    # Computing the regional mutation frequency in the held out sample using the bins with the largest differential mutation density
    if (names(types[j]) == "cg2at") {
      pos.dens.lo <- sum(assays(se.lo)[[v.pyR1]][pos.bins,]) / sum(assays(se.lo)[[e.pyR1]][pos.bins,])
      neg.dens.lo <- sum(assays(se.lo)[[v.pyR1]][neg.bins,]) / sum(assays(se.lo)[[e.pyR1]][neg.bins,])
    } else {
      pos.dens.lo <- (sum(assays(se.lo)[[v.pyR1]][pos.bins,]) + sum(assays(se.lo)[[v.puR1]][pos.bins,])) / (sum(assays(se.lo)[[e.pyR1]][pos.bins,]) + sum(assays(se.lo)[[e.puR1]][pos.bins,]))
      neg.dens.lo <- (sum(assays(se.lo)[[v.pyR1]][neg.bins,]) + sum(assays(se.lo)[[v.puR1]][neg.bins,])) / (sum(assays(se.lo)[[e.pyR1]][neg.bins,]) + sum(assays(se.lo)[[e.puR1]][neg.bins,]))
    }

    multimf.lo <- (pos.dens.lo - neg.dens.lo) * 1e6

    feature.name <- paste0("multimf_", names(types[j]))

    # Saving the regional mutation frequency for the held out sample
    colData(se)[i,feature.name] <- multimf.lo

    # Saving the bin sets
    metadata(se)[[i]][[names(types[j])]]$pos.bins <- pos.bins
    metadata(se)[[i]][[names(types[j])]]$neg.bins <- neg.bins

  } # End of j loop

} # End of i loop
#----------------------------------------------------------------------------------------------------------------------


#-----------------------------------------------------
# Computing the GEMINI score [C>A] for each patient
#--------------------------------------------------------------------------------------------------------------------------------------------------
gemini_fit_cg2at <- glm(factor(patient_type, levels = c("Control", "Case")) ~ multimf_cg2at, data = colData(se), family = "binomial")
se$gemini_score_cg2at <- gemini_fit_cg2at$fitted.values
print(se$gemini_score_cg2at)
#--------------------------------------------------------------------------------------------------------------------------------------------------



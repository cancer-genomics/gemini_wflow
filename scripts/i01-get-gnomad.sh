#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=100G
#$ -l h_vmem=100G
#$ -l cancergen

module load conda_R/4.0.x

#--------
# Input
#-------------------------
tmpDir=../temp
outDir=../outDir/gnomad
#-------------------------

mkdir -p $tmpDir $outDir

# The gnomAD vcf file was initially available as follows:
# wget https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.vcf.bgz -O $tmpDir/gnomad.genomes.r3.0.sites.vcf.gz
# wget https://storage.googleapis.com/gnomad-public/release/3.0/vcf/genomes/gnomad.genomes.r3.0.sites.vcf.bgz.tbi -O $tmpDir/gnomad.genomes.r3.0.sites.vcf.gz.tbi

# They now appear to be hosted here (the md5sum values match those from the files above):
wget http://hgdownload-euro.soe.ucsc.edu/gbdb/hg38/gnomAD/vcf/gnomad.genomes.r3.0.sites.vcf.gz -O $tmpDir/gnomad.genomes.r3.0.sites.vcf.gz
wget http://hgdownload-euro.soe.ucsc.edu/gbdb/hg38/gnomAD/vcf/gnomad.genomes.r3.0.sites.vcf.gz.tbi -O $tmpDir/gnomad.genomes.r3.0.sites.vcf.gz.tbi

Rscript ./i01-get-gnomad.R $tmpDir $outDir

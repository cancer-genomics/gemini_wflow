#!/bin/bash

#$ -cwd
#$ -j y
#$ -l mem_free=12G
#$ -l h_vmem=12G
#$ -pe local 16
#$ -t 1-163
#$ -l cancergen

#-------
# Input
#---------------------------------------------------------------------------------------------------------------
bamDir=../bams
outDir=../outDir
tmpDir=../temp
nProcesses=16
#---------------------------------------------------------------------------------------------------------------

module load conda_R/4.0.x

bamFile=$(ls -1v $bamDir/*bam | head -n $SGE_TASK_ID | tail -n 1)

mkdir -p $outDir

Rscript ./i05-count-sbs.R $bamFile $outDir $tmpDir $nProcesses

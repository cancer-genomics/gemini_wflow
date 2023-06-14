outDir=../outDir/chain

mkdir -p $outDir

cd $outDir

wget https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/hg19ToHg38.over.chain.gz
gzip -d hg19ToHg38.over.chain.gz


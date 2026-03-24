#Software dependency
python
bowtie2
cutadapt (This is cutadapt 3.7 with Python 3.6.9)
samtools
bedtools

conda install bioconda::multiqc

#python package dependency
sys
gzip
pysam
click
os
pickle
numpy
pandas
progress.bar
collections

#R package dependency
ggplot2
dplyr
gridExtra
plyr
labeling


library(ggplot2)
library(dplyr)
library(gridExtra)
library(plyr)
library(labeling)




#STEP 0

mkdir Genome
cd Genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz 
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
wget https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/hg38.full.blacklist.bed
gunzip hg38.fa.gz

#STEP 0

bedtools maskfasta  -fi hg38.fa -bed hg38.full.blacklist.bed -fo hg38.mitoMask.fa
bowtie2-build hg38.mitoMask.fa hg38.mitoMask  ## This takes time and will generate following, if on WI server ln -s /lab/solexa_weissman/cweng/Genomes/GRCH38/GRCH38_Bowtie2_MitoMask/*.bt2 ./

#     color: 0
#     reverse: 1
# Total time for backward call to driver() for mirror index: 00:22:35

#/disk1/cai016/redeemV-master



## Restart From here
## assign path

REDEEM_V=/disk1/cai016/redeemV-master #The loacation where the REDEEM-V is downloaded to (eg. /lab/solexa_weissman/cweng/Packages/REDEEM-V/)
MyMultiome=$REDEEM_V/MyMultiome/
mitoConsensus=$REDEEM_V/mitoConsensus/

#STEP 1
ln -s $REDEEM_V/source/Example.R1.fq.gz ./
ln -s $REDEEM_V/source/Example.R2.fq.gz ./
ln -s $REDEEM_V/source/Example.i5.fq.gz ./

## install fastqc in ubuntu
#sudo apt-get -y install fastqc
#conda install multiqc #failed
#pip install multiqc

#mkdir fastqc
#fastqc -o fastqc Example.R1.fq.gz Example.i5.fq.gz Example.R2.fq.gz
#multiqc fastqc

#BiocManager :: install ( "qvalue") 
##
chmod u+x *.sh ##change the property
chmod u+x *.R ##change the property
chmod u+x *.py ##change the property

cutadapt --cores=$CORE -a CTGTCTCTTATA -A CTGTCTCTTATA -o $Read1.trim -p $Read2.trim $Read1 $Read2



## STEP 1 
cutadapt -a CTGTCTCTTATA -A CTGTCTCTTATA -o Example1.trim -p Example.trim Example.R1.fq.gz  Example.R2.fq.gz 
## STEP 1 

ln -s /disk1/cai016/Genome/GRCH38/GRCH38_Bowtie2_MitoMask/*.bt2 ./

#bowtie2Index=/lab/solexa_weissman/cweng/Genomes/GRCH38/GRCH38_Bowtie2_MitoMask/hg38.mitoMask
bowtie2Index=/disk1/cai016/Genomes/GRCH38/GRCH38_Bowtie2_MitoMask/hg38.mitoMask
bowtie2 -X 1200  --very-sensitive -p 1 -x $bowtie2Index -1 $Read1.trim  -2 $Read2.trim | samtools view -@ $CORE -bS - > $name.tmp.bam

Example.R1.fq.gz.trim

bowtie2 -X 1200  --very-sensitive -p 1 -x $bowtie2Index -1 Example.R1.fq.gz.trim  -2 Example.R2.fq.gz.trim | samtools view -@ 1 -bS - > $name.tmp.bam


#bowtie2Index=/lab/solexa_weissman/cweng/Genomes/GRCH38/GRCH38_Bowtie2_MitoMask/hg38.mitoMask
bowtie2 -X 1200  --very-sensitive -p $CORE -x $bowtie2Index -1 $Read1.trim  -2 $Read2.trim | samtools view -@ $CORE -bS - > $name.tmp.bam


CORE = 10
bowtie2Index=/disk1/cai016/Genomes/GRCH38/GRCH38_Bowtie2_MitoMask/hg38.mitoMask

## STEP 2
bowtie2 -X 1200  --very-sensitive -p 2 -x hg38.mitoMask -1 Example.R1.fq.gz.trim  -2 Example.R2.fq.gz.trim | samtools view -@ 2 -bS - > Tutorial_Mito.tmp.bam
samtools sort -@ 2 Tutorial_Mito.tmp.bam > Tutorial_Mito.bam
samtools index -@ 2 Tutorial_Mito.bam
## STEP 2


$MyMultiome/MultiomeATAC_mito.sh -n Tutorial_Mito -1 Example.R1.fq.gz -2 Example.R2.fq.gz -i Example.i5.fq.gz -c 0 -t 8 -m $MyMultiome -b Genome/hg38.mitoMask  # add -q for quick, which skip the QC steps
$MyMultiome/MultiomeATAC_mito.sh -h


$MyMultiome/MultiomeATAC_mito.sh -n Tutorial_Mito -1 Example.R1.fq.gz -2 Example.R2.fq.gz -i Example.i5.fq.gz -c 0 -t 8 -m $MyMultiome -b /disk1/cai016/Genome/hg38.mitoMask  # add -q for quick, which skip the QC steps






ln -s $REDEEM_V/source/barcodes.tsv.gz ./
Rscript $MyMultiome/Helpers/RNAbc2ATAC.R $REDEEM_V barcodes.tsv.gz Tutorial_atac.barcodes.tsv




Threads=24
WD=`pwd`/Out_mitoConsensus
python $mitoConsensus/Preprocess.py -i Tutorial_Mito.uniqmapped.bam  -c $Threads -b Tutorial_atac.barcodes.tsv -o $WD -g rCRS -bt BC -sd $mitoConsensus


sudo apt-get install git-lfs

echo $WD
for ((i=1;i<=$Threads;i++))
do
bsub python $mitoConsensus/mitoConsensus.py  barcodes.$i $WD BC 30
done


$mitoConsensus/Finalize.sh $WD $Threads $mitoConsensus






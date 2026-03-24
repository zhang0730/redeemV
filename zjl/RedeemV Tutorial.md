Updated on 2022-10-25
# Introduction:
This tutorial use a small example to go through ReDeeM-V pipeline, which starts with mtDNA library fastq files and end with consensus variant calling 

#本教程使用一个小示例来介绍ReDeeM-V管道，该管道从mtDNA库fastq文件开始，以共识变量调用结束 

- **Inputs**: R1.fq (150nt), R2.fq (150nt), i5.fq (24nt)
- **Key Output:** QCplot.png (QC plot assesing mtDNA mapping as well as library complexity, etc)
- **Key Output:** QualifiedTotalCts, RawGenotypes.Total.StrandBalance, RawGenotypes.VerySensitive.StrandBalance, RawGenotypes.Sensitive.StrandBalance, RawGenotypes.Specific.StrandBalance
            

# Prerequisite 

## Download the package
```console
git clone https://github.com/chenweng1991/REDEEM-V.git
```

#这里要先把他的文件下载下来 路径后面会用到

然后下面的这些包也很重要 建议先下齐  不然会有摸不着头脑的Error

### Software dependency 

- python
- bowtie2
- cutadapt (This is cutadapt 3.7 with Python 3.6.9)
- samtools
- bedtools

### python package dependency
- sys
- gzip
- pysam
- click
- os
- pickle
- numpy
- pandas
- progress.bar
- collections

### R package dependency
- ggplot2
- dplyr
- gridExtra
- plyr
- labeling
- qvalue

# Main pipeline
## Assign paths
```console
REDEEM_V=ThePathToREDEEM-V #The location where the REDEEM-V is downloaded to (eg. /lab/solexa_weissman/cweng/Packages/REDEEM-V/)

#我的是REDEEM_V=/disk1/cai026/redeem/REDEEM-V/

MyMultiome=$REDEEM_V/MyMultiome/
mitoConsensus=$REDEEM_V/mitoConsensus/
#目录变量 之前的git clone下载到哪里了？
```
## Step 0: Make mitoMask bowtie2 index
Download hg38 fasta file and mitochondrial blacklist file 

下载人类基因组 hg38 的 FASTA 文件、染色体大小文件以及一个线粒体黑名单 文件。接着解压 FASTA 文件。

```console
mkdir Genome
cd Genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz 
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes
wget https://raw.githubusercontent.com/caleblareau/mitoblacklist/master/combinedBlacklist/hg38.full.blacklist.bed
gunzip hg38.fa.gz
```
Make hg38 mito masked bowtie2 index

#对基因组的线粒体黑名单进行掩码处理，并随后构建Bowtie2索引

```console
bedtools maskfasta  -fi hg38.fa -bed hg38.full.blacklist.bed -fo hg38.mitoMask.fa 

#使用 bedtools maskfasta 来对基因组的线粒体黑名单序列进行掩码处理。这个命令会将基因组中的某些区域（由BED文件定义）标记为“N”（即掩码），从而使得这些区域在后续的分析中被忽略。

bowtie2-build hg38.mitoMask.fa hg38.mitoMask  

## This takes time and will generate following, if on WI server ln -s /lab/solexa_weissman/cweng/Genomes/GRCH38/GRCH38_Bowtie2_MitoMask/*.bt2 ./
#用于使用Bowtie2构建索引，以便后续可以快速进行比对。hg38.mitoMask.fa: 输入的FASTA文件，这是经过掩码处理的基因组。
#hg38.mitoMask: 输出的Bowtie2索引的前缀，索引文件将以 hg38.mitoMask.*.bt2 的形式生成。
#这一步将花费大量时间
```

## Step1: Get sequencing result and QC 
fastqc and multiqc need to be installed to run Step 1, but theyr'e not required for REDEEM-V </br>
For this tutorial, Link the example fq files

#创建软链接（symbolic links），以便在当前目录下创建指向原始文件的链接。这些命令主要用于方便地在当前工作目录中使用这些文件，而不需要复制它们。

所以上面的路径正确就显得很重要了

```console
ln -s $REDEEM_V/source/Example.R1.fq.gz ./
ln -s $REDEEM_V/source/Example.R2.fq.gz ./
ln -s $REDEEM_V/source/Example.i5.fq.gz ./
```
Fastqc and multiqc

#使用FastQC工具对多个测序文件进行质量控制，并使用MultiQC工具汇总这些质量控制报告。

```console
mkdir fastqc
fastqc -o fastqc Example.R1.fq.gz Example.i5.fq.gz Example.R2.fq.gz #安装fastqc
multiqc fastqc #安装multiqc
```
## Step2: mapping and basic QC
```console
$MyMultiome/MultiomeATAC_mito.sh -n Tutorial_Mito -1 Example.R1.fq.gz -2 Example.R2.fq.gz -i Example.i5.fq.gz -c 0 -t 8 -m $MyMultiome -b Genome/hg38.mitoMask  # add -q for quick, which skip the QC steps
#上面两行是一条命令
#用于执行一个脚本 MultiomeATAC_mito.sh  仔细学习这个脚本
首先第一步使用cutadapt修剪fastq文件中的适配器序列
第二步是mapping
第三步是提取细胞barcode
第四步是使用samtools过滤BAM文件，获取唯一映射 后生成bigwig文件并且调用峰值
第五步使用samtools从唯一映射的BAM中提取线粒体基因组 创建索引 输出包含线粒体读取的BAM
第六步使用samtools将BAM转换成bed文件 并添加cell barcode
第七步使用python脚本去除单细胞水平上的重复事件  生成一系列文件：
         .uniqmapped.fragment.tsv
         .uniqmapped.fragment.1000cut.tsv
第八步使用grep命令从去重后的单细胞片段tsv中提取与线粒体基因组（chrM）相关的数据
第九步将峰值与单细胞数据关联起来 生成矩阵文件
第十步使用python脚本汇总去重后的单细胞片段，包括非线粒体和线粒体
第十一步计算Readcounts
第十二步生成质量控制图表
清理中间文件
rm -rf tmp
rm *trim
rm -rf *tmp.bam
rm -rf *uniqmapped.RawBed
rm -rf *uniqmapped.fragment.0.cut.tsv
rm -rf *uniqmapped.fragment.0.cut.summary
```
After this step a QCplot is generated
![Old1_HSPC_Mito QCplot](https://user-images.githubusercontent.com/43254272/198318936-1c7f1f4b-c203-4b93-8997-b1d82adc3b62.png)

## Step3: Parse the Cellranger qualified cell barcode
In a ReDeeM experiment, the ATAC and RNA data is analyzed by Cellranger arc. The barcodes.tsv.gz from that will be parsed here.
Becasue the barcode in barcodes.tsv.gz is RNA barcode, we need to translate into ATAC barcode (via RNAbc2ATAC.R) to be able to match mtDNA data

使用Cellranger arc对ATAC和RNA数据进行分析。 这里将解析其中的barcode .tsv.gz。  

由于barcode .tsv.gz中的条形码是RNA条形码，我们需要翻译成ATAC条形码(通过RNAbc2ATAC.R)才能匹配mtDNA数据 

```console
ln -s $REDEEM_V/source/barcodes.tsv.gz ./
Rscript $MyMultiome/Helpers/RNAbc2ATAC.R $REDEEM_V barcodes.tsv.gz Tutorial_atac.barcodes.tsv
```
## Step4, prepare for variant calling

- 为变量调用做准备

```console
Threads=24
WD=`pwd`/Out_mitoConsensus

python $mitoConsensus/Preprocess.py -i Tutorial_Mito.uniqmapped.bam  -c $Threads -b Tutorial_atac.barcodes.tsv -o $WD -g rCRS -bt BC -sd $mitoConsensus
#上面两个是一行代码
#在Out_mitoConsensus文件夹里面有3个文件了
```
## Step5, run consensus variant calling
Run variant calling in parallel
```console
echo $WD
for ((i=1;i<=$Threads;i++))
do
bsub python $mitoConsensus/mitoConsensus.py  barcodes.$i $WD BC 30
done
#做循环 
bsub命令是在LSF包内
其实也可以不用bsub  直接跑循环
```
Concat together
```console
$mitoConsensus/Finalize.sh $WD $Threads $mitoConsensus
```

# Expected result 重点
- **Location: The results are saved in $WD/final**

- **Major 5 files:** 
  
  - QualifiedTotalCts, 
  - RawGenotypes.Total.StrandBalance(Least stringent), 
  - RawGenotypes.VerySensitive.StrandBalance(Less stringent), 
  - RawGenotypes.Sensitive.StrandBalance(Stringent), 
  - RawGenotypes.Specific.StrandBalance (Most Stringent)
  
- **QualifiedTotalCts** is a table with 6 columnes that show mtDNA coverage per position per cell

  **QualifiedTotalCts**是一个有6列的表，显示每个细胞每个位置的mtDNA覆盖率 

| Cellbarcode| coordinates on mt genome|# unique frag(total)|# unique frag(less stringent)|# unique frag(stringent)|# unique frag(very stringent)|
| ------------- |----------------------|--------------------|-----------------------------|------------------------|-----------------------------|

- **RawGenotypes** is a table with 14 columnes that show the consensus variant calling.  Each row is a molecule with a potential variant

  RawGenotypes是一个有14列的表，显示了共识变量。 每一行都是一个有潜在变异的分子 

|MoleculeID | CellBC | Pos | Variant | V | Ref | FamSize | V-counts | CSS | DB_Cts | SG_Cts | Is+ | Is- | TotalDepth|
|-----------|--------|-----|--------|---------|----------|---------|---------|-----|--------|---------|-----|-----|-----------|

1. MoleculeID: Cellbarcode+start+end which is the identifier to **define a molecule**
2. CellBC: Cell barcode
3. Pos: The coordinate of the variant
4. Variant: A description of the variant
5. V: The variant base called on Pos
6. Ref: The reference base on Pos
7. FamSize: The consensus family size, or the total number of PCR copies for the given molecule
8. V-count: Number of PCR copies that support the variant
9. CSS: Consensus score, which is the proportion of PCR copies that support the variant
10. DB_**Ct**s: Number of double cover copies, which are positions that sequenced by both Read1 and Read2
11. SG_Cts: Number of single cover copies, which are positions that sequenced by only Read1 or Read2
12. Is+: If the variant is discovered on plus strand
13. Is- : If the variant is discovered on minus strand
14. TotalDepth: On this given position, total number of unique fragment in the given cell

- Strand biased variant has been removed

   链偏变体被移除
  ![StrandBiase](https://user-images.githubusercontent.com/43254272/198328824-40977739-6fdf-4813-9461-9c5bee18d53a.png)

- **QualifiedTotalCts and RawGenotypes**. are the inputs for REDEEM-R for downstream mutation filtering and phylogenetic tree reconstruction

  QualifiedTotalCts和RawGenotypes。 *为REDEEM-R用于下游突变过滤和系统发育树重建的输入 


# (Optional) Cell Hashing

- First, generate barcodes.tsv.16bp from the Cell Ranger result
```
zcat ../CellRanger/OBC/outs/filtered_feature_bc_matrix/barcodes.tsv.gz | awk '{print substr($1,1,16)}' > barcodes.tsv.16bp
```
- Second, Create HashDic.csv like below
```
## Create a HashDic.csv
#TGTCTTTCCTGCCAG,SBM1037
#CTCCTCTGCAATTAC,SBM1253
#CAGTAGTCACGGTCA,SBM1324
#ATTGACCCGCGTTAG,SBM1218
```
- Third run the script below, which will generate a plot.pdf and a folder pymulti that contains a result table
```
python3 $REDEEMV/JohnnyCellHash/RunCellHash.py $REDEEMV OBC OBC_Hash_S1_L001_R1_001.fastq.gz.trim OBC_Hash_S1_L001_R3_001.fastq.gz.trim
```

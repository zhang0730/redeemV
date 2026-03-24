## Introduction



ReDeeM: single-cell **Re**gulatory multi-omics with **Dee**p **M**itochondrial mutation profiling. ReDeeM is a single-cell multiomics platform featuring ultra-sensitive mitochondrial DNA (mtDNA) variant calling and joint RNA+ATAC profiling. ReDeeM enables fine-scale lineage tracing at single cell level, allowing for subclonal and phylogenetic analyses, with simultaneous integrative analyses of cell-state and gene regulatory circuits.

The analytical pipelines for ReDeeM analysis includes two parts: **两部分：python在上游  R在下游**

- [redeemV](https://github.com/chenweng1991/redeemV) is set of in-house Bash pipeline and python scripts for mapping and deep mitochondrial mutation calling. **(Input from rawdata)**
- [redeemR](https://github.com/chenweng1991/redeemR) is an in-house R package for downstream lineage tracing and single cell integrative analysis. (**This page**, **Input from redeemV)**



**redeemV** is a streamlined pipeline taking advantage of endogenous unique molecular identifier (eUMI) for consensus-based error correction in single-cell mitochondrial DNA mutation detection.

单细胞线粒体DNA突变检测

![image-20240908173938090](C:\Users\13623\AppData\Roaming\Typora\typora-user-images\image-20240908173938090.png)

​	

**redeemR** is designed to refine the somatic mitochondiral mutations, build the single-cell phylogenetic tree and facilitate genealogical and clonal/subclonal-resolved single-cell multiomics analysis.

构建单细胞系统发育树，促进谱系和克隆/亚克隆分辨单细胞多组学分析。 

![image-20240908173948782](C:\Users\13623\AppData\Roaming\Typora\typora-user-images\image-20240908173948782.png)



### 网址：

https://github.com/chenweng1991/redeemV?tab=readme-ov-file

https://github.com/sankaranlab/redeemR?tab=readme-ov-file



## **redeemV**

教程网址：

https://github.com/chenweng1991/redeemV/blob/master/Tutorial_20221025.md
# scHuman_organoids_COVID19

# Title: Integration of human organoids single-cell transcriptomic profiles and human genetics repurposes critical cell type-specific drug targets for severe COVID-19 
Authors: Yunlong Ma1,2,*, Yijun Zhou1,2,*, Dingping Jiang1,*, Wei Dai3*, Cheng Chen1 , Jingjing Li1, Chunyu Deng4, Gongwei Zheng1, Yaru Zhang1,2, Fei Qiu1, Haojun Sun1, Yinghao Yao2, Haijun Han5, Jia Qu1, Nan Wu6, Jianzhong Su1,2,#

1 School of Biomedical Engineering, School of Ophthalmology & Optometry and Eye Hospital, Wenzhou Medical University, Wenzhou, 325027, China.
2 Oujiang Laboratory, Zhejiang Lab for Regenerative Medicine, Vision and Brain Health, Wenzhou 325101, Zhejiang, China.
3 Wenzhou Institute, University of Chinese Academy of Sciences, Wenzhou, 325011, China.
4 School of Life Science and Technology, Harbin Institute of Technology, Harbin, 150080, China.
5 School of Medicine, Hangzhou City University, Hangzhou 310015, China.
6 Beijing Key Laboratory for Genetic Research of Skeletal Deformity, Key laboratory of big data for spinal deformities, Department of Orthopedic Surgery, Peking Union Medical College Hospital, Peking Union Medical College and Chinese Academy of Medical Sciences, Beijing 100730, China.
*These authors contributed equally to this work.
#Corresponding author

# Keywords: 
```
Human organoids, Single-cell RNA sequencing, GWAS, COVID-19, Drug targets
```

# Summary
Human organoids recapitulate the cell type diversity and function of their primary organs holding tremendous potentials for basic and translational research. Advances in single-cell RNA sequencing (scRNA-seq) technology and genome-wide association study (GWAS) have accelerated the biological and therapeutic interpretation of trait-relevant cell types or states. Here, we constructed a computational framework to integrate atlas-level organoid scRNA-seq data, GWAS summary statistics, expression quantitative trait loci, and gene-drug interaction data for distinguishing critical cell populations and drug targets relevant to COVID-19 severity. We found that 39 cell types across eight kinds of organoids were significantly associated with COVID-19 outcomes. Notably, subset of lung mesenchymal stem cells (MSCs) increased proximity with fibroblasts predisposed to repair COVID-19-damaged lung tissue. We repurposed 33 druggable genes, including IFNAR2, TYK2, and VIPR2, and their interacting drugs for COVID-19 in a cell-type-specific manner. Overall, our results showcase that host genetic determinants have cellular specific contribution to COVID-19 severity, and identification of cell type-specific drug targets may facilitate to develop effective therapeutics for treating severe COVID-19 and its complications.


# Computational Framework
To facilitate the data integration and minimize the batch effects, we have built a unified pipeline to conduct re-alignment, quality control, and standard analysis of all human organoids (n = 1,159,206 cells) and fetal scRNA-seq datasets (n =223,334 cells, Supplementary Figure S1 and Table S1). To distinguish critical cell types/subpopulations and repurpose potential drugs and interacting targets for the treatment of severe COVID-19, we devise a computational framework to incorporate these organoids and fetal scRNA-seq data and large-scale meta-GWAS summary statistics on three COVID-19 phenotypes (i.e., very severe, hospitalized, and susceptible COVID-19; Figure 1, and Supplementary Table S2 and Figures S2-S3). There are three main sections: (1) integrating GWAS summary statistics with human organoids scRNA-seq datasets to genetically map trait-relevant single-cell landscapes for three COVID-19 outcomes (Figure 1A); (2) combining GWAS summary statistics with eQTL data in the GTEx database to identify putative risk genes and critical pathways associated with COVID-19 severities (Figure 1B); and (3) prioritization of cell type-specific gene-drug interaction pairs for treating severe COVID-19 and related complications at a fine-grained resolution (Figure 1C). 

![Figure 1](C:/Users/MYL/Desktop/images/Figure_1.png "Figure_1.png")

# Scripts
In this study, we have built a unified pipeline and devised a computational framework, which contains a series of bioinformatic tools to fulfil the integrative genomic analyses. For example:
```
library(scPagwas)
 #1.Leveraging scPagwas to integrating GWAS summary statistics on three COVID-19 outcomes and scRNA-seq data on human lung, intestinal, and brain organoids.
 Pagwas_data<-scPagwas_main(Pagwas = NULL,
                     gwas_data =system.file("extdata", "GWAS_summ_example.txt", package = "scPagwas"), # The GWAS Summary statistics files 
                     Single_data =system.file("extdata", "scRNAexample.rds", package = "scPagwas"),# scRNA-seq data in seruat format with "RNA" assays and normalized.
                     output.prefix="test", # the prefix name for output files
                     output.dirs="scPagwastest_output",# the directory file's name for output
                     block_annotation = block_annotation,# gene position in chromosome is provided by package.
                     assay="RNA", # the assays for scRNA-seq data to use.
                     Pathway_list=Genes_by_pathway_kegg,# pathway list is provided by package, including gene symbols.
                     n.cores=1,
                     iters_singlecell = 100,
                     chrom_ld = chrom_ld,# The LD data is provided by package.
                     singlecell=T, # Whether to run the singlecell process.
                     celltype=T# Whether to run the celltype process.)
 ```

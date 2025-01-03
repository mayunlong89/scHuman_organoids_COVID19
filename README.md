# scHuman_organoids_COVID19

This project have been done and published in [Ma et al. Cell Proliferation, 2024](https://onlinelibrary.wiley.com/doi/full/10.1111/cpr.13558).

All the single-cell data on multiple human organoids has been documented in the [scHOB](https://schob.su-lab.org/) database.

This analytic results and visualization are shown in the Website, see [Link](https://schob.su-lab.org/function/).

We applied our recently developed method scPagwas to integrate human organoids scRNA-seq data with GWAS summary data on three COVID-19 phenotypes. (see the method code in [Github link](https://github.com/mayunlong89/scPagwas_main), and the paper in [Ma et al. Cell Genomics, 2023](https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00180-5)).

# All COVID-19-related projects in our group:
1) Meta-analysis of large-scale GWAS data to uncover novel loci for COVID-19. see [Ma et al. Human Molecular Genetics, 2021](https://academic.oup.com/hmg/article/30/13/1247/6265026), and see related [Github codes](https://github.com/mayunlong89/Host_genetics_for_COVID-19).
2) COVID-19 Quarantine Reveals That Behavioral Changes Have an Effect on Myopia Progression. see [Xu, Ma et al. Ophthalmology, 2021](https://www.sciencedirect.com/science/article/pii/S0161642021002578), see related [Github codes](https://github.com/mayunlong89/MIES).
3) Identification of genetics-influenced immune cell sub-populations relevant to severe COVID-19. see [Ma et al. Genome Medicine, 2022](https://link.springer.com/article/10.1186/s13073-022-01021-1), and see related [Github codes](https://github.com/mayunlong89/COVID19_scRNA).
4) Development of novel polygenic regression method scPagwas for integrating scRNA-seq data with GWAS on complex diseases. see [Ma et al. Cell Genomics, 2023](https://www.cell.com/cell-genomics/fulltext/S2666-979X(23)00180-5), and see related [Github codes](https://github.com/mayunlong89/scPagwas_main).
5) Repurposing cell type-specific durg targets for severe COVID-19 based on human organoids scRNA-seq atlas. see [Ma et al. Cell Proliferation, 2024](https://onlinelibrary.wiley.com/doi/full/10.1111/cpr.13558), and see related [Github codes](https://github.com/mayunlong89/scHuman_organoids_COVID19).

# Keywords: 
```
Human organoids, Single-cell RNA sequencing, GWAS, COVID-19, Drug targets
```

# Summary
Human organoids recapitulate the cell type diversity and function of their primary organs holding tremendous potentials for basic and translational research. Advances in single-cell RNA sequencing (scRNA-seq) technology and genome-wide association study (GWAS) have accelerated the biological and therapeutic interpretation of trait-relevant cell types or states. Here, we constructed a computational framework to integrate atlas-level organoid scRNA-seq data, GWAS summary statistics, expression quantitative trait loci, and gene-drug interaction data for distinguishing critical cell populations and drug targets relevant to COVID-19 severity. We found that 39 cell types across eight kinds of organoids were significantly associated with COVID-19 outcomes. Notably, subset of lung mesenchymal stem cells (MSCs) increased proximity with fibroblasts predisposed to repair COVID-19-damaged lung tissue. We repurposed 33 druggable genes, including IFNAR2, TYK2, and VIPR2, and their interacting drugs for COVID-19 in a cell-type-specific manner. Overall, our results showcase that host genetic determinants have cellular specific contribution to COVID-19 severity, and identification of cell type-specific drug targets may facilitate to develop effective therapeutics for treating severe COVID-19 and its complications.


# Computational Framework
To facilitate the data integration and minimize the batch effects, we have built a unified pipeline to conduct re-alignment, quality control, and standard analysis of all human organoids (n = 1,159,206 cells) and fetal scRNA-seq datasets (n =223,334 cells, Supplementary Figure S1 and Table S1). To distinguish critical cell types/subpopulations and repurpose potential drugs and interacting targets for the treatment of severe COVID-19, we devise a computational framework to incorporate these organoids and fetal scRNA-seq data and large-scale meta-GWAS summary statistics on three COVID-19 phenotypes (i.e., very severe, hospitalized, and susceptible COVID-19; Figure 1, and Supplementary Table S2 and Figures S2-S3). There are three main sections: (1) integrating GWAS summary statistics with human organoids scRNA-seq datasets to genetically map trait-relevant single-cell landscapes for three COVID-19 outcomes (Figure 1A); (2) combining GWAS summary statistics with eQTL data in the GTEx database to identify putative risk genes and critical pathways associated with COVID-19 severities (Figure 1B); and (3) prioritization of cell type-specific gene-drug interaction pairs for treating severe COVID-19 and related complications at a fine-grained resolution (Figure 1C). 

![Figure 1](https://github.com/mayunlong89/scHuman_organoids_COVID19/blob/main/figures/Figure%201.png)



# Citation
Ma et al. Integration of human organoids single-cell transcriptomic profiles and human genetics repurposes critical cell type-specific drug targets for severe COVID-19, [Cell Proliferation, 2023,e13558](https://onlinelibrary.wiley.com/doi/full/10.1111/cpr.13558).


# Scripts
In this study, we have built a unified pipeline and devised a computational framework, which contains a series of bioinformatic tools to fulfil the integrative genomic analyses. For example:
```r
# 1.Leveraging scPagwas to integrating GWAS summary statistics on three COVID-19 outcomes and scRNA-seq data on human lung, intestinal, and brain organoids.
library(scPagwas)
library(readr)
library(dplyr)
library(Seurat)
library(tidyverse)
library(rhdf5)
library(ggplot2)
library(grDevices)
library(stats)
library(FactoMineR)
library(scales)
library(reshape2)
library(ggdendro)
library(grImport2)
library(gridExtra)
library(grid)
library(sisal)
require("RColorBrewer")
require("ggsci")
require("ggpubr")
source(system.file("extdata", "plot_scpathway_contri_dot.R", package = "scPagwas"))

Args <- commandArgs(T)
memo = print(Args[1])
gwas = print(Args[2])
stage = print(Args[3])

data_path <-  file.path(glue::glue("analysis/{memo}/drived_data"))
seu_path <-  file.path(glue::glue("analysis/{memo}/drived_data/{memo}.rds"))
figure_path <- file.path(glue::glue("analysis/{memo}/figures/{stage}"))

if(!dir.exists(figure_path)){
dir.create(figure_path)
}else{
    print('Dir exists')
}

seu <- readRDS(seu_path)

rf_10x <- Read10X('/share2/pub/zhouyj/zhouyj/organoid/data/brain/GSE117512/CR_S1/outs/filtered_feature_bc_matrix')
rf_seu <- CreateSeuratObject(counts = rf_10x)

seu@assays$RNA@data@Dimnames[[1]] <- rf_seu@assays$RNA@data@Dimnames[[1]]

Idents(seu) <- seu@meta.data$annotation2
table(Idents(seu))

Pagwas <- scPagwas_main(Pagwas = NULL,
                             gwas_data = gwas,
                             Single_data = seu,
                             output.prefix=glue::glue("{memo}_{stage}"),
                             output.dirs=data_path,
                             singlecell=T,
                             celltype=T,
                             Pathway_list=Genes_by_pathway_kegg,
                             assay="RNA",
                             block_annotation = block_annotation,
                             seruat_return=T,
                             ncores = 1,
                             chrom_ld=chrom_ld
)

#saveRDS(Pagwas,file = glue::glue("analysis/{memo}/drived_data/{stage}_pagwas.rds"))

color26 <- c("#D9DD6B","#ECEFA4","#D54C4C","#8D2828","#FDD2BF","#E98580","#DF5E5E","#492F10","#334257","#476072","#548CA8",
"#00A19D","#ECD662","#5D8233","#284E78","#3E215D","#835151","#F08FC0","#C6B4CE","#BB8760","#FFDADA","#3C5186",
"#558776","#E99497","#FFBD9B","#0A1D37")

png(glue::glue("analysis/{memo}/figures/{stage}/plot_gpas_trs_pval.png"),width = 600, height = 2000)
scPagwas_Visualization(Single_data=Pagwas,
                        p_thre = 0.05,
                        FigureType = "umap",
                        width = 7,
                        height = 7,
                        lowColor = "white", 
                        highColor = "red",
                        output.dirs="figure",
                        size = 0.5,
                        do_plot = T)
dev.off()

png(glue::glue("analysis/{memo}/figures/{stage}/plot_bar_positie_nagtive1.png"),width = 600, height = 600)
plot_bar_positie_nagtive(seurat_obj=Pagwas,
                         var_ident="annotation2",
                         var_group="positiveCells",
                         vec_group_colors=c("#E8D0B3","#7EB5A6"),
                         do_plot = T)
dev.off()

png(glue::glue("analysis/{memo}/figures/{stage}/plot_bar_positie_nagtive2.png"),width = 600, height = 600)
plot_bar_positie_nagtive(seurat_obj=Pagwas,
                              var_ident="positiveCells",
                              var_group="annotation2",
                              p_thre = 0.01,
                              vec_group_colors=NULL,
                              f_color=colorRampPalette(brewer.pal(n=10, name="RdYlBu")),
                              do_plot = T)
dev.off()


ct <-  as.vector(unique(Pagwas@meta.data$annotation2))
png(glue::glue("analysis/{memo}/figures/{stage}/plot_scpathway_dot.png"),width = 1000, height = 600)
plot_scpathway_dot(Pagwas=Pagwas,
                   celltypes=ct,
                   topn_path_celltype=5,
                   filter_p=0.05,
                   max_logp=15,
                   display_max_sizes=F,
                   size_var ="logrankPvalue" ,
                   col_var="proportion",
                   shape.scale = 8,
                   cols.use=c("lightgrey", "#E45826"),
                   dend_x_var = "logrankPvalue",
                   dist_method="euclidean",
                   hclust_method="ward.D",
                   do_plot = T,
                   #figurenames = "Pathway_plot.pdf",
                   width = 7,
                   height = 7)
dev.off()

pdf(glue::glue("analysis/{memo}/figures/{stage}/heritability_cor_scatterplot.pdf"))
heritability_cor_scatterplot(gene_heri_cor=Pagwas@misc$gene_heritability_correlation,
                             topn_genes_label=10,
                             color_low="#035397",
                             color_high ="#F32424",
                             color_mid = "white",
                             text_size=2,
                             do_plot=T,
                             max.overlaps =20,
                             width = 7,
                             height = 7)
dev.off()

top5genes<-rownames(Pagwas@misc$gene_heritability_correlation)[order(Pagwas@misc$gene_heritability_correlation,decreasing = T)[1:5]]
pdf(glue::glue("analysis/{memo}/figures/{stage}/plot_vln_Corgenes.pdf"),width = 6, height =7)
plot_vln_Corgenes(seurat_obj=Pagwas,
             assay="RNA", slot="data",
             var_group="annotation2",
             vec_features=top5genes,
             vec_group_colors= color26,
             do_plot = T
             )
dev.off()

pdf(glue::glue("analysis/{memo}/figures/{stage}/bootstrap_p_Plot.pdf"))
Bootstrap_P_Barplot(p_results=Pagwas@misc$bootstrap_results$bp_value[-1],
                    p_names=rownames(Pagwas@misc$bootstrap_results)[-1],
                    figurenames = "Bootstrap_P_Barplot.pdf",
                    width = 5,
                    height = 7,
                    do_plot=T,
                    title = glue::glue("{stage}_celltype")
                    )
dev.off()

pdf(glue::glue("analysis/{memo}/figures/{stage}/estimate_Plot.pdf"))
Bootstrap_estimate_Plot(Pagwas=Pagwas,
                        width = 9,
                        height = 7,
                        do_plot=T)
dev.off()

 ```

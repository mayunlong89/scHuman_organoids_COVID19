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

library(Seurat)
library(ggplot2)
library(dplyr)

# load data
seu = readRDS("all.clean.raw.rds")
sue = NormalizeData(seu)

# read meta
meta = readRDS("all.clean.meta.rds")

# add meta
seu@meta.data = meta[colnames(seu),]

# marker genes of CD8T
# Figure 1f
if(T){
  # subset
  seu_subset = subset( seu, cells  = rownames(seu@meta.data)[seu$major_cluster == "CD8T"] )
  seu_subset$minor_cluster = factor(seu_subset$minor_cluster, 
                                    levels = seu_subset$minor_cluster %>% droplevels() %>% levels %>% rev)
  
  DotPlot(seu_subset, 
          dot.scale = 8,
    features = c(
    "SELL","CCR7",#"TCF7","LEF1", "S1PR1", # Tn/Tcm
    "CX3CR1","NKG7", # Teff
    "GZMK","CD44", # Tem
    "MYADM","CD69",# Trm
    "CD160","KLRC2", # IEL
    "IL17A","IL22","IL26","PDCD1", # Tc17
    "SLC4A10",'ZBTB16', # MAIT
    "ZNF683","IKZF2" # act T
  ), group.by = "minor_cluster") +
    xlab("") + ylab("") +  
    scale_color_gradient2() + 
    theme_classic2() + 
    theme(
      legend.position = "right",
      axis.text.x = element_text(family = "ArialMT", face = "italic",size = 16,  color = "black", angle = 90, hjust = 1),
      axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "italic", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    )

  # save plot
  # ggsave("CD8T_marker_genes.dotplot.pdf", height = 3.8, width = 8.5)
}

# marker genes of myeloid 
# Figure 1i
if(T){
  # subset
  seu_subset = subset( seu, cells  = rownames(seu@meta.data)[seu$major_cluster == "Myeloid"] )
  seu_subset$minor_cluster = factor(seu_subset$minor_cluster, 
                                    levels = seu_subset$minor_cluster %>% droplevels() %>% levels %>% rev)
  
  genes = c('LYZ', 'FCGR3A',"FCN1","S100A8",
            "CXCL2","IL1B",
            'C1QA', 'C1QB', 
            'APOE', 'APOC1',
            'LYVE1', 'CCL24',
            'AREG', 'CPM',
            'MKI67', 'TOP2A',
            "CLEC9A","CD1C","CLEC4C",
            'LAMP3', 'CCR7',
            "TPSAB1","PPBP")
  
  DotPlot(seu_subset, dot.scale = 8,
          features = genes, group.by = "minor_cluster") +
    xlab("") + ylab("") +  
    scale_color_gradient2() + 
    theme_classic2() + 
    theme(
      legend.position = "right",
      axis.text.x = element_text(family = "ArialMT", face = "italic",size = 16,  color = "black", angle = 90, hjust = 1),
      #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "italic", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    )

  # save plot
  # ggsave("Myeloid_marker_genes.dotplot.pdf", height = 4.5, width = 10)
}

# Gene expression of TFs in myeloid cells
# TFs: MAF, NR1H3, ETV5, KLF9
# Extended Data Figure 5d
if(T){
  # subset
  seu_subset = subset( seu, cells  = rownames(seu@meta.data)[seu$major_cluster == "Myeloid"] )
  seu_subset$minor_cluster = factor(seu_subset$minor_cluster, 
                                    levels = seu_subset$minor_cluster %>% droplevels() %>% levels %>% rev)
  # plot
  DotPlot(seu_subset, 
          dot.scale = 8,
          features = c("MAF","NR1H3","ETV5","KLF9"), group.by = "minor_cluster") +
    xlab("") + ylab("") +  
    scale_color_gradient2() + 
    theme_classic2() + 
    theme(
      legend.position = "right",
      axis.text.x = element_text(family = "ArialMT", face = "italic",size = 16,  color = "black", angle = 90, hjust = 1),
      #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "italic", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    )
  
  # ggsave("myeloid_regulon.dotplot.pdf", width = 6.5, height = 4.5)
}

# Selected marker genes of non-immune cells
# Extended Data Figure 6e
if(T){
  # subset
  seu_subset = seu
  seu_subset$minor_cluster = factor(seu_subset$minor_cluster, 
                                    levels = seu_subset$minor_cluster %>% droplevels() %>% levels %>% rev)
  
  genes = c("CCL23","CCL20","ICAM2",
            "DUOXA2","DUOX2","PI3",
            "IL11","IL24","IL13RA2","CXCL1","CXCL6",
            "CCL19","CCL21","GREM1","C3",
            "RGS5","RGS16","MEF2C",
            "SELP","RAB3C","HLA-DQA1")
  
  # plot
  Idents(seu_subset) = seu_subset$minor_cluster
  DotPlot(seu_subset, dot.scale = 8, idents = c("DUOX2+ epithelial","M-like cell","Inflammatory fibroblast","Reticular fibroblast","Immature pericyte","Adult venous EC (SELE+)","Cycling EC"),
          features = genes, group.by = "minor_cluster") +
    xlab("") + ylab("") +  
    scale_color_gradient2() + 
    theme_classic2() + 
    theme(
      legend.position = "right",
      axis.text.x = element_text(family = "ArialMT", face = "italic",size = 16,  color = "black", angle = 90, hjust = 1),
      #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "italic", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    )
  
  # save plot
  ggsave("selected_subset.marker_genes.dotplot.pdf", height = 4, width = 9.5)
}

# Gene expression of MHC-II in epithelial cells
# Extended Data Figure 10b
if(T){
  # subset
  seu_subset = subset( seu, cells  = rownames(seu@meta.data)[seu$major_cluster == "Epithelial"] )
  seu_subset$minor_cluster = factor(seu_subset$minor_cluster, 
                                    levels = seu_subset$minor_cluster %>% droplevels() %>% levels %>% rev)
  
  # genes
  hla2genes = rownames(seu)[ grepl("HLA-",rownames(seu))] %>% sort
  hla2genes = hla2genes[4:17] %>% sort(decreasing = T)
  
  # plot
  DotPlot(seu_subset, 
          dot.scale = 8,
          features = hla2genes, group.by = "minor_cluster") +
    xlab("") + ylab("") +  
    scale_color_gradient2() + 
    theme_classic2() + 
    theme(
      legend.position = "right",
      axis.text.x = element_text(family = "ArialMT", face = "italic",size = 16,  color = "black", angle = 90, hjust = 1),
      #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "italic", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    )
  
  # save
  # ggsave("hla2_epi.dotplot.pdf", height = 5, width = 10)
}
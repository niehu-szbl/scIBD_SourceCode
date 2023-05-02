library(Seurat)
library(dplyr)
library(ggplot2)
library(ggsci)
library(ggpubr)

# prepare data
if(T){
  # load data
  seu = readRDS("./AvivRegev_2019cell_uc.raw.rds")
  seu@meta.data %>% colnames()
  colnames(seu@meta.data)[4] = "cell_id"
  
  # check data
  seu@assays$RNA@data[1,] %>% sum
  
  # Normalize Data
  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
  
  # Cluster
  result = seu@meta.data %>% group_by(Cluster) %>% summarise(count=n())
  result
  
  # Disease status
  result = seu@meta.data %>% group_by(Health) %>% summarise(count=n())
  result
}

# violin plot
# Expression of LAMP3, CCR7
# group by cluster
# Extended Data Figure 4d
if(T){
  
  # load myeloid
  imm = read.table("/home/niehu/niehu/ibd_public_data_20210821/AvivRegev_2019Cell_uc/download/SCP259/metadata/all.meta2.txt", header = T, stringsAsFactors = F, sep = "\t")
  myeloid_cells = imm[imm$Cluster %in% c("Macrophages","DC1","DC2","CD69+ Mast","CD69- Mast","Inflammatory Monocytes","Cycling Monocytes"),]$NAME
  
  genes = c("LAMP3","CCR7")
  cells = myeloid_cells
  
  title = "Gene expression"
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% t %>% as.data.frame()
  
  rownames(gex) = rownames(seu@meta.data[cells, ])
  gex$cell_id = seu@meta.data[cells, ]$cell_id
  gex = reshape2::melt(gex)
  gex = left_join(gex, seu@meta.data[,c("cell_id","Health",'Cluster')])
  gex %>% head
  
  colnames(gex)[2:3] = c("gene","exp")
  
  ggplot(gex, aes(x = Cluster, y = log2(exp+0.1), fill = Cluster, color = Cluster)) +
    facet_wrap(~gene, scales = "free_y", ncol = 4) + 
    geom_violin(trim = TRUE, scale = "width") + 
    scale_fill_d3(palette = 'category10')+ 
    scale_color_d3(palette = 'category10') + 
    stat_summary(fun = "mean", geom = "crossbar", width = 0.3, colour = "lightblue4") + 
    theme_classic() +
    labs(x="",y="Gene expression\nlog2(Norm + 0.1)", title =title) + 
    theme(legend.position = "right",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    )
  #ggsave("LAMP3_CCR7_in_Myeloid.SCP259.pdf", width = 8, height = 3)
}

# violin plot
# Expression of DUOX2, DUOXA2
# group by disease status
# Extended Data Figure 4e
if(T){

  epi = read.table("/home/niehu/niehu/ibd_public_data_20210821/AvivRegev_2019Cell_uc/download/SCP259/cluster/Epi.tsne.txt", header = F, stringsAsFactors = F, sep = "\t")
  epi_cells = epi$V1

  genes = c("DUOX2","DUOXA2")
  cells = epi_cells
  title = ""
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% t %>% as.data.frame()
  rownames(gex) = rownames(seu@meta.data[cells, ])
  gex$cell_id = seu@meta.data[cells, ]$cell_id
  gex = reshape2::melt(gex)
  gex = left_join(gex, seu@meta.data[,c("cell_id","Health",'Cluster')])
  gex %>% head
  
  colnames(gex)[2:3] = c("gene","exp")
  
  gex = gex %>% 
    filter(Cluster == 'Enterocytes')
  
  my_comparisons <- list( c("Inflamed", "Healthy"), c("Non-inflamed", "Healthy"),
                          c("Inflamed","Non-inflamed"))
  
  ggplot(gex, aes(x = Health, y = log2(exp+0.1), fill = Health, color = Health)) +
    facet_wrap(~gene, scales = "free_y", ncol = 4) + 
    geom_violin(trim = TRUE, scale = "width") + 
    scale_fill_d3(palette = 'category10')+ 
    scale_color_d3(palette = 'category10') + 
    stat_summary(fun = "mean", geom = "crossbar", width = 0.3, colour = "lightblue4") + 
    theme_classic() +
    labs(x="",y="Gene expression\nlog2(Norm + 0.1)", title =title) + 
    theme(legend.position = "right",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")#, label = "p.signif") 
  # ggsave("DUOX2_DUOXA2_in_Enterocytes.SCP259.pdf", width = 7, height = 4)
  
  # summary
  gex %>% colnames()
  gex$Health %>% table
  gex %>% group_by(Health) %>% count()
  
  # statistic test
  wilcox.test(gex[gex$gene == "DUOX2" & gex$Health == "Healthy", ]$exp,
              gex[gex$gene == "DUOX2" & gex$Health == "Inflamed", ]$exp)
}

# violin plot
# Expression of DUOX2, DUOXA2
# group by cluster
# Extended Data Figure 4f
if(T){
  
  epi = read.table("/home/niehu/niehu/ibd_public_data_20210821/AvivRegev_2019Cell_uc/download/SCP259/cluster/Epi.tsne.txt", header = F, stringsAsFactors = F, sep = "\t")
  epi_cells = epi$V1

  genes = c("DUOX2","DUOXA2")
  cells = epi_cells
  title = "Gene expression"
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% t %>% as.data.frame()
  rownames(gex) = rownames(seu@meta.data[cells, ])
  gex$cell_id = seu@meta.data[cells, ]$cell_id
  gex = reshape2::melt(gex)
  gex = left_join(gex, seu@meta.data[,c("cell_id","Health",'Cluster')])
  gex %>% head
  
  colnames(gex)[2:3] = c("gene","exp")
  
  ggplot(gex, aes(x = Cluster, y = log2(exp+0.1), fill = Cluster, color = Cluster)) +
    facet_wrap(~gene, scales = "free_y", ncol = 4) + 
    geom_violin(trim = TRUE, scale = "width") + 
    scale_fill_d3(palette = 'category20b')+ 
    scale_color_d3(palette = 'category20b') + 
    stat_summary(fun = "mean", geom = "crossbar", width = 0.3, colour = "lightblue4") + 
    theme_classic() +
    labs(x="",y="Gene expression\nlog2(Norm + 0.1)", title =title) + 
    theme(legend.position = "right",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    )
  # ggsave("DUOX2_DUOXA2_in_Epi.SCP259.pdf", width = 10, height = 4.3)
}

# violin plot
# Tc17 markers
# group by cluster
# Extended Data Figure 4g
if(T){
  
  # load myeloid
  imm = read.table("/home/niehu/niehu/ibd_public_data_20210821/AvivRegev_2019Cell_uc/download/SCP259/metadata/all.meta2.txt", header = T, stringsAsFactors = F, sep = "\t")
  cd8t = c("CD8+ IELs","CD8+ IL17+","CD8+ LP")
  
  genes = c("IL17A","IL22","IL26","RORC")
  cells = imm[imm$Cluster %in% cd8t,]$NAME
  
  title = "Gene expression"
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% t %>% as.data.frame()
  rownames(gex) = rownames(seu@meta.data[cells, ])
  gex$cell_id = seu@meta.data[cells, ]$cell_id
  gex = reshape2::melt(gex)
  gex = left_join(gex, seu@meta.data[,c("cell_id","Health",'Cluster')])
  gex %>% head
  
  colnames(gex)[2:3] = c("gene","exp")
  
  ggplot(gex, aes(x = Cluster, y = log2(exp+0.1), fill = Cluster, color = Cluster)) +
    facet_wrap(~gene, scales = "free_y", ncol = 4) + 
    geom_violin(trim = TRUE, scale = "width") + 
    scale_fill_d3(palette = 'category10')+ 
    scale_color_d3(palette = 'category10') + 
    stat_summary(fun = "mean", geom = "crossbar", width = 0.3, colour = "lightblue4") + 
    theme_classic() +
    labs(x="",y="Gene expression\nlog2(Norm + 0.1)", title =title) + 
    theme(legend.position = "bottom",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    )
  # ggsave("Tc17_markers_in_CD8T.SCP259.pdf", width = 6, height = 3)
}

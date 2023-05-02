library(dplyr)
library(lisi)
library(patchwork)

setwd("/data1/niehu/ibd_public_data_20210821/analysis_20220111/02.integrate/combine_clean_V3/integrate_without_batch_correction/")

# load color
color.minor_cluster = read.table("../color_set/major_cluster_color.txt", header = F, stringsAsFactors = F, sep = "\t", comment.char = "!")
colnames(color.minor_cluster) = c("minor_cluster","color")
color.minor_cluster %>% head
rownames(color.minor_cluster) = color.minor_cluster$minor_cluster
color.minor_cluster$minor_cluster = NULL

# load meta
meta = readRDS("../update_cell_type/all.clean.update_metadata.rds")
colnames(meta) 
minor_clusters = levels(meta$minor_cluster)
studies = levels(meta$study)

# plot umap
source("/data1/niehu/soft/script/color_palette_89.R")
color.study = color_palette_89[ seq( from = 1,to = length(color_palette_89), 
                                 length.out = seu$study %>% droplevels %>% levels %>% length) %>% as.integer ]
names(color.study) = seu$study %>% levels()

# load seu
#seu = readRDS("../update_cell_type/all.clean.raw.rds")

# CD8T
if(T){
  
  # set cluster
  cluster = 'CD8T'
  
  # load data
  obs= read.table("./CD8T.no_batch_correction.obs.csv", header = T, sep = ",", stringsAsFactors = F, row.names = 1)
  head(obs)
  
  cd8t_minor_clusters = minor_clusters[minor_clusters %in% obs$minor_cluster]
  obs$minor_cluster = factor(obs$minor_cluster, levels = cd8t_minor_clusters)
  obs$study = factor(obs$study, levels = studies)
  mycolor = color.minor_cluster[cd8t_minor_clusters,]
  
  p1 = ggplot(obs, aes(x = UMAP_1, y = UMAP_2, color = study)) + 
    geom_point(size= 0.1, alpha = 0.7) + 
    scale_color_manual(values = color.study[obs$study %>% levels]) + 
    #labs(x="UMAP_1",y="UMAP_2", title =paste0(cluster," cells")) + 
    guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))+
    theme_classic2() + 
    theme(
      legend.position = "none",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "plain", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) + coord_fixed()
  
  p2 = ggplot(meta[obs %>% rownames, ], aes(x = UMAP_1, y = UMAP_2, color = study)) + 
    geom_point(size= 0.1, alpha = 0.7) + 
    scale_color_manual(values = color.study[obs$study %>% levels]) +
    #labs(x="UMAP_1",y="UMAP_2", title =paste0(cluster," cells")) + 
    guides(colour = guide_legend(override.aes = list(size=1, alpha = 1)))+
    theme_classic2() + 
    theme(
      legend.position = "none",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "plain", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) + coord_fixed()
  
  p3 = ggplot(obs, aes(x = UMAP_1, y = UMAP_2, color = minor_cluster)) + 
    geom_point(size= 0.1, alpha = 0.7) + 
    scale_color_manual(values = mycolor) + 
    #labs(x="UMAP_1",y="UMAP_2", title =paste0(cluster," cells")) + 
    guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))+
    theme_classic2() + 
    theme(
      legend.position = "none",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "plain", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) + coord_fixed()
  
  
  p4 = ggplot(meta[obs %>% rownames, ], aes(x = UMAP_1, y = UMAP_2, color = minor_cluster)) + 
    geom_point(size= 0.1, alpha = 0.7) + 
    scale_color_manual(values = mycolor) + 
    #labs(x="UMAP_1",y="UMAP_2", title =paste0(cluster," cells")) + 
    guides(colour = guide_legend(override.aes = list(size=1, alpha = 1)))+
    theme_classic2() + 
    theme(
      legend.position = "none",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "plain", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) + coord_fixed()
  #ggsave(paste0(cluster,"_umap.png"), width = 5, height = 5, device = "png",)
  p = (p1 | p2) / (p3 | p4)
  ggsave(filename = "cd8t.png", width = 8.5, height = 4.5, dpi = 600, plot = p, device = "png")
  
  
  p5 = ggplot(obs, aes(x = UMAP_1, y = UMAP_2, color = study)) + 
    geom_point(size= 0.1, alpha = 0) + 
    scale_color_manual(values = color.study[obs$study %>% levels]) + 
    #labs(x="UMAP_1",y="UMAP_2", title =paste0(cluster," cells")) + 
    guides(colour = guide_legend(override.aes = list(size=2, alpha = 1)))+
    theme_classic2() + 
    theme(
      legend.position = "right",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "plain", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) #+ coord_fixed()
  
  p6 = ggplot(meta[obs %>% rownames, ], aes(x = UMAP_1, y = UMAP_2, color = minor_cluster)) + 
    geom_point(size= 0.1, alpha = 0) + 
    scale_color_manual(values = mycolor) + 
    #labs(x="UMAP_1",y="UMAP_2", title =paste0(cluster," cells")) + 
    guides(colour = guide_legend(override.aes = list(size=2, alpha = 1)))+
    theme_classic2() + 
    theme(
      legend.position = "right",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "plain", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) #+ coord_fixed()
  
  p = p5 | p6
  ggsave(filename = "cd8t_blank.pdf", width = 8, height = 3.5, dpi = 600, plot = p, device = "pdf")
}

# Epi
if(T){
  
  # set cluster
  cluster = 'Epithelial'
  
  # load data
  obs= read.table("./Epi.no_batch_correction.obs.csv", header = T, sep = ",", stringsAsFactors = F, row.names = 1)
  head(obs)
  
  epi_minor_clusters = minor_clusters[minor_clusters %in% obs$minor_cluster]
  obs$minor_cluster = factor(obs$minor_cluster, levels = epi_minor_clusters)
  
  obs$study = factor(obs$study, levels = studies)
  mycolor = color.minor_cluster[epi_minor_clusters,]
  
  p1 = ggplot(obs, aes(x = UMAP_1, y = UMAP_2, color = study)) + 
    geom_point(size= 0.1, alpha = 0.7) + 
    scale_color_manual(values = color.study[obs$study %>% levels]) + 
    #labs(x="UMAP_1",y="UMAP_2", title =paste0(cluster," cells")) + 
    guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))+
    theme_classic2() + 
    theme(
      legend.position = "none",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "plain", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) + coord_fixed()
  
  p2 = ggplot(meta[obs %>% rownames, ], aes(x = UMAP_1, y = UMAP_2, color = study)) + 
    geom_point(size= 0.1, alpha = 0.7) + 
    scale_color_manual(values = color.study[obs$study %>% levels]) +
    #labs(x="UMAP_1",y="UMAP_2", title =paste0(cluster," cells")) + 
    guides(colour = guide_legend(override.aes = list(size=1, alpha = 1)))+
    theme_classic2() + 
    theme(
      legend.position = "none",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "plain", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) + coord_fixed()
  
  p3 = ggplot(obs, aes(x = UMAP_1, y = UMAP_2, color = minor_cluster)) + 
    geom_point(size= 0.1, alpha = 0.7) + 
    scale_color_manual(values = mycolor) + 
    #labs(x="UMAP_1",y="UMAP_2", title =paste0(cluster," cells")) + 
    guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))+
    theme_classic2() + 
    theme(
      legend.position = "none",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "plain", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) + coord_fixed()
  
  
  p4 = ggplot(meta[obs %>% rownames, ], aes(x = UMAP_1, y = UMAP_2, color = minor_cluster)) + 
    geom_point(size= 0.1, alpha = 0.7) + 
    scale_color_manual(values = mycolor) + 
    #labs(x="UMAP_1",y="UMAP_2", title =paste0(cluster," cells")) + 
    guides(colour = guide_legend(override.aes = list(size=1, alpha = 1)))+
    theme_classic2() + 
    theme(
      legend.position = "none",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "plain", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) + coord_fixed()
  #ggsave(paste0(cluster,"_umap.png"), width = 5, height = 5, device = "png",)
  p = (p1 | p2) / (p3 | p4)
  ggsave(filename = "epi.png", width = 8.5, height = 4.5, dpi = 600, plot = p, device = "png")
  
  
  p5 = ggplot(obs, aes(x = UMAP_1, y = UMAP_2, color = study)) + 
    geom_point(size= 0.1, alpha = 0) + 
    scale_color_manual(values = color.study[obs$study %>% levels]) + 
    #labs(x="UMAP_1",y="UMAP_2", title =paste0(cluster," cells")) + 
    guides(colour = guide_legend(override.aes = list(size=2, alpha = 1)))+
    theme_classic2() + 
    theme(
      legend.position = "right",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "plain", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) #+ coord_fixed()
  
  p6 = ggplot(meta[obs %>% rownames, ], aes(x = UMAP_1, y = UMAP_2, color = minor_cluster)) + 
    geom_point(size= 0.1, alpha = 0) + 
    scale_color_manual(values = mycolor) + 
    #labs(x="UMAP_1",y="UMAP_2", title =paste0(cluster," cells")) + 
    guides(colour = guide_legend(override.aes = list(size=2, alpha = 1)))+
    theme_classic2() + 
    theme(
      legend.position = "right",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "plain", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) #+ coord_fixed()
  
  p = p5 | p6
  ggsave(filename = "epi_blank.pdf", width = 8, height = 3.5, dpi = 600, plot = p, device = "pdf")
}

# LISI
if(T){
  obs= read.table("./CD8T.no_batch_correction.obs.csv", header = T, sep = ",", stringsAsFactors = F, row.names = 1)
  res1 <- compute_lisi(obs[,c("UMAP_1","UMAP_2")], obs, c('study'))
  
  after = meta[obs %>% rownames(),c("UMAP_1","UMAP_2","study")]
  res2 <- compute_lisi(after[,c("UMAP_1","UMAP_2")], after, c('study'))
  
  res3 = cbind(res1, res2)
  colnames(res3) = c("before","after")
  res3$group = "CD8T"
  
  obs= read.table("./Epi.no_batch_correction.obs.csv", header = T, sep = ",", stringsAsFactors = F, row.names = 1)
  res4 <- compute_lisi(obs[,c("UMAP_1","UMAP_2")], obs, c('study'))
  
  after = meta[obs %>% rownames(),c("UMAP_1","UMAP_2","study")]
  res5 <- compute_lisi(after[,c("UMAP_1","UMAP_2")], after, c('study'))
  
  res6 = cbind(res4, res5)
  colnames(res6) = c("before","after")
  res6$group = "Epithelial"
  
  res = rbind(res3, res6)
  res = reshape2::melt(res)
  head(res)
  
  # summary number of cells per group
  res %>% group_by(group, variable) %>% count() # 184220, 195735
  
  # plot
  my_comparisons  = list( c("before","after"))
  ggplot(res, aes(x = variable, y = value, fill = variable)) + 
    facet_grid(~group) +
    geom_boxplot(outlier.size = 0.01) +
    labs(x="", y= "LISI", title = "") + 
    scale_fill_d3(alpha = 0.6, palette = "category20") + 
    theme_classic2() + 
    theme(
      legend.position = "none",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black", angle = 45, hjust = 1),
      #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "plain", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) + 
    stat_compare_means(comparisons = my_comparisons, 
                       method = "wilcox.test") #label = "p.signif")
  ggsave("batch.pdf", width = 3, height = 4)
  
  wilcox.test(res[res$group == "Epithelial" & res$variable == "before",]$value, 
              res[res$group == "Epithelial" & res$variable == "after",]$value)
  
  wilcox.test(res[res$group == "CD8T" & res$variable == "before",]$value, 
              res[res$group == "CD8T" & res$variable == "after",]$value)
}
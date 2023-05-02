library(ggplot2)
library(dplyr)

# prepare data
if(T){
  # load data
  meta = readRDS("all.clean.meta.rds")
  colnames(meta)
  head(meta)
  
  # random order row
  set.seed(100)
  meta = meta[sample(1:nrow(meta), size = nrow(meta)), ]
  
  # set color
  major_color = c("#1C4E97","#416EB0","#6899CB","#44A052","#7EB97B",
                  "#BAD399","#E37E27","#D44A20","#C31820")
  names(major_color) = meta$major_cluster %>% levels
}

# umap plot
# major cluster
# Figure 1a
if(T){
  
  # input data
  umap_data = meta[,c("gUMAP_1","gUMAP_2","major_cluster")]
  colnames(umap_data)
  head(umap_data)
  
  colnames(umap_data) = c("UMAP_1","UMAP_2","label")
  ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = label)) +
    geom_point(size = 0.01) +
    scale_color_manual("Label", values = major_color[umap_data$label %>% levels ]) + 
    xlab("UMAP_1") +
    ylab("UMAP_2") + 
    xlim( min(umap_data$UMAP_1), max(umap_data$UMAP_1) ) + 
    ylim( min(umap_data$UMAP_2), max(umap_data$UMAP_2) ) + 
    guides(color = guide_legend(override.aes = list(size = 3))) + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
          axis.text.y = element_text(size=10, family="Times", color = "black"),
          axis.title.x = element_text(size=12, family="Times", color = "black"),
          axis.title.y = element_text(size=12, family="Times", color = "black"),     
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.position = "right") +
    theme(axis.line = element_line(color = "black"))
  # ggsave( "scIBD_major_cluster_umap.png", height = 6, width = 8, dpi = 300)
}

# umap plot
# tissue
# Figure 1b
if(T){
  
  # input data
  umap_data = meta[,c("gUMAP_1","gUMAP_2","tissue.sub")]
  colnames(umap_data)
  head(umap_data)
  
  # set color
  mycolor = c("#D92220","#496D9D","#3FA35E","#7B6C82","#C86549","#F5B222","#D8C131","#B25C47","#DC83B1","#999999")
  names(mycolor) = c("appendix","blood","cecum","colon","duodenum","hindgut","ileum","jejunum","lymph node","rectum")
  
  colnames(umap_data) = c("UMAP_1","UMAP_2","label")
  umap_data$label = factor(umap_data$label, levels = names(mycolor))
  
  # plot
  ggplot(umap_data, aes(x = UMAP_1, y = UMAP_2, color = label)) +
    geom_point(size = 0.01) +
    scale_color_manual("Label", values = mycolor) + 
    xlab("UMAP_1") +
    ylab("UMAP_2") + 
    xlim( min(umap_data$UMAP_1), max(umap_data$UMAP_1) ) + 
    ylim( min(umap_data$UMAP_2), max(umap_data$UMAP_2) ) + 
    guides(color = guide_legend(override.aes = list(size = 3))) + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
          axis.text.y = element_text(size=10, family="Times", color = "black"),
          axis.title.x = element_text(size=12, family="Times", color = "black"),
          axis.title.y = element_text(size=12, family="Times", color = "black"),     
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.position = "right") +
    theme(axis.line = element_line(color = "black"))
  # ggsave( "scIBD_tissue.sub_umap.png", height = 6, width = 8, dpi = 300)
}

# umap plot
# disease
# Figure 1c
if(T){
  
  umap_data = meta[, c("gUMAP_1","gUMAP_2","disease")]
  colnames(umap_data)
  head(umap_data)

  #set color
  mycolor = c("#7CC17C","#AFB3C0","#DCB6A4","#F6D08C",
              "#E2E89C","#4B78A4","#A62786","#CD2D4D","#A65B29","#656665")
  names(mycolor) = c("CD_PBMC","CD_inflamed","CD_inflamed/CD_non_inflamed",
                     "CD_non_inflamed","Colitis_inflamed","Healthy","Healthy_PBMC",
                     "UC_PBMC","UC_inflamed","UC_non_inflamed") 
  
  colnames(umap_data) = c("gUMAP_1","gUMAP_2","label")
  umap_data$label = factor(umap_data$label, levels = names(mycolor))
  
  # plot
  ggplot(umap_data, aes(x = gUMAP_1, y = gUMAP_2, color = label)) +
    geom_point(size = 0.05) +
    scale_color_manual("Label", values = mycolor) +
    xlab("UMAP_1") +
    ylab("UMAP_2") + 
    xlim( min(seu$gUMAP_1), max(seu$gUMAP_1) ) + 
    ylim( min(seu$gUMAP_2), max(seu$gUMAP_2) ) + 
    guides(color = guide_legend(override.aes = list(size = 3))) + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
          axis.text.y = element_text(size=10, family="Times", color = "black"),
          axis.title.x = element_text(size=12, family="Times", color = "black"),
          axis.title.y = element_text(size=12, family="Times", color = "black"),     
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.position = "right") +  # show legend on right
    theme(axis.line = element_line(color = "black"))
  # ggsave( paste0("landscape_umpa_by_disease.png"), device = "png", height = 5, width = 6)
}

# Percentage of cells in each major cluster in scIBD
# Extended Data Figure 2a
if(T){
  
  # summary
  result = meta %>% dplyr::group_by(major_cluster) %>% 
    dplyr::summarise(count = n()) %>% 
    dplyr::mutate(pct = count/sum(count))
  head(result)
  
  # plot
  ggplot(result, aes(x = major_cluster, y = pct, fill = major_cluster)) + 
    geom_bar(stat = "identity") + 
    xlab("") + 
    ylab("Percentage of cells") + 
    theme_classic() + 
    scale_fill_manual(values = major_color) + 
    theme(axis.text.x = element_text(size = 10, color = "black", angle = 45, vjust = 0.5, hjust = 1))
  # ggsave( "Cell_composition_of_major_cluster.pdf", height = 4, width = 6)
}

# Cell composition by stage
# Extended Data Figure 2b
if(T){
  
  # summary
  result = meta %>% dplyr::group_by(major_cluster, stage) %>% 
    dplyr::summarise(count = n()) %>% group_by(stage) %>% 
    dplyr::mutate(pct = count/sum(count))
  head(result)
  
  # plot
  ggplot(result, aes(x = stage, y = pct, fill = major_cluster)) + 
    geom_bar(stat = "identity", position = position_stack()) + 
    xlab("") + 
    ylab("Percentage of cells") + 
    theme_classic() + 
    scale_fill_manual(values = major_color) + 
    theme(axis.text.x = element_text(size = 10, color = "black", angle = 45, vjust = 1, hjust = 1))
  # ggsave( "Cell_composition_of_major_cluster_by_stage.pdf", height = 4, width = 4)
}

# umap plot
# development stage
# Extended Data Figure 2c
if(T){
  
  umap_data = meta[, c("gUMAP_1","gUMAP_2","stage")]
  colnames(umap_data)
  head(umap_data)
  
  #set color
  mycolor = c("#7AC07B","#99B2A2","#666666")
  names(mycolor) = c("fetal","pediatric","adult") 
  
  colnames(umap_data) = c("gUMAP_1","gUMAP_2","label")
  umap_data$label = factor(umap_data$label, levels = names(mycolor))
  
  # plot
  ggplot(umap_data, aes(x = gUMAP_1, y = gUMAP_2, color = label)) +
    geom_point(size = 0.05) +
    scale_color_manual("Label", values = mycolor) +
    xlab("UMAP_1") +
    ylab("UMAP_2") + 
    xlim( min(seu$gUMAP_1), max(seu$gUMAP_1) ) + 
    ylim( min(seu$gUMAP_2), max(seu$gUMAP_2) ) + 
    guides(color = guide_legend(override.aes = list(size = 3))) + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
          axis.text.y = element_text(size=10, family="Times", color = "black"),
          axis.title.x = element_text(size=12, family="Times", color = "black"),
          axis.title.y = element_text(size=12, family="Times", color = "black"),     
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.position = "right") +  # show legend on right
    theme(axis.line = element_line(color = "black"))
  # ggsave( paste0("landscape_umpa_by_stage.png"), device = "png", height = 5, width = 6)
}

# umap plot
# study or dataset
# Extended Data Figure 2d
if(T){
  
  umap_data = meta[, c("gUMAP_1","gUMAP_2","study")]
  colnames(umap_data)
  head(umap_data)
  
  #set color
  mycolor = c("#1D2E61","#2773AB","#62AD6F","#706EA6",
              "#A06EA7","#AF488C","#C74A4A","#DAE5AE",
              "#E2A444","#E05070","#F5DCD2","#F8EF9A")
  names(mycolor) = c("K. Parikh et al., Nature, 2019",
                     "Kinchen et al., Cell, 2018",
                     "D. Corridoni et al., Nat Med, 2020",
                     "D. Fawkner−Corbett et al., Cell, 2021",
                     "C. S. Smillie et al., Cell, 2019",
                     "R. Elmentaite et al., Nature, 2021",
                     "R. Elmentaite et al., Dev Cell, 2020",
                     "B. Huang et al., Cell, 2019",
                     "J. C. Martin et al., Cell, 2019",
                     "N. Jaeger et al., Nat Commun, 2021",
                     "B. S. Boland et al., Sci Immunol, 2020",
                     "M. Friedrich et al., Nat Med, 2021")
  
  colnames(umap_data) = c("gUMAP_1","gUMAP_2","label")
  umap_data$label = factor(umap_data$label, levels = names(mycolor))
  
  # plot
  ggplot(umap_data, aes(x = gUMAP_1, y = gUMAP_2, color = label)) +
    geom_point(size = 0.05) +
    scale_color_manual("Label", values = mycolor) +
    xlab("UMAP_1") +
    ylab("UMAP_2") + 
    xlim( min(seu$gUMAP_1), max(seu$gUMAP_1) ) + 
    ylim( min(seu$gUMAP_2), max(seu$gUMAP_2) ) + 
    guides(color = guide_legend(override.aes = list(size = 3))) + 
    theme_bw() +
    theme(axis.text.x = element_text(size = 10, family="Times", color = "black"),
          axis.text.y = element_text(size=10, family="Times", color = "black"),
          axis.title.x = element_text(size=12, family="Times", color = "black"),
          axis.title.y = element_text(size=12, family="Times", color = "black"),     
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_blank(),
          legend.position = "right") +  # show legend on right
    theme(axis.line = element_line(color = "black"))
  # ggsave( paste0("landscape_umpa_by_stdy.png"), device = "png", height = 5, width = 6)
}

# umap plot
# all cell subtypes of each major cluter
# Extended Data Figure 3d
if(T){
  # get major cluster
  major_clusters = meta$major_cluster %>% levels
  major_clusters

  plots1 = list()
  for(i in 1:length(major_clusters)){
    
    # subset
    meta.sub = meta %>% dplyr::filter(major_cluster == major_clusters[i])
    meta.sub$minor_cluster = meta.sub$minor_cluster %>% droplevels()
    
    # set color
    clusters = meta.sub$minor_cluster %>% levels()
    mycolor = color.minor_cluster[clusters,]
    
    p1 = ggplot(meta.sub, aes(x = UMAP_1, y = UMAP_2, color = minor_cluster)) + 
      geom_point(size= 0.1, alpha = 0.7) + 
      scale_color_manual(values = mycolor) + 
      labs(x="UMAP_1",y="UMAP_2", title = major_clusters[i]) + 
      guides(colour = guide_legend(override.aes = list(size=3, alpha = 1)))+
      theme_classic2() + 
      theme(
        legend.position = "none",
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),#element_text(family = "ArialMT", size = 0,  color = "black"),
        axis.text.y = element_blank(),#element_text(family = "", size = 0, color = 'black'),
        axis.title.x = element_blank(),#element_text(family = "ArialMT",size = 0, color = 'black'),
        axis.title.y = element_blank(),#element_text(family = "ArialMT",size = 0, color = 'black'),
        legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
        legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
        strip.text = element_text(family = "ArialMT",face = "plain", size = 16, color ='black'),
        strip.background = element_blank(),
        plot.title = element_text(family = "ArialMT",size = 14, color ='black', hjust = 0.5)
      ) #+ coord_fixed()
    ggsave(paste0(major_clusters[i],"_umap.png"), width = 5, height = 5, device = "png", dpi = 600, plot = p1)
    plots1[[i]] = p1

  }
  names(plots1) = major_clusters

  png("nine_major_cluster.png", width = 700, height = 700, units = "px")
  (plots1[[1]]  | plots1[[2]] | plots1[[3]]) / (plots1[[4]]  | plots1[[5]] | plots1[[6]]) / (plots1[[7]]  | plots1[[8]] | plots1[[9]]) 
  dev.off()
}


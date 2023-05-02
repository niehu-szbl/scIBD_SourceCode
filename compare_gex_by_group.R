library(dplyr)
library(coin)
library(ggplot)
library(paletteer)

# prepare data
if(T){
  # load data
  seu = readRDS("all.clean.raw.rds")

  # load meta data
  meta = readRDS("all.clean.meta.rds")
  seu@meta.data = meta
}

# Extended Data Figure 7b
# NOS2 in DUOX2+ epithelial (HC VS UC)
if(T){
  genes = c("NOS2")
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% c("DUOX2+ epithelial")) %>% 
    dplyr::filter(disease %in% c("Healthy","UC_inflamed")) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = "NOS2 in DUOX2+ epithelial"
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  my_comparisons = list( c("UC_inflamed","Healthy"))
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin(scale = "width") +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif")
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  
  # save plot
  # ggsave("NOS2_duox2_epi.pdf", width = 2, height = 3)  
}

# FERMT1 in M−like cell (HC VS UC)
if(T){
  genes = c("FERMT1")
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% c("M-like cell")) %>% 
    dplyr::filter(disease %in% c("Healthy","UC_inflamed")) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = "FERMT1 in M−like cell"
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  my_comparisons = list( c("UC_inflamed","Healthy"))
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin() +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") 
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  
  # save plot
  # ggsave("FERMT1_M_cell.pdf", width = 2, height = 3) 
}

# PLCG2 in Tuft (HC VS UC)
if(T){
  genes = c("PLCG2")
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% c("Tuft")) %>% 
    dplyr::filter(disease %in% c("Healthy","UC_inflamed")) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = "PLCG2 in Tuft"
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  my_comparisons = list( c("UC_inflamed","Healthy"))
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin() +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") 
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  # save plot
  # ggsave("PLCG2_tuft.pdf", width = 2, height = 3)
}

# ITLN1 in Goblet (HC VS UC)
if(T){
  genes = c("ITLN1")
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% c("Goblet")) %>% 
    dplyr::filter(disease %in% c("Healthy","UC_inflamed","CD_inflamed")) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = "ITLN1 in Goblet"
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  my_comparisons = list( c("UC_inflamed","Healthy"),c("CD_inflamed","Healthy"),c("UC_inflamed","CD_inflamed"))
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin() +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1, 2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1, 2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") 
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "CD_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "CD_inflamed",]$gex)$p.value
  
  # save plot
  # ggsave("ITLN1_goblet.pdf", width = 2, height = 3) 
}

# PROX1 in LEC (HC vs CD)
if(T){
  genes = c("PROX1")
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% c("LEC")) %>% 
    dplyr::filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = "PROX1 in LEC"
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  my_comparisons = list( c("CD_inflamed","Healthy"),
                         c("UC_inflamed","Healthy"),
                         c("CD_inflamed","UC_inflamed"))
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin() +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") 
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "CD_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "CD_inflamed",]$gex)$p.value
  
  # save plot
  # ggsave("PROX1_LEC.pdf", width = 2.5, height = 3)
}

# LY75 in LAMP3+ DC cells (CD sig)
if(T){
  genes = c("LY75")
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% c("LAMP3+ DC")) %>% 
    dplyr::filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = "LY75 in LAMP3+ DC cells"
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  my_comparisons = list( c("CD_inflamed","Healthy"),
                         c("Healthy","UC_inflamed"),
                         c("CD_inflamed","UC_inflamed"))
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin(scale = "width") +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + 
    stat_compare_means(comparisons = my_comparisons, 
                       method = "wilcox.test", label = "p.signif") 
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "CD_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "CD_inflamed",]$gex)$p.value
  
  # save plot
  # ggsave("LY75_lamp3_dc.pdf", width = 2.5, height = 3)
}

# IL2RA in Treg
if(T){
  genes = c("IL2RA")
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% c("CD4+ Treg")) %>% 
    dplyr::filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = "IL2RA in Treg"
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  my_comparisons = list( c("CD_inflamed","Healthy"),
                         c("UC_inflamed","Healthy"),
                         c("CD_inflamed","UC_inflamed"))
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin(scale = "width") +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") 
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "CD_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "CD_inflamed",]$gex)$p.value
  
  # save plot
  # ggsave("IL2RA_treg.pdf", width = 2.5, height = 3)
}

# EDNRB in Stromal 2
if(T){
  genes = c("EDNRB")
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% c("Stromal 2")) %>% 
    dplyr::filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = "EDNRB in Stromal 2"
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  my_comparisons = list( c("CD_inflamed","Healthy"),
                         c("UC_inflamed","Healthy"),
                         c("CD_inflamed","UC_inflamed"))
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin(scale = "width") +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") 
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "CD_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "CD_inflamed",]$gex)$p.value
  
  # save plot
  # ggsave("EDNRB_Stromal2.pdf", width = 2.5, height = 3)
}

# Extended Data Figure 8b
# PTGS2 in Inflammatory monocyte
if(T){
  genes = c("PTGS2")
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% c("Inflammatory monocyte")) %>% 
    dplyr::filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = "PTGS2 in Inflammatory monocyte"
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  my_comparisons = list( c("CD_inflamed","Healthy"),
                         c("UC_inflamed","Healthy"),
                         c("CD_inflamed","UC_inflamed"))
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin(scale = "width") +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif")
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "CD_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "CD_inflamed",]$gex)$p.value
  
  # save plot
  ggsave("PTGS2_IAF.pdf", width = 2.5, height = 3)
}

# CSF2RA in LAMP3+ DC
if(T){
  genes = c("CSF2RA")
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% c("LAMP3+ DC")) %>% 
    dplyr::filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = "CSF2RA in LAMP3+ DC"
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  my_comparisons = list( c("CD_inflamed","Healthy"),
                         c("UC_inflamed","Healthy"),
                         c("CD_inflamed","UC_inflamed"))
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin(scale = "width") +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif")
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "CD_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "CD_inflamed",]$gex)$p.value
  
  # plot
  ggsave("CSF2RA_gex.pdf", width = 2.5, height = 3)
}

# MMP9 in APOE+ macrophage
if(T){
  genes = c("MMP9")
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% c("APOE+ macrophage")) %>% 
    dplyr::filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = "MMP9 in APOE+ macrophage"
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  my_comparisons = list( c("CD_inflamed","Healthy"),
                         c("UC_inflamed","Healthy"),
                         c("CD_inflamed","UC_inflamed"))
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin(scale = "width") +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif")
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "CD_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "CD_inflamed",]$gex)$p.value
  
  # save plot
  ggsave("MMP9_gex.pdf", width = 2.5, height = 3)
}

# CCL11 in Stromal 1
if(T){
  genes = c("CCL11")
  cluster = c("Stromal 1")
  disease_status = c("CD_inflamed","Healthy","UC_inflamed")
  mycombn = combn(disease_status, m = 2) 
  my_comparisons = c()
  for(i in 1:3){
    my_comparisons[[i]] = c(mycombn[1,i], mycombn[2,i])
  }
  
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% cluster) %>% 
    dplyr::filter(disease %in% disease_status) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = paste0(genes, " in ", cluster)
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin(scale = "width") +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif")
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "CD_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "CD_inflamed",]$gex)$p.value
  
  # save plot
  ggsave(paste0(genes, '_gex.pdf'), width = 2.5, height = 3)
}

# SCN7A in adult glia
if(T){
  genes = c("SCN7A")
  cluster = c("Adult glia")
  disease_status = c("CD_inflamed","Healthy","UC_inflamed")
  mycombn = combn(disease_status, m = 2) 
  my_comparisons = c()
  for(i in 1:3){
    my_comparisons[[i]] = c(mycombn[1,i], mycombn[2,i])
  }
  
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% cluster) %>% 
    dplyr::filter(disease %in% disease_status) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = paste0(genes, " in ", cluster)
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin(scale = "width") +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif")
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "CD_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "CD_inflamed",]$gex)$p.value
  
  # save plot
  ggsave(paste0(genes, '_gex.pdf'), width = 2.5, height = 3)
}

# TPH1 in enteroendocrine
if(T){
  genes = c("TPH1")
  cluster = c("Enteroendocrine")
  disease_status = c("Healthy","UC_inflamed")
  mycombn = combn(disease_status, m = 2) 
  my_comparisons = list(c("Healthy","UC_inflamed"))
  
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% cluster) %>% 
    dplyr::filter(disease %in% disease_status) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = paste0(genes, " in ", cluster)
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin(scale = "width") +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif")
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  
  # save plot
  ggsave(paste0(genes, '_gex.pdf'), width = 2.5, height = 3)
}

# PDE4C in DUOX2+ epithelial
if(T){
  genes = c("PDE4C")
  cluster = c("DUOX2+ epithelial")
  disease_status = c("Healthy","UC_inflamed")
  mycombn = combn(disease_status, m = 2) 
  my_comparisons = list(c("Healthy","UC_inflamed"))
  
  cells = meta %>% 
    dplyr::filter(stage == "adult") %>% 
    dplyr::filter(minor_cluster %in% cluster) %>% 
    dplyr::filter(disease %in% disease_status) %>% 
    dplyr::filter(tissue %in% c("largeInt","smallInt")) %>% rownames()
  cells %>% length
  
  title = paste0(genes, " in ", cluster)
  gex = data.frame(seu@assays$RNA@data[genes,cells]) %>% as.data.frame()
  meta.sub = meta[cells, ]
  meta.sub$gex = gex %>% unlist
  meta.sub$disease %>% droplevels() %>% table
  
  ggplot(meta.sub, aes(x = disease, y = gex, fill = disease, color = disease)) + 
    geom_violin(scale = "width") +
    stat_summary(fun = "mean", geom = "crossbar", width = 0.2, colour = "lightblue4") + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) +
    theme_classic() + 
    labs(x="",y="Relative gene expression", title =title) + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 0,  color = "black", angle = 45, hjust = 1),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black', face = 'italic'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 10, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif")
  
  # summary p value
  colnames(meta.sub)
  meta.sub %>% group_by(disease) %>% count()
  wilcox.test(meta.sub[meta.sub$disease == "UC_inflamed",]$gex, 
              meta.sub[meta.sub$disease == "Healthy",]$gex)$p.value
  
  # save plot
  ggsave(paste0(genes, '_gex.pdf'), width = 2.5, height = 3)
}

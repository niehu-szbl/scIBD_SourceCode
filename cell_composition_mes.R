library("tictoc")
library("ggplot2")
library("ggsci")
library("ggpubr")
library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
library("dplyr")

# load meta data
meta = readRDS("all.clean.meta.rds")

# set meta
seu@meta.data = meta[rownames(seu@meta.data), ]

# Mesenchymal cells
# disease distribution
# Ro/e, heatmap plot
if(T){
  # select data
  data.sub = meta %>%
    filter(major_cluster %in% c("Mesenchymal")) %>%
    filter(stage == "adult") %>% 
    filter(tissue %in% c("smallInt","largeInt"))
  data.sub$minor_cluster = data.sub$minor_cluster %>% droplevels()
  
  # count
  N.o <- table(data.sub[["minor_cluster"]],data.sub[["disease"]])
  N.o = N.o[,colSums(N.o) > 0]
  N.o
  N.o = N.o[!grepl("Fetal", N.o %>% rownames()),]
  N.o = N.o[!grepl("Pediatric", N.o %>% rownames()),]
  
  # order disease
  N.o = N.o[,c("UC_inflamed","Healthy","CD_inflamed","UC_non_inflamed","CD_non_inflamed")]
  
  # filter cluster
  N.o = N.o[which(rowSums(N.o > 20) > 0) %>% names(),]
  N.o %>% rownames()
  
  # order cluster
  subset_order = seu$minor_cluster %>% levels
  N.o = N.o[subset_order[subset_order %in% rownames(N.o)],]
  
  # calculate expected value
  res.chisq <- chisq.test(N.o)
  
  # calculate Ro/e
  R.oe <- (res.chisq$observed)/(res.chisq$expected)
  
  mark = R.oe
  mark[R.oe == 0] = '-'
  mark[R.oe > 0 & R.oe < 0.5] = '+/-'
  mark[R.oe >= 0.5 & R.oe <= 1] = '+'
  mark[R.oe > 1 & R.oe <= 2] = '++'
  mark[R.oe > 2] = '+++'
  
  cell_function = function(j, i, x, y, width, height, fill){
    grid.text(mark[i, j], x, y, gp = gpar(fontsize = 10))
  }
  
  getPalette = colorRampPalette(brewer.pal(11, "Blues")[1:6])
  R.oe.v2 = R.oe
  R.oe.v2[R.oe.v2 > 3] = 3
  pdf("scIBD.cell_type_compare.between_disease.endo.pdf", width = 4, height = 4.5)
  Heatmap(R.oe.v2,
          column_title = "",
          cluster_rows = FALSE, cluster_columns = FALSE,
          #row_split = c(rep("1-Myeloid", myeloid$V1[myeloid$V1 %in% rownames(R.oe.v2)] %>% length ), 
          #              rep("2-CD4+ T", cd4t$V1[cd4t$V1 %in% rownames(R.oe.v2)] %>% length), 
          #              rep("3-CD8+ T", cd8t$V1[cd8t$V1 %in% rownames(R.oe.v2)] %>% length),
          #              rep("4-ILC", ilc$V1[ilc$V1 %in% rownames(R.oe.v2)] %>% length), 
          #              rep("5-B/Plasma",B_Plasma$V1[B_Plasma$V1 %in% rownames(R.oe.v2)] %>% length)),
          cell_fun = cell_function,
          col = getPalette(5),
          heatmap_legend_param = list(
            title = "Ro/e",
            break_dist = 1,
            #col_fun = colorRamp2(c(0, 1, 2, 4, max(R.oe.v2)), getPalette(5)),
            #at = c(0, 1, 2, 4, max(R.oe.v2)),
            col_fun = colorRamp2(c(0, 0.5, 1, 2, 4), getPalette(5)),
            at = c(0, 1, 2, 3, 4),
            labels = c('0','0.5','1','3','max')
            #labels = c('0','0.5','1','4',sprintf("%.1f", max(R.oe.v2)))
          )
  )
  dev.off()
}

# Inflammatory fibroblast
if(T){
  tmp = meta %>% filter(stage == "adult") %>% 
    filter(major_cluster == "Mesenchymal") %>%
    filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>%
    filter(!(tissue.sub %in% c("blood","lymph node"))) %>% 
    group_by(sample, minor_cluster) %>% 
    summarise(count = n()) %>% reshape2::dcast(sample ~ minor_cluster)
  tmp[is.na(tmp)] = 0
  rownames(tmp) = tmp$sample
  tmp$sample = NULL
  table(rowSums(tmp) > 100)
  
  tmp.filter = tmp[rowSums(tmp) > 100, ]
  tmp.filter = tmp.filter/rowSums(tmp.filter)
  #View(tmp.filter)
  
  tmp.meta = meta[,c("sample","disease")] %>% unique
  tmp.meta = tmp.meta[!is.na(tmp.meta$sample),]
  
  tmp.filter$sample = rownames(tmp.filter)
  tmp.filter = left_join(tmp.filter, tmp.meta, by = "sample", kep)
  
  cluster = "Inflammatory fibroblast"
  tmp.filter.sub = tmp.filter[,c(cluster,"disease")]
  colnames(tmp.filter.sub) = c("pct", "group")
  tmp.filter.sub$group = tmp.filter.sub$group %>% droplevels()
  
  # summary sample size in each group
  tmp.filter %>% group_by(disease) %>% count() # CD_inflamed 2, Healthy 39, UC_inflamed 29
  
  my_comparisons = list( c("CD_inflamed","Healthy"),
                         c("UC_inflamed","Healthy"),
                         c("CD_inflamed","UC_inflamed"))
  ggplot(tmp.filter.sub, aes(x=group, y = pct, fill = group, color = group)) + 
    labs(x = "",y = "Percentage of cells", title = "") + 
    geom_boxplot( outlier.size = -1, linewidth =0.5) +
    geom_point( position = position_jitter(width = 0.1), size = 1) + 
    scale_fill_manual(values = alpha(paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)], 0.6) ) + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) + 
    scale_x_discrete(labels = c("CD","HC","UC")) + 
    theme_classic2() + 
    theme(
      legend.position = "none",
      #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black", angle = 45, hjust = 1),
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "italic", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, 
                           method = "t.test", label = "p.signif", 
                           method.args = list(alternative = "greater")) 
  #ggsave("IAFs.pct.boxplot.pdf", width = 2, height = 3)
  
  t.test(tmp.filter.sub[tmp.filter.sub$group == "Healthy",]$pct, 
                tmp.filter.sub[tmp.filter.sub$group == "CD_inflamed",]$pct)$p.value/2
  t.test(tmp.filter.sub[tmp.filter.sub$group == "Healthy",]$pct, 
         tmp.filter.sub[tmp.filter.sub$group == "UC_inflamed",]$pct)$p.value/2
  t.test(tmp.filter.sub[tmp.filter.sub$group == "CD_inflamed",]$pct, 
         tmp.filter.sub[tmp.filter.sub$group == "UC_inflamed",]$pct)$p.value/2
}

# Reticular fibroblast
if(T){
  tmp = meta %>% filter(stage == "adult") %>% 
    filter(major_cluster == "Mesenchymal") %>%
    filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>%
    filter(!(tissue.sub %in% c("blood","lymph node"))) %>% 
    group_by(sample, minor_cluster) %>% 
    summarise(count = n()) %>% reshape2::dcast(sample ~ minor_cluster)
  tmp[is.na(tmp)] = 0
  rownames(tmp) = tmp$sample
  tmp$sample = NULL
  table(rowSums(tmp) > 100)
  
  tmp.filter = tmp[rowSums(tmp) > 100, ]
  tmp.filter = tmp.filter/rowSums(tmp.filter)
  #View(tmp.filter)
  
  tmp.meta = meta[,c("sample","disease")] %>% unique
  tmp.meta = tmp.meta[!is.na(tmp.meta$sample),]
  
  tmp.filter$sample = rownames(tmp.filter)
  tmp.filter = left_join(tmp.filter, tmp.meta, by = "sample", kep)
  
  cluster = "Reticular fibroblast"
  tmp.filter.sub = tmp.filter[,c(cluster,"disease")]
  colnames(tmp.filter.sub) = c("pct", "group")
  tmp.filter.sub$group = tmp.filter.sub$group %>% droplevels()
  
  my_comparisons = list( c("CD_inflamed","Healthy"),
                         c("Healthy","UC_inflamed"),
                         c("CD_inflamed","UC_inflamed"))
  ggplot(tmp.filter.sub, aes(x=group, y = pct, fill = group, color = group)) + 
    labs(x = "",y = "Percentage of cells", title = "") + 
    geom_boxplot( outlier.size = -1, linewidth =0.5) +
    geom_point( position = position_jitter(width = 0.1), size = 1) + 
    scale_fill_manual(values = alpha(paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)], 0.6) ) + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) + 
    scale_x_discrete(labels = c("CD","HC","UC")) + 
    theme_classic2() + 
    theme(
      legend.position = "none",
      #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black", angle = 45, hjust = 1),
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "italic", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, 
                           method = "t.test", label = "p.signif", 
                           method.args = list(alternative = "less")) 
  #ggsave("Reticular_fibroblast.pct.boxplot.pdf", width = 2, height = 3)
  
  t.test(tmp.filter.sub[tmp.filter.sub$group == "Healthy",]$pct, 
         tmp.filter.sub[tmp.filter.sub$group == "CD_inflamed",]$pct)$p.value/2
  t.test(tmp.filter.sub[tmp.filter.sub$group == "Healthy",]$pct, 
         tmp.filter.sub[tmp.filter.sub$group == "UC_inflamed",]$pct)$p.value/2
  t.test(tmp.filter.sub[tmp.filter.sub$group == "CD_inflamed",]$pct, 
         tmp.filter.sub[tmp.filter.sub$group == "UC_inflamed",]$pct)$p.value/2
}

# Cycling stromal
# too low 
if(T){
  tmp = meta %>% filter(stage == "adult") %>% 
    filter(major_cluster == "Mesenchymal") %>%
    filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>%
    filter(!(tissue.sub %in% c("blood","lymph node"))) %>% 
    group_by(sample, minor_cluster) %>% 
    summarise(count = n()) %>% reshape2::dcast(sample ~ minor_cluster)
  tmp[is.na(tmp)] = 0
  rownames(tmp) = tmp$sample
  tmp$sample = NULL
  table(rowSums(tmp) > 100)
  
  tmp.filter = tmp[rowSums(tmp) > 100, ]
  tmp.filter = tmp.filter/rowSums(tmp.filter)
  #View(tmp.filter)
  
  tmp.meta = meta[,c("sample","disease")] %>% unique
  tmp.meta = tmp.meta[!is.na(tmp.meta$sample),]
  
  tmp.filter$sample = rownames(tmp.filter)
  tmp.filter = left_join(tmp.filter, tmp.meta, by = "sample", kep)
  
  cluster = "Cycling stromal"
  tmp.filter.sub = tmp.filter[,c(cluster,"disease")]
  colnames(tmp.filter.sub) = c("pct", "group")
  
  my_comparisons = list( c("Healthy","CD_inflamed"),
                         c("Healthy","UC_inflamed"),
                         c("CD_inflamed","UC_inflamed"))
  ggplot(tmp.filter.sub, aes(x=group, y = pct, fill = group, color = group)) + 
    labs(x = "",y = "Percentage of cells", title = "") + 
    geom_boxplot( outlier.size = -1, linewidth =0.5) +
    geom_point( position = position_jitter(width = 0.1), size = 1) + 
    scale_fill_manual(values = alpha(paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)], 0.6) ) + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) + 
    scale_x_discrete(labels = c("CD","HC","UC")) + 
    theme_classic2() + 
    theme(
      legend.position = "none",
      #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black", angle = 45, hjust = 1),
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "italic", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, 
                           method = "t.test", label = "p.signif", 
                           method.args = list(alternative = "less")) 
  ggsave("Cycling_stromal.pct.boxplot.pdf", width = 2, height = 3)
}

# Immature pericyte
if(T){
  tmp = meta %>% filter(stage == "adult") %>% 
    filter(major_cluster == "Mesenchymal") %>%
    filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>%
    filter(!(tissue.sub %in% c("blood","lymph node"))) %>% 
    group_by(sample, minor_cluster) %>% 
    summarise(count = n()) %>% reshape2::dcast(sample ~ minor_cluster)
  tmp[is.na(tmp)] = 0
  rownames(tmp) = tmp$sample
  tmp$sample = NULL
  table(rowSums(tmp) > 100)
  
  tmp.filter = tmp[rowSums(tmp) > 100, ]
  tmp.filter = tmp.filter/rowSums(tmp.filter)
  #View(tmp.filter)
  
  tmp.meta = meta[,c("sample","disease")] %>% unique
  tmp.meta = tmp.meta[!is.na(tmp.meta$sample),]
  
  tmp.filter$sample = rownames(tmp.filter)
  tmp.filter = left_join(tmp.filter, tmp.meta, by = "sample", kep)
  
  cluster = "Immature pericyte"
  tmp.filter.sub = tmp.filter[,c(cluster,"disease")]
  colnames(tmp.filter.sub) = c("pct", "group")
  
  my_comparisons = list( c("Healthy","CD_inflamed"),
                         c("Healthy","UC_inflamed"),
                         c("CD_inflamed","UC_inflamed"))
  ggplot(tmp.filter.sub, aes(x=group, y = pct, fill = group, color = group)) + 
    labs(x = "",y = "Percentage of cells", title = "") + 
    geom_boxplot( outlier.size = -1, linewidth =0.5) +
    geom_point( position = position_jitter(width = 0.1), size = 1) + 
    scale_fill_manual(values = alpha(paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)], 0.6) ) + 
    scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,2,3)] ) + 
    scale_x_discrete(labels = c("CD","HC","UC")) + 
    theme_classic2() + 
    theme(
      legend.position = "none",
      #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black", angle = 45, hjust = 1),
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "italic", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    ) + stat_compare_means(comparisons = my_comparisons, 
                           method = "t.test", label = "p.signif", 
                           method.args = list(alternative = "less")) 
  #ggsave("Immature_pericyte.pct.boxplot.pdf", width = 2, height = 3)
  
  t.test(tmp.filter.sub[tmp.filter.sub$group == "Healthy",]$pct, 
         tmp.filter.sub[tmp.filter.sub$group == "CD_inflamed",]$pct)$p.value/2
  t.test(tmp.filter.sub[tmp.filter.sub$group == "Healthy",]$pct, 
         tmp.filter.sub[tmp.filter.sub$group == "UC_inflamed",]$pct)$p.value/2
  t.test(tmp.filter.sub[tmp.filter.sub$group == "CD_inflamed",]$pct, 
         tmp.filter.sub[tmp.filter.sub$group == "UC_inflamed",]$pct)$p.value/2
}

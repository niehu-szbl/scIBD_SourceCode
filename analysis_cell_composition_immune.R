library(ComplexHeatmap)
library(RColorBrewer)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(circlize)


# load meta data
meta = readRDS("all.clean.meta.rds")

# order subset
subset_order = c(meta[meta$major_cluster == "Myeloid",]$minor_cluster %>% droplevels() %>% levels,
  meta[meta$major_cluster == "CD4T",]$minor_cluster %>% droplevels() %>% levels,
  meta[meta$major_cluster == "CD8T",]$minor_cluster %>% droplevels() %>% levels,
  meta[meta$major_cluster == "ILC",]$minor_cluster %>% droplevels() %>% levels,
  meta[meta$major_cluster == "B_Plasma",]$minor_cluster %>% droplevels() %>% levels)

# disease distribution
# Ro/e, heatmap
# Figure 1d
# Extended Data Figure 5b
if(T){
  
  # select data
  data.sub = meta %>%
    filter(major_cluster %in% c("Myeloid","CD4T","CD8T","ILC","B_Plasma")) %>%
    filter(stage == "adult") %>%
    filter(tissue %in% c("smallInt","largeInt","blood"))
  
  # count
  N.o <- table(data.sub[["minor_cluster"]],data.sub[["disease"]])
  
  # order disease
  N.o = N.o[,c("CD_PBMC","UC_PBMC","Healthy_PBMC","CD_inflamed","Healthy","UC_inflamed",
               "CD_non_inflamed","UC_non_inflamed")]
  
  # filter cluster
  N.o = N.o[which(rowSums(N.o > 20) > 0) %>% names(),]
  
  # order cluster
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
    grid.text(mark[i, j], x, y, gp = gpar(fontsize = 8))
  }
  
  #getPalette = colorRampPalette(brewer.pal(11, "Spectral")[2:6] %>% rev)
  getPalette = colorRampPalette(brewer.pal(9, "Blues")[1:6])
  R.oe.v2 = R.oe
  R.oe.v2[R.oe.v2 > 3] = 3
  pdf("scIBD.cell_type_compare.between_disease_statues.immune.pdf", width = 5, height = 7)
  Heatmap(R.oe.v2,
          column_title = "",
          cluster_rows = FALSE, cluster_columns = FALSE,
          row_split = c(rep("5-Myeloid", myeloid$V1[myeloid$V1 %in% rownames(R.oe.v2)] %>% length ), 
                        rep("1-CD4+ T", cd4t$V1[cd4t$V1 %in% rownames(R.oe.v2)] %>% length), 
                        rep("2-CD8+ T", cd8t$V1[cd8t$V1 %in% rownames(R.oe.v2)] %>% length),
                        rep("3-ILC", ilc$V1[ilc$V1 %in% rownames(R.oe.v2)] %>% length), 
                        rep("4-B/Plasma",B_Plasma$V1[B_Plasma$V1 %in% rownames(R.oe.v2)] %>% length)),
          column_split = c(rep("PBMC",3), rep("Intestine-1",3), rep("Intestine-2",2)),
          cell_fun = cell_function,
          col = getPalette(5),
          heatmap_legend_param = list(
            title = "Ro/e",
            break_dist = 1,
            #col_fun = colorRamp2(c(0, 1, 2, 4, max(R.oe.v2)), getPalette(5)),
            #at = c(0, 1, 2, 4, max(R.oe.v2)),
            col_fun = colorRamp2(c(0, 0.5, 1, 2, 3), getPalette(5)),
            at = c(0, 1, 2, 3, 4),
            labels = c('0','0.5','1','2','max')
            #labels = c('0','0.5','1','4',sprintf("%.1f", max(R.oe.v2)))
          )
  )
  dev.off()
}

# tissue distribution
# Ro/e, heatmap plot
# Extended Data Figure 5a
if(T){
  # select data
  data.sub = meta %>%
    filter(major_cluster %in% c("Myeloid","CD4T","CD8T","ILC","B_Plasma")) %>%
    filter(stage == "adult")
  
  # count
  N.o <- table(data.sub[["minor_cluster"]],data.sub[["tissue.sub"]])
  
  # order tissue
  N.o = N.o[,c("blood","cecum","colon","rectum","ileum","lymph node")]
  
  # filter cluster
  N.o = N.o[which(rowSums(N.o > 20) > 0) %>% names(),]
  
  # order cluster
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
  
  getPalette = colorRampPalette(brewer.pal(9, "Blues")[1:6])
  R.oe.v2 = R.oe
  R.oe.v2[R.oe.v2 > 3] = 3
  pdf("scIBD.cell_type_compare.between_tissue.immune.pdf", width = 4.5, height = 7)
  Heatmap(R.oe.v2,
          column_title = "",
          cluster_rows = FALSE, cluster_columns = FALSE,
          row_split = c(rep("5-Myeloid", myeloid$V1[myeloid$V1 %in% rownames(R.oe.v2)] %>% length ), 
                        rep("1-CD4+ T", cd4t$V1[cd4t$V1 %in% rownames(R.oe.v2)] %>% length), 
                        rep("2-CD8+ T", cd8t$V1[cd8t$V1 %in% rownames(R.oe.v2)] %>% length),
                        rep("3-ILC", ilc$V1[ilc$V1 %in% rownames(R.oe.v2)] %>% length), 
                        rep("4-B/Plasma",B_Plasma$V1[B_Plasma$V1 %in% rownames(R.oe.v2)] %>% length)),
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

# disease distribution
# Box plot, by sample
# Figure 1e
# The final p values are divided by 2 to calculate on-tailed p value
if(T){
  # collect all to plot
  all_res = data.frame()
  
  # by sample
  # boxplot
  # CD4T, CD4+ Treg
  if(T){
    tmp = meta %>% filter(stage == "adult") %>%
      filter(major_cluster == "CD4T") %>%
      filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>%
      filter( tissue %in% c("smallInt", "largeInt")) %>% 
      group_by(sample, minor_cluster) %>% 
      summarise(count = n()) %>% reshape2::dcast(sample ~ minor_cluster)
    tmp[is.na(tmp)] = 0
    rownames(tmp) = tmp$sample
    tmp$sample = NULL
    table(rowSums(tmp) > 100)
    
    tmp.filter = tmp[rowSums(tmp) > 100, ]
    tmp.filter = tmp.filter/rowSums(tmp.filter)
    
    tmp.meta = meta[,c("sample","disease")] %>% unique
    tmp.meta = tmp.meta[!is.na(tmp.meta$sample),]
    
    tmp.filter$sample = rownames(tmp.filter)
    tmp.filter = left_join(tmp.filter, tmp.meta, by = "sample", kep)
    
    cluster = "CD4+ Treg"
    tmp.filter.sub = tmp.filter[,c(cluster,"disease")]
    colnames(tmp.filter.sub) = c("pct", "group")
    tmp.filter.sub$cluster = cluster
    all_res = rbind(all_res, tmp.filter.sub)
    
    # plot
    if(F){
      my_comparisons = list( c("Healthy","CD_inflamed"),
                             c("Healthy","UC_inflamed"),
                             c("CD_inflamed","UC_inflamed"))
      
      ggplot(tmp.filter.sub, aes(x=group, y = pct, fill = group, color = group)) + 
        geom_boxplot(outlier.size = -1, linewidth =1) +
        geom_point( position = position_jitter(width = 0.1), size = 2) + 
        labs(x = "",y = "Percentage of cells", title = cluster) + 
        scale_x_discrete(labels = c("CD", "HC", "UC")) + 
        scale_fill_d3(alpha = 0.6) + 
        scale_color_d3() + 
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
        ) + stat_compare_means(comparisons = my_comparisons, method = "t.test")
    }
  }
  
  # CD8T, CD8+ Tc17
  if(T){
    tmp = meta %>% filter(stage == "adult") %>%
      filter(major_cluster == "CD8T") %>%
      filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>%
      filter( tissue %in% c("smallInt", "largeInt")) %>% 
      group_by(sample, minor_cluster) %>% 
      summarise(count = n()) %>% reshape2::dcast(sample ~ minor_cluster)
    tmp[is.na(tmp)] = 0
    rownames(tmp) = tmp$sample
    tmp$sample = NULL
    table(rowSums(tmp) > 100)
    
    tmp.filter = tmp[rowSums(tmp) > 100, ]
    tmp.filter = tmp.filter/rowSums(tmp.filter)
    
    tmp.meta = meta[,c("sample","disease")] %>% unique
    tmp.meta = tmp.meta[!is.na(tmp.meta$sample),]
    
    tmp.filter$sample = rownames(tmp.filter)
    tmp.filter = left_join(tmp.filter, tmp.meta, by = "sample", kep)
    
    cluster = "CD8+ Tc17"
    tmp.filter.sub = tmp.filter[,c(cluster,"disease")]
    colnames(tmp.filter.sub) = c("pct", "group")
    tmp.filter.sub$cluster = cluster
    all_res = rbind(all_res, tmp.filter.sub)
    
    # plot
    if(F){
      my_comparisons = list( c("Healthy","CD_inflamed"),
                             c("Healthy","UC_inflamed"),
                             c("CD_inflamed","UC_inflamed"))
      
      ggplot(tmp.filter.sub, aes(x=group, y = pct, fill = group, color = group)) + 
        geom_boxplot(outlier.size = -1, linewidth =1) +
        geom_point( position = position_jitter(width = 0.1), size = 2) + 
        labs(x = "",y = "Percentage of cells", title = cluster) + 
        scale_x_discrete(labels = c("CD", "HC", "UC")) + 
        scale_fill_d3(alpha = 0.6) + 
        scale_color_d3() + 
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
        ) + stat_compare_means(comparisons = my_comparisons, method = "t.test")  
    }
  }
  
  # IgG plasma
  if(T){
    tmp = meta %>% filter(stage == "adult") %>%
      filter(major_cluster == "B_Plasma") %>%
      filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>%
      filter( tissue %in% c("smallInt", "largeInt")) %>% 
      group_by(sample, minor_cluster) %>% 
      summarise(count = n()) %>% reshape2::dcast(sample ~ minor_cluster)
    tmp[is.na(tmp)] = 0
    rownames(tmp) = tmp$sample
    tmp$sample = NULL
    table(rowSums(tmp) > 100)
    
    tmp.filter = tmp[rowSums(tmp) > 100, ]
    tmp.filter = tmp.filter/rowSums(tmp.filter)
    
    tmp.meta = meta[,c("sample","disease")] %>% unique
    tmp.meta = tmp.meta[!is.na(tmp.meta$sample),]
    
    tmp.filter$sample = rownames(tmp.filter)
    tmp.filter = left_join(tmp.filter, tmp.meta, by = "sample", kep)
    
    cluster = "IgG plasma"
    tmp.filter.sub = tmp.filter[,c(cluster,"disease")]
    colnames(tmp.filter.sub) = c("pct", "group")
    tmp.filter.sub$cluster = cluster
    all_res = rbind(all_res, tmp.filter.sub)
    
    # plot
    if(F){
      my_comparisons = list( c("Healthy","CD_inflamed"),
                             c("Healthy","UC_inflamed"),
                             c("UC_inflamed","CD_inflamed"))
      ggplot(tmp.filter.sub, aes(x=group, y = pct, fill = group, color = group)) + 
        geom_boxplot(outlier.size = -1, linewidth =1) +
        geom_point( position = position_jitter(width = 0.1), size = 2) + 
        labs(x = "",y = "Percentage of cells", title = cluster) + 
        scale_x_discrete(labels = c("CD", "HC", "UC")) + 
        scale_fill_d3(alpha = 0.6) + 
        scale_color_d3() + 
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
                               method = "t.test")  
    }
  }
  
  # GC B
  if(T){
    tmp = meta %>% filter(stage == "adult") %>%
      filter(major_cluster == "B_Plasma") %>%
      filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>%
      filter( tissue %in% c("smallInt", "largeInt")) %>% 
      group_by(sample, minor_cluster) %>% 
      summarise(count = n()) %>% reshape2::dcast(sample ~ minor_cluster)
    tmp[is.na(tmp)] = 0
    rownames(tmp) = tmp$sample
    tmp$sample = NULL
    table(rowSums(tmp) > 100)
    
    tmp.filter = tmp[rowSums(tmp) > 100, ]
    tmp.filter = tmp.filter/rowSums(tmp.filter)
    
    tmp.meta = meta[,c("sample","disease")] %>% unique
    tmp.meta = tmp.meta[!is.na(tmp.meta$sample),]
    
    tmp.filter$sample = rownames(tmp.filter)
    tmp.filter = left_join(tmp.filter, tmp.meta, by = "sample", kep)
    
    cluster = "GC B"
    tmp.filter.sub = tmp.filter[,c(cluster,"disease")]
    colnames(tmp.filter.sub) = c("pct", "group")
    tmp.filter.sub$cluster = cluster
    all_res = rbind(all_res, tmp.filter.sub)
    
    # plot
    if(F){
      my_comparisons = list( c("Healthy","CD_inflamed"),
                             c("Healthy","UC_inflamed"),
                             c("CD_inflamed","UC_inflamed"))
      ggplot(tmp.filter.sub, aes(x=group, y = pct, fill = group, color = group)) + 
        geom_boxplot(outlier.size = -1, linewidth =1) +
        geom_point( position = position_jitter(width = 0.1), size = 2) + 
        labs(x = "",y = "Percentage of cells", title = cluster) + 
        scale_x_discrete(labels = c("CD", "HC", "UC")) + 
        scale_fill_d3(alpha = 0.6) + 
        scale_color_d3() + 
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
        ) + stat_compare_means(comparisons = my_comparisons, method = "t.test")  
    }
  }
  
  # Cycling GC B
  if(T){
    tmp = meta %>% filter(stage == "adult") %>%
      filter(major_cluster == "B_Plasma") %>%
      filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>%
      filter( tissue %in% c("smallInt", "largeInt")) %>% 
      group_by(sample, minor_cluster) %>% 
      summarise(count = n()) %>% reshape2::dcast(sample ~ minor_cluster)
    tmp[is.na(tmp)] = 0
    rownames(tmp) = tmp$sample
    tmp$sample = NULL
    table(rowSums(tmp) > 100)
    
    tmp.filter = tmp[rowSums(tmp) > 100, ]
    tmp.filter = tmp.filter/rowSums(tmp.filter)
    
    tmp.meta = meta[,c("sample","disease")] %>% unique
    tmp.meta = tmp.meta[!is.na(tmp.meta$sample),]
    
    tmp.filter$sample = rownames(tmp.filter)
    tmp.filter = left_join(tmp.filter, tmp.meta, by = "sample", kep)
    
    cluster = "Cycling GC B"
    tmp.filter.sub = tmp.filter[,c(cluster,"disease")]
    colnames(tmp.filter.sub) = c("pct", "group")
    tmp.filter.sub$cluster = cluster
    all_res = rbind(all_res, tmp.filter.sub)
    
    # plot
    if(F){
      my_comparisons = list( c("Healthy","CD_inflamed"),
                             c("Healthy","UC_inflamed"),
                             c("UC_inflamed","CD_inflamed"))
      ggplot(tmp.filter.sub, aes(x=group, y = pct, fill = group, color = group)) + 
        geom_boxplot(outlier.alpha = 0.5, outlier.size = 2, linewidth =1) +
        geom_point( position = position_jitter(width = 0.1), size = 2) + 
        labs(x = "",y = "Percentage of cells", title = cluster) + 
        scale_x_discrete(labels = c("CD", "HC", "UC")) + 
        scale_fill_d3(alpha = 0.6) + 
        scale_color_d3() + 
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
                               method = "t.test")  
    }
  }
  
  # Cycling plasma
  if(T){
    tmp = meta %>% filter(stage == "adult") %>%
      filter(major_cluster == "B_Plasma") %>%
      filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>%
      filter( tissue %in% c("smallInt", "largeInt")) %>% 
      group_by(sample, minor_cluster) %>% 
      summarise(count = n()) %>% reshape2::dcast(sample ~ minor_cluster)
    tmp[is.na(tmp)] = 0
    rownames(tmp) = tmp$sample
    tmp$sample = NULL
    table(rowSums(tmp) > 100)
    
    tmp.filter = tmp[rowSums(tmp) > 100, ]
    tmp.filter = tmp.filter/rowSums(tmp.filter)
    
    tmp.meta = meta[,c("sample","disease")] %>% unique
    tmp.meta = tmp.meta[!is.na(tmp.meta$sample),]
    
    tmp.filter$sample = rownames(tmp.filter)
    tmp.filter = left_join(tmp.filter, tmp.meta, by = "sample", kep)
    
    cluster = "Cycling plasma"
    tmp.filter.sub = tmp.filter[,c(cluster,"disease")]
    colnames(tmp.filter.sub) = c("pct", "group")
    tmp.filter.sub$cluster = cluster
    all_res = rbind(all_res, tmp.filter.sub)
    
    # plot
    if(F){
      my_comparisons = list( c("Healthy","CD_inflamed"),
                             c("Healthy","UC_inflamed"),
                             c("UC_inflamed","CD_inflamed"))
      ggplot(tmp.filter.sub, aes(x=group, y = pct, fill = group, color = group)) + 
        geom_boxplot(outlier.size = -1, linewidth =1) +
        geom_point( position = position_jitter(width = 0.1), size = 2) + 
        labs(x = "",y = "Percentage of cells", title = cluster) + 
        scale_x_discrete(labels = c("CD", "HC", "UC")) + 
        scale_fill_d3(alpha = 0.6) + 
        scale_color_d3() + 
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
                               method = "t.test")  
    }
  }
  
  # Myeloid, Inflammatory monocyte
  if(T){
    tmp = meta %>% filter(stage == "adult") %>%
      filter(major_cluster == "Myeloid") %>%
      filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>%
      filter( tissue %in% c("smallInt", "largeInt")) %>% 
      group_by(sample, minor_cluster) %>% 
      summarise(count = n()) %>% reshape2::dcast(sample ~ minor_cluster)
    tmp[is.na(tmp)] = 0
    rownames(tmp) = tmp$sample
    tmp$sample = NULL
    table(rowSums(tmp) > 100)
    
    tmp.filter = tmp[rowSums(tmp) > 100, ]
    tmp.filter = tmp.filter/rowSums(tmp.filter)
    
    tmp.meta = meta[,c("sample","disease")] %>% unique
    tmp.meta = tmp.meta[!is.na(tmp.meta$sample),]
    
    tmp.filter$sample = rownames(tmp.filter)
    tmp.filter = left_join(tmp.filter, tmp.meta, by = "sample", kep)
    
    cluster = "Inflammatory monocyte"
    tmp.filter.sub = tmp.filter[,c(cluster,"disease")]
    colnames(tmp.filter.sub) = c("pct", "group")
    tmp.filter.sub$cluster = cluster
    all_res = rbind(all_res, tmp.filter.sub)
    
    # plot
    if(F){
      my_comparisons = list( c("Healthy","CD_inflamed"),
                             c("Healthy","UC_inflamed"),
                             c("UC_inflamed","CD_inflamed"))
      ggplot(tmp.filter.sub, aes(x=group, y = pct, fill = group, color = group)) + 
        geom_boxplot(outlier.alpha = 0.5, outlier.size = 2, linewidth =1) +
        geom_point( position = position_jitter(width = 0.1), size = 2) + 
        labs(x = "",y = "Percentage of cells", title = cluster) + 
        scale_x_discrete(labels = c("CD", "HC", "UC")) + 
        scale_fill_d3(alpha = 0.6) + 
        scale_color_d3() + 
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
        ) + stat_compare_means(comparisons = my_comparisons, method = "t.test")  
    }
  }
  
  # Myeloid, APOE macrophage
  if(T){
    tmp = meta %>% filter(stage == "adult") %>%
      filter(major_cluster == "Myeloid") %>%
      filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>%
      filter( tissue %in% c("smallInt", "largeInt")) %>% 
      group_by(sample, minor_cluster) %>% 
      summarise(count = n()) %>% reshape2::dcast(sample ~ minor_cluster)
    tmp[is.na(tmp)] = 0
    rownames(tmp) = tmp$sample
    tmp$sample = NULL
    table(rowSums(tmp) > 100)
    
    tmp.filter = tmp[rowSums(tmp) > 100, ]
    tmp.filter = tmp.filter/rowSums(tmp.filter)
    
    tmp.meta = meta[,c("sample","disease")] %>% unique
    tmp.meta = tmp.meta[!is.na(tmp.meta$sample),]
    
    tmp.filter$sample = rownames(tmp.filter)
    tmp.filter = left_join(tmp.filter, tmp.meta, by = "sample", kep)
    
    cluster = "APOE+ macrophage"
    tmp.filter.sub = tmp.filter[,c(cluster,"disease")]
    colnames(tmp.filter.sub) = c("pct", "group")
    tmp.filter.sub$cluster = cluster
    all_res = rbind(all_res, tmp.filter.sub)
    
    # plot
    if(F){
      my_comparisons = list( c("Healthy","CD_inflamed"),
                             c("Healthy","UC_inflamed"),
                             c("UC_inflamed","CD_inflamed"))
      ggplot(tmp.filter.sub, aes(x=group, y = pct, fill = group, color = group)) + 
        geom_boxplot(outlier.size = -1, linewidth =1) +
        geom_point( position = position_jitter(width = 0.1), size = 2) + 
        labs(x = "",y = "Percentage of cells", title = cluster) + 
        scale_x_discrete(labels = c("CD", "HC", "UC")) + 
        scale_fill_d3(alpha = 0.6) + 
        scale_color_d3() + 
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
        ) + stat_compare_means(comparisons = my_comparisons, method = "t.test")  
    }
  }
  
  # Myeloid, AREG macrophage
  if(T){
    tmp = meta %>% filter(stage == "adult") %>%
      filter(major_cluster == "Myeloid") %>%
      filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>%
      filter( tissue %in% c("smallInt", "largeInt")) %>% 
      group_by(sample, minor_cluster) %>% 
      summarise(count = n()) %>% reshape2::dcast(sample ~ minor_cluster)
    tmp[is.na(tmp)] = 0
    rownames(tmp) = tmp$sample
    tmp$sample = NULL
    table(rowSums(tmp) > 100)
    
    tmp.filter = tmp[rowSums(tmp) > 100, ]
    tmp.filter = tmp.filter/rowSums(tmp.filter)
    
    tmp.meta = meta[,c("sample","disease")] %>% unique
    tmp.meta = tmp.meta[!is.na(tmp.meta$sample),]
    
    tmp.filter$sample = rownames(tmp.filter)
    tmp.filter = left_join(tmp.filter, tmp.meta, by = "sample", kep)
    
    cluster = "AREG+ macrophage"
    tmp.filter.sub = tmp.filter[,c(cluster,"disease")]
    colnames(tmp.filter.sub) = c("pct", "group")
    tmp.filter.sub$cluster = cluster
    all_res = rbind(all_res, tmp.filter.sub)
    
    # plot
    if(F){
      my_comparisons = list( c("Healthy","CD_inflamed"),
                             c("Healthy","UC_inflamed"),
                             c("UC_inflamed","CD_inflamed"))
      ggplot(tmp.filter.sub, aes(x=group, y = pct, fill = group, color = group)) + 
        geom_boxplot(outlier.size = -1, linewidth =1) +
        geom_point( position = position_jitter(width = 0.1), size = 2) + 
        labs(x = "",y = "Percentage of cells", title = cluster) + 
        scale_x_discrete(labels = c("CD", "HC", "UC")) + 
        scale_fill_d3(alpha = 0.6) + 
        scale_color_d3() + 
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
        ) + stat_compare_means(comparisons = my_comparisons, method = "t.test")  
    }
  }
  
  # Myeloid, LAMP3+ DC
  if(T){
    tmp = meta %>% filter(stage == "adult") %>%
      filter(major_cluster == "Myeloid") %>%
      filter(disease %in% c("CD_inflamed","Healthy","UC_inflamed")) %>%
      filter( tissue %in% c("smallInt", "largeInt")) %>% 
      group_by(sample, minor_cluster) %>% 
      summarise(count = n()) %>% reshape2::dcast(sample ~ minor_cluster)
    tmp[is.na(tmp)] = 0
    rownames(tmp) = tmp$sample
    tmp$sample = NULL
    table(rowSums(tmp) > 100)
    
    tmp.filter = tmp[rowSums(tmp) > 100, ]
    tmp.filter = tmp.filter/rowSums(tmp.filter)
    
    tmp.meta = meta[,c("sample","disease")] %>% unique
    tmp.meta = tmp.meta[!is.na(tmp.meta$sample),]
    
    tmp.filter$sample = rownames(tmp.filter)
    tmp.filter = left_join(tmp.filter, tmp.meta, by = "sample", kep)
    
    cluster = "LAMP3+ DC"
    tmp.filter.sub = tmp.filter[,c(cluster,"disease")]
    colnames(tmp.filter.sub) = c("pct", "group")
    tmp.filter.sub$cluster = cluster
    all_res = rbind(all_res, tmp.filter.sub)
    
    # plot
    if(F){
      my_comparisons = list( c("Healthy","CD_inflamed"),
                             c("Healthy","UC_inflamed"),
                             c("UC_inflamed","CD_inflamed"))
      ggplot(tmp.filter.sub, aes(x=group, y = pct, fill = group, color = group)) + 
        geom_boxplot(outlier.size = -1, linewidth =1) +
        geom_point( position = position_jitter(width = 0.1), size = 2) + 
        labs(x = "",y = "Percentage of cells", title = cluster) + 
        scale_x_discrete(labels = c("CD", "HC", "UC")) + 
        scale_fill_d3(alpha = 0.6) + 
        scale_color_d3() + 
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
        ) + stat_compare_means(comparisons = my_comparisons, method = "t.test")  
    }
  }
  
  # plot together
  if(T){
    my_comparisons = list( c("Healthy","CD_inflamed"),
                           c("Healthy","UC_inflamed"),
                           c("UC_inflamed","CD_inflamed"))
    #tmp = levels(meta$minor_cluster)
    #cluster_order = tmp[tmp %in% all_res$cluster]
    all_res$cluster = factor(all_res$cluster, 
                             levels = c("CD4+ Treg","CD8+ Tc17","IgG plasma","GC B","Cycling GC B","Cycling plasma",
                                        "Inflammatory monocyte","APOE+ macrophage","AREG+ macrophage","LAMP3+ DC"))
    
    colnames(all_res)
    # check sample size
    all_res %>% group_by(group, cluster) %>% count() %>% reshape2::dcast(cluster ~ group, value.var = "n")
    
    # plot
    ggplot(all_res, aes(x=group, y = pct, fill = group, color = group)) + 
      facet_wrap(~ cluster, ncol = 5, scales = "free_y") + 
      geom_boxplot(linewidth = 1, outlier.size = -1) +
      geom_point( position = position_jitter(width = 0.1), size = 1.5) + 
      labs(x = "",y = "Percentage of cells", title = "") + 
      scale_x_discrete(labels = c("CD", "HC", "UC")) + 
      scale_fill_d3(alpha = 0.6) + 
      scale_color_d3() + 
      theme_classic2() + 
      theme(
        legend.position = "right",
        #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black", angle = 45, hjust = 1),
        #axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
        axis.text.x = element_text(size = 0),
        axis.text.y = element_text(family = "ArialMT", size = 12, color = 'black'),
        axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
        axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
        legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
        legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
        strip.text = element_text(family = "ArialMT",size = 12, color ='black'),
        strip.background = element_blank(),
        plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
      ) + stat_compare_means(comparisons = my_comparisons, 
                             method = "t.test")  
    #ggsave("cell_composition.scIBD.pdf",width = 8.5, height = 4)
  }
}


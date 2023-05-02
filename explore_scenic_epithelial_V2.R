library(SCopeLoomR)
library(dplyr)
library(ggrepel)
library(RColorBrewer)
library(DT)

# load data
scenicLoomPath = "./Epithelial.pyscenic_output.updated.loom"
loom <- open_loom(scenicLoomPath)
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")
regulonAUC <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
cell_annotation = get_cell_annotation(loom)
loom$close()

# load meta_data
meta_data = readRDS("all.clean.meta.rds")

meta_data.myeloid = meta_data[cell_annotation %>% rownames, ]
cell_annotation = cbind(cell_annotation, meta_data.myeloid)
cell_annotation$Label = NULL
cell_annotation$Percent_mito = NULL
cell_annotation$nGene = NULL
cell_annotation$nUMI = NULL

# compare regulon between healthy and UC
# Extended Data Figure 10e
if(T){
  g1_cells = cell_annotation %>%
    filter(stage == "adult") %>%
    filter(disease == "UC_inflamed") %>% rownames
  g1_cells %>% length # 19339
  
  g2_cells = cell_annotation %>% 
    filter(stage == "adult") %>%
    filter(disease == "Healthy") %>% rownames
  g2_cells %>% length #
  
  # select cells
  group <- list(g1_cells = g1_cells, g2_cells = g2_cells)
  regulonActivity_mean_byDisease = sapply(group,
                                          function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
  tfs = rownames(regulonAUC)
  
  # get result
  result <- sapply(tfs, function(tf){
    t_test_res = t.test( getAUC(regulonAUC)[tf, g1_cells] %>% as.vector(),
                         getAUC(regulonAUC)[tf, g2_cells] %>% as.vector())
    res = c(t_test_res$estimate[1], t_test_res$estimate[2],
            t_test_res$estimate[1]/(t_test_res$estimate[2]), t_test_res$p.value)
    names(res) <- c("g1_avg","g2_avg","fold_change","p_value")
    # return result
    return(res)
  }) %>% t %>% as.data.frame
  
  # remove NaN, the mean value of two group are both 0
  #result = result[!is.nan(result$p_value),]
  
  # adjust p-value and round
  result$p.adjust = p.adjust(result$p_value, "BH")
  result$log2FoldChange = log2(result$fold_change)
  #result$log2FoldChange = round(result$log2FoldChange, 6)
  #result$g1_avg <- round(result$g1_avg, 6)
  #result$g2_avg <- round(result$g2_avg, 6)
  
  # write to file
  colnames(result) = c("UC_inflamed", "Healthy", "FoldChange", "p-value", "p.adjust", "log2FoldChange")
  #write.table(result, file = "Epithelial.diff_regulon.UC_vs_Healthy.txt", 
  #            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  # set p.adjust to non-zero value
  #result[result$p.adjust == 0,]$p.adjust = runif(sum(result$p.adjust == 0), min = 0.000001, max = 0.00001)
  
  # set color
  fc_cutoff = 1.25
  result$Diff <- "NO"
  result$Diff[result$log2FoldChange > log2(fc_cutoff) &
                result$p.adjust < 0.01] <- "UP"
  result$Diff[result$log2FoldChange < -1*log2(fc_cutoff) & 
                result$p.adjust < 0.01] <- "DOWN"
  result$delabel <- NA
  result$delabel[result$Diff != "NO"] <- rownames(result)[result$Diff != "NO"]
  
  mycolors <- c("#E41A1C", "#377EB8", "grey50")
  names(mycolors) <- c("UP", "DOWN", "NO")
  
  # limit log2foldchange
  result$x = result$log2FoldChange
  thres_pos = 2
  if( sum(result$x > thres_pos) > 0){
    result[result$log2FoldChange > thres_pos, ]$x = 
      thres_pos + ( result[result$log2FoldChange > thres_pos, ]$log2FoldChange - thres_pos)/10
  }
  thres_neg = -1*thres_pos
  if(sum(result$log2FoldChange < thres_neg) > 0){
    result[result$log2FoldChange < thres_neg, ]$x = 
      thres_neg + (result[result$log2FoldChange < thres_neg, ]$log2FoldChange - thres_neg)/10
  }
  # In case of log2FoldChange equal to Inf or -Inf
  if(sum(result$log2FoldChange == -Inf) > 0){
    result[result$log2FoldChange == -Inf, ]$x = runif(sum(result$log2FoldChange == -Inf), -3, -2)
  }
  if(sum(result$log2FoldChange == Inf) > 0){
    result[result$log2FoldChange == -Inf, ]$x = runif(sum(result$log2FoldChange == -Inf), 2, 3)
  }
  
  # limit -log10(adjusted p-value)
  result$y = -log10(result$p.adjust)
  result$y2 = result$y
  
  y_cutoff = 5
  if( sum(result$y > y_cutoff) > 0){
    result[result$y > y_cutoff, ]$y2 = 
      y_cutoff + ( result[result$y > y_cutoff, ]$y - y_cutoff )/100
  }
  result[result$y == Inf, ]$y2 = runif(sum(result$y == Inf), 8, 10)
  
  ggplot(data=result, aes(x=x, y=y2, col=Diff, label=delabel)) + 
    geom_point(size = 1, alpha = 0.5) +
    xlim(-3,3) +
    ylim(0, 10) + 
    scale_color_manual(values=mycolors)+
    geom_text_repel(max.overlaps = 10, label.size = 0.1) +
    geom_vline(xintercept=c(-log2(fc_cutoff), log2(fc_cutoff)), col="black", linetype="dotted") +
    geom_hline(yintercept=2, col="black", linetype="dotted") +
    theme_classic() + 
    xlab("log2FoldChange") + 
    ylab("-log10(adjustd p-value") + 
    ggtitle('Epithelial') + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = 'Epithelial.regulon_diff.UC_VS_healthy.pdf', width = 8, height = 8, dpi = 600)
}

# compare regulon between healthy and CD
# Extended Data Figure 10f
if(T){
  g1_cells = cell_annotation %>%
    filter(stage == "adult") %>%
    filter(disease == "CD_inflamed") %>% rownames
  g1_cells %>% length # 19339
  
  g2_cells = cell_annotation %>% 
    filter(stage == "adult") %>%
    filter(disease == "Healthy") %>% rownames
  g2_cells %>% length #
  
  # select cells
  group <- list(g1_cells = g1_cells, g2_cells = g2_cells)
  regulonActivity_mean_byDisease = sapply(group,
                                          function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
  tfs = rownames(regulonAUC)
  
  # get result
  result <- sapply(tfs, function(tf){
    t_test_res = t.test( getAUC(regulonAUC)[tf, g1_cells] %>% as.vector(),
                         getAUC(regulonAUC)[tf, g2_cells] %>% as.vector())
    res = c(t_test_res$estimate[1], t_test_res$estimate[2],
            t_test_res$estimate[1]/(t_test_res$estimate[2]), t_test_res$p.value)
    names(res) <- c("g1_avg","g2_avg","fold_change","p_value")
    # return result
    return(res)
  }) %>% t %>% as.data.frame
  
  # remove NaN, the mean value of two group are both 0
  #result = result[!is.nan(result$p_value),]
  
  # adjust p-value and round
  result$p.adjust = p.adjust(result$p_value, "BH")
  result$log2FoldChange = log2(result$fold_change)
  #result$log2FoldChange = round(result$log2FoldChange, 6)
  #result$g1_avg <- round(result$g1_avg, 6)
  #result$g2_avg <- round(result$g2_avg, 6)
  
  # write to file
  colnames(result) = c("CD_inflamed", "Healthy", "FoldChange", "p-value", "p.adjust", "log2FoldChange")
  #write.table(result, file = "Epithelial.diff_regulon.UC_vs_Healthy.txt", 
  #            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  # set p.adjust to non-zero value
  #result[result$p.adjust == 0,]$p.adjust = runif(sum(result$p.adjust == 0), min = 0.000001, max = 0.00001)
  
  # set color
  fc_cutoff = 1.25
  result$Diff <- "NO"
  result$Diff[result$log2FoldChange > log2(fc_cutoff) &
                result$p.adjust < 0.01] <- "UP"
  result$Diff[result$log2FoldChange < -1*log2(fc_cutoff) & 
                result$p.adjust < 0.01] <- "DOWN"
  result$delabel <- NA
  result$delabel[result$Diff != "NO"] <- rownames(result)[result$Diff != "NO"]
  
  mycolors <- c("#E41A1C", "#377EB8", "grey50")
  names(mycolors) <- c("UP", "DOWN", "NO")
  
  # limit log2foldchange
  result$x = result$log2FoldChange
  thres_pos = 2
  if( sum(result$x > thres_pos) > 0){
    result[result$log2FoldChange > thres_pos, ]$x = 
      thres_pos + ( result[result$log2FoldChange > thres_pos, ]$log2FoldChange - thres_pos)/10
  }
  thres_neg = -1*thres_pos
  if(sum(result$log2FoldChange < thres_neg) > 0){
    result[result$log2FoldChange < thres_neg, ]$x = 
      thres_neg + (result[result$log2FoldChange < thres_neg, ]$log2FoldChange - thres_neg)/10
  }
  # In case of log2FoldChange equal to Inf or -Inf
  if(sum(result$log2FoldChange == -Inf) > 0){
    result[result$log2FoldChange == -Inf, ]$x = runif(sum(result$log2FoldChange == -Inf), -3, -2)
  }
  if(sum(result$log2FoldChange == Inf) > 0){
    result[result$log2FoldChange == -Inf, ]$x = runif(sum(result$log2FoldChange == -Inf), 2, 3)
  }
  
  # limit -log10(adjusted p-value)
  result$y = -log10(result$p.adjust)
  result$y2 = result$y
  
  y_cutoff = 5
  if( sum(result$y > y_cutoff) > 0){
    result[result$y > y_cutoff, ]$y2 = 
      y_cutoff + ( result[result$y > y_cutoff, ]$y - y_cutoff )/100
  }
  result[result$y == Inf, ]$y2 = runif(sum(result$y == Inf), 8, 10)
  
  ggplot(data=result, aes(x=x, y=y2, col=Diff, label=delabel)) + 
    geom_point(size = 1, alpha = 0.5) +
    xlim(-3,3) +
    ylim(0, 10) + 
    scale_color_manual(values=mycolors)+
    geom_text_repel(max.overlaps = 15, label.size = 0.1) +
    geom_vline(xintercept=c(-log2(fc_cutoff), log2(fc_cutoff)), col="black", linetype="dotted") +
    geom_hline(yintercept=2, col="black", linetype="dotted") +
    theme_classic() + 
    xlab("log2FoldChange") + 
    ylab("-log10(adjustd p-value") + 
    ggtitle('Epithelial') + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  
  ggsave(filename = 'Epithelial.regulon_diff.CD_VS_healthy.pdf', width = 8, height = 8, dpi = 600)
}








# compare regulon between healthy and CD
if(F){
  g1_cells = cell_annotation %>%
    filter(stage == "adult") %>%
    filter(disease == "CD_inflamed") %>% rownames
  g1_cells %>% length # 1123
  
  g2_cells = cell_annotation %>% 
    filter(stage == "adult") %>%
    filter(disease == "Healthy") %>% rownames
  g2_cells %>% length # 63411
  
  # select cells
  group <- list(g1_cells = g1_cells, g2_cells = g2_cells)
  regulonActivity_mean_byDisease = sapply(group,
                                          function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
  tfs = rownames(regulonAUC)
  
  # get result
  result <- sapply(tfs, function(tf){
    t_test_res = t.test( getAUC(regulonAUC)[tf,g1_cells] %>% as.vector(),
                         getAUC(regulonAUC)[tf,g2_cells] %>% as.vector())
    res = c(t_test_res$estimate[1], t_test_res$estimate[2],
            t_test_res$estimate[1]/(t_test_res$estimate[2]+0.001), t_test_res$p.value)
    names(res) <- c("g1_avg","g2_avg","fold_change","p_value")
    # return result
    return(res)
  }) %>% t %>% as.data.frame
  
  # remove NaN, the mean value of two group are both 0
  result = result[!is.nan(result$p_value),]
  
  # remove regulons with too low activity
  result = result[ result$fold_change != 0, ]
  
  # adjust p-value and round
  result$p.adjust = round(p.adjust(result$p_value, "BH"), 6)
  result$log2FoldChange = log2(result$fold_change)
  result$log2FoldChange = round(result$log2FoldChange, 6)
  result$g1_avg <- round(result$g1_avg, 6)
  result$g2_avg <- round(result$g2_avg, 6)
  
  # write to file
  colnames(result) = c("UC_inflamed", "Healthy", "FoldChange", "p-value", "p.adjust", "log2FoldChange")
  #write.table(result, file = "Epithelial.diff_regulon.CD_vs_Healthy.txt", 
  #            sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  # set p.adjust to non-zero value
  result[result$p.adjust == 0,]$p.adjust = runif(sum(result$p.adjust == 0), 
                                                 min = 0.000001, max = 0.00001)
  
  # limit log2foldchange
  thres_pos = 1.2
  if( sum(result$log2FoldChange > thres_pos) > 0){
    result[result$log2FoldChange > thres_pos, ]$log2FoldChange = 
      thres_pos + ( result[result$log2FoldChange > thres_pos, ]$log2FoldChange - thres_pos)/10
  }
  thres_neg = -1*thres_pos
  if(sum(result$log2FoldChange < thres_neg) > 0){
    result[result$log2FoldChange < thres_neg, ]$log2FoldChange = 
      thres_neg + (result[result$log2FoldChange < thres_neg, ]$log2FoldChange - thres_neg)/10
  }
  
  result$Diff <- "NO"
  result$Diff[result$log2FoldChange > log2(thres_pos) &
                result$p.adjust < 0.01] <- "UP"
  result$Diff[result$log2FoldChange < -log2(thres_pos) & 
                result$p.adjust < 0.01] <- "DOWN"
  result$delabel <- NA
  #result[rownames(result) == 'PITX1',"delabel"] = 'PITX1'
  result$delabel[result$Diff != "NO"] <- rownames(result)[result$Diff != "NO"]
  
  mycolors <- c("#E41A1C", "#377EB8", "grey50")
  names(mycolors) <- c("UP", "DOWN", "NO")
  
  ggplot(data=result, aes(x=log2FoldChange, y=-log10(p.adjust), col=Diff, label=delabel)) + 
    geom_point(size = 1, alpha = 0.5) +
    xlim(-3,3) +
    ylim(0, max(-log10(result$p.adjust))) + 
    scale_color_manual(values=mycolors)+
    geom_text_repel(max.overlaps = 20, label.size = 0.1) +
    geom_vline(xintercept=c(-1*log2(thres_pos), log2(thres_pos)), col="black", linetype="dotted") +
    geom_hline(yintercept=2, col="black", linetype="dotted") +
    theme_classic() + 
    ggtitle('Epithelial') + 
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5))
  outfile = paste0("Epithelial.CD_VS_Healthy.regulon_diff.pdf")
  ggsave(filename = outfile, width = 5, height = 5, dpi = 600)
}


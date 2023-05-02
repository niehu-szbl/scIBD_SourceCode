library(dplyr)
library(Seurat)
library(paletteer)
library(ggpubr)
library(Seurat)
library(ggplot2)
library(EnhancedVolcano)

# prepare data
if(T){
  # load data
  seu = readRDS("all.clean.raw.rds")
  
  # update meta
  meta = readRDS("all.clean.meta.rds")
  
  # set meta
  seu@meta.data = meta[rownames(seu@meta.data), ]
  
  # filter genes
  if(T){
    # check
    # grep("^IG[HKL]",rownames(seu), value = T)
    # grep("^IGH",rownames(seu), value = T)
    # substr(grep("^IGH",rownames(seu), value = T), 4,4) %>% unique # "M" "A" "E" "G" "D" "J" "V"
    # grep("^IGK",rownames(seu), value = T)
    # substr(grep("^IGK",rownames(seu), value = T), 4,4) %>% unique # "C" "J" "V"
    # grep("^IGL",rownames(seu), value = T)
    # substr(grep("^IGL",rownames(seu), value = T), 4,4) %>% unique # "O" "V" "J" "C" "L"
    
    exclude = c("IGHM","IGHGA1","IGHGA2","IGHE","IGHG1","IGHG2","IGHG3","IGHG4","IGHD")
    exclude = c(exclude, grep( "^IG[HKL][VJC]", rownames(seu), perl = T, value = T))
    exclude = c(exclude, grep( "^TR[AB][VJC]", rownames(seu), perl = T, value = T))
    exclude = c(exclude, grep("^TR[GD][VJC]", rownames(seu), perl = T, value =T))
    exclude = c(exclude, grep( "^A[LPC]\\d{5,}", rownames(seu), perl = T, value = T))
    exclude = c(exclude, grep( "^LINC", rownames(seu), perl = T, value = T))
    exclude = c(exclude, grep( "^RP[LS]", rownames(seu), perl = T, value = T))
    exclude = c(exclude, grep( "^MT-", rownames(seu), perl = T, value = T))
    exclude = c(exclude, grep( "^MTRNR", rownames(seu), perl = T, value = T))
    #exclude = c(exclude, grep( "^HSP", rownames(seu), perl = T, value = T))
    exclude = exclude %>% as.vector() %>% unique()
    exclude %>% length
    
    
    # subset
    genes = rownames(seu)[!(rownames(seu) %in% exclude)]
    seu = seu[genes]
  }
  
  # normalize data
  seu <- NormalizeData(seu)
}

# DEG analysis
if(T){
  
  # get data
  cells = seu@meta.data %>% filter(minor_cluster %in% c("CD8+ Tc17")) %>%
    filter(disease %in% c("CD_inflamed","UC_inflamed")) %>% 
    filter(stage == "adult") %>% 
    filter(tissue %in% c("smallInt","largeInt"))
  cells %>% group_by(disease) %>% summarise(count =n())
  
  # subset
  seu.sub = subset(seu, cells = cells %>% rownames)
  
  # filter low expressed genes
  nrow(seu.sub) # 15503
  sum(rowSums(seu.sub@assays$RNA@data > 0) >= 3) # 11499
  sum(rowSums(seu.sub@assays$RNA@data > 0) >= 10) # 10157
  seu.sub <- seu.sub[rowSums(seu.sub@assays$RNA@data > 0) >= 10, ]
  seu.sub
  
  # identify DEG
  Idents(seu.sub) <- 'disease'
  degs <- FindMarkers(seu.sub, ident.1 = "UC_inflamed", ident.2 = "CD_inflamed",
                      logfc.threshold = 0, only.pos = FALSE)
  write.table(x = degs, file = "DEG_analysis_Tc17_vs_Th17.degs.txt", sep = "\t", col.names = T, row.names = T, quote = F)
  
  fc = 2
  p_val_adj_cutoff = 0.05
  
  deg = degs
  # explore up-regulated and down-regulated genes
  # up_genes = deg %>% filter(p_val_adj < p_val_adj_cutoff & avg_log2FC >= log2(fc)) %>% arrange(desc(avg_log2FC)) %>% rownames()
  # write.table(up_genes, file = "up_genes.txt", quote = F, sep = "\t", col.names = F, row.names = F)
  
  # down_genes = deg %>% filter(p_val_adj < p_val_adj_cutoff & avg_log2FC <= -1*log2(fc)) %>% arrange(avg_log2FC) %>% rownames()
  # write.table(down_genes, file = "down_genes.txt", quote = F, sep = "\t", col.names = F, row.names = F)
  
  title = "UC_inflamed VS CD_inflamed"
  
  ncell = seu.sub$disease %>% droplevels() %>% table %>% as_tibble() 
  colnames(ncell) = c("group",'count')
  ncell
  
  keyvals <- ifelse(
    degs$avg_log2FC < -log2(fc) & degs$p_val_adj < p_val_adj_cutoff, 'royalblue',
    ifelse(degs$avg_log2FC > log2(fc) & degs$p_val_adj < p_val_adj_cutoff, 'red2','grey90'))
  names(keyvals)[keyvals == 'red2'] <- 'Up'
  names(keyvals)[keyvals == 'black'] <- 'NotSig'
  names(keyvals)[keyvals == 'royalblue'] <- 'Down'
  
  up_regulated_genes = degs %>% filter(p_val_adj < 0.05 & avg_log2FC > log2(fc)) %>% top_n(8, wt=avg_log2FC) %>% rownames()
  down_regulated_genes = degs %>% filter(p_val_adj < 0.05 & avg_log2FC < -1*log2(fc)) %>% top_n(8, wt=-1*avg_log2FC) %>% rownames()
  other = c("LAYN","PDCD1","TIGIT","IL17A","MYADM","CTLA4","IL21R","DUSP4","NR4A2")
  EnhancedVolcano(degs, lab = rownames(degs), x = 'avg_log2FC', y = 'p_val_adj', 
                  pointSize = 2,
                  title = title,
                  subtitle = paste0('number of cells ',ncell[2,2], ' vs ',ncell[1,2]), 
                  #ylim = c(0,5), 
                  #xlim = c(-4,4),
                  colCustom = keyvals,
                  drawConnectors = TRUE,
                  arrowheads = F,
                  colAlpha = 1,
                  # widthConnectors = 12,
                  selectLab = c(up_regulated_genes, down_regulated_genes, other) %>% unique,
                  pCutoff = p_val_adj_cutoff,
                  FCcutoff= log2(fc),
                  labSize = 5,
                  labFace = "italic",
                  labCol = 'black'
  ) + theme_classic() + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "Helvetica", size = 16,  color = "black"),
          axis.text.y = element_text(family = "Helvetica", size = 16, color = 'black'),
          axis.title.y = element_text(family = "Helvetica",size = 16, color = 'black'),
          legend.text = element_text(family = "Helvetica",size = 16, color ='black'),
          legend.title  = element_text(family = "Helvetica",size = 16, color ='black'),
          strip.text = element_text(family = "Helvetica",size = 16, color ='black'),
          plot.title = element_text(family = "Helvetica",size = 18, color ='black', hjust = 0.5)
    )

  # save plot
  # ggsave("Tc17.UC_vs_CD.pdf", width = 5, height = 5)
}

# calculate AUC score
# Exhuastion
if(T){
  
  library(GSEABase)
  library(AUCell)
  
  # rank
  seu.sub = subset(seu, cells = seu@meta.data[seu$minor_cluster == "CD8+ Tc17",] %>% rownames)
  cells_rankings <- AUCell_buildRankings(seu.sub@assays$RNA@data, plotStats=T)
  
  # gene set
  geneSets = c('HAVCR2','CXCL13','PDCD1','LAYN','TOX','IFNG','GZMB','TNFRSF9','LAG3','ENTPD1','TIGIT')
  geneSets <- GeneSet(geneSets, setName="Exhaustion")
  
  # score
  cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings)
  
  # plot
  #cells_assignment <- AUCell_exploreThresholds(cells_AUC, plotHist=TRUE, assign=TRUE)
  
  # Tex score, 2 disease status
  # group by disease and sample
  if(T){
      seu.sub$Tex_score = as.vector(cells_AUC@assays@data$AUC)
      
      # select samples
      cts = meta %>% filter(major_cluster == "CD8T") %>% 
                      filter(stage == "adult") %>%
                        filter(tissue %in% c("smallInt","largeInt","blood")) %>% 
                          filter(disease %in% c("CD_inflamed","UC_inflamed")) %>% 
                          group_by(sample) %>% summarise(count = n())
      samples = cts[cts$count > 100, "sample"] %>% unlist %>% as.vector
      samples %>% length
      samples %>% head
      
      # check number of samples per group
      colnames(meta)
      meta %>% dplyr::filter(sample %in% samples) %>% dplyr::select(sample, disease)  %>% dplyr::distinct() %>% 
        group_by(disease) %>% count()
      
      seu.sub.filter = subset(seu.sub, cells = seu.sub@meta.data %>% filter(sample %in% samples) %>% rownames )
      seu.sub.filter
      
      tmp = seu.sub@meta.data[,c("sample", "disease", "Tex_score")]
      tmp = tmp %>% group_by(sample, disease) %>% summarise(avg = mean(Tex_score))
      tmp = tmp %>% filter(disease %in% c("UC_inflamed","CD_inflamed"))
      
      # save result
      write.table(tmp, file = "Tex_score_of_Tc17_in_UC_CD.txt", sep = "\t",
                  quote = F, col.names = T, row.names = F)
      
      my_comparisons = list(c("UC_inflamed","CD_inflamed"))
      ggplot(tmp, aes(x=disease, y = avg, fill = disease, color = disease)) + 
        geom_boxplot( outlier.size = -1, linewidth =1) +
        geom_point( position = position_jitter(width = 0.1), size = 1.5) + 
        scale_fill_manual(values = alpha(paletteer::paletteer_d("ggsci::category10_d3")[c(1,3)], 0.6) ) + 
        scale_color_manual(values = paletteer::paletteer_d("ggsci::category10_d3")[c(1,3)] ) + 
        scale_x_discrete(labels = c("CD","UC")) + 
        labs(x="", y= "Exhuastion score", title = "By sample") + 
        theme_classic2() + 
        theme(
          legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
          axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 12, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 12, color ='black'),
          strip.text = element_text(family = "ArialMT",face = "italic", size = 16, color ='black'),
          strip.background = element_blank(),
          plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
        ) + stat_compare_means(comparisons = my_comparisons, 
                               method = "t.test", label = "p.signif", 
                               method.args = list(alternative = "greater")) 
      # save plot
      # ggsave("Tex_AUC_score.two_disease_status.boxplot.pdf", width = 1.8, height = 3)}
}


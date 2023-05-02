library(Seurat)
library(ggplot2)
library(dplyr)
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
  cells = seu@meta.data %>% dplyr::filter(minor_cluster %in% c("CD8+ Tc17","CD4+ Th17"))
  cells %>% dplyr::group_by(disease) %>% dplyr::summarise(count =n())
  cells %>% dplyr::group_by(minor_cluster) %>% dplyr::summarise(count =n())
  
  # subset
  seu.sub = subset(seu, cells = cells %>% rownames)
  
  # filter low expressed genes
  nrow(seu.sub) # 15503
  sum(rowSums(seu.sub@assays$RNA@data > 0) >= 3) # 13258
  sum(rowSums(seu.sub@assays$RNA@data > 0) >= 10) # 11823
  seu.sub <- seu.sub[rowSums(seu.sub@assays$RNA@data > 0) >= 10, ]
  seu.sub
  
  # identify DEG
  Idents(seu.sub) <- 'minor_cluster'
  degs <- FindMarkers(seu.sub, ident.1 = "CD8+ Tc17", ident.2 = "CD4+ Th17", 
                      logfc.threshold = 0, only.pos = FALSE)
  nrow(degs) # 3964
  write.table(x = degs, file = "DEG_analysis_Tc17_vs_Th17.degs.txt", sep = "\t", col.names = T, row.names = T, quote = F)
  
  degs = read.table("./DEG_analysis_Tc17_vs_Th17.degs.txt", header = T, sep = "\t")
  deg = degs
  fc = 2
  p_val_adj_cutoff = 0.05
  
  # explore up-regulated and down-regulated genes
  # up_genes = deg %>% filter(p_val_adj < p_val_adj_cutoff & avg_log2FC >= log2(fc)) %>% arrange(desc(avg_log2FC)) %>% rownames()
  # write.table(up_genes, file = "up_genes.txt", quote = F, sep = "\t", col.names = F, row.names = F)
  
  # down_genes = deg %>% filter(p_val_adj < p_val_adj_cutoff & avg_log2FC <= -1*log2(fc)) %>% arrange(avg_log2FC) %>% rownames()
  # write.table(down_genes, file = "down_genes.txt", quote = F, sep = "\t", col.names = F, row.names = F)
  
  title = "CD8+ Tc17 vs CD4+ Th17"
  ncell = seu.sub$minor_cluster %>% droplevels() %>% table %>% as_tibble() 
  colnames(ncell) = c("group",'count')
  ncell

  keyvals <- ifelse(
    degs$avg_log2FC < -log2(fc) & degs$p_val_adj < p_val_adj_cutoff, 'royalblue',
    ifelse(degs$avg_log2FC > log2(fc) & degs$p_val_adj < p_val_adj_cutoff, 'red2','grey90'))
  names(keyvals)[keyvals == 'red2'] <- 'Up'
  names(keyvals)[keyvals == 'black'] <- 'NotSig'
  names(keyvals)[keyvals == 'royalblue'] <- 'Down'
  
  up_regulated_genes = degs %>% filter(p_val_adj < p_val_adj_cutoff & avg_log2FC > log2(fc)) %>% top_n(20, wt=avg_log2FC) %>% rownames()
  down_regulated_genes = degs %>% filter(p_val_adj < p_val_adj_cutoff & avg_log2FC < -1*log2(fc)) %>% top_n(20, wt=-1*avg_log2FC) %>% rownames()
  other = c()
  EnhancedVolcano(degs, lab = rownames(degs), x = 'avg_log2FC', y = 'p_val_adj', 
                  pointSize = 2,
                  title = title,
                  subtitle = paste0('number of cells ',ncell[2,2], ' vs ',ncell[1,2]), 
                  #ylim = c(0,5), 
                  #xlim = c(-2,2),
                  colCustom = keyvals,
                  drawConnectors = TRUE,
                  arrowheads = F,
                  colAlpha = 1,
                  # widthConnectors = 12,
                  selectLab = c(up_regulated_genes, down_regulated_genes, other) %>% unique,
                  pCutoff = p_val_adj_cutoff,
                  FCcutoff= log2(fc),
                  labSize = 4,
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
  # ggsave("Tc17_vs_Th17.deg.pdf", width = 6, height = 4.5)
  
}
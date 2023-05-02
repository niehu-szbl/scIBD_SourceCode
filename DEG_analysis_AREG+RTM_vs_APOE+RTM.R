library(clusterProfiler)
library(org.Hs.eg.db)

# set working directory
setwd("/data1/niehu/ibd_public_data_20210821/analysis_20220111/02.integrate/combine_clean_V3/myeloid_v3")

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
  cells = seu@meta.data %>% filter(minor_cluster %in% c("AREG+ macrophage","APOE+ macrophage"))
  cells %>% group_by(disease) %>% summarise(count =n())
  cells %>% group_by(minor_cluster) %>% summarise(count =n())
  
  # subset
  seu.sub = subset(seu, cells = cells %>% rownames)
  
  # filter low expressed genes
  nrow(seu.sub) # 15503
  sum(rowSums(seu.sub@assays$RNA@data > 0) >= 3) # 13658
  sum(rowSums(seu.sub@assays$RNA@data > 0) >= 10) # 12221
  seu.sub <- seu.sub[rowSums(seu.sub@assays$RNA@data > 0) >= 10, ]
  seu.sub
  
  # identify DEG
  Idents(seu.sub) <- 'minor_cluster'
  degs <- FindMarkers(seu.sub, ident.1 = "AREG+ macrophage", ident.2 = "APOE+ macrophage", 
                      logfc.threshold = 0, only.pos = FALSE)
  write.table(x = degs, file = "DEG_analysis_AREG+RTM_vs_APOE+RTM.degs.txt", sep = "\t", col.names = T, row.names = T, quote = F)
  
  degs = read.table("./DEG_analysis_AREG+RTM_vs_APOE+RTM.degs.txt", header = T, sep = "\t")
  fc = 2
  p_val_adj_cutoff = 0.05
  
  # explore up-regulated and down-regulated genes
  # up_genes = deg %>% filter(p_val_adj < p_val_adj_cutoff & avg_log2FC >= log2(fc)) %>% arrange(desc(avg_log2FC)) %>% rownames()
  # write.table(up_genes, file = "up_genes.txt", quote = F, sep = "\t", col.names = F, row.names = F)
  
  # down_genes = deg %>% filter(p_val_adj < p_val_adj_cutoff & avg_log2FC <= -1*log2(fc)) %>% arrange(avg_log2FC) %>% rownames()
  # write.table(down_genes, file = "down_genes.txt", quote = F, sep = "\t", col.names = F, row.names = F)
  
  title = "AREG+ macrophage vs APOE+ macrophage"
  
  ncell = seu.sub$minor_cluster %>% droplevels() %>% table %>% as_tibble() 
  colnames(ncell) = c("group",'count')
  ncell
  
  library(EnhancedVolcano)
  keyvals <- ifelse(
    degs$avg_log2FC < -log2(fc) & degs$p_val_adj < p_val_adj_cutoff, 'royalblue',
    ifelse(degs$avg_log2FC > log2(fc) & degs$p_val_adj < p_val_adj_cutoff, 'red2','grey90'))
  names(keyvals)[keyvals == 'red2'] <- 'Up'
  names(keyvals)[keyvals == 'black'] <- 'NotSig'
  names(keyvals)[keyvals == 'royalblue'] <- 'Down'
  
  up_regulated_genes = degs %>% filter(p_val_adj < 0.05 & avg_log2FC > log2(fc)) %>% top_n(10, wt=avg_log2FC) %>% rownames()
  down_regulated_genes = degs %>% filter(p_val_adj < 0.05 & avg_log2FC < -1*log2(fc)) %>% top_n(10, wt=-1*avg_log2FC) %>% rownames()
  other = c("THBS1","NLRP3","CXCR4","CCL4L2","CPM","FCN1","CLEC10A","CD55")
  selected_genes = c(up_regulated_genes, down_regulated_genes, other) %>% unique
  selected_genes = selected_genes[!(selected_genes %in% c("MALAT1","TSC22D3"))]
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
                  selectLab = selected_genes,
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
  # ggsave("AREG_RTM_VS_APOE_RTM.deg.pdf", width = 5, height = 4.5)
  
}

# GO enrichment analysis by GSEA
if(T){
  
  ego_result = list()
  ego_result_names = c()
  
  # Generate input for GSEA
  deg = read.table("./DEG_analysis_AREG+RTM_vs_APOE+RTM.degs.txt", header = T, sep = "\t")
  geneList = deg$avg_log2FC
  names(geneList) = rownames(deg)
  geneList = geneList %>% sort(decreasing = T)
  
  for(onto in c("BP","CC","MF")){
    ego <- gseGO(geneList     = geneList,
                 OrgDb        = org.Hs.eg.db,
                 ont          = onto,
                 keyType = 'SYMBOL',
                 minGSSize    = 10,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 verbose      = FALSE)
    saveRDS(ego, file = paste0("deg.gsea.",onto,".rds"))
    write.table(ego, file = paste0("deg.gsea.",onto,".txt"), 
                sep = "\t", quote = F, col.names = T)
    ego_result = c(ego_result, ego)
    ego_result_names = c(ego_result_names, onto)
  }
  
  names(ego_result) = ego_result_names
  saveRDS(ego_result, file = "deg.all.gsea_result.rds")
}

# plot, GSEA
if(T){
  
  # no need to simplify
  # just select go terms in specific level
  ego = readRDS("./deg.gsea.BP.rds")
  
  # simplify
  ego <- enrichplot::pairwise_termsim(ego)
  ego@result %>% nrow
  ego <- simplify(ego, cutoff=0.7, by="p.adjust", select_fun=min)
  ego@result %>% nrow
  
  # check
  ego@result %>% nrow
  ego@result$ID %>% unique %>% length
  colnames(ego@result)
  
  # get go levels
  go_levels = GOxploreR::GOTermBPOnLevel(ego@result$ID)
  go_levels %>% nrow
  go_levels$Term %>% unique %>% length
  colnames(go_levels)
  colnames(go_levels)[1] = 'ID'
  
  # combine
  result = left_join(ego@result, go_levels, by = "ID")
  result %>% head
  #result = result %>% filter(setSize < 350)
  result$Level %>% table
  #result = result %>% filter(Level %in% 4:7)
  result = result %>% filter(!(ID %in% c("GO:0045934","GO:0010558","GO:0045935")))
  result.up = result %>% filter(NES > 0)
  result.up %>% nrow
  
  result.down = result %>% filter(NES < 0)
  result.down %>% nrow
  
  result.up.show = result.up %>% filter(p.adjust < 0.01) %>% arrange(desc(NES)) %>% head(10)
  result.down.show = result.down %>% filter(p.adjust < 0.01) %>% arrange(NES) %>% head(10) %>% arrange(desc(NES))
  result.show = rbind(result.up.show, result.down.show)
  result.show = result.show[nrow(result.show):1,]
  result.show$Description = factor(result.show$Description, 
                                   levels = result.show$Description)
  result.show$group = 'UP'
  result.show[result.show$NES < 0, ]$group = 'DOWN'
  
  # save result.show
  write.table(result.show, file = "GSEA_on_GO_terms.AREG_RTM_vs_APOE_RTM.txt", 
              sep = "\t", col.names = T, row.names = F, quote = F)
  
  ggplot(result.show, aes(x = Description, y = NES, fill = group)) + 
    geom_col() +
    scale_fill_manual(values = paletteer::paletteer_d("ggsci::category20c_d3")[c(6,7)]) + 
    labs(y="NES", x= "", title ="GSEA") + 
    coord_flip() + 
    theme_classic() + 
    theme(legend.position = "none",
          axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
          axis.text.y = element_text(family = "ArialMT", size = 12, color = 'black'),
          axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
          legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
          legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
          strip.text = element_text(family = "ArialMT",size = 16, color ='black'),
          plot.title = element_text(family = "ArialMT",size = 18, color ='black', hjust = 0.5)
    )

  # save plot
  # ggsave("trm_gsea.pdf", width = 6.5, height = 4)
}
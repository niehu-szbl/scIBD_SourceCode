library(ggtree)
library(aplot)
library(GSVA)
library(Seurat)
library(reshape2)
library(ggplot2)
library(aplot)
library(RColorBrewer)
library(scales)
library(paletteer)
library(dplyr)

# prepare data
if(T){
  # load data
  seu = readRDS("all.clean.raw.rds")

  # normalize
  seu = NormalizeData(seu, normalization.method = "RC", scale.factor = 10000)

  # check data
  seu@assays$RNA@data %>% colSums %>% head # CP10K

  # update meta
  meta = readRDS("all.clean.meta.rds")

  # set meta
  seu@meta.data = meta[rownames(seu@meta.data), ]
}

# heatmap plot
# subtype-specific GWAS risk genes
# Figure 2c
if(T){
  
  genes = read.table("./GWAS_risk_genes.hit_gene.minor.pick.txt", header = F, stringsAsFactors = F)
  genes = genes$V1
  subset_seu = seu[genes]
  
  # get expression
  expression = data.frame(t(as.matrix(subset_seu@assays$RNA@data)))
  colnames(expression) = rownames(subset_seu)
  expression$label = subset_seu$minor_cluster
  expression$sample = subset_seu$sample
  expression = melt(expression)
  colnames(expression) <- c("minor_cluster","sample","gene","expression")
  
  # expression, long format
  avg_exp = expression %>%
    dplyr::group_by(gene, minor_cluster, sample) %>% 
    dplyr::summarise(avg_exp = mean(expression)) %>% 
    dplyr::group_by(gene, minor_cluster) %>%
    dplyr::summarise(avg_exp = mean(avg_exp))
  head(avg_exp)
  
  # expression, matrix
  mtx = dcast(avg_exp,  gene ~ minor_cluster, value.var = "avg_exp")
  rownames(mtx) = mtx[,1]
  mtx = mtx[,-1]
  #mtx = mtx %>% t %>% scale %>% t %>% as.data.frame() %>% mutate(gene = rownames(.))
  mtx = apply(mtx, MARGIN = 1, function(x) rescale(x, to=c(0,1)) ) %>% t %>% as.data.frame() %>% dplyr::mutate(gene = rownames(.))
  
  # expression, long format, scaled
  avg_exp_scaled = mtx %>% melt
  colnames(avg_exp_scaled) <- c("gene","minor_cluster","avg_exp")
  
  major_minor_mapping = subset_seu@meta.data %>%
    dplyr::select(minor_cluster, major_cluster) %>% 
    unique %>% 
    dplyr::mutate(group = major_cluster, p = "group")
  
  p = ggplot(avg_exp_scaled, aes(x = minor_cluster, y = gene)) +
    geom_tile(aes(fill = avg_exp), color = "white", lwd = 0.25, linetype = 0)+
    #scale_fill_viridis_c("Expression", option = "viridis") +
    scale_fill_paletteer_c("grDevices::Blues 3", direction = -1) + 
    xlab("")+
    ylab("")+
    theme_bw() + 
    theme(axis.text.x = element_text(size = 8,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 8,  color = "black", face = "italic"),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size=12, color = "black"),
          axis.title.y = element_blank(),  
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent', color='white'),
          panel.border = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"),
          legend.position = "right") + 
    geom_vline(xintercept = 0.5 + major_minor_mapping$major_cluster %>% table %>% cumsum,  size = 0.5, color = "white")
  
  p2 = ggplot(major_minor_mapping, aes(y = p, x = minor_cluster, fill = major_cluster)) + 
    geom_tile(aes(fill = major_cluster)) +
    scale_fill_manual(values = color_palette_9) + 
    scale_y_discrete(position = "right")+ 
    xlab(NULL)+
    ylab(NULL) + 
    theme_minimal() + 
    labs(fill = "Cell type") + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  mtx = mtx[,-ncol(mtx)]
  #phc <- hclust(dist(t(mtx))) %>% 
  #  ggtree() + layout_dendrogram()
  
  phr <- hclust(dist(mtx)) %>% 
    ggtree(layout="rectangular",branch.length="none")
  
  p %>%
    insert_top(p2, height = 0.05) %>%
    #insert_top(phc, height = 0.1) %>%
    insert_left(phr, width = 0.1)

  # save plot
  # ggsave("heatmap.risk_genes.minor_deg.pdf", height = 6, width = 12)
}

# heatmap plot
# gene expression of genes
# FDA-arrpoval target genes
# Figure 2d
if(T){
  
  # genes
  genes = ["TNF","ITGA4","ITGB7","S1PR1","S1PR5","IL23R","IL12RB1","IL12RB2","JAK1","JAK3"]
  subset_seu = seu[genes]
  
  # get expression
  expression = data.frame(t(as.matrix(subset_seu@assays$RNA@data)))
  colnames(expression) = rownames(subset_seu)
  expression$label = subset_seu$minor_cluster
  expression$sample = subset_seu$sample
  expression = melt(expression)
  colnames(expression) <- c("minor_cluster","sample","gene","expression")
  
  # expression, long format
  avg_exp = expression %>% 
    group_by(gene, minor_cluster, sample) %>% 
    summarise(avg_exp = mean(expression)) %>% 
    group_by(gene, minor_cluster) %>%
    summarise(avg_exp = mean(avg_exp))
  
  # expression, matrix
  mtx = dcast(avg_exp,  gene ~ minor_cluster, value.var = "avg_exp")
  rownames(mtx) = mtx[,1]
  mtx = mtx[,-1]
  #mtx = mtx %>% t %>% scale %>% t %>% as.data.frame() %>% mutate(gene = rownames(.))
  mtx = apply(mtx, MARGIN = 1, function(x) rescale(x, to=c(0,1)) ) %>% t %>% as.data.frame() %>% mutate(gene = rownames(.))
  
  # expression, long format, scaled
  avg_exp_scaled = mtx %>% melt
  colnames(avg_exp_scaled) <- c("gene","minor_cluster","avg_exp")
  
  major_minor_mapping = subset_seu@meta.data %>%
    dplyr::select(minor_cluster, major_cluster) %>% 
    unique %>%
    mutate(group = major_cluster, p = "group")
  
  avg_exp_scaled$gene = factor(avg_exp_scaled$gene, levels = genes %>% rev)
  p = ggplot(avg_exp_scaled, aes(x = minor_cluster, y = gene)) +
    geom_tile(aes(fill = avg_exp), color = "white", lwd = 0.25, linetype = 0)+
    #scale_fill_viridis_c("Expression", option = "viridis") +
    scale_fill_paletteer_c("grDevices::Blues 3", direction = -1) + 
    xlab("")+
    ylab("")+
    theme_bw() + 
    theme(axis.text.x = element_text(size = 8,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 8,  color = "black", face = "italic"),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size=12, color = "black"),
          axis.title.y = element_blank(),  
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent', color='white'),
          panel.border = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"),
          legend.position = "right") + 
    geom_vline(xintercept = 0.5 + major_minor_mapping$major_cluster %>% table %>% cumsum,  size = 0.5, color = "white")
  
  p2 = ggplot(major_minor_mapping, aes(y = p, x = minor_cluster, fill = major_cluster)) + 
    geom_tile(aes(fill = major_cluster)) +
    scale_fill_manual(values = color_palette_9) + 
    scale_y_discrete(position = "right")+ 
    xlab(NULL)+
    ylab(NULL) + 
    theme_minimal() + 
    labs(fill = "Cell type") + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  mtx = mtx[,-ncol(mtx)]
  # phc <- hclust(dist(t(mtx))) %>% 
  #   ggtree() + layout_dendrogram()
  
  phr <- hclust(dist(mtx)) %>% 
    ggtree(layout="rectangular",branch.length="none")
  
  p %>%
    insert_top(p2, height = 0.05) #%>%
    #insert_top(phc, height = 0.1) %>%
    #insert_left(phr, width = 0.1)

  # save plot
  ggsave("heatmap.drug_targets.minor_deg.pdf", height = 3.5, width = 12)
}

# heatmap plot
# major cluster-specific GWAS risk genes
# Extended Data Figure 7a
if(T){
  
  genes = read.table("./GWAS_risk_genes.hit_gene.major.pick.txt", header = F, stringsAsFactors = F)
  genes = genes$V1
  subset_seu = seu[genes]
  
  # get expression
  expression = data.frame(t(as.matrix(subset_seu@assays$RNA@data)))
  colnames(expression) = rownames(subset_seu)
  expression$label = subset_seu$minor_cluster
  expression$sample = subset_seu$sample
  expression = melt(expression)
  colnames(expression) <- c("minor_cluster","sample","gene","expression")
  
  # expression, long format
  avg_exp = expression %>% 
    group_by(gene, minor_cluster, sample) %>% 
    summarise(avg_exp = mean(expression)) %>% 
    group_by(gene, minor_cluster) %>%
    summarise(avg_exp = mean(avg_exp))
  
  # expression, matrix
  mtx = dcast(avg_exp,  gene ~ minor_cluster, value.var = "avg_exp")
  rownames(mtx) = mtx[,1]
  mtx = mtx[,-1]
  #mtx = mtx %>% t %>% scale %>% t %>% as.data.frame() %>% mutate(gene = rownames(.))
  mtx = apply(mtx, MARGIN = 1, function(x) rescale(x, to=c(0,1)) ) %>% t %>% as.data.frame() %>% mutate(gene = rownames(.))
  
  # expression, long format, scaled
  avg_exp_scaled = mtx %>% melt
  colnames(avg_exp_scaled) <- c("gene","minor_cluster","avg_exp")
  
  major_minor_mapping = subset_seu@meta.data %>%
    dplyr::select(minor_cluster, major_cluster) %>% 
    dplyr::distinct() %>% 
    dplyr::mutate(group = major_cluster, p = "group")
  
  head(avg_exp_scaled)
  write.table(avg_exp_scaled, "Subtype_specific_risk_genes.heatmap.txt",
              sep = "\t", col.names = T, row.names = F, quote = F)
  
  p = ggplot(avg_exp_scaled, aes(x = minor_cluster, y = gene)) +
    geom_tile(aes(fill = avg_exp), color = "white", lwd = 0.25, linetype = 0)+
    #scale_fill_viridis_c("Expression", option = "viridis") +
    scale_fill_paletteer_c("grDevices::Blues 3", direction = -1) + 
    xlab("")+
    ylab("")+
    theme_bw() + 
    theme(axis.text.x = element_text(size = 8,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 8,  color = "black", face = "italic"),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size=12, color = "black"),
          axis.title.y = element_blank(),  
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent', color='white'),
          panel.border = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"),
          legend.position = "right") + 
    geom_vline(xintercept = 0.5 + major_minor_mapping$major_cluster %>% table %>% cumsum,  size = 0.5, color = "white")
  
  p2 = ggplot(major_minor_mapping, aes(y = p, x = minor_cluster, fill = major_cluster)) + 
    geom_tile(aes(fill = major_cluster)) +
    scale_fill_manual(values = color_palette_9) + 
    scale_y_discrete(position = "right")+ 
    xlab(NULL)+
    ylab(NULL) + 
    theme_minimal() + 
    labs(fill = "Cell type") + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  mtx = mtx[,-ncol(mtx)]
  #phc <- hclust(dist(t(mtx))) %>% 
  #  ggtree() + layout_dendrogram()
  
  phr <- hclust(dist(mtx)) %>% 
    ggtree(layout="rectangular",branch.length="none")
  
  p %>%
    insert_top(p2, height = 0.05) %>%
    #insert_top(phc, height = 0.1) %>%
    insert_left(phr, width = 0.1)
  # ggsave("heatmap.risk_genes.major_deg.pdf", height = 3.5, width = 12)
}

# heatmap, minor cluster
# subtype-specific therapy targets
# Extended Data Figure 8a
if(T){
  
  genes = read.table("./Therapy_targets.hit_gene.minor.pick.txt", header = F, stringsAsFactors = F)
  genes = genes$V1
  length(genes)
  subset_seu = seu[genes]
  
  # get expression
  expression = data.frame(t(as.matrix(subset_seu@assays$RNA@data)))
  colnames(expression) = rownames(subset_seu)
  expression$label = subset_seu$minor_cluster
  expression$sample = subset_seu$sample
  expression = melt(expression)
  colnames(expression) <- c("minor_cluster","sample","gene","expression")
  
  # expression, long format
  avg_exp = expression %>% 
    group_by(gene, minor_cluster, sample) %>% 
    summarise(avg_exp = mean(expression)) %>% 
    group_by(gene, minor_cluster) %>%
    summarise(avg_exp = mean(avg_exp))
  
  # expression, matrix
  mtx = dcast(avg_exp,  gene ~ minor_cluster, value.var = "avg_exp")
  rownames(mtx) = mtx[,1]
  mtx = mtx[,-1]
  #mtx = mtx %>% t %>% scale %>% t %>% as.data.frame() %>% mutate(gene = rownames(.))
  mtx = apply(mtx, MARGIN = 1, function(x) rescale(x, to=c(0,1)) ) %>% t %>% as.data.frame() %>% mutate(gene = rownames(.))
  
  # expression, long format, scaled
  avg_exp_scaled = mtx %>% melt
  colnames(avg_exp_scaled) <- c("gene","minor_cluster","avg_exp")
  
  major_minor_mapping = subset_seu@meta.data %>%
    select(minor_cluster, major_cluster) %>% 
    unique %>% 
    mutate(group = major_cluster, p = "group")
  
  p = ggplot(avg_exp_scaled, aes(x = minor_cluster, y = gene)) +
    geom_tile(aes(fill = avg_exp), color = "white", lwd = 0.25, linetype = 0)+
    #scale_fill_viridis_c("Expression", option = "viridis") +
    scale_fill_paletteer_c("grDevices::Blues 3", direction = -1) + 
    xlab("")+
    ylab("")+
    theme_bw() + 
    theme(axis.text.x = element_text(size = 8,  color = "black", angle = 90, vjust = 0.5, hjust = 1),
          axis.text.y = element_text(size = 8,  color = "black", face = "italic"),
          axis.ticks.y = element_blank(),
          axis.title.x = element_text(size=12, color = "black"),
          axis.title.y = element_blank(),  
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_rect(fill='transparent', color='white'),
          panel.border = element_rect(fill = NA, colour = "black", size = 1, linetype = "solid"),
          legend.position = "right") + 
    geom_vline(xintercept = 0.5 + major_minor_mapping$major_cluster %>% table %>% cumsum,  size = 0.5, color = "white")
  
  p2 = ggplot(major_minor_mapping, aes(y = p, x = minor_cluster, fill = major_cluster)) + 
    geom_tile(aes(fill = major_cluster)) +
    scale_fill_manual(values = color_palette_9) + 
    scale_y_discrete(position = "right")+ 
    xlab(NULL)+
    ylab(NULL) + 
    theme_minimal() + 
    labs(fill = "Cell type") + 
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  mtx = mtx[,-ncol(mtx)]
  #phc <- hclust(dist(t(mtx))) %>% 
  #  ggtree() + layout_dendrogram()
  
  phr <- hclust(dist(mtx)) %>% 
    ggtree(layout="rectangular",branch.length="none")
  
  p %>%
    insert_top(p2, height = 0.05) %>%
    #insert_top(phc, height = 0.1) %>%
    insert_left(phr, width = 0.1)
  ggsave("heatmap.therapy_targets_more.minor_deg.pdf", height = 6, width = 12)
}

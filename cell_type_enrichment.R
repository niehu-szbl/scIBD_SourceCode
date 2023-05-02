library(ggtree)
library(aplot)
library(GSVA)
library(Seurat)
library(reshape2)
library(dplyr)
library(ggplot2)
library(aplot)
library(RColorBrewer)
library(reshape2)
library(GeneOverlap)

# load data
seu = readRDS("all.clean.raw.rds")

# update meta
meta = readRDS("all.clean.meta.rds")

# set meta
seu@meta.data = meta[rownames(seu@meta.data), ]

# get expressed genes
if(T){
  select = rowSums(seu@assays$RNA@data>0) > 200
  expressed_genes = rownames(seu)[ select  ]
  expressed_genes %>% length
  write.table(expressed_genes, "scIBD.expressed_genes.txt", row.names = F, col.names = F, sep = "\t", quote = F)
}

# load expressed genes
if(T){
  expressed_genes = read.table("scIBD.expressed_genes.txt", header = F, stringsAsFactors = F, sep = "\t")
  expressed_genes = expressed_genes$V1
  expressed_genes %>% length
}

all_minor_clusters = seu$minor_cluster %>% levels()
mapping = seu@meta.data[,c("minor_cluster","major_cluster")] %>% distinct
mapping = mapping[!is.na(mapping$minor_cluster),]
mapping %>% nrow
rownames(mapping) = mapping$minor_cluster

# input
markers_df = read.table("./scIBD_markers.txt", header = T, stringsAsFactors = F, sep = "\t")
markers_df %>% nrow
markers_df = markers_df %>% filter(p_val_adj < 0.05 & avg_log2FC > log2(1.5))
markers_df %>% nrow
expressed_genes = read.table("./scIBD.expressed_genes.txt", header = F, stringsAsFactors = F, sep = "\t")
expressed_genes = expressed_genes$V1

run_gsea = function(seuObj, expressed_genes, input_geneset, markers_df){
  
  all_pval = c()
  all_ods = c()
  
  # test
  for( i in 1:length(all_minor_clusters) ){
    cs_deg = markers_df[markers_df$cluster == all_minor_clusters[i], ]$gene
    go.obj <- newGeneOverlap(cs_deg, input_geneset, genome.size=length(expressed_genes))
    go.obj <- testGeneOverlap(go.obj)
    pval = getPval(go.obj)
    ods = getOddsRatio(go.obj)
    all_pval = c(all_pval, pval)
    all_ods = c(all_ods, ods)
  }
  
  # generate result
  res = cbind("p-value"=all_pval, "odds ratio"= all_ods)
  rownames(res) = all_minor_clusters
  return(res)
}

# cell type enrichment of CD risk genes
if(T){
  
  # prepare gene set
  input_geneset = read.table("./gsct.cd.gene.txt", header = F, stringsAsFactors = F, sep = "\t")
  input_geneset = input_geneset$V1
  input_geneset = input_geneset[input_geneset %in% rownames(seu)]
  input_geneset %>% length # ibd 184,  cd 132, uc 139
  
  # run gsea
  if(T){
    res = run_gsea(seu, expressed_genes, input_geneset, markers_df)
    res = data.frame(res)
    res$cell_type = rownames(res)
    res %>% colnames
    res$major_cluster = mapping[res %>% rownames, "major_cluster"]
    res = res[all_minor_clusters, ]
    res$cell_type = factor(res$cell_type, levels = all_minor_clusters)
    res$major_cluster = factor(res$major_cluster, levels = res$major_cluster %>% unique)
    res$padj = round(p.adjust(res$p.value, "BH"), 5)
    res = res[,c("padj","odds.ratio","cell_type","major_cluster")]
    
    data = res
    colnames(data) = c("padj","value","label","group")
    data$value = data$value + 0.05
  }
  
  # Complex circular bar plot
  if(T){
    # Set a number of 'empty bar' to add at the end of each group
    empty_bar <- 3
    to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
    colnames(to_add) <- colnames(data)
    to_add$group <- rep(levels(data$group), each=empty_bar)
    data <- rbind(data, to_add)
    data <- data %>% arrange(group)
    data$id <- seq(1, nrow(data))
    
    # prepare a data frame for base lines
    base_data <- data %>% 
      group_by(group) %>% 
      summarize(start=min(id), end=max(id) - empty_bar) %>% 
      rowwise() %>% 
      mutate(title=mean(c(start, end)))
    
    # Get the name and the y position of each label
    label_data <- data
    number_of_bar <- nrow(label_data)
    angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar
    label_data$hjust <- ifelse( angle < -90, 1, 0)
    label_data$angle <- ifelse(angle < -90, angle+180, angle)
    label_data[which(label_data$padj >= 0.01), "label"] = NA
    
    # load minor color
    major_color = read.table("/data1/niehu/ibd_public_data_20210821/analysis_20220111/02.integrate/combine_clean_V3/scIBD_figures/color/major_cluster_color.txt",
                             header = F, stringsAsFactors = F, sep = "\t", comment.char = ">")
    name = major_color$V1
    col = major_color$V2
    major_color = col
    names(major_color) = name
    
    # save output
    write.table(data, file = "cell_type_enrichment_of_CD_risk_genes.txt",
                sep = "\t", col.names = T, row.names = T, quote = F)
    
    # plot
    ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +
      geom_bar(stat="identity") +
      #ylim(-1, max(data$value)) +
      ylim(-0.7, max(data$value, na.rm = T)*1.2) + 
      scale_fill_manual(values = major_color) +
      theme_minimal() +
      theme(
        legend.position = "right",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-2,4), "cm")
      ) +
      coord_polar() +
      geom_text(data=label_data, aes(x=id, y=value + 0.1, label=label, hjust=hjust), 
                color="black", fontface="bold",alpha=0.8, size=2, 
                angle= label_data$angle, inherit.aes = FALSE ) +
      geom_segment(data = base_data, aes(x = start, y = -0.1, xend = end, yend = -0.1), 
                   colour = "black", alpha = 0.8, size = 0.25 , inherit.aes = FALSE ) #+
    #geom_text( data = base_data, aes(x = title, y = -0.5, label = group), 
    #           colour = "black", alpha=0.8, size=4, fontface="bold", inherit.aes = FALSE)
    # ggsave("GWAS_risk_gene.CD.circular_barplot.pdf", height = 7, width = 7)
  }
}

# cell type enrichment of UC risk genes
if(T){
  
  # prepare gene set
  input_geneset = read.table("./gsct.uc.gene.txt", header = F, stringsAsFactors = F, sep = "\t")
  input_geneset = input_geneset$V1
  input_geneset = input_geneset[input_geneset %in% rownames(seu)]
  input_geneset %>% length # ibd 184,  cd 132, uc 139
  
  # run gsea
  if(T){
    res = run_gsea(seu, expressed_genes, input_geneset, markers_df)
    res = data.frame(res)
    res$cell_type = rownames(res)
    res %>% colnames
    res$major_cluster = mapping[res %>% rownames, "major_cluster"]
    res = res[all_minor_clusters, ]
    res$cell_type = factor(res$cell_type, levels = all_minor_clusters)
    res$major_cluster = factor(res$major_cluster, levels = res$major_cluster %>% unique)
    res$padj = round(p.adjust(res$p.value, "BH"), 5)
    res = res[,c("padj","odds.ratio","cell_type","major_cluster")]
    
    data = res
    colnames(data) = c("padj","value","label","group")
    data$value = data$value + 0.05
  }
  
  # Complex circular bar plot
  if(T){
    # Set a number of 'empty bar' to add at the end of each group
    empty_bar <- 3
    to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group), ncol(data)) )
    colnames(to_add) <- colnames(data)
    to_add$group <- rep(levels(data$group), each=empty_bar)
    data <- rbind(data, to_add)
    data <- data %>% arrange(group)
    data$id <- seq(1, nrow(data))
    
    # prepare a data frame for base lines
    base_data <- data %>% 
      group_by(group) %>% 
      summarize(start=min(id), end=max(id) - empty_bar) %>% 
      rowwise() %>% 
      mutate(title=mean(c(start, end)))
    
    # Get the name and the y position of each label
    label_data <- data
    number_of_bar <- nrow(label_data)
    angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar
    label_data$hjust <- ifelse( angle < -90, 1, 0)
    label_data$angle <- ifelse(angle < -90, angle+180, angle)
    label_data[which(label_data$padj >= 0.01), "label"] = NA
    
    # load minor color
    major_color = read.table("/data1/niehu/ibd_public_data_20210821/analysis_20220111/02.integrate/combine_clean_V3/scIBD_figures/color/major_cluster_color.txt",
                             header = F, stringsAsFactors = F, sep = "\t", comment.char = ">")
    name = major_color$V1
    col = major_color$V2
    major_color = col
    names(major_color) = name
    
    # save output
    write.table(data, file = "cell_type_enrichment_of_UC_risk_genes.txt",
                sep = "\t", col.names = T, row.names = T, quote = F)
    
    # plot
    ggplot(data, aes(x=as.factor(id), y=value, fill=group)) +
      geom_bar(stat="identity") +
      #ylim(-1, max(data$value)) +
      ylim(-0.7, max(data$value, na.rm = T)*1.2) + 
      scale_fill_manual(values = major_color) +
      theme_minimal() +
      theme(
        legend.position = "right",
        axis.text = element_blank(),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(rep(-2,4), "cm")
      ) +
      coord_polar() +
      geom_text(data=label_data, aes(x=id, y=value + 0.1, label=label, hjust=hjust), 
                color="black", fontface="bold",alpha=0.8, size=2, 
                angle= label_data$angle, inherit.aes = FALSE ) +
      geom_segment(data = base_data, aes(x = start, y = -0.1, xend = end, yend = -0.1), 
                   colour = "black", alpha = 0.8, size = 0.25 , inherit.aes = FALSE )
    # ggsave("GWAS_risk_gene.CD.circular_barplot.pdf", height = 7, width = 7)
  }
}


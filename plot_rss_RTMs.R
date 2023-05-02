library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggpubr)

# plot rss
# APOE+ macrophage
# Figure 1i
if(T){
  # load regulon
  regulon = read.table("./Myeloid.rss.txt", header = T, sep = "\t")
  head(regulon)
  
  tmp = data.frame(tf = rownames(regulon), rss = regulon$APOE..macrophage) %>% arrange(desc(rss))
  tmp$rank = 1:nrow(tmp)
  tmp$color = 'A'
  tmp[1:10,]$color = 'B'
  tmp2 = tmp %>% head(10)
  ggplot(tmp, aes(x = rank, y = rss, color = color)) + 
    geom_point(size = 0.5) +
    scale_color_manual(values = c("#EE9572","#1874CD") %>% rev)+
    labs(x="Rank",y="Regulon specificity score",title="APOE+ RTM") + 
    geom_text_repel(data = tmp2, aes(x=rank, y = rss, label=tf), max.overlaps = 20)+
    theme_classic2() + 
    theme(
      legend.position = "none",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "italic", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 16, color ='black', hjust = 0.5)
    )

  # save plot
  ggsave("APOE_regulon.pdf", width = 3, height = 4.5)
}

# plot rss
# AREG+ macrophage
# Figure 1m
if(T){
  # load regulon
  regulon = read.table("./Myeloid.rss.txt", header = T, sep = "\t")
  head(regulon)
  
  tmp = data.frame(tf = rownames(regulon), rss = regulon$AREG..macrophage) %>% arrange(desc(rss))
  tmp$rank = 1:nrow(tmp)
  tmp$color = 'A'
  tmp[1:10,]$color = 'B'
  tmp2 = tmp %>% head(10)
  ggplot(tmp, aes(x = rank, y = rss, color = color)) + 
    geom_point(size = 0.5) +
    scale_color_manual(values = c("#EE9572","#1874CD") %>% rev)+
    labs(x="Rank",y="Regulon specificity score",title="APOE+ RTM") + 
    geom_text_repel(data = tmp2, aes(x=rank, y = rss, label=tf), max.overlaps = 20)+
    theme_classic2() + 
    theme(
      legend.position = "none",
      axis.text.x = element_text(family = "ArialMT", size = 16,  color = "black"),
      axis.text.y = element_text(family = "ArialMT", size = 16, color = 'black'),
      axis.title.y = element_text(family = "ArialMT",size = 16, color = 'black'),
      axis.title.x = element_text(family = "ArialMT",size = 16, color = 'black'),
      legend.text = element_text(family = "ArialMT",size = 16, color ='black'),
      legend.title  = element_text(family = "ArialMT",size = 16, color ='black'),
      strip.text = element_text(family = "ArialMT",face = "italic", size = 16, color ='black'),
      strip.background = element_blank(),
      plot.title = element_text(family = "ArialMT",size = 16, color ='black', hjust = 0.5)
    )

  # save plot
  ggsave("AREG_regulon.pdf", width = 3, height = 4.5)
}

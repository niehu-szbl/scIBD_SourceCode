library(ggplot2)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(dplyr)

# prepare data
if(T){
  # update meta
  meta = readRDS("all.clean.meta.rds")

  # load color
  minor_cluster.color = read.table("minor_cluster_color.txt", header = F, sep = "\t", comment.char = "@")
  tmp = minor_cluster.color$V2
  names(tmp) = minor_cluster.color$V1
  minor_cluster.color = tmp
}

# myeloid
# selected dataset
if(T){
  studies = meta$study %>% levels
  
  keep =  c("D. Fawkner-Corbett et al., Cell, 2021","C. S. Smillie et al., Cell, 2019","R. Elmentaite et al., Nature, 2021",
            "R. Elmentaite et al., Dev Cell, 2020","B. Huang et al., Cell, 2019","J. C. Martin et al., Cell, 2019",
            "B. S. Boland et al., Sci Immunol, 2020")
  tmp = meta %>% dplyr::filter(major_cluster == 'Myeloid') %>% 
    dplyr::filter(study %in% keep) %>% 
    dplyr::filter(stage %in% c("pediatric","adult")) %>% 
    dplyr::group_by(study, minor_cluster) %>% summarise(count = n())
  
  tmp$study = factor(tmp$study, levels = studies[tmp$study %in% studies] %>% rev)
  ggplot(tmp, aes(x = study, y = count, fill = minor_cluster)) + 
    geom_bar(position="fill", stat="identity") + 
    scale_fill_manual(values = minor_cluster.color[tmp$minor_cluster %>% droplevels() %>% levels]) + 
    #geom_text(aes(label = count), hjust = -0.1)+
    labs(x = "", y = "Percentage of cells", title = "") +
    coord_flip() +
    theme_classic() + 
    theme(
      axis.text.x = element_text(size = 14, family="ArialMT", color = "black"),
      axis.text.y = element_text(size=14, family="ArialMT", color = "black"),
      axis.title.x = element_text(size=14, family="ArialMT", color = "black"),
      axis.title.y = element_text(size=14, family="ArialMT", color = "black"),     
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      legend.position = "bottom") +  # show legend on right
    theme(axis.line = element_line(color = "black"))
  ggsave( "Cell_composition.myeloid.by_study.pdf", width = 7, height  = 3.5)
}

# CD4+
# selected dataset
if(T){
  studies = meta$study %>% levels
  
  keep =  c("D. Fawkner-Corbett et al., Cell, 2021","C. S. Smillie et al., Cell, 2019","R. Elmentaite et al., Nature, 2021",
            "R. Elmentaite et al., Dev Cell, 2020","B. Huang et al., Cell, 2019","J. C. Martin et al., Cell, 2019",
            "B. S. Boland et al., Sci Immunol, 2020","N. Jaeger et al., Nat Commun, 2021")
  tmp = meta %>% dplyr::filter(major_cluster == 'CD4T') %>% 
    dplyr::filter(study %in% keep) %>% 
    dplyr::filter(stage %in% c("pediatric","adult")) %>% 
    dplyr::group_by(study, minor_cluster) %>% summarise(count = n())
  
  tmp$study = factor(tmp$study, levels = studies[tmp$study %in% studies] %>% rev)
  ggplot(tmp, aes(x = study, y = count, fill = minor_cluster)) + 
    geom_bar(position="fill", stat="identity") + 
    scale_fill_manual(values = minor_cluster.color[tmp$minor_cluster %>% droplevels() %>% levels]) + 
    #geom_text(aes(label = count), hjust = -0.1)+
    labs(x = "", y = "Percentage of cells", title = "") +
    coord_flip() +
    theme_classic() + 
    theme(
      axis.text.x = element_text(size = 14, family="ArialMT", color = "black"),
      axis.text.y = element_text(size=14, family="ArialMT", color = "black"),
      axis.title.x = element_text(size=14, family="ArialMT", color = "black"),
      axis.title.y = element_text(size=14, family="ArialMT", color = "black"),     
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      legend.position = "bottom") +  # show legend on right
    theme(axis.line = element_line(color = "black"))
  ggsave( "Cell_composition.CD4T.by_study.pdf", width = 7, height  = 3.8)
}

# CD8+
# selected dataset
if(T){
  studies = meta$study %>% levels
  
  keep =  c("D. Fawkner-Corbett et al., Cell, 2021","C. S. Smillie et al., Cell, 2019","R. Elmentaite et al., Nature, 2021",
            "R. Elmentaite et al., Dev Cell, 2020","B. Huang et al., Cell, 2019","J. C. Martin et al., Cell, 2019",
            "B. S. Boland et al., Sci Immunol, 2020","D. Corridoni et al., Nat Med, 2020","N. Jaeger et al., Nat Commun, 2021")
  tmp = meta %>% dplyr::filter(major_cluster == 'CD8T') %>% 
    dplyr::filter(study %in% keep) %>% 
    dplyr::filter(stage %in% c("pediatric","adult")) %>% 
    dplyr::group_by(study, minor_cluster) %>% summarise(count = n())
  
  tmp$study = factor(tmp$study, levels = studies[tmp$study %in% studies] %>% rev)
  ggplot(tmp, aes(x = study, y = count, fill = minor_cluster)) + 
    geom_bar(position="fill", stat="identity") + 
    scale_fill_manual(values = minor_cluster.color[tmp$minor_cluster %>% droplevels() %>% levels]) + 
    #geom_text(aes(label = count), hjust = -0.1)+
    labs(x = "", y = "Percentage of cells", title = "") +
    coord_flip() +
    theme_classic() + 
    theme(
      axis.text.x = element_text(size = 14, family="ArialMT", color = "black"),
      axis.text.y = element_text(size=14, family="ArialMT", color = "black"),
      axis.title.x = element_text(size=14, family="ArialMT", color = "black"),
      axis.title.y = element_text(size=14, family="ArialMT", color = "black"),     
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      legend.position = "bottom") +  # show legend on right
    theme(axis.line = element_line(color = "black"))
  ggsave( "Cell_composition.CD8T.by_study.pdf", width = 7, height  = 4)
}

# ILC
# selected dataset
if(T){
  studies = meta$study %>% levels
  
  keep =  c("D. Fawkner-Corbett et al., Cell, 2021","C. S. Smillie et al., Cell, 2019","R. Elmentaite et al., Nature, 2021",
            "R. Elmentaite et al., Dev Cell, 2020","B. Huang et al., Cell, 2019","J. C. Martin et al., Cell, 2019",
            "B. S. Boland et al., Sci Immunol, 2020")
  tmp = meta %>% dplyr::filter(major_cluster == 'ILC') %>% 
    dplyr::filter(study %in% keep) %>% 
    dplyr::filter(stage %in% c("pediatric","adult")) %>% 
    dplyr::group_by(study, minor_cluster) %>% summarise(count = n())
  
  tmp$study = factor(tmp$study, levels = studies[tmp$study %in% studies] %>% rev)
  ggplot(tmp, aes(x = study, y = count, fill = minor_cluster)) + 
    geom_bar(position="fill", stat="identity") + 
    scale_fill_manual(values = minor_cluster.color[tmp$minor_cluster %>% droplevels() %>% levels]) + 
    #geom_text(aes(label = count), hjust = -0.1)+
    labs(x = "", y = "Percentage of cells", title = "") +
    coord_flip() +
    theme_classic() + 
    theme(
      axis.text.x = element_text(size = 14, family="ArialMT", color = "black"),
      axis.text.y = element_text(size=14, family="ArialMT", color = "black"),
      axis.title.x = element_text(size=14, family="ArialMT", color = "black"),
      axis.title.y = element_text(size=14, family="ArialMT", color = "black"),     
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      legend.position = "bottom") +  # show legend on right
    theme(axis.line = element_line(color = "black"))
  ggsave( "Cell_composition.ILC.by_study.pdf", width = 9, height  = 3)
}

# B/Plasma
# selected dataset
if(T){
  studies = meta$study %>% levels
  
  keep =  c("D. Fawkner-Corbett et al., Cell, 2021","C. S. Smillie et al., Cell, 2019","R. Elmentaite et al., Nature, 2021",
            "R. Elmentaite et al., Dev Cell, 2020","B. Huang et al., Cell, 2019","J. C. Martin et al., Cell, 2019",
            "B. S. Boland et al., Sci Immunol, 2020")
  tmp = meta %>% dplyr::filter(major_cluster == 'B_Plasma') %>% 
    dplyr::filter(study %in% keep) %>% 
    dplyr::filter(stage %in% c("pediatric","adult")) %>% 
    dplyr::group_by(study, minor_cluster) %>% summarise(count = n())
  
  tmp$study = factor(tmp$study, levels = studies[tmp$study %in% studies] %>% rev)
  ggplot(tmp, aes(x = study, y = count, fill = minor_cluster)) + 
    geom_bar(position="fill", stat="identity") + 
    scale_fill_manual(values = minor_cluster.color[tmp$minor_cluster %>% droplevels() %>% levels]) + 
    #geom_text(aes(label = count), hjust = -0.1)+
    labs(x = "", y = "Percentage of cells", title = "") +
    coord_flip() +
    theme_classic() + 
    theme(
      axis.text.x = element_text(size = 14, family="ArialMT", color = "black"),
      axis.text.y = element_text(size=14, family="ArialMT", color = "black"),
      axis.title.x = element_text(size=14, family="ArialMT", color = "black"),
      axis.title.y = element_text(size=14, family="ArialMT", color = "black"),     
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      legend.position = "bottom") +  # show legend on right
    theme(axis.line = element_line(color = "black"))
  ggsave( "Cell_composition.B_Plasma.by_study.pdf", width = 7, height  = 4)
}

# Epithelial
# selected dataset
if(T){
  studies = meta$study %>% levels
  
  keep =  c("K. Parikh et al., Nature, 2019","D. Fawkner-Corbett et al., Cell, 2021","C. S. Smillie et al., Cell, 2019",
            "R. Elmentaite et al., Nature, 2021","R. Elmentaite et al., Dev Cell, 2020","B. Huang et al., Cell, 2019")
  tmp = meta %>% dplyr::filter(major_cluster == 'Epithelial') %>% 
    dplyr::filter(study %in% keep) %>% 
    dplyr::filter(stage %in% c("pediatric","adult")) %>% 
    dplyr::group_by(study, minor_cluster) %>% summarise(count = n())
  
  tmp$study = factor(tmp$study, levels = studies[tmp$study %in% studies] %>% rev)
  ggplot(tmp, aes(x = study, y = count, fill = minor_cluster)) + 
    geom_bar(position="fill", stat="identity") + 
    scale_fill_manual(values = minor_cluster.color[tmp$minor_cluster %>% droplevels() %>% levels]) + 
    #geom_text(aes(label = count), hjust = -0.1)+
    labs(x = "", y = "Percentage of cells", title = "") +
    coord_flip() +
    theme_classic() + 
    theme(
      axis.text.x = element_text(size = 14, family="ArialMT", color = "black"),
      axis.text.y = element_text(size=14, family="ArialMT", color = "black"),
      axis.title.x = element_text(size=14, family="ArialMT", color = "black"),
      axis.title.y = element_text(size=14, family="ArialMT", color = "black"),     
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      legend.position = "bottom") +  # show legend on right
    theme(axis.line = element_line(color = "black"))
  ggsave( "Cell_composition.Epithelial.by_study.pdf", width = 10, height  = 4)
}

# Mesenchymal
# selected dataset
if(T){
  studies = meta$study %>% levels
  
  keep =  c("Kinchen et al., Cell, 2018","D. Fawkner-Corbett et al., Cell, 2021","C. S. Smillie et al., Cell, 2019",
            "R. Elmentaite et al., Nature, 2021","R. Elmentaite et al., Dev Cell, 2020","B. Huang et al., Cell, 2019",
            "J. C. Martin et al., Cell, 2019","M. Friedrich et al., Nat Med, 2021")
  tmp = meta %>% dplyr::filter(major_cluster == 'Mesenchymal') %>% 
    dplyr::filter(study %in% keep) %>% 
    dplyr::filter(stage %in% c("pediatric","adult")) %>% 
    dplyr::group_by(study, minor_cluster) %>% summarise(count = n())
  
  tmp$study = factor(tmp$study, levels = studies[tmp$study %in% studies] %>% rev)
  ggplot(tmp, aes(x = study, y = count, fill = minor_cluster)) + 
    geom_bar(position="fill", stat="identity") + 
    scale_fill_manual(values = minor_cluster.color[tmp$minor_cluster %>% droplevels() %>% levels]) + 
    #geom_text(aes(label = count), hjust = -0.1)+
    labs(x = "", y = "Percentage of cells", title = "") +
    coord_flip() +
    theme_classic() + 
    theme(
      axis.text.x = element_text(size = 14, family="ArialMT", color = "black"),
      axis.text.y = element_text(size=14, family="ArialMT", color = "black"),
      axis.title.x = element_text(size=14, family="ArialMT", color = "black"),
      axis.title.y = element_text(size=14, family="ArialMT", color = "black"),     
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_blank(),
      legend.position = "bottom") +  # show legend on right
    theme(axis.line = element_line(color = "black"))
  ggsave( "Cell_composition.Mesenchymal.by_study.pdf", width = 8, height  = 4)
}

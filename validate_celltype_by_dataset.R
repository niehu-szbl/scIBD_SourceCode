library(dplyr)
library(ggplot2)

# prepare data
if(T){
  # update meta
  meta = readRDS("all.clean.meta.rds")

  # load color
  color.dataset = read.table("./dataset_color.txt", header = T, sep = "\t", stringsAsFactors = F, comment.char = "/")
  tmp = color.dataset$color
  names(tmp) = color.dataset$dataset
  color.dataset = tmp
  rm(tmp)
}

# Tc17
# selected dataset
if(T){
  
  studies = meta.data$study %>% levels
  
  drop =  c("K. Parikh et al., Nature, 2019","Kinchen et al., Cell, 2018", "M. Friedrich et al., Nat Med, 2021")
  tmp =meta.data %>% filter(minor_cluster == 'CD8+ Tc17') %>% 
    group_by(study, .drop = F) %>% summarise(count = n()) %>% filter(!(study %in% drop))
  
  tmp$study = factor(tmp$study, levels = studies[tmp$study %in% studies])
  ggplot(tmp, aes(x = study, y = count, fill = study)) + 
    geom_bar(stat = 'identity') + 
    scale_fill_manual(values = color.dataset[tmp$study]) + 
    geom_text(aes(label = count), hjust = -0.1)+
    labs(x = "", y = "Number of cells", title = "Number of CD8+ Tc17 in individual studies") +
    ylim(0, 1600)+
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
      legend.position = "none") +  # show legend on right
    theme(axis.line = element_line(color = "black"))
  ggsave( "Number_of_Tc17_by_study.pdf", width = 8, height  = 4)
}

# LAMP3+ DC
# selected dataset
if(T){
  
  studies = meta.data$study %>% levels
  
  keep =  c("D. Fawkner-Corbett et al., Cell, 2021","C. S. Smillie et al., Cell, 2019","R. Elmentaite et al., Nature, 2021",
            "R. Elmentaite et al., Dev Cell, 2020","B. Huang et al., Cell, 2019","J. C. Martin et al., Cell, 2019","B. S. Boland et al., Sci Immunol, 2020")
  tmp =meta.data %>% filter(minor_cluster == 'LAMP3+ DC') %>% 
    group_by(study, .drop = F) %>% summarise(count = n()) %>% filter(study %in% keep)
  
  tmp$study = factor(tmp$study, levels = studies[tmp$study %in% studies])
  ggplot(tmp, aes(x = study, y = count, fill = study)) + 
    geom_bar(stat = 'identity') + 
    scale_fill_manual(values = color.dataset[tmp$study]) + 
    geom_text(aes(label = count), hjust = -0.1)+
    labs(x = "", y = "Number of cells", title = "Number of LAMP3+ DC cells in individual studies") +
    ylim(0, 300)+
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
      legend.position = "none") +  # show legend on right
    theme(axis.line = element_line(color = "black"))
  ggsave( "Number_of_lamp3_dc_by_study.pdf", width = 8, height  = 4)
}

# DUOX2+ epithelial
# selected dataset
if(T){
  
  studies = meta.data$study %>% levels
  
  keep =  c("K. Parikh et al., Nature, 2019",
            "D. Fawkner-Corbett et al., Cell, 2021","C. S. Smillie et al., Cell, 2019",
            "R. Elmentaite et al., Nature, 2021","R. Elmentaite et al., Dev Cell, 2020","B. Huang et al., Cell, 2019",
            "M. Friedrich et al., Nat Med, 2021")
  tmp = meta.data %>% filter(minor_cluster == 'DUOX2+ epithelial') %>% 
    group_by(study, .drop = F) %>% summarise(count = n()) %>% filter(study %in% keep)
  
  tmp$study = factor(tmp$study, levels = studies[tmp$study %in% studies])
  ggplot(tmp, aes(x = study, y = count, fill = study)) + 
    geom_bar(stat = 'identity') + 
    scale_fill_manual(values = color.dataset[tmp$study]) + 
    geom_text(aes(label = count), hjust = -0.1)+
    labs(x = "", y = "Number of cells", title = "Number of DUOX2+ epithelial cells in individual studies") +
    ylim(0, 3500)+
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
      legend.position = "none") +  # show legend on right
    theme(axis.line = element_line(color = "black"))
  ggsave( "Number_of_duox2_epi_by_study.pdf", width = 8, height  = 4)
}



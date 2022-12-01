library(flowCore)
library(flowWorkspace)
library(ggcyto)
library(flowAI)
library(gridExtra)
library(tidyverse)
library(flowStats)
library(CytoML)
library(Rtsne)
library(FlowSOM)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(lme4)
library(multcomp)
library(dplyr)

######### STACKED BAR CHART
######### STACKED BAR CHART
######### STACKED BAR CHART

####Prepare data

gs_get_pop_paths(gs)[c(11,15,18,22,23,27,30,31,32,35)]
#counts_table <-gs_pop_get_count_fast(gs[1])
counts_table <- gs_pop_get_count_fast(gs, format = "long", subpopulations = gs_get_pop_paths(gs)[c(11,15,18,22,23,27,30,31,32,35)])
counts_table
counts_table <- counts_table %>% pivot_wider(id_cols = name, 
                                             names_from = Population, 
                                             values_from = c("Count", "ParentCount"))
write.csv(counts_table, "C:/Users/edmondsonef/Desktop/counts_table.csv")


counts_table <- read.csv("C:/Users/edmondsonef/Desktop/counts_table.csv")
counts_table_t <- counts_table[-1] %>% t() %>% as.data.frame() %>% setNames(counts_table[,1])
props_table <- t(t(counts_table_t) / colSums(counts_table_t[])) * 100
props_table_t <- props_table[] %>% t() %>% as.data.frame() %>% setNames(row.names(props_table))
counts <- as.data.frame.matrix(counts_table_t)
props <- as.data.frame.matrix(props_table)
write.csv(props, "C:/Users/edmondsonef/Desktop/props.csv")


props <- read.csv("C:/Users/edmondsonef/Desktop/props.csv", header = T, stringsAsFactors = F)
props <- data.frame(props[,-1], row.names = props[,1])
ggdf <- reshape2::melt(data.frame(cluster = rownames(props), props), 
                       id.vars = "cluster", value.name = "proportion", variable.name = "sample_id")
write.csv(ggdf, "C:/Users/edmondsonef/Desktop/ggdf.csv")
#ggdf <- read.csv("C:/Users/edmondsonef/Desktop/ggdf.csv", header = T, stringsAsFactors = F)




setwd("C:/Users/edmondsonef/Desktop/R-plots/")
color_clusters <- c("#7BAFDE", "#7570B3", "#882E72", "#B17BA6", "#DC050C",
                    "#33A02C", "#B2DF8A", "#A6761D","#FB8072","#d4b7b7")
#"#BEAED4", "#aeae5c", "#1e90ff", "#00bfff", 
#"#999999", "#FF7F00","#56ff0d", "#ffff00"
#"#aa8282", "#8600bf", "#ba5ce3", "#808000")







ggdf <- read_excel("C:/Users/edmondsonef/Desktop/Humanized/agg (for stacked barplots).xlsx", sheet = "Final.123")
ggdf <- dplyr::filter(ggdf, ggdf$`Remove` != "1")


# ggdf <- ggdf %>% group_by(clusters, Group) %>% summarise(prop = mean(proportion), SD= sd(proportion), GroupX = first(GroupX), n = n()) 
ggdf <- ggdf %>% group_by(clusters, Group, Tissue) %>% summarise(prop = mean(proportion), SD= sd(proportion), n = n()) 

plot <- ggplot(ggdf, aes(x = Group, y = prop, fill = clusters))+
  geom_bar(stat = "identity")+
  labs(title="Humanized NSG, MDA-MB-231T") +
  theme_bw() +
  theme(axis.title.x=element_blank(), text = element_text(size = 12))+
  facet_wrap(~ Tissue, scales = "free_x") +
  geom_text(data = subset(ggdf, prop > 0.1), 
            aes(label = paste(round(prop,0),"%")), 
            size = 3, color = "white", #fontface = "bold", family = "Fira Sans",
            fill = "white", label.size = 0, 
            position = position_stack(vjust =  0.5))+
  #geom_errorbar(aes(ymin = prop-SD, ymax = prop+SD), width = 0.3, position = "identity")+
  scale_fill_manual(values = color_clusters) 


plot


tiff("plots_123.tiff", units="in", width=10, height=7, res=300)
plot
dev.off()









ggdf <- read_excel("C:/Users/edmondsonef/Desktop/Humanized/agg (for stacked barplots).xlsx", sheet = "counts.1")
ggdf <- ggdf %>% group_by(weeks, Groups, Tissue) %>% summarise(CD45_count = mean(`CD45 Cells`), Percent_human = mean(`% Human`), n = n()) 

ggdf_blood <- dplyr::filter(ggdf, Tissue == "Blood")
ggdf_other <- dplyr::filter(ggdf, Tissue == c("Marrow", "Spleen"))


##Blood: Percentage Human Over Time
ggplot(ggdf_blood, aes(weeks, Percent_human, color = Group_sim)) +
  geom_count(aes(size = Percent_human)) +
  theme_bw() +
  geom_line()



ggplot(ggdf, aes(Group, CD45_count, color = Tissue)) +
  geom_point(aes(size = Percent_human)) +
  theme_bw()


ggplot(ggdf_other, aes(x=Group, y=`% Human`, fill=Tissue)) +
  geom_jitter(aes(color = `Tissue`), width = 0.2, height = 0.001, size = 3) +
  theme_bw()







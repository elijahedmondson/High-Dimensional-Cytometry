
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
#manual
#Load data
#myfiles <- list.files(path="C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15701 02Feb2022 Simone/", pattern = ".FCS", ignore.case = TRUE)
#fs <- read.flowSet(myfiles, path="C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15701 02Feb2022 Simone/")#, truncate_max_range = FALSE)

# data_dir <- "C:/Users/edmondsonef/Desktop/Humanized/Flow/Aggregate/FCS/" 
# ws <- open_flowjo_xml("C:/Users/edmondsonef/Desktop/Humanized/Flow/Aggregate/FCS/19-331-122 Full Dataset WSP 220414.wsp")
# #ws <- open_flowjo_xml(paste0(data_dir,"15761 12Apr2022 LASP.wsp"))
# ws
# 
# head(fj_ws_get_samples(ws, group_id = c(1)))
# tail(fj_ws_get_sample_groups(ws))
# fj_ws_get_keywords(ws, 28)[1:8]
# 
# myfiles <- list.files(path=data_dir, pattern = "Samples", ignore.case = TRUE)
# 
# path_comp <- "C:/Users/edmondsonef/Desktop/0-Autospill/6-10Mar2022/autospill_compensation.csv"
# compensation_matrix <- read.csv(path_comp, check.names = FALSE, row.names = 1)
# 






#####
##### Creating manual plots




data_dir <- "C:/Users/edmondsonef/Desktop/Humanized/Flow/Aggregate/FCS/marrow/19-331-122/" 
#myfiles <- list.files(path=data_dir, pattern = "Samples", ignore.case = TRUE)
#fs <- read.flowSet(myfiles, path=data_dir, truncate_max_range = FALSE)



ws <- open_flowjo_xml("C:/Users/edmondsonef/Desktop/Humanized/Flow/Aggregate/FCS/19-331-122 Full Dataset WSP 220414.wsp")
gs <- flowjo_to_gatingset(ws, name = 1, path=data_dir)
plot(gs)

pData(gs) %>% head(15)
pData(gs)$animal <-  gsub(".*_.*_...(.*)_.*.fcs_.*","\\1",sampleNames(gs))
pData(gs) %>% head(15)

out <- pData(gs)
write.csv(out, "C:/Users/edmondsonef/Desktop/Human1.csv")
out <- read.csv("C:/Users/edmondsonef/Desktop/Human1.csv", header = T, stringsAsFactors = F)
pData(gs)$strain <- out$strain


colnames(gs)
colnames(gs)[colnames(gs)=="Comp-BB515-A"] <- "CD8"
colnames(gs)[colnames(gs)=="Comp-BB700-P-A"] <- "CD4"
colnames(gs)[colnames(gs)=="Comp-APC-A"] <- "CD11b"
colnames(gs)[colnames(gs)=="Comp-APC-Cy7-A"] <- "CD19"
colnames(gs)[colnames(gs)=="Comp-BV421-A"] <- "CD3"
colnames(gs)[colnames(gs)=="Comp-BV786-A"] <- "CD33"
colnames(gs)[colnames(gs)=="Comp-BUV395-A"] <- "mCD45"
colnames(gs)[colnames(gs)=="Comp-BUV805-A"] <- "huCD45"
colnames(gs)[colnames(gs)=="Comp-PE-A"] <- "CD56"
colnames(gs)[colnames(gs)=="Comp-PE-CF594-A"] <- "CD66b"
colnames(gs)[colnames(gs)=="Comp-PE-Cy7-A"] <- "CD25"
colnames(gs)

ggcyto(gs[1],aes(x="FSC-A",y="FSC-H"),subset="root")+geom_gate("sing")+
  geom_hex(bins = 200)+
  geom_stats(adjust = 0.8)+ggcyto_par_set(limits = "instrument") +
  facet_wrap(~animal,ncol = 10)

colnames(gs)
plot(gs)
gs_get_pop_paths(gs, path = 2)
nodelist <- gs_get_pop_paths(gs, path = "auto")[c(11,15,18,22,23,27,30,31,32,35)]
nodelist
node <- nodelist[3]
g <- gs_pop_get_gate(gs, node)
g
gs_pop_get_stats(gs)[1:10,]


autoplot(gs, nodelist[1], digits = 2)+
  facet_wrap(~animal) + geom_hex(bins = 300)

autoplot(gs, nodelist[2], digits = 2)+
  facet_wrap(~animal)
autoplot(gs, nodelist[3], digits = 2)+
  facet_wrap(~animal)
autoplot(gs, nodelist[4], digits = 2)+
  facet_wrap(~animal)
autoplot(gs, nodelist[5], digits = 2)+
  facet_wrap(~animal)
autoplot(gs, nodelist[6], digits = 2)+
  facet_wrap(~animal)

ps <- gs_pop_get_count_with_meta(gs)
ps <- ps %>% mutate(percent_of_parent=Count/ParentCount)
head(ps)




######### % Human
######### % Human
######### % Human
human_counts <- gs_pop_get_count_fast(gs, format = "long", subpopulations = gs_get_pop_paths(gs)[1])
human_counts
human_counts <- human_counts %>% pivot_wider(id_cols = name, 
                                             names_from = Population, 
                                             values_from = c("Count", "ParentCount"))
write.csv(human_counts, "C:/Users/edmondsonef/Desktop/Humanized.count.csv")


######### STACKED BAR CHART
######### STACKED BAR CHART
######### STACKED BAR CHART
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




data <- read_excel("C:/Users/edmondsonef/Desktop/agg.xlsx", sheet = "Final.blood")
ggdf <- data


color_clusters <- c("#7BAFDE", "#7570B3", "#882E72", "#B17BA6", "#FF7F00", 
                    "#DC050C", "#FB8072", "#FDB462", #"#E7298A", "#E78AC3", 
                    "#33A02C", "#B2DF8A", "#55A1B1", "#8DD3C7", "#A6761D", 
                    "#E6AB02", "#1965B0", "#BEAED4", "#666666", "#999999", 
                    "#aa8282", "#d4b7b7", "#8600bf", "#ba5ce3", "#808000", 
                    "#aeae5c", "#1e90ff", "#00bfff", "#56ff0d", "#ffff00")

plot <- ggplot(ggdf, aes(x = ID, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ GroupX, scales = "free_x") +
  theme_bw() +
  theme(axis.title.x=element_blank(), text = element_text(size = 10))+
  labs(title="Blood") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values = color_clusters) 
plot
setwd("C:/Users/edmondsonef/Desktop/R-plots/")
tiff("blood_plots.tiff", units="in", width=10, height=4, res=600)
plot
dev.off()
######### STACKED BAR CHART
######### STACKED BAR CHART
######### STACKED BAR CHART

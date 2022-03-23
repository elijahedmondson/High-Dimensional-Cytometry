
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(flowAI)
library(gridExtra)
library(tidyverse)
library(flowStats)
library(flowWorkspace)
library(CytoML)
library(Rtsne)
library(FlowSOM)

myfiles <- list.files(path="C:/Users/edmondsonef/Desktop/Humanized/Flow/28Feb2022/", pattern = ".fcs", ignore.case = TRUE)
wsp <- open_flowjo_xml("C:/Users/edmondsonef/Desktop/Humanized/Flow/28Feb2022/15716 28Feb2022 Simone.wsp")
#fcs_file <- "C:/Users/edmondsonef/Desktop/Humanized/Flow/02Feb2022/Samples_Tube_018 Animal 120 BMC_018.fcs"
#fs <- read.flowSet(myfiles, path="C:/Users/edmondsonef/Desktop/Humanized/Flow/02Feb2022", truncate_max_range = FALSE)
#fs_comp <-compensate(fs, spillover(fs[[1]])$SPILL)

#tail(fj_ws_get_sample_groups(wsp))
fj_ws_get_samples(wsp, group_id = 1)




#Removing stuff
#Removing stuff
#Removing stuff
gs <- flowjo_to_gatingset(wsp, name = 1, path ="C:/Users/edmondsonef/Desktop/Humanized/Flow/28Feb2022/")
plot(gs)
#autoplot(gs[[1]])

gs_get_pop_paths(gs)
recompute(gs)


sampStats <- gs_pop_get_stats(gs)
write.csv(sampStats, "C:/Users/edmondsonef/Desktop/sampStats.csv")
sampStats <- read.csv("C:/Users/edmondsonef/Desktop/sampStats.csv")

popMFI <- gs_pop_get_stats(gs, type = pop.MFI)
write.csv(popMFI, "C:/Users/edmondsonef/Desktop/popMFI.csv")
popMFI <- read.csv("C:/Users/edmondsonef/Desktop/popMFI.csv")

FULL1 <- dplyr::right_join(sampStats, popMFI, by = "X")
write.csv(FULL1, "C:/Users/edmondsonef/Desktop/28Feb2022.csv")

### Get Stats on Manual Gates or Nodes
### Get Stats on Manual Gates or Nodes
### Get Stats on Manual Gates or Nodes
### Get Stats on Manual Gates or Nodes
files = myfiles

gs <- gs[1]
for(file in files){}

nodes <- c("/scatter/sing/", "/scatter/sing/hCD45+",
           "Q6: CD3+ , CD4 [PCP55]+",
           "Q10: CD3+ , CD8 [FITC]+",
           "Q13: CD3- , CD19 [AFire750]+",
           "Q17: CD3- , CD56+",
           "Q18: CD3+ , CD56+",
           "Q29: CD66b [PEDazz]- , CD11b [AF647]+",
           "Q30: CD66b [PEDazz]+ , CD11b [AF647]+",
           "Q31: CD66b [PEDazz]+ , CD11b [AF647]-",
           "Q33: CD33- , CD11b [AF647]+",
           "Q38: CD25+ , CD3+",
           "Q35: CD33+ , CD11b [AF647]-",
           "Q39: CD25+ , CD3-")
gs_pop_get_stats(gs, nodes, "percent")
nodeCount <- gs_pop_get_stats(gs, nodes, "count")
nodeCount
#write.csv(nodeCount, "C:/Users/edmondsonef/Desktop/nodeCount.csv")
### Pass a build-in function
nodePopMFI <- gs_pop_get_stats(gs, nodes, type = pop.MFI)
nodePopMFI
# compute the stats based on the raw data scale
nodePopMFI.inv <- gs_pop_get_stats(gs, nodes, type = pop.MFI, inverse.transform = TRUE)
nodePopMFI.inv
# supply user-defined stats fun
pop.quantiles <- function(fr){
  chnls <- colnames(fr)
  res <- matrixStats::colQuantiles(exprs(fr), probs = 0.75)
  names(res) <- chnls
  res
}
quants <- gs_pop_get_stats(gs, nodes, type = pop.quantiles)


#nodeCount <- as.data.frame(nodeCount)
#nodePopMFI <- as.data.frame(nodePopMFI)
#write.csv(nodeCount, "C:/Users/edmondsonef/Desktop/nodeCount.csv")
#nodeCount <- read.csv("C:/Users/edmondsonef/Desktop/nodeCount.csv")
#write.csv(nodePopMFI, "C:/Users/edmondsonef/Desktop/nodePopMFI.csv")
#nodePopMFI <- read.csv("C:/Users/edmondsonef/Desktop/nodePopMFI.csv")

FULL <- dplyr::right_join(nodeCount, nodePopMFI, by = "X")
write.csv(FULL, "C:/Users/edmondsonef/Desktop/FULL.csv")
### Get Stats on Manual Gates or Nodes
### Get Stats on Manual Gates or Nodes
### Get Stats on Manual Gates or Nodes
### Get Stats on Manual Gates or Nodes
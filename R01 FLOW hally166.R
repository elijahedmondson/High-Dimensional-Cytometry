
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(flowAI)
library(gridExtra)
library(tidyverse)
library(flowStats)

#manual
#Load data
#myfiles <- list.files(path="C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15701 02Feb2022 Simone/", pattern = ".FCS", ignore.case = TRUE)
#fs <- read.flowSet(myfiles, path="C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15701 02Feb2022 Simone/")#, truncate_max_range = FALSE)
myfiles <- list.files(path="C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/", pattern = ".FCS", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path="C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022", truncate_max_range = FALSE)
fs_comp <-compensate(fs, spillover(fs[[1]])$SPILL)


###https://jchellmuth.com/posts/FACS-with-R/
###https://jchellmuth.com/posts/FACS-with-R/
###https://jchellmuth.com/posts/FACS-with-R/
###https://jchellmuth.com/posts/FACS-with-R/
###https://jchellmuth.com/posts/FACS-with-R/



pData(fs) %>% head(3)
pData(fs)$well <- gsub(".*_.*_(.*)_.*.fcs","\\1",sampleNames(fs)) # extract well from name and add new 'well' column
pData(fs) %>% head(3)
colnames(fs)

colnames(fs)[colnames(fs)=="BB515-A"] <- "CD8"
colnames(fs)[colnames(fs)=="BB700-P-A"] <- "CD4"
colnames(fs)[colnames(fs)=="APC-A"] <- "CD11b"
colnames(fs)[colnames(fs)=="APC-Cy7-A"] <- "CD19"
colnames(fs)[colnames(fs)=="BV421-A"] <- "CD3"
colnames(fs)[colnames(fs)=="BV786-A"] <- "CD33"
colnames(fs)[colnames(fs)=="BUV395-A"] <- "mCD45"
colnames(fs)[colnames(fs)=="BUV805-A"] <- "huCD45"
colnames(fs)[colnames(fs)=="PE-A"] <- "CD56"
colnames(fs)[colnames(fs)=="PE-CF594-A"] <- "CD66b"
colnames(fs)[colnames(fs)=="PE-Cy7-A"] <- "CD25"
#fs <- fsApply(fs, function(x) {transform(x, estimateLogicle(x, c("CD8", "CD4", "CD11b","CD19", "CD3", "CD33",
#                                                                       "mCD45", "huCD45", "CD56", "CD66b")))})
#fs <- fs_trans
gs <- GatingSet(fs)
g.singlets <- polygonGate(filterId = "Singlets","FSC-A"=c(2e4,26e4,26e4,2e4),"FSC-H"=c(0e4,16e4,24e4,6e4)) # define gate
gs_pop_add(gs,g.singlets) # add gate to GatingSet

ggcyto(fs[[1]], aes(x="FSC-A",y="FSC-H"),subset="root") + geom_hex(bins = 500) + 
  labs(title = "huCD45 vs FSC", y = "FSC-H", x = "FSC-A") + ggcyto_par_set(limits = "instrument") +
  geom_gate(g.singlets) +
  theme_bw() + theme(legend.position = "none", aspect.ratio = 1) + facet_wrap(~name, ncol = 4) 
gs_pop_add(gs,g.singlets) # add gate to GatingSet


recompute(gs)

ggcyto(gs[[1]],aes(x="FSC-A",y="SSC-A"),subset="root")+geom_hex(bins = 200)+ggcyto_par_set(limits = "instrument")
ggcyto(gs[[1]],aes(x="FSC-A",y="SSC-A"),subset="Singlets")+geom_hex(bins = 200)+ggcyto_par_set(limits = "instrument")

ggcyto(gs,aes(x="FSC-A",y="FSC-H"),subset="root")+geom_hex(bins = 100)+geom_gate("Singlets")+
  geom_stats(adjust = 0.8)+ggcyto_par_set(limits = "instrument")+
  facet_wrap(~well,ncol = 10)


g.live <- polygonGate(filterId = "Live","FSC-A"=c(8e4,28e4,28e4,8e4),"SSC-A"=c(3e4,4e4,28e4,28e4)) # define gate
ggcyto(gs[[1]],aes(x="FSC-A",y="SSC-A"),subset="Singlets")+geom_hex(bins = 200)+geom_gate(g.live)+ggcyto_par_set(limits = "instrument") # check gate
gs_pop_add(gs,g.live,parent="Singlets") # add gate to GatingSet
recompute(gs)


ggcyto(gs,aes(x="FSC-A",y="SSC-A"),subset="Singlets")+geom_hex(bins = 100)+geom_gate("Live")+
  geom_stats(adjust = 0.8)+ggcyto_par_set(limits = "instrument")+
  facet_wrap(~well,ncol = 10)

g.huCD45 <- rectangleGate(filterId="huCD45 positive","huCD45"=c(1000, Inf)) # set gate
ggcyto(gs[[1]],aes(x=huCD45),subset="Live")+geom_density(fill="forestgreen")+geom_gate(g.huCD45)+ggcyto_par_set(limits = "instrument")+scale_x_flowJo_biexp() # check gate
gs_pop_add(gs,g.huCD45,parent="Live")
recompute(gs)


ggcyto(gs,aes(x=huCD45),subset="Live",)+geom_density(fill="forestgreen")+geom_gate("huCD45 positive")+
  geom_stats()+
  ggcyto_par_set(limits = "instrument")+scale_x_flowJo_biexp()+
  facet_wrap(~well,ncol = 10)

gs_pop_get_count_fast(gs) %>% head
ps <- data.frame(gs_pop_get_count_fast(gs))

ps$percent_of_parent <- ps$Count/ps$ParentCount*100
psm <- merge(ps,pData(fs),by="name")

###https://jchellmuth.com/posts/FACS-with-R/
###https://jchellmuth.com/posts/FACS-with-R/
###https://jchellmuth.com/posts/FACS-with-R/
###https://jchellmuth.com/posts/FACS-with-R/
###https://jchellmuth.com/posts/FACS-with-R/
## Not run: 



# fr is a flowFrame
sg <- gate_singlet(fs, area = "FSC-A", height = "FSC-H")
sg
# plot the gate 
xyplot(`FSC-H` ~ `FSC-A`, fr, filter = sg)

## End(Not run)







fs_comp_clean <- flow_auto_qc(fs_comp)
trans <- estimateLogicle(fs_comp_clean[[1]], colnames(fs_comp_clean[,3:15]))
fs_comp_clean_trans <- transform(fs_comp_clean, trans)

#Visualise file
fs_comp_clean_trans[[1]]
autoplot(fs_comp_clean_trans[[1]])

#Basic gating
ggcyto(fs_comp_clean_trans[[1]], aes(x="FSC-A", y="SSC-A"))+geom_hex(bins=512)
ggcyto(fs_comp_clean_trans[[2]], aes(x="BUV805-A", y="BUV395-A"))+geom_hex(bins=512)
ggcyto(fs_comp_clean_trans[[2]], aes(x="BV421-A", y="APC-Cy7-A"))+geom_hex(bins=512)
ggcyto(fs_comp_clean_trans[[2]], aes(x="BUV805-A", y="SSC-A"))+geom_hex(bins=512)


names(fs_comp_clean_trans[[2]])
#exprs(fs_comp_clean_trans[[2]])
#each_col(fs_comp_clean_trans[[2]], median)


#create the empty gating set
gs<-GatingSet(fs_comp_clean_trans)


#Cells - FSC SSC
#CREATES RECTANGLE FOR GATING
rg1<-rectangleGate("BUV805-A"=c(2, Inf), filterId = "huCD45")
gs_pop_add(gs, rg1, parent="root")
recompute(gs)
gs_get_pop_paths(gs)
ggcyto(fs_comp_clean_trans[[1]], aes(x="BUV805-A", y="BUV395-A"))+geom_hex(bins=256)+geom_gate(gs_pop_get_gate(gs, "huCD45"))
gs_pop_get_stats(gs)

#Singlet gating
ggcyto(fs_comp_clean_trans[[1]], aes(x = "FSC-H", y = 'FSC-W'))+ geom_hex(bins = 256)
rg2 <- rectangleGate("FSC-H"=c(3.6, 4.2),"FSC-W"=c(50, 240))
gs_pop_add(gs, rg2, parent = "NoneDebris", name = "singlets")
gs_get_pop_paths(gs)
recompute(gs)
ggcyto(fs_comp_clean_trans, aes(x = "FSC-H", y = 'FSC-W'))+ geom_hex(bins = 256)+ geom_gate(gs_pop_get_gate(gs, "singlets"))

#exploring the gatingset
plot(gs)
gs_pop_get_stats(gs)
gs_pop_get_stats(gs, "NoneDebris", "percent")


#automatic
#Load data
myfiles <- list.files(path="C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZV_files", pattern = ".FCS", ignore.case = TRUE)
fs <- read.flowSet(myfiles[1:2], path="C:/Users/chall/Downloads/FlowRepository_FR-FCM-ZZZV_files", alter.names=TRUE)
fs_comp <-compensate(fs,spillover(fs[[1]])$SPILL)
#fix compensation matrix
matrix<-spillover(fs[[1]])$SPILL
matrix
colnames(matrix)<-c("X.FITC.A.", "X.Pacific.Blue.A.", "X.Alexa.680.A.", "X.APC.A.", "X.PE.Cy7.A.", "X.PE.Cy55.A.", "X.PE.Tx.RD.A.", "X.PE.Green.laser.A.")
fs_comp <-compensate(fs,matrix)
#continue
fs_comp_clean <- flow_auto_qc(fs_comp)
trans <- estimateLogicle(fs_comp_clean[[1]], colnames(fs_comp_clean[,c(4,6:12)]))
fs_comp_clean_trans <- transform(fs_comp_clean, trans)
autoplot(fs_comp_clean_trans[[1]])

#create the empty gating set
auto_gs<-GatingSet(fs_comp_clean_trans)

#cell gate
fs_data<- gs_pop_get_data(auto_gs)
noneDebris_gate<- fsApply(fs_data, function(fr) openCyto:::.flowClust.2d(fr, channels= c("FSC.A","SSC.A")))
gs_pop_add(auto_gs, noneDebris_gate, parent = "root", name="noneDebris_gate")
recompute(auto_gs)
autoplot(auto_gs, x="FSC.A", y="SSC.A", "noneDebris_gate", bins=256)

#Singlet gate
fs_data <- gs_pop_get_data(auto_gs, "noneDebris_gate") #get parent data
singlet_gate <- fsApply(fs_data, function(fr) openCyto:::.singletGate(fr, channels =c("FSC.A", "FSC.H")))
gs_pop_add(auto_gs, singlet_gate, parent = "noneDebris_gate", name = "singlets")
recompute(auto_gs)
autoplot(auto_gs, x = 'FSC.A', y = 'FSC.H', "singlets", bins = 256)

#Quad gate
fs_data <- gs_pop_get_data(auto_gs, "singlets") #get parent data
BGquad_gate <- fsApply(fs_data, function(fr) openCyto:::.quadGate.seq(fr, gFunc="mindensity", min=c(3,3), channels =c('X.FITC.A.', 'X.PE.Tx.RD.A.')))
gs_pop_add(auto_gs, BGquad_gate, parent = "singlets", names = c("1", "2", "3", "4"))
recompute(auto_gs)
gs_get_pop_paths(auto_gs[[1]])
plot(auto_gs)
autoplot(auto_gs, x = 'X.FITC.A.', y = 'X.PE.Tx.RD.A.', gs_get_pop_paths(auto_gs)[4:7], bins = 256)

#fix plot
p<-ggcyto(auto_gs[1:2],aes(x = 'X.FITC.A.', y = 'X.PE.Tx.RD.A.'), subset="singlets", arrange = FALSE)
p<- p + geom_hex(bins=256)
p<- p + geom_gate(gs_get_pop_paths(auto_gs)[4:7]) 
p<- p + geom_stats(gs_get_pop_paths(auto_gs)[4:7])
p<- p + theme(strip.text = element_text(size = 7))
myPars <- ggcyto_par_set(limits = list(y = c(3,5), x = c(3,5)))
p<- p  + myPars
p

#Removing stuff
gs_pop_remove(auto_gs, "singlets")

#statistics
gs_pop_get_stats(auto_gs)
gs_pop_get_stats(auto_gs, "noneDebris_gate", "percent")
gs_pop_get_stats(auto_gs, "noneDebris_gate", type = pop.MFI)

pop.quantiles <- function(fr){
  chnls <- colnames(fr)
  res <- matrixStats::colQuantiles(exprs(fr), probs = 0.75)
  names(res) <- chnls
  res
}
gs_pop_get_stats(auto_gs, gs_get_pop_paths(auto_gs), type = pop.quantiles)

pop.mean <- function(fr){
  chnls <- colnames(fr)
  res <- colMeans(exprs(fr))
  names(res) <- chnls
  res
}
gs_pop_get_stats(auto_gs, gs_get_pop_paths(auto_gs), type = pop.mean)


















markers_of_interest <- c("SSC-A", "CD8", "CD4", "CD3", "CD33",
                         "CD11b", "CD25", "mCD45", "huCD45", "CD56", 
                         "CD19", "CD66b")







PlotDimRed(fsom, colsToUse = fsom$map$colsUsed, colorBy = "metaclusters")


#https://rpubs.com/artur_matysik/flow
###https://jchellmuth.com/posts/FACS-with-R/

library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(flowViz)
library(flowStats)
library(scales)
library(dplyr)
library(stringr)
library(viridis)
library(flowCore)
library(flowWorkspace)
library(openCyto)
library(ggcyto)
library(flowAI)
library(gridExtra)
library(tidyverse)
library(flowStats)

#manual
#Load data
myfiles <- list.files(path="C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15701 02Feb2022 Simone/", pattern = ".FCS", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path="C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15701 02Feb2022 Simone/")#, truncate_max_range = FALSE)
myfiles <- list.files(path="C:/Users/edmondsonef/Desktop/samp/", pattern = ".FCS", ignore.case = TRUE)
fs <- read.flowSet(myfiles, path="C:/Users/edmondsonef/Desktop/samp/", truncate_max_range = FALSE)
fs_comp <-compensate(fs, spillover(fs[[1]])$SPILL)





pData(fs) %>% head(3)
pData(fs)$well <- gsub(".*_.*_(.*)_.*.fcs","\\1",sampleNames(fs)) # extract well from name and add new 'well' column
pData(fs) %>% head(3)
colnames(fs)

colnames(fs)[colnames(fs)=="BB515-A"] <- "CD8"
colnames(fs)[colnames(fs)=="BB700-P-A"] <- "CD4"
colnames(fs)[colnames(fs)=="APC-A"] <- "CD11b"
colnames(fs)[colnames(fs)=="APC-Cy7-A"] <- "CD19"
colnames(fs)[colnames(fs)=="BV421-A"] <- "CD3"
colnames(fs)[colnames(fs)=="BV786-A"] <- "CD33"
colnames(fs)[colnames(fs)=="BUV395-A"] <- "mCD45"
colnames(fs)[colnames(fs)=="BUV805-A"] <- "huCD45"
colnames(fs)[colnames(fs)=="PE-A"] <- "CD56"
colnames(fs)[colnames(fs)=="PE-CF594-A"] <- "CD66b"
colnames(fs)[colnames(fs)=="PE-Cy7-A"] <- "CD25"
#fs <- fsApply(fs, function(x) {transform(x, estimateLogicle(x, c("CD8", "CD4", "CD11b","CD19", "CD3", "CD33",
#                                                                       "mCD45", "huCD45", "CD56", "CD66b")))})
#fs <- fs_trans
gs <- GatingSet(fs)
g.singlets <- polygonGate(filterId = "Singlets","FSC-A"=c(2e4,26e4,26e4,2e4),"FSC-H"=c(0e4,16e4,24e4,6e4)) # define gate
gs_pop_add(gs,g.singlets) # add gate to GatingSet

ggcyto(fs[[1]], aes(x="FSC-A",y="FSC-H"),subset="root") + geom_hex(bins = 500) + 
  labs(title = "huCD45 vs FSC", y = "FSC-H", x = "FSC-A") + ggcyto_par_set(limits = "instrument") +
  geom_gate(g.singlets) +
  theme_bw() + theme(legend.position = "none", aspect.ratio = 1) + facet_wrap(~name, ncol = 4) 
gs_pop_add(gs,g.singlets) # add gate to GatingSet


recompute(gs)

ggcyto(gs[[1]],aes(x="FSC-A",y="SSC-A"),subset="root")+geom_hex(bins = 200)+ggcyto_par_set(limits = "instrument")
ggcyto(gs[[1]],aes(x="FSC-A",y="SSC-A"),subset="Singlets")+geom_hex(bins = 200)+ggcyto_par_set(limits = "instrument")

ggcyto(gs,aes(x="FSC-A",y="FSC-H"),subset="root")+geom_hex(bins = 100)+geom_gate("Singlets")+
  geom_stats(adjust = 0.8)+ggcyto_par_set(limits = "instrument")+
  facet_wrap(~well,ncol = 10)


g.live <- polygonGate(filterId = "Live","FSC-A"=c(8e4,28e4,28e4,8e4),"SSC-A"=c(3e4,4e4,28e4,28e4)) # define gate
ggcyto(gs[[1]],aes(x="FSC-A",y="SSC-A"),subset="Singlets")+geom_hex(bins = 200)+geom_gate(g.live)+ggcyto_par_set(limits = "instrument") # check gate
gs_pop_add(gs,g.live,parent="Singlets") # add gate to GatingSet
recompute(gs)


ggcyto(gs,aes(x="FSC-A",y="SSC-A"),subset="Singlets")+geom_hex(bins = 100)+geom_gate("Live")+
  geom_stats(adjust = 0.8)+ggcyto_par_set(limits = "instrument")+
  facet_wrap(~well,ncol = 10)

g.huCD45 <- rectangleGate(filterId="huCD45 positive","huCD45"=c(1000, Inf)) # set gate
ggcyto(gs[[1]],aes(x=huCD45),subset="Live")+geom_density(fill="forestgreen")+geom_gate(g.huCD45)+ggcyto_par_set(limits = "instrument")+scale_x_flowJo_biexp() # check gate
gs_pop_add(gs,g.huCD45,parent="Live")
recompute(gs)


ggcyto(gs,aes(x=huCD45),subset="Live",)+geom_density(fill="forestgreen")+geom_gate("huCD45 positive")+
  geom_stats()+
  ggcyto_par_set(limits = "instrument")+scale_x_flowJo_biexp()+
  facet_wrap(~well,ncol = 10)

gs_pop_get_count_fast(gs) %>% head
ps <- data.frame(gs_pop_get_count_fast(gs))

ps$percent_of_parent <- ps$Count/ps$ParentCount*100
psm <- merge(ps,pData(fs),by="name")




library(flowCore)
library(FlowSOM)
library(ggplot2)

# 18. Compute the FlowSOM object

SOM_x <- 10
SOM_y <- 10
n_meta <- 8
seed <- 2020
scaling <- FALSE


names(fs)
exprs(fs)
markers_of_interest <- c("SSC-A", "CD8", "CD4", "CD3", "CD33",
                         "CD11b", "CD25", "mCD45", "huCD45", "CD56", 
                         "CD19", "CD66b")
markers_of_interest <- c("SSC-A", "CD4", "CD3", "CD33",
                         "CD11b", "CD25", "mCD45", "huCD45", "CD56", 
                         "CD19", "CD66b")


fsom <- FlowSOM(input = fs,
                scale = scaling,
                colsToUse = c(4,7:16),
                seed = seed,
                nClus = n_meta,
                xdim = SOM_x, ydim = SOM_y)
saveRDS(fsom, paste(dir_RDS, "fsom.rds"))

PlotDimRed(fsom, colsToUse = fsom$map$colsUsed, colorBy = "metaclusters")



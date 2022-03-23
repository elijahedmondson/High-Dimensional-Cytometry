#' # Identify the files
#' 
#' "C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15701 02Feb2022 Simone/"
#'
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

ws <- open_flowjo_xml("C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/15726 10Mar2022 Simone.wsp")
tail(fj_ws_get_sample_groups(ws))
fj_ws_get_samples(ws, group_id = 1)
fj_ws_get_keywords(ws, 205)[1:30]

gs <- flowjo_to_gatingset(ws, name = 1, path ="C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/")
plot(gs)
gs_pop_get_gate(gs, "hCD45+")
gs_get_compensations(gs)[1]
head(gs_pop_get_stats(gs, xml = TRUE))
gs_pop_get_data(gs)
sampleNames(gs)

fj_ws_get_samples(ws)
fj_ws_get_sample_groups(ws)$groupName


gs_pop_get_count_fast(gs)
ggcyto::autoplot(gs[[1]], c("/scatter/sing"))
ggcyto::autoplot(gs[[1]], c("/scatter/sing/hCD45+"), bins=200)
ggcyto::autoplot(gs[[1]], c("/scatter/sing/hCD45+/Q6: CD3+ , CD4 [PCP55]+"), bins=200)
ggcyto::autoplot(gs[[1]], c("/scatter/sing/hCD45+"), bins=200)

gatingResult <- GetFlowJoLabels(fcs_file, ws,
                                cell_types = cell_types,
                                getData = TRUE)

fsom <- FlowSOM(gatingResult$flowFrame, 
                #compensate = TRUE, 
                transform = TRUE,
                toTransform = c(7:8,11), 
                colsToUse = c(7:8,11),
                nClus = 3,
                seed = 1)





myfiles <- list.files(path="C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/", pattern = ".FCS", ignore.case = TRUE)
#fs <- read.flowSet(myfiles, path="C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/")#, truncate_max_range = FALSE)

wsp_file <- "C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/15726 10Mar2022 Simone.wsp"
fcs_file <- "C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/Samples_Tube_015 Animal 137 blood_015.fcs"


# Specify the cell types of interest for assigning one label per cell
cell_types <- c("/scatter/sing/hCD45+")

# Parse the FlowJo workspace   


gatingResult <- GetFlowJoLabels(fcs_file, wsp_file,
                                cell_types = cell_types,
                                getData = TRUE)

# Check the number of cells assigned to each gate
colSums(gatingResult$matrix)

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
colnames(fs)
# Build a FlowSOM tree
fsom <- FlowSOM(gatingResult$flowFrame, 
                       #compensate = TRUE, 
                       #transform = TRUE,
                       #toTransform = c(7:8,11), 
                       colsToUse = c(7:17),
                       nClus = 20,
                       seed = 1)
PlotStars(fsom)
PlotFlowSOM(fsom, equalNodeSize = F)

PlotPies(fsom, gatingResult$manual,
         backgroundValues = fsom$metaclustering)

PlotManualBars(fsom, manualVector = gatingResult$manual,
               manualOrder = c(cellTypes = "/scatter/sing/hCD45+/Q6: CD3+ , CD4 [PCP55]+"))


dups <- duplicated(fsom)
fsom <- fsom[!dups, ]


PlotDimRed(fsom, cTotal = 10000, colsToUse = fsom$map$colsUsed, colorBy = "metaclusters", check_duplicates = FALSE)

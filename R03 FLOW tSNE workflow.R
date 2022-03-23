# load packages

library(flowCore)
library(FlowSOM)
library(Rtsne)
library(ggplot2)

#############################
### LOAD AND PREPARE DATA ###
#############################
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

ws <- open_flowjo_xml("C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/15726 10Mar2022 Simone.wsp")
gs <- flowjo_to_gatingset(ws, name = 1, path ="C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/")
plot(gs)
fj_ws_get_samples(ws)
fj_ws_get_sample_groups(ws)$groupName
gs_pop_get_count_fast(gs)

ggcyto::autoplot(gs[[1]], c("/scatter/sing"))
ggcyto::autoplot(gs[[1]], c("/scatter/sing/hCD45+"), bins=200)
ggcyto::autoplot(gs[[1]], c("/scatter/sing/hCD45+/Q6: CD3+ , CD4 [PCP55]+"), bins=200)
ggcyto::autoplot(gs[[1]], c("/scatter/sing/hCD45+"), bins=200)


# load data (download file from link above or GitHub repository)

file <- "C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/Samples_Tube_015 Animal 137 blood_015.fcs"
fs <- flowCore::exprs(flowCore::read.FCS(file, transformation = FALSE, truncate_max_range = FALSE))

head(fs)
dim(fs)
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
head(fs)
dim(fs)

# select protein marker columns to use for clustering

marker_cols <- 7:17

# apply arcsinh transformation
# (with standard scale factor of 5 for CyTOF data; alternatively 150 for flow 
# cytometry data; see Bendall et al. 2011, Science, Supplementary Figure S2)

asinh_scale <- 5
fs[, marker_cols] <- asinh(fs[, marker_cols] / asinh_scale)

summary(fs)

# create flowFrame object (required input format for FlowSOM)

data_FlowSOM <- flowCore::flowFrame(fs)

###################
### RUN FLOWSOM ###
###################

# set seed for reproducibility

set.seed(1234)

# run FlowSOM (initial steps prior to meta-clustering)

out <- FlowSOM::ReadInput(data_FlowSOM, transform = FALSE, scale = FALSE)
out <- FlowSOM::BuildSOM(out, colsToUse = marker_cols)
out <- FlowSOM::BuildMST(out)

# optional visualization

FlowSOM::PlotStars(out)

# extract cluster labels (pre meta-clustering) from output object

labels_pre <- out$map$mapping[, 1]

# specify final number of clusters for meta-clustering (can also be selected 
# automatically, but this often does not perform well)

k <- 40

# run meta-clustering

# note: In the current version of FlowSOM, the meta-clustering function 
# FlowSOM::metaClustering_consensus() does not pass along the seed argument 
# correctly, so results are not reproducible. We use the internal function 
# ConsensusClusterPlus::ConsensusClusterPlus() to get around this. However, this
# will be fixed in the next update of FlowSOM (version 1.5); then the following 
# (simpler) code can be used instead:
#seed <- 1234
#out <- FlowSOM::metaClustering_consensus(out$map$codes, k = k, seed = seed)

seed <- 1234
out <- ConsensusClusterPlus::ConsensusClusterPlus(t(out$map$codes), maxK = k, seed = seed)
out <- out[[k]]$consensusClass


# extract cluster labels from output object

labels <- out[labels_pre]

# summary of cluster sizes and number of clusters

table(labels)
length(table(labels))

# save cluster labels

res <- data.frame(cluster = labels)

write.table(res, file = "C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/cluster_labels_FlowSOM.txt", 
            row.names = FALSE, quote = FALSE, sep = "\t")




#################
### RUN RTSNE ###
#################

# subsampling (required due to runtime)

n_sub <- 10000

set.seed(1234)
ix <- sample(1:length(labels), n_sub)

# prepare data for Rtsne (matrix format required)

data_Rtsne <- data[ix, marker_cols]
data_Rtsne <- as.matrix(data_Rtsne)

head(data_Rtsne)
dim(data_Rtsne)

# remove any near-duplicate rows (required by Rtsne)

dups <- duplicated(data_Rtsne)
data_Rtsne <- data_Rtsne[!dups, ]

dim(data_Rtsne)


# run Rtsne (Barnes-Hut-SNE algorithm; runtime: 2-3 min)

# note initial PCA is not required, since we do not have too many dimensions
# (i.e. not thousands, which may be the case in other domains)

set.seed(1234)
out_Rtsne <- Rtsne(data_Rtsne, pca = FALSE, verbose = TRUE)




###################
### CREATE PLOT ###
###################

# load cluster labels (if not still loaded)

file_labels <- "C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/cluster_labels_FlowSOM.txt"
data_labels <- read.table(file_labels, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
labels <- data_labels[, "cluster"]

# select points used by Rtsne

labels_plot <- labels[ix][!dups]
length(labels_plot)  ## should be same as number of rows in data_Rtsne

# prepare Rtsne output data for plot

data_plot <- as.data.frame(out_Rtsne$Y)
colnames(data_plot) <- c("tSNE_1", "tSNE_2")

head(data_plot)
dim(data_plot)  ## should match length of labels_plot (otherwise labels will not match up correctly)

data_plot[, "cluster"] <- as.factor(labels_plot)

head(data_plot)


# plot 2-dimensional t-SNE projection

ggplot(data_plot, aes(x = tSNE_1, y = tSNE_2, color = cluster)) + 
  geom_point(size = 0.2) + 
  coord_fixed(ratio = 1) + 
  ggtitle("t-SNE projection with FlowSOM clustering") + 
  theme_bw()

ggsave("C:/Users/edmondsonef/Desktop/samp/15726 10Mar2022/FlowSOM_Rtsne_plot.pdf", height = 6, width = 7)

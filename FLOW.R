
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("cytofkit")



library(FlowSOM)
browseVignettes("FlowSOM")

library(Rtsne) # Load package



library(flowCore)
browseVignettes("flowCore")

library(ggcyto)
browseVignettes("ggcyto")

library(M3C)
browseVignettes("M3C")

library(flowClust)
browseVignettes("flowClust")

library(flowStats)
browseVignettes("flowStats")


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
# load packages
library(flowCore)
library(FlowSOM)
library(Rtsne)
library(ggplot2)


#############################
### LOAD AND PREPARE DATA ###
#############################

FCS_Path = "C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15709 15Feb2022 Simone/"
list.files(FCS_Path)

FF_2.15.22 <- read.flowSet(path = FCS_Path)

PlotFileScatters(FF_2.15.22)

##TRAINING A FLOWSOM MODEL##
FlowSOM(FF_2.15.22, scale=F, colsToUse=1:3, xdim=10, ydim=10, nClus=10, seed=NULL)






# load data (download file from link above or GitHub repository)

file <- "C:/Users/edmondsonef/Desktop/Humanized Mouse Models/Flow/15709 15Feb2022 Simone/Samples_Tube_014 Animal 141_026.fcs"
data <- flowCore::exprs(flowCore::read.FCS(file, transformation = FALSE, truncate_max_range = FALSE))

head(data)
dim(data)

# select protein marker columns to use for clustering

marker_cols <- 1:18

# apply arcsinh transformation
# (with standard scale factor of 5 for CyTOF data; alternatively 150 for flow 
# cytometry data; see Bendall et al. 2011, Science, Supplementary Figure S2)

asinh_scale <- 5
data[, marker_cols] <- asinh(data[, marker_cols] / asinh_scale)

summary(data)

# create flowFrame object (required input format for FlowSOM)

data_FlowSOM <- flowCore::flowFrame(data)




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

write.table(res, file = "results/cluster_labels_FlowSOM.txt", 
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

file_labels <- "results/cluster_labels_FlowSOM.txt"
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

ggsave("plots/FlowSOM_Rtsne_plot.pdf", height = 6, width = 7)

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


########### 1. Generate Counts From Manual Gating ########### 
########### 
plot_dir <- "C:/Users/edmondsonef/Desktop/Humanized/Flow/Plots/"
data_dir <- "C:/Users/edmondsonef/Desktop/Humanized/Flow/" 
study_dir1 <- "1-05Jan2022/"
study_dir2 <- "2-02Feb2022/"
study_dir3 <- "3-15Feb2022/"
study_dir4 <- "4-28Feb2022/"
study_dir5 <- "5-02Mar2022/"
study_dir6 <- "6-10Mar2022/"
ws <- open_flowjo_xml("C:/Users/edmondsonef/Desktop/Humanized/Flow/5-02Mar2022/15719 02Mar2022 Simone.wsp")
ws
head(fj_ws_get_samples(ws, group_id = 2))
gs <- flowjo_to_gatingset(ws, name = 2, path=paste0(data_dir,study_dir6))
gs_get_pop_paths(gs)
recompute(gs)

png(paste0(plot_dir,"6-10Mar2022-01_Gates.png"), width = 3000, height =1800,res = 165)
plot(gs)
dev.off()
png(paste0(plot_dir,"6-10Mar2022-02_scatter.png"), width = 3000, height =1800,res = 165)
autoplot(gs, "/scatter")
dev.off()
png(paste0(plot_dir,"6-10Mar2022-03_scatter-sing.png"), width = 3000, height =1800,res = 165)
autoplot(gs, "/scatter/sing")
dev.off()
png(paste0(plot_dir,"6-10Mar2022-04_hCD45.png"), width = 3000, height =2400,res = 165)
autoplot(gs, "hCD45+")
dev.off()
png(paste0(plot_dir,"6-10Mar2022-05_CD3-CD4.png"), width = 3000, height =2400,res = 165)
autoplot(gs, "Q6: CD3+ , CD4 [PCP55]+")
dev.off()
png(paste0(plot_dir,"6-10Mar2022-06_CD19.png"), width = 3000, height =2400,res = 165)
autoplot(gs, "Q13: CD3- , CD19 [AFire750]+")
dev.off()
png(paste0(plot_dir,"6-10Mar2022-07_CD33-CD11b.png"), width = 3000, height =2400,res = 165)
autoplot(gs, "Q34: CD33+ , CD11b [AF647]+")
dev.off()


########### 2. Preprocessing and QC ###########
###########  ###########
# 2. Define the general and preprocessing variables

#study_dir <- "1-05Jan2022/"
#study_dir <- "2-02Feb2022/"
#study_dir <- "3-15Feb2022/"
#study_dir <- "4-28Feb2022/"
#study_dir <- "5-02Mar2022/"
study_dir <- "6-10Mar2022/"

data_dir <- "C:/Users/edmondsonef/Desktop/Humanized/Flow/" 
file_pattern <- "\\d.fcs"
reference_file <- read.FCS(paste0(data_dir,study_dir,'Samples_Tube_020 Animal 120 spleen_020.fcs'), truncate_max_range = FALSE)
reference_marker <- "PE-A" # Scatter values will be scaled to have the same range

markers_of_interest <- c("BB515-A",
                         "BB700-P-A",
                         "APC-A",
                         "APC-Cy7-A",
                         "BV421-A",
                         "BV786-A",
                         "BUV395-A",
                         "BUV805-A",
                         "PE-A",
                         "PE-CF594-A",
                         "PE-Cy7-A")
live_gate <- flowCore::polygonGate(filterId = "Live",
                                   .gate = matrix(data = c(60000, 100000, 150000, 
                                                           250000, 250000, 60000, 
                                                           60000, 1.6, 1.9, 2.5,
                                                           2.5, -0.3, -0.3, 1.6),
                                                  ncol = 2,
                                                  dimnames = list(c(), 
                                                                  c("FSC-A", 
                                                                    "APC-Cy7-A"))))

# 3. Define and create the directories
dir_prepr <- paste0(data_dir,study_dir,"1-Preprocessed/") #where the preprocessed data will be stored
dir_QC <- paste0(data_dir,study_dir,"2-QC/") #where the data QC results will be stored
dir_RDS <- paste0(data_dir,study_dir,"3-RDS/") #where the R objects will be stored
dir_results <- paste0(data_dir,study_dir,"4-Results/") #where the results will be stored
dir_raw <- paste0(data_dir,study_dir) #where the raw data is located
path_comp <- "C:/Users/edmondsonef/Desktop/0-Autospill/6-10Mar2022/autospill_compensation.csv" #where comp matrix is located

for (path in c(dir_prepr, dir_QC, dir_RDS, dir_results)){
  dir.create(path)
}

# 4. Prepare some additional information for preprocessing the files 
# given the variable choices of step 2.
files <- list.files(path = dir_raw, pattern = "Samples")
files
channels_of_interest <- GetChannels(object = reference_file,
                                    markers = markers_of_interest, 
                                    exact = FALSE)
compensation_matrix <- read.csv(path_comp, 
                                check.names = FALSE, row.names = 1)

colnames(compensation_matrix) <- sub(" :: .*", "",         
                                     colnames(compensation_matrix))

# Compute transformation list
ff_m <- PeacoQC::RemoveMargins(reference_file, channels_of_interest)
names(ff_m)
exprs(ff_m)
each_col(ff_m, median)

ff_c <- flowCore::compensate(ff_m, compensation_matrix)
translist <- estimateLogicle(ff_c, colnames(compensation_matrix))
ff_t <- flowCore::transform(ff_c, translist)
q5_goal <- quantile(exprs(ff_t)[,reference_marker], 0.05)
q95_goal <- quantile(exprs(ff_t)[,reference_marker], 0.95)
q5_SSCA <- quantile(exprs(ff_t)[,"SSC-A"], 0.05)
q95_SSCA <- quantile(exprs(ff_t)[,"SSC-A"], 0.95)
SSCA_a <- (q95_goal - q5_goal) / (q95_SSCA - q5_SSCA)
SSCA_b <- q5_goal - q5_SSCA * (q95_goal - q5_goal) / (q95_SSCA - q5_SSCA)
translist <- c(translist, 
               transformList("SSC-A", flowCore::linearTransform(a = SSCA_a,
                                                                b = SSCA_b)))
translist

# 5. Read the first fcs file into a flowframe
ff <- read.FCS(paste0(dir_raw, files[7]), truncate_max_range = FALSE)
# 6. Remove margin events
ff_m <- PeacoQC::RemoveMargins(ff, channels_of_interest)
# 7. Compensate
ff_c <- flowCore::compensate(ff_m, compensation_matrix)
# 8. Transform, logicle for marker channels, linear for scatter channel
ff_t <- flowCore::transform(ff_c, translist)
# 9. Remove doublets and filter live cells
ff_s <- PeacoQC::RemoveDoublets(ff_t)
selected_live <- flowCore::filter(ff_s, live_gate)
ff_l <- ff_s[selected_live@subSet, ]
# 10. QC with PeacoQC
PQC <- PeacoQC::PeacoQC(ff = ff_s,
                        channels = channels_of_interest,
                        plot = TRUE, save_fcs = FALSE,
                        output_directory = dir_QC)

# 11. Save the preprocessed data
write.FCS(PQC$FinalFF,
          file = paste0(dir_prepr, files[1]))

# 12. Visualize the preprocessing
filter_plot <- function(ff_pre, ff_post, title, channel_x, channel_y){
  df <- data.frame(x = exprs(ff_pre)[,channel_x],
                   y = exprs(ff_pre)[,channel_y])
  i <- sample(nrow(df), 10000)
  if (!"Original_ID" %in% colnames(exprs(ff_pre))) {
    ff_pre@exprs <- cbind(ff_pre@exprs,
                          Original_ID = seq_len(nrow(ff_pre@exprs)))
  }
  p <- ggplot(df[i,], aes(x = x, y = y)) +
    geom_point(size = 0.5,
               color = ifelse(exprs(ff_pre)[i,"Original_ID"] %in%
                                exprs(ff_post)[,"Original_ID"], 'blue', 'red')) +
    xlab(GetMarkers(ff_pre, channel_x)) + 
    ylab(GetMarkers(ff_pre, channel_y)) +
    theme_minimal() + theme(legend.position = "none") +
    ggtitle(title)
  return(p)
}
to_plot <- list(list(ff_pre = ff,
                     ff_post = ff_m,
                     title = "Removed margin events",
                     channel_x = "BUV395-A",
                     channel_y = "BUV805-A"),
                list(ff_pre = ff_t,
                     ff_post = ff_s,
                     title = "Removed doublets",
                     channel_x = "FSC-A",
                     channel_y = "FSC-H"),
                list(ff_pre = ff_s,
                     ff_post = ff_l,
                     title = "Removed debris and dead cells",
                     channel_x = "FSC-A",
                     channel_y = "BUV395-A"),
                list(ff_pre = ff_l,
                     ff_post = PQC$FinalFF,
                     title = "Removed low quality events",
                     channel_x = "Time",
                     channel_y = "BUV395-A"))

plot_list <- list()
for (plot in to_plot) {
  plot_list[[length(plot_list) + 1]] <- filter_plot(ff_pre = plot$ff_pre,
                                                    ff_post = plot$ff_post,
                                                    title = plot$title,
                                                    channel_x = plot$channel_x,
                                                    channel_y = plot$channel_y)
}

png(paste0(dir_QC, sub("fcs", "png", files[1])), width = 1920)
print(ggpubr::ggarrange(plotlist = plot_list, nrow = 1))
dev.off()

# 13. Run the preprocessing pipeline for all the files
for (file in files){
  ff <- read.FCS(paste0(dir_raw, file), truncate_max_range = FALSE)
  ff_m <- PeacoQC::RemoveMargins(ff, channels_of_interest)
  ff_c <- flowCore::compensate(ff_m, compensation_matrix)
  ff_t <- flowCore::transform(ff_c, translist)
  ff_s <- PeacoQC::RemoveDoublets(ff_t)
  selected_live <- flowCore::filter(ff_s, live_gate)
  ff_l <- ff_s[selected_live@subSet, ]
  
  PQC <- PeacoQC::PeacoQC(ff = ff_l,
                          channels = channels_of_interest,
                          plot = TRUE, save_fcs = FALSE,
                          output_directory = dir_QC)
  
  write.FCS(PQC$FinalFF,
            file = paste0(dir_prepr, file))
  
  to_plot <- list(list(ff_pre = ff,
                       ff_post = ff_m,
                       title = "Removed margin events",
                       channel_x = "BUV395-A",
                       channel_y = "BUV805-A"),
                  list(ff_pre = ff_t,
                       ff_post = ff_s,
                       title = "Removed doublets",
                       channel_x = "FSC-A",
                       channel_y = "FSC-H"),
                  list(ff_pre = ff_s,
                       ff_post = ff_l,
                       title = "Removed debris and dead cells",
                       channel_x = "FSC-A",
                       channel_y = "BUV395-A"),
                  list(ff_pre = ff_l,
                       ff_post = PQC$FinalFF,
                       title = "Removed low quality events",
                       channel_x = "Time",
                       channel_y = "BUV395-A"))

  plot_list <- list()
  for (plot in to_plot) {
    plot_list[[length(plot_list) + 1]] <- filter_plot(ff_pre = plot$ff_pre,
                                                      ff_post = plot$ff_post,
                                                      title = plot$title,
                                                      channel_x = plot$channel_x,
                                                      channel_y = plot$channel_y)
  }
  
  png(paste0(dir_QC, sub("fcs", "png", file)), width = 1920)
  print(ggpubr::ggarrange(plotlist = plot_list, nrow = 1))
  dev.off()
}
  

  ########### 3. Additional QC ###########
  ###########  ###########
  # 14. Perform quality control between all files
  # 14.(A) Plot the signal per channel and per file
  # 14.(A)(i) Define the variables
  #file_names <- sub(".*15_(.*).fcs", "\\1", files)
  write.csv(files, paste0(dir_raw,"names.csv"))
  
  data <- read_excel("C:/Users/edmondsonef/Desktop/Humanized/Flow/Group.Names.xlsx", 
                     sheet = "6-10Mar")
  file_groups <- data$Group
  # 14.(A)(ii) Make the overview plot
  PlotFileScatters(input = paste0(dir_prepr, files),
                   channels = channels_of_interest,
                   #names = file_groups, 
                   legend = T,
                   groups = file_groups, nrow = 4,
                   plotFile = paste0(dir_QC, "file_scatters.png"))
  
  # 14.(B) Perform principal commponent analysis (PCA)
  # 14.(B)(i) Retrieve the median marker expression values per file
  medians <- matrix(data = NA,
                    nrow = length(files), ncol = length(channels_of_interest),
                    dimnames = list(files, channels_of_interest))
  
  for (file in files){
    ff <- read.FCS(paste0(dir_prepr, file))
    medians[file,] <- apply(exprs(ff)[,channels_of_interest], 2, median)
  }
  
  # 14.(B)(ii) Calculate the PCs
  pc <- prcomp(medians, scale. = TRUE)
  
  # 14.(B)(iii) Visualize the PCs
  ggplot(data.frame(pc$x[,1:2], file_groups)) + 
    geom_point(aes(x= PC1, y = PC2, col = file_groups)) +
    theme_minimal()
  ggsave(paste0(dir_QC, "file_PCA.png"), width = 10)
  


  ########### 4. Create Aggegregate Files ###########
  ###########  ###########
  data_dir <- "C:/Users/edmondsonef/Desktop/Humanized/Flow/" 
  #study_dir <- "1-05Jan2022/"
  #study_dir <- "2-02Feb2022/"
  #study_dir <- "3-15Feb2022/"
  #study_dir <- "4-28Feb2022/"
  #study_dir <- "5-02Mar2022/"
  study_dir <- "6-10Mar2022/"
  dir_prepr <- paste0(data_dir,study_dir,"1-Preprocessed/") #where the preprocessed data will be stored
  dir_QC <- paste0(data_dir,study_dir,"2-QC/") #where the data QC results will be stored
  dir_RDS <- paste0(data_dir,study_dir,"3-RDS/") #where the R objects will be stored
  dir_results <- paste0(data_dir,study_dir,"4-Results/") #where the results will be stored
  dir_raw <- paste0(data_dir,study_dir) #where the raw data is located
  
  #dir_group <- paste0(dir_raw,"/1-Preprocessed/NSG/")
  #dir_group <- paste0(dir_raw,"/1-Preprocessed/NSG-IL15/")
  dir_group <- paste0(dir_raw,"/1-Preprocessed/spleen/")
  #dir_group <- paste0(dir_raw,"/1-Preprocessed/")
  
  # 15. Choose the number of cells to include in the aggregate file
  n <- 700000
  
  # 16. Make an aggregate file
  set.seed(2020)
  
  
  files <- list.files(path = dir_group, pattern = "Samples")
  files
  paste0(dir_group, files)
  agg <- AggregateFlowFrames(paste0(dir_group, files),
                             cTotal = n,
                             writeOutput = TRUE,
                             outputFile = paste0(dir_group, "aggregate.fcs"))
  


  ########### 5. Train FlowSOM model ###########
  ###########  ###########
  agg <- read.FCS(paste0(dir_group,'aggregate.fcs'), truncate_max_range = FALSE)
  
  
  #Level 1 - create model to separate human cells from mouse
  #SOM_x <- 10
  #SOM_y <- 10
  n_meta <- 10
  seed <- 2020
  scaling <- FALSE
  #scaling <- TRUE
  
  ###All...###
  ###All...###
  ###All...###
  fsom <- FlowSOM(input = agg,
                  scale = F, 
                  transform = T,
                  toTransform = c(7:17),
                  colsToUse = c(7:17),
                  seed = seed,
                  nClus = n_meta)
  PlotStars(fsom = fsom,
            backgroundValues = fsom$metaclustering)
  
  tsne <- Rtsne::Rtsne(fsom$map$codes, perplexity = 6)
  PlotStars(fsom = fsom, view = tsne$Y,
            backgroundValues = fsom$metaclustering)
  
  FlowSOMmary(fsom = fsom,
              plotFile = paste0(dir_group, "ALL_fsom.trans_summary.pdf"))
  saveRDS(fsom, paste(dir_RDS, "fsom.rds"))
  fsom$prettyColnames
  
  QueryStarPlot(fsom, query="BUV805-A")
  ManualVector()
  QueryMultiple(fsom = fsom_level1, cellTypes = )
    ###All...###
  ###All...###
  ###All...###
  ###All...###
  
  
  
  fsom_level1 <- FlowSOM(input = agg,
                  scale = F,
                  transform = T,
                  toTransform = c(7:17),
                  colsToUse = c("BUV395-A","BUV805-A"),
                  seed = seed,
                  nClus = 5)

  PlotLabels(fsom_level1, labels = fsom_level1$metaclustering)
  p <- PlotMarker(fsom_level1, "BUV805-A")
  p <- PlotMarker(fsom_level1, "BUV395-A")
  print(p, newpage = T)
  
  
  PlotStars(fsom = fsom_level1,backgroundValues = fsom_level1$metaclustering)
  Plot2DScatters(fsom = fsom_level1, 
                 channelpairs = list(c("BUV805-A", "BUV395-A")),
                 metaclusters = 1:5, 
                 plotFile = paste0(dir_group, "mCD45-hCD45_level1.png"))
  FlowSOMmary(fsom = fsom_level1,
              plotFile = paste0(dir_group, "L1_fsom_summary.pdf"))
   
  saveRDS(fsom_level1, paste(dir_RDS, "fsom_level1.rds"))
  #MC 1 and 5 are the huCD45+
  #MC 2, 3, and 4 are the mCD45+
  
  # Subset the original fcs file
  fsom_tmp <- NewData(fsom_level1, agg)
  clustering <- GetMetaclusters(fsom_tmp)
  agg_tmp_hCD45 <- agg[clustering %in% c(4),]
  agg_tmp_mCD45 <- agg[clustering %in% c(2,3),]
  
  #HUMAN: Create hCD45 subset
  fsom_level2_hCD45 <- FlowSOM(input = agg_tmp_hCD45,
                         scale = F,
                         colsToUse = c(4,7:17),
                         #colsToUse = c(7:12,15:17),
                         seed = 2020)
  PlotStars(fsom = fsom_level2_hCD45,
            backgroundValues = fsom_level2_hCD45$metaclustering)
  FlowSOMmary(fsom = fsom_level2_hCD45,
              plotFile = paste0(dir_group, "L2_hCD45_fsom_summary.pdf"))
  
  saveRDS(fsom_level2_hCD45, paste(dir_RDS, "fsom_level2_hCD45.rds"))
  
  #MOUSE: Create mCD45 subset
  fsom_level2_mCD45 <- FlowSOM(input = agg_tmp_mCD45,
                               scale = F,
                               colsToUse = c(7:17),
                               seed = 2020)
  saveRDS(fsom_level2_mCD45, paste(dir_RDS, "fsom_level2_mCD45.rds"))
  FlowSOMmary(fsom = fsom_level2_mCD45,
              plotFile = paste0(dir_group, "L2_mCD45_fsom_summary.pdf")) 
  
  
  
  
  ########### 6. Test Quality ###########
  ###########  ###########
  
  fsom <- fsom_level2_hCD45
  agg <- agg_tmp_hCD45
  
  
  # 20. Check the FlowSOM quality
  # 20.(A) Make 2D scatter plots
  # 20.(A)(i) Specify the parameters
  fsom$prettyColnames
  channel_pairs = list(c("FSC-A", "SSC-A"),
                       c("BB700-P-A", "BB515-A"),
                       c("BV421-A", "APC-Cy7-A"),
                       c("APC-A", "BV786-A"),
                       c("PE-CF594-A", "BV786-A"),
                       c("BV786-A", "PE-A"),
                       c("BUV805-A", "BUV395-A"))
  metaclusters_of_interest <- seq_len(n_meta)
  clusters_of_interest <- NULL
  
  # 20.(A)(ii) Make the 2D scatter plots
  Plot2DScatters(fsom = fsom,
                 channelpairs = channel_pairs,
                 metaclusters = metaclusters_of_interest,
                 clusters = clusters_of_interest,
                 plotFile = paste0(dir_group, "fsom_2D_scatters.png"))
  

  ########### 7. Test with Manual Gating ###########
  ########### 
  fsom <- fsom_level2_hCD45
  agg <- agg_tmp_hCD45
  fsom$prettyColnames
  #wspFile = "C:/Users/edmondsonef/Desktop/Humanized/Flow/5-02Mar2022/15719 02Mar2022 Simone.wsp"
  wspFile = "C:/Users/edmondsonef/Desktop/Humanized/Flow/2-02Feb2022/15701 02Feb2022 Simone.wsp"
  wspFile = "C:/Users/edmondsonef/Desktop/Humanized/Flow/6-10Mar2022/15726 10Mar2022 Simone.wsp"
  ws <- open_flowjo_xml(wspFile)
  ws
  head(fj_ws_get_samples(ws, group_id = 2))
  files <- list.files(path = dir_group, pattern = "Samples")
  head(files)
  
  # 20.(B) Check the consistency with manual labeling
  # 20.(B)(i) Extract the gating information from the wsp file
  gating <- GetFlowJoLabels(files = files, cellTypes = "hCD45+",
                            wspFile = wspFile, group =2,
                            path = dir_raw)
  ####*EFE TIP: PHYSICALLY REMOVE CONTROLS, etc####
  
  # 20.(B)(ii) Get an overview of the gatenames and define the cell types of interest
  print(levels(gating[[1]][["manual"]]))
  cell_types_of_interest <- c(#"Unlabeled", 
                              #"hCD45+",
                              "Q6: CD3+ , CD4 [PCP55]+",
                              "Q7: CD3+ , CD4 [PCP55]-",
                              "Q8: CD3- , CD4 [PCP55]-",
                              "Q9: CD3- , CD8 [FITC]+",
                              "Q10: CD3+ , CD8 [FITC]+",
                              "Q11: CD3+ , CD8 [FITC]-",
                              "Q12: CD3- , CD8 [FITC]-",
                              "Q13: CD3- , CD19 [AFire750]+",
                              "Q14: CD3+ , CD19 [AFire750]+",
                              "Q15: CD3+ , CD19 [AFire750]-",
                              "Q16: CD3- , CD19 [AFire750]-",
                              "Q17: CD3- , CD56+",
                              "Q18: CD3+ , CD56+",
                              "Q19: CD3+ , CD56-",
                              "Q20: CD3- , CD56-",
                              "Q29: CD66b [PEDazz]- , CD11b [AF647]+",
                              "Q30: CD66b [PEDazz]+ , CD11b [AF647]+",
                              "Q31: CD66b [PEDazz]+ , CD11b [AF647]-",
                              "Q32: CD66b [PEDazz]- , CD11b [AF647]-",
                              "Q33: CD33- , CD11b [AF647]+",
                              "Q34: CD33+ , CD11b [AF647]+",
                              "Q35: CD33+ , CD11b [AF647]-",
                              "Q36: CD33- , CD11b [AF647]-",
                              "Q37: CD25- , CD3+",
                              "Q38: CD25+ , CD3+",
                              "Q39: CD25+ , CD3-",
                              "Q40: CD25- , CD3-")
  cell_types_of_interest <- c("hCD45+", "Q6: CD3+ , CD4 [PCP55]+",
                              "Q10: CD3+ , CD8 [FITC]+",
                              "Q13: CD3- , CD19 [AFire750]+",
                              "Q30: CD66b [PEDazz]+ , CD11b [AF647]+",
                              "Q34: CD33+ , CD11b [AF647]+")
  cell_types_of_interest <- c("Q15: CD3+ , CD19 [AFire750]-",
                              "Q31: CD66b [PEDazz]+ , CD11b [AF647]-",
                              "Q10: CD3+ , CD8 [FITC]+",
                              "Q17: CD3- , CD56+",
                              "Q34: CD33+ , CD11b [AF647]+","Q6: CD3+ , CD4 [PCP55]+",
                              "Q13: CD3- , CD19 [AFire750]+",
                              "Q29: CD66b [PEDazz]- , CD11b [AF647]+")
  
  # 20.(B)(iii) Compile the labels of the aggregate file
  aggregate_labels <- c()
  for (file in unique(exprs(agg)[, "File"])) {
    aggregate_labels <- c(aggregate_labels, 
                          as.character(ManualVector(gating[[file]][["matrix"]],
                                                    cell_types_of_interest)
                                       [exprs(agg)[, "Original_ID"]
                                         [exprs(agg)[, "File"] == file]]))
  }
  # 20.(B)(iv) Show the manual labeling on the FlowSOM tree
  PlotPies(fsom = fsom,
           cellTypes = factor(aggregate_labels, levels = c("Unlabeled",
                                                           cell_types_of_interest)))

  ggsave(paste0(dir_group, "Human_fsom_Labeled.pdf"))
  
  
  
  # 19.(B)(v) Calculate the purity of the FlowSOM clustering
  Purity(realClusters = aggregate_labels,
         predictedClusters = GetClusters(fsom))
  
  # 20.(C) Inspect the file contribution per cluster
  # 20.(C)(i) Specify a color vector (optional)
  file_colors <- c("#990000", "#cc0000", "#ff0000", #Different shades within the groups
                   "#1d1d77", "#2b3b92", "#3859ac", "#4677c7") 
  file_colors <- c("#990000", "#4677c7") 
  # 20.(C)(ii) Show the file contribution
  p <- PlotPies(fsom = fsom,
                cellTypes = factor(files[fsom$data[,"File"]]),
                #cellTypes = cell_types_of_interest,
                colorPalette = file_colors)
  AddStarsPies(p = p, # Legend to show how it should be
               arcs = data.frame(
                 x0 = rep(0, length(files)),
                 y0 = rep(0, length(files)),
                 start = seq(0, 2 * pi, length.out = 8)[-8],
                 end = seq(0, 2 * pi, length.out = 8)[-1],
                 value = rep(1, length(files)),
                 Markers = files),
               colorPalette = file_colors)
  ggsave(paste0(dir_group, "fsom_filecontribution.pdf"))
  
PlotManualBars(fsom, list_insteadof_plots = T, 
               manualVector = factor(aggregate_labels, levels = c(cell_types_of_interest)))
ggsave(paste0(dir_group, "fsom_filecontribution.pdf"))  


  ########## 8. Discovery and downstream analysis ###########
  ########## 
  fsom <- readRDS("C:/Users/edmondsonef/Desktop/Humanized/Flow/6-10Mar2022/3-RDS/ fsom_level2_hCD45.rds")
                 
  # 21. Explore the FlowSOM result
  
  # 21.(B) Look for nodes with a specific pattern
  # 21.(B)(i) Specify the query
  query <- list("B cells" = c("CD19 [APC-F750]" = "high", "CD3" = "low"),
                #"Activated T cells" = c("CD3" = "high", "CD25"="high"),
                "Neutrophils" = c("CD33"="high","CD66b [PE-Dazz]" = "high"),
                "CD33/CD11b" = c("CD33"="high","CD11b [AF647]" = "high"),
                "CD33" = c("CD33"="high"),
                "B cells Activated" = c("CD19 [APC-F750]" = "high", "CD25"="high"),
                #"Mouse CD45+" = c("mCD45" = "high"),
                #"Human CD45+" = c("huCD45" = "high"),
                "NK Cell" = c("CD56" = "high", "CD3" = "low"),
                "NK T Cell" = c("CD56" = "high", "CD3" = "high"),
                "CD4 T Cell Activated" = c("CD4 [PerCPCy55]" = "high", "CD3"="high", "CD25"="high"),
                "CD4 T Cell" = c("CD4 [PerCPCy55]" = "high", "CD3"="high"),
                "CD8 T Cell" = c("CD8 [FITC]" = "high", "CD3"="high"),
                "CD8 T Cell Activated" = c("CD8 [FITC]" = "high", "CD3"="high", "CD25"="high"))#,
                #"T cells" = c("CD3" = "high"))
  
  # 21.(B)(ii) Retrieve the cluster labels based on the query
  labels <- QueryMultiple(fsom = fsom,
                          cellTypes = query,
                          plotFile = paste0(dir_group, "fsom_QueryStarPlot.pdf"))
  
  # 21.(B)(iii) Show the retrieved labels on the FlowSOM tree
  PlotVariable(fsom = fsom,
               variable = labels)
  ggsave(paste0(dir_results, "fsom_query.pdf"))
  
  # 22. Get features per fcs file
  # Specify the variables of interest
  
  types <- c("counts", "percentages", "MFIs")
  
  MFIs <- c("CD3", "CD19 [APC-F750]",
            "CD8 [FITC]","CD56","CD25","CD33",
            "CD66b [PE-Dazz]","CD11b [AF647]","CD4 [PerCPCy55]")
  
  MFIs <- c("CD3", "CD19 [APC-F750]",
            "CD56","CD4 [PerCPCy55]")  
  # Get the features
  features <- GetFeatures(fsom = fsom,
                          files = paste0(dir_group, files),
                          filenames = files,
                          type = types,
                          MFI = MFIs)
  

  ########## 9. Compare Groups ###########
  ########## 
  file_names = paste0(dir_group, files)
  # 23. Define the groups and feature you would want to compare.
  feature <- "MFIs"
  
  stat <- "fold changes"
  
  # 24. Compare the 2 groups of interest
  stats <- GroupStats(features$cluster_MFIs,                     
                      groups = list("NSG-IL15" = files[1], "NSG" = files[2]))
  
  # 25. Show the findings of step 24 on the trees
  # Define the plotting variables
  stat_levels <- c(paste0(names(grouplist)[2], " underrepresented compared to ",
                          names(grouplist)[1]),
                   paste0(names(grouplist)[1], " underrepresented compared to ",
                          names(grouplist)[2]),
                   "--")
  colors <- c("blue", "red", "white")
  
  # Show statistical findings on FlowSOM trees
  cluster_stat <- stats[stat,]
  cluster_stat <- factor(ifelse(cluster_stat < -2.5, stat_levels[1],
                                ifelse(cluster_stat > 2.5, stat_levels[2],
                                       stat_levels[3])), 
                         levels = stat_levels)
  cluster_stat[is.na(cluster_stat)] <- stat_levels[3]
  gr_1 <- PlotStars(fsom = fsom, title = names(grouplist)[1], 
                    nodeSizes = stats[paste0("medians ", names(grouplist)[1]),], 
                    backgroundValues = cluster_stat,
                    backgroundColors = colors, 
                    list_insteadof_ggarrange = TRUE)
  gr_2 <- PlotStars(fsom = fsom, title = names(grouplist)[2], 
                    nodeSizes = stats[paste0("medians ", names(grouplist)[2]),],
                    backgroundValues = cluster_stat,
                    backgroundColors = colors,
                    list_insteadof_ggarrange = TRUE)
  ggpubr::ggarrange(plotlist = list(gr_1$tree, gr_2$tree, gr_2$starLegend, 
                                    gr_2$backgroundLegend), 
                    heights = c(3,1))
  ggsave(paste0(dir_results, "NSG_vs_IL15_fsom_groups.pdf"), width = 10, height = 7.5)
  
  
  
  p <- PlotVariable(fsom, title = "Fold change group 1 vs. group 2",
                    variable = C_stats["fold changes", ])
  print(p, newpage = FALSE)

  ########## 10. Map new data onto FlowSOM object ###########
  ########## 
  
  # 26. Map new data on the FlowSOM object
  for (file in files){
    ff_prepr <- read.FCS(paste0(dir_prepr, file))
    ff_raw <- read.FCS(paste0(dir_raw, file))
    fsom_tmp <- NewData(fsom = fsom,
                        input = ff_prepr)
    clustering <- GetClusters(fsom_tmp)
    clustering_raw <- matrix(data = rep(0, nrow(exprs(ff_raw))),
                             ncol = 1, dimnames = list(c(), "FlowSOM"))
    clustering_raw[exprs(ff_prepr)[,"Original_ID"]] <- clustering
    ff_tmp <- flowCore::fr_append_cols(ff_raw, clustering_raw)
    write.FCS(ff_tmp, paste0(dir_prepr, "FlowSOM_", file))
  }
  
  
  

  ########## 11. Additional ###########
  ########## Applying FlowSOM to files or groups separately and then meta-cluster on all ####
    # Compute separate FlowSOM objects
  fsom_KO <- FlowSOM(input = paste0(dir_prepr, files[1:3]),
                     scale = FALSE, colsToUse = channels_of_interest,
                     seed = 2020)
  
  fsom_WT <- FlowSOM(input = paste0(dir_prepr, files[4:7]),
                     scale = FALSE, colsToUse = channels_of_interest,
                     seed = 2020)
  
  # Extract the cluster median fluorescence intensity values (MFIs)
  MFI_KO <- GetClusterMFIs(fsom = fsom_KO, prettyColnames = TRUE, colsUsed = TRUE)
  rownames(MFI_KO) <- paste0("KO", rownames(MFI_KO))
  MFI_WT <- GetClusterMFIs(fsom = fsom_WT, prettyColnames = TRUE, colsUsed = TRUE)
  rownames(MFI_WT) <- paste0("WT", rownames(MFI_WT))
  
  # Obtain the meta-clusters by hierarchical clustering
  all_clusters <- rbind(MFI_KO, MFI_WT)
  hclust <- hclust(dist(all_clusters))
  metaclustering <- cutree(hclust, 15) #MC 14 corresponds to the NK cells
  
  # Generate one clustering heatmap from all clusters
  ann <- data.frame(cohort = rep(c("KO", "WT"), each = 100), 
                    row.names = rownames(all_clusters))
  p <- pheatmap::pheatmap(t(all_clusters), cluster_rows = F, cutree_cols = 15, 
                          cellwidth = 5, fontsize_col = 3, annotation_col = ann,
                          cluster_cols = hclust)
  ggsave(p, filename = paste0(dir_results, "Higher_level_clustering.pdf"), width = 17)
  
  # Generate the meta-cluster percentages boxplots
  fsom_KO$metaclustering <- factor(unname(metaclustering[1:100]), levels = 1:15)
  fsom_WT$metaclustering <- factor(unname(metaclustering[101:200]), levels =  1:15)
  perc_KO <- GetFeatures(fsom = fsom_KO, 
                         files = paste0(dir_prepr, files[1:3]),
                         level = "metaclusters", type = "percentages", 
                         filenames = files[1:3])
  perc_WT <- GetFeatures(fsom = fsom_WT, 
                         files = paste0(dir_prepr, files[4:7]),
                         level = "metaclusters", type = "percentages", 
                         filenames = files[4:7])
  
  df <- data.frame(rbind(perc_KO[[1]], perc_WT[[1]])*100, 
                   cohort = rep(c("KO", "WT"), c(3, 4)), check.names = FALSE)
  df_g <- tidyr::gather(df, "MC", "percentage", -cohort)
  ggplot(df_g, aes(x = cohort, y = percentage)) +
    geom_boxplot() +
    facet_wrap(~MC, scales = "free") +
    theme_minimal()
  ggsave(filename = "Results/FlowSOM_boxplot.pdf", width = 10, height = 10)
  

  


  ########## 12. Subset and hierarchical approach ###########
  wspFile = "C:/Users/edmondsonef/Desktop/Humanized/Flow/5-02Mar2022/15719 02Mar2022 Simone.wsp"
  ws <- open_flowjo_xml("C:/Users/edmondsonef/Desktop/Humanized/Flow/5-02Mar2022/15719 02Mar2022 Simone.wsp")
  ws
  head(fj_ws_get_samples(ws, group_id = 5))
  files <- list.files(path = dir_raw, pattern = "Samples")
  files
  
  # 20.(B) Check the consistency with manual labeling
  # 20.(B)(i) Extract the gating information from the wsp file
  gating <- GetFlowJoLabels(files = files,
                            wspFile = wspFile, group =5,
                            path = dir_raw)
  
  # 20.(B)(ii) Get an overview of the gatenames and define the cell types of interest
  print(levels(gating[[1]][["manual"]]))
  cell_types_of_interest <- c("hCD45+","Q6: CD3+ , CD4 [PCP55]+",
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
  
  
  # 20.(B)(iii) Compile the labels of the aggregate file
  aggregate_labels <- c()
  for (file in unique(exprs(agg)[, "File"])) {
    aggregate_labels <- c(aggregate_labels, 
                          as.character(ManualVector(gating[[file]][["matrix"]],
                                                    cell_types_of_interest)
                                       [exprs(agg)[, "Original_ID"]
                                         [exprs(agg)[, "File"] == file]]))
  }
    
    
    
    
    
    
  # Read in preprocessed fcs file, lymphocyte panel
  ff <- read.FCS(paste0(data_dir,study_dir3,'1-Preprocessed/NSG/Samples_Tube_015 Animal 142_027.fcs'), truncate_max_range = FALSE)
  #manual_labels <- readRDS(paste0(dir_raw, "attachments/lympho_labels.rds"))
  
  # Perform a first level clustering to isolate the lymphocytes
  fsom_level1 <- FlowSOM(input = ff,
                         scale = FALSE, nClus = 3,
                         colsToUse = c(13:14),
                         seed = 2020)
  fsom_level1$prettyColnames
  PlotStars(fsom = fsom_level1,
            backgroundValues = fsom_level1$metaclustering)
  
  # Inspect the 2D scatter plots to identify the meta-clusters of interest
  Plot2DScatters(fsom = fsom_level1, 
                 channelpairs = list(c("BUV805-A", "BUV395-A")),
                 metaclusters = 1:2, 
                 plotFile = paste0(data_dir, "hierarchy_level1.png")) 
  
  #MC 1, 4 and 5 are the lymphocytes (CD3+, CD161-)
  
  # Subset the original fcs file
  fsom_tmp <- NewData(fsom_level1, ff)
  clustering <- GetMetaclusters(fsom_tmp)
  ff_tmp <- ff[clustering %in% c(2),]
  
  # Perform a second level clustering to characterize the lymphocytes
  fsom_level2 <- FlowSOM(input = ff_tmp,
                         scale = FALSE,
                         colsToUse = c(7:17),
                         seed = 2020)
  
  # Plot the lymphocytes FlowSOM tree
  PlotStars(fsom = fsom_level2,
            backgroundValues = fsom_level2$metaclustering)
  
  # Show the manual labels on the FlowSOM trees
  PlotPies(fsom = fsom_level1, cellTypes = manual_labels,
           title = "First level clustering")
  PlotPies(fsom = fsom_level2, cellTypes = factor(manual_labels[clustering %in% c(1, 4, 5)]),
           backgroundValues = fsom_level2$metaclustering, title = "Second level clustering")

  
  
   
  

  ########## 13. PlotDimRed ###########
  
  PlotDimRed(fsom, cTotal = 500,
    colsToUse = fsom$map$colsUsed,
    colorBy = "metaclusters",
    dimred = Rtsne::Rtsne)
  
  
  
  
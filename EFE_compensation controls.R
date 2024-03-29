

library(flowCore)

######################
# Parameter settings #
######################

# Files to use
files = list.files(path="C:/Users/edmondsonef/Desktop/Humanized/Flow/6-10Mar2022/", pattern = "Samples", ignore.case = TRUE)

# Compensation matrix
compensationFile = file.path("C:/Users/edmondsonef/Desktop/10Mar2022/table_compensation/autospill_compensation.csv")
colsToCompensate = c(7:17)

# GatingML file
#gatingFile <- file.path(dataID,"attachments/151030_21-10-15_Tube_028.fcs_gates.xml")

wsp_file <- "C:/Users/edmondsonef/Desktop/Humanized/Flow/6-10Mar2022/15726 10Mar2022 Simone.wsp"
gatingFile <- open_flowjo_xml(wsp_file)

# Indices of the gates of interest
tail(fj_ws_get_sample_groups(gatingFile))
fj_ws_get_samples(gatingFile, group_id = 1)
gs <- flowjo_to_gatingset(gatingFile, name = 1, path ="C:/Users/edmondsonef/Desktop/Humanized/Flow/6-10Mar2022/")
plot(gs)
autoplot(gs[[1]])
gate_ids <- c("/scatter/sing/", "/scatter/sing/hCD45+",
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

# Cell types to use as final labels
#cellTypes <- c("B cells","T cells","NK cells","NK T cells",
#               "Macrophages","DCs","Neutrophils")
cellTypes <- c("Q6: CD3+ , CD4 [PCP55]+",
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
# Columns to use for analysis
colsToCluster = c(7:17)

###################
# Helper function #
###################

ProcessGatingML <- function(flowFrame,gatingFile,gateIDs,
                            cellTypes,silent=FALSE){
  gating_xml <- XML::xmlToList(XML::xmlParse(gatingFile))
  flowEnv <- new.env()
  flowUtils::read.gatingML(gatingFile, flowEnv)
  #  A. Read the gates from xml
  filterList <- list()
  for(cellType in names(gateIDs)){
    filterList[[cellType]] <-  flowEnv[[
      as.character(gating_xml[[gateIDs[cellType]]]$.attrs["id"])
    ]]
  }
  #  B. Process the fcs file for all specified gates
  results <- matrix(NA,nrow=nrow(flowFrame),ncol=length(gateIDs),
                    dimnames = list(NULL,names(gateIDs)))
  for(cellType in names(gateIDs)){
    if(!silent){message(paste0("Processing ",cellType))}
    results[,cellType] <- flowCore::filter(flowFrame,
                                           filterList[[cellType]])@subSet
  }
  #  C. Assign one celltype to each cell
  manual <- rep("Unknown",nrow(flowFrame))
  for(celltype in cellTypes){
    manual[results[,celltype]] <- celltype
  }
  manual <- factor(manual,levels = c("Unknown",cellTypes))
  
  list("matrix"=results,"manual"=manual)
}

########################
# Preprocess the data  #
########################

# Read the compensation matrix
comp <- read.csv(compensationFile,row.names=1,check.names = FALSE)
colnames(comp) <- rownames(comp) <- gsub(" ::.*","",colnames(comp))

for(file in files){
  message(paste0("Processing ",file))
  # Load the raw data
  ff <- read.FCS(file)
  # and compensate
  ff <- compensate(ff,comp)
  ff@description$SPILL <- comp
  colnames(ff)[colsToCompensate] <- paste("Comp-",
                                          colnames(ff)[colsToCompensate],sep="")
  
  # Extract manual gating
  gatingRes <- ProcessGatingML(ff, gatingFile, gate_ids, cellTypes) 
  gatingMatrix <- gatingRes$matrix
  manual <- gatingRes$manual
  
  # Cells to use as input for the algorithms
  selected <- gatingMatrix[,"Live single cells"]
  
  # Save selected cells to fcs file for algorithms which can only take
  # a file as input. File is already compensated, so identity matrix as
  # compensation matrix
  new_comp <- diag(length(colnames(ff)[7:19]))
  colnames(new_comp) <- colnames(ff)[7:19]
  ff@description$SPILL <- new_comp
  write.FCS(ff[selected,],
            file=gsub(".fcs","_selected.fcs",file))
  # Save also a subset of only the first 10.000 cells for faster processing
  write.FCS(ff[selected,][1:10000,],
            file=gsub(".fcs","_selected_10000.fcs",file))
  
  # Transform the data
  ff_t <- transform(ff,transformList(colnames(ff)[7:19],logicleTransform()))
  
  # Save results so this step can be skipped next time
  save(ff,ff_t,selected,manual,gatingMatrix,colsToCluster,
       file=gsub(".fcs",".Rdata",file))
}

## Load multiple packages simultaneously
Packages <- c("dplyr","ggplot2","plotly","reshape2","stringr","tidyr")
lapply(Packages, library, character.only = TRUE)

#################################################################################

## General Functions
zf_qPCR_Analysis <- function(Chemical_DF, Chemical, TranscriptList,
                             Direction, TimePt, NewFile) {
  
  ## Subset beta-actin Ct values from DF
  Bactin <- Chemical_DF[grep("Bactin", Chemical_DF$X),]
  
  ## Iterate through each transcript
  for (x in 1:length(TranscriptList)) {
    
    Target = TranscriptList[x]
    
    ## Subset Ct values for transcript of interest
    TargetDF <- Chemical_DF[grep(Target, Chemical_DF$X),]
    ## Merge beta-actin and transcript of interest Ct values into single DF
    Ct_Ref_Target <- rbind(Bactin, TargetDF)
    
    DF_long <- melt(Ct_Ref_Target[,2:(ncol(Ct_Ref_Target))], id.vars=c("X.1"))
    
    ## Remove rows with "NA" Ct values 
    DF_long <- DF_long[complete.cases (DF_long$value),]
    colnames(DF_long) <- c("Condition", "SampleInfo", "CtValue")
    
    ## Index biological replicates for transcript
    ExptNames_Index <- unique(DF_long$SampleInfo)
    ExptNames_Index <- data.frame(SampleInfo = ExptNames_Index)
    
    ExptNames_Index$BRIndex <- with(ExptNames_Index, 
                                      as.numeric(factor(SampleInfo,
                                                        levels=unique(SampleInfo))))
    
    DF_long$BRIndex <- ExptNames_Index[match(DF_long$SampleInfo,
                                               ExptNames_Index$SampleInfo),2]
    
    PValueDF <- data.frame(matrix(ncol = 12, nrow = 0))
    PValue_nRow = 0
      
    FoldChange <- data.frame(matrix(ncol = 8, nrow = 0))
    colnames(FoldChange) <- c("Condition", "BRIndex", "Summary", 
                              "FoldChange", "TechRep", "Transcript", "TimePt",
                              "Chemical")
    a=0
      
    BioRepCount = as.numeric(length(unique(DF_long$BRIndex)))
    
    ## Iterate through each biological replicate
    for (j in 1:BioRepCount) {
      
      ## Subset by biological replicate of interest
      BioRep <- DF_long[grep(paste("^",j,"$", sep=""), DF_long$BRIndex),]
      
      ControlReference <- BioRep[grep("Control_Reference", BioRep$Condition),]
      ControlReference <- subset(ControlReference, select = c("CtValue"))
      ControlReference$CtValue <- as.numeric(ControlReference$CtValue)
      
      ControlTarget <- BioRep[grep("Control_Target", BioRep$Condition),]
      ControlTarget <- subset(ControlTarget, select = c("CtValue"))
      ControlTarget$CtValue <- as.numeric(ControlTarget$CtValue)
      TechnicalRep_Control <- as.numeric(nrow(ControlTarget))
      
      TreatmentReference <- BioRep[grep("Treatment_Reference", BioRep$Condition),]
      TreatmentReference <- subset(TreatmentReference, select = c("CtValue"))
      TreatmentReference$CtValue <- as.numeric(TreatmentReference$CtValue)
      
      TreatmentTarget <- BioRep[grep("Treatment_Target", BioRep$Condition),]
      TreatmentTarget <- subset(TreatmentTarget, select = c("CtValue"))
      TreatmentTarget$CtValue <- as.numeric(TreatmentTarget$CtValue)
      TechnicalRep_Treatment <- as.numeric(nrow(TreatmentTarget))
      
      ## Compute possible delta Ct values 
      deltaCt_Control <- data.frame()
      for (x in 1:length(ControlReference$CtValue)) {
        for (i in 1:length((ControlTarget)$CtValue)) {
          deltaCt_Control[x,i] <-
            (ControlTarget$CtValue[[i]] - ControlReference$CtValue[[x]])
        }}
      
      deltaCt_Treatment <- data.frame()
      for (x in 1:length(TreatmentReference$CtValue)) {
        for (i in 1:length((TreatmentTarget)$CtValue)) {
          deltaCt_Treatment[x,i] <-
            (TreatmentTarget$CtValue[[i]] - TreatmentReference$CtValue[[x]])
        }}
      
      ## Compute statistics 
      PValueDF[PValue_nRow + 1, 1] <- "Treatment"
      PValueDF[PValue_nRow + 1, 2] <- as.numeric(j)
      
      deltaCt_Control_Mean <- data.frame(matrix
                                         (ncol = ncol(deltaCt_Control), nrow = 0))
      for (x in 1:ncol(deltaCt_Control)) {
        deltaCt_Control_Mean[1,x] <- mean(unlist(deltaCt_Control[,x]))
      }
      
      deltaCt_Treatment_Mean <- data.frame(matrix
                                           (ncol = ncol(deltaCt_Treatment), nrow = 0))
      for (x in 1:ncol(deltaCt_Treatment)) {
        deltaCt_Treatment_Mean[1,x] <- mean(unlist(deltaCt_Treatment[,x]))
      }
      
      ## Run t-test
      qPCR_Test <- t.test(unlist(deltaCt_Control_Mean), unlist(deltaCt_Treatment_Mean),
                          ## unequal variance between two samples
                          var.equal = FALSE,
                          alternative = Direction)
      
      PValueDF[PValue_nRow + 1, 3] <- qPCR_Test$estimate[[1]]
      PValueDF[PValue_nRow + 1, 4] <- sd(unlist(deltaCt_Control_Mean))
      
      PValueDF[PValue_nRow + 1, 5] <- qPCR_Test$estimate[[2]]
      PValueDF[PValue_nRow + 1, 6] <- sd(unlist(deltaCt_Treatment_Mean))
      
      PValueDF[PValue_nRow + 1, 7] <- qPCR_Test$p.value
      PValueDF[PValue_nRow + 1, 8] <- qPCR_Test$alternative
      PValueDF[PValue_nRow + 1, 9] <- qPCR_Test$parameter
      PValueDF[PValue_nRow + 1, 10] <- qPCR_Test$method
      PValueDF[PValue_nRow + 1, 11] <- Target
      PValueDF[PValue_nRow + 1, 12] <- paste("(", qPCR_Test$conf.int[1], " - ", 
                                             qPCR_Test$conf.int[2], ")",
                                             sep = "")
      
      PValue_nRow = PValue_nRow + 1
      
      ## Compute possible delta delta Ct values 
      deltadeltaCt_Control <- data.frame()
      for (x in 1:length(unlist(deltaCt_Control))) {
        for (i in 1:length(unlist(deltaCt_Control))) {
          deltadeltaCt_Control[x,i] <-
            (unlist(deltaCt_Control)[[i]] - unlist(deltaCt_Control)[[x]])
        }}
      
      deltadeltaCt_Treatment <- data.frame()
      for (x in 1:length(unlist(deltaCt_Control))) {
        for (i in 1:length(unlist(deltaCt_Treatment))) {
          deltadeltaCt_Treatment[x,i] <-
            (unlist(deltaCt_Treatment)[[i]] - unlist(deltaCt_Control)[[x]])
        }}
      
      a = nrow(FoldChange)
      for (x in 1:length(unlist(deltadeltaCt_Control))) {
        FoldChange[a+x,1] <- "Control"
        FoldChange[a+x,2] <- as.numeric(j)
        FoldChange[a+x,3] <- "Summary"
        FoldChange[a+x,4] <- 2 ** (-1*(unlist(deltadeltaCt_Control)[[x]]))
        FoldChange[a+x,5] <- TechnicalRep_Control
        FoldChange[a+x,6] <- Target
        FoldChange[a+x,7] <- TimePt
        FoldChange[a+x,8] <- Chemical
      }
      
      a = nrow(FoldChange)
      for (x in 1:length(unlist(deltadeltaCt_Treatment))) {
        FoldChange[a+x,1] <- "Treatment"
        FoldChange[a+x,2] <- as.numeric(j)
        FoldChange[a+x,3] <- "Summary"
        FoldChange[a+x,4] <- 2 ** (-1*(unlist(deltadeltaCt_Treatment)[[x]]))
        FoldChange[a+x,5] <- TechnicalRep_Treatment
        FoldChange[a+x,6] <- Target
        FoldChange[a+x,7] <- TimePt
        FoldChange[a+x,8] <- Chemical
      }}
    
    FoldChange$ExptName <- ExptNames_Index[match(FoldChange$BRIndex,
                                                 ExptNames_Index$BRIndex),1]
    
    colnames(PValueDF) <- c("Condition", "BRIndex", 
                            "Control_Mean", "Control_StdDev", 
                            "Treatment_Mean", "Treatment_StdDev",
                            "pValue", "Direction", "DF", "Method",
                            "Transcript", "Percent95_CI")
    
    ## Summarize fold-changes by transcript & biological replicate 
    FoldChange_Stats_BioRep <- 
      FoldChange %>% group_by(Condition, Transcript, BRIndex, 
                              ExptName, TimePt,
                              Chemical) %>% summarise(Mean = mean(FoldChange),
                                                      Std = sd(FoldChange),
                                                      TechRep = mean(TechRep))
    
    FoldChange_Stats_Pre <- full_join(FoldChange_Stats_BioRep, PValueDF, 
                                      by = c('Condition' , 'BRIndex', 'Transcript'))
    
    FoldChange_Stats_All <- 
      FoldChange %>% group_by(Condition, Transcript,
                              Summary, TimePt,
                              Chemical) %>% summarise(Mean = mean(FoldChange), 
                                                      Std = sd(FoldChange), 
                                                      TechRep = mean(TechRep))
    
    FoldChange_Stats <- full_join(FoldChange_Stats_Pre, FoldChange_Stats_All,
                                  by = c("Transcript", "Condition", 
                                         "Mean", "Std", "TechRep", "TimePt",
                                         "Chemical"))
    
    FoldChange_Stats <- 
      FoldChange_Stats %>% select(ExptName, TimePt, Chemical, BRIndex,
                                  Condition, Mean, Std, TechRep,
                                  pValue, Direction,
                                  Control_Mean, Control_StdDev,
                                  Treatment_Mean, Treatment_StdDev,
                                  DF, Method, Transcript, Percent95_CI)
    
    a <- nrow(ExptNames_Index) +1
    print(Target)
    print(FoldChange_Stats[a:(2*nrow(ExptNames_Index)), 1:4])
    print(FoldChange_Stats [a:(2*nrow(ExptNames_Index)),8:14])
    
    FileName = paste("TEMP/", Target, "_", Chemical, 
                     "_", TimePt, ".csv", 
                     sep = "")
    
    ## Append statistics to csv file
    write.table(FoldChange_Stats, file = FileName,
                sep=",", row.names = FALSE, col.names = TRUE,
                append = FALSE)
  }
  
  ## Merge individual stats file
  compileDF <- data.frame()
  
  for (x in 1:length(TranscriptList)) {
      FileName = paste("TEMP/", TranscriptList[x], 
                       "_", Chemical, "_", TimePt,".csv", sep = "")
      
      Transcript_DF <- read.csv(file = FileName,
                                header = TRUE, sep = ",", skipNul = TRUE, 
                                stringsAsFactors = FALSE)
      
      compileDF <- rbind(compileDF, Transcript_DF)
      file.remove(FileName)
    }
  
  StatsFile = paste("TEMP/", NewFile, "_Stats.csv", sep = "")
  write.table(compileDF, file = StatsFile,
                sep=",", row.names = FALSE, col.names = TRUE,
                append = FALSE)
  
  ## Generate CSV file with only summary stats for FIGURES
  ## Subset summary statistic 
  DataGraph_Summary <- dplyr::filter(compileDF,is.na(ExptName)) 
  
  ## Remove empty columns 
  DataGraph_Summary <- DataGraph_Summary[colSums(!is.na(DataGraph_Summary)) > 0]
  DataGraph_Summary <- subset(DataGraph_Summary, select = -c(Chemical, TechRep) )
  
  for (x in 1:nrow(DataGraph_Summary)) {
    if (isTRUE(DataGraph_Summary$Condition[x] == "Control")) {
      DataGraph_Summary$Condition[x] <- "DMSO"
    } else {
      DataGraph_Summary$Condition[x] <- Chemical
    }}
  
  DataGraph_Summary$Condition = factor(DataGraph_Summary$Condition,
                                       levels = c("DMSO", Chemical))
  
  DataGraph_Summary$Name <- paste(DataGraph_Summary$Condition, 
                                  TimePt,
                                  sep = "_")

  SummaryFile = paste("TEMP/", NewFile, "_Summary.csv", sep = "")
  write.table(DataGraph_Summary, file = SummaryFile,
              sep=",", row.names = FALSE, col.names = TRUE,
              append = FALSE)
}

plot_annot <- function(NewFile, TranscriptList) {
  
  ## Calculate y-axis of annotation 
  DataGraph_Summary$FCmax <- NA
  
  for (x in 1:nrow(DataGraph_Summary)) {
    if (isTRUE(DataGraph_Summary$Condition[x] == "DMSO")) {
      if(isTRUE(DataGraph_Summary$Mean[x] >= DataGraph_Summary$Mean[x+1])) {
        DataGraph_Summary$FCmax[x]  <- round(DataGraph_Summary$Mean[x] + 
                                               DataGraph_Summary$Std[x], digits = 2)
      } else {
        DataGraph_Summary$FCmax[x+1]  <- round(DataGraph_Summary$Mean[x+1] + 
                                               DataGraph_Summary$Std[x+1], digits = 2)
      }}}
  
  print("y-axis of annotation layer")
  cat("\n")
  print(DataGraph_Summary[,c("Condition", "Transcript", "FCmax")])
  cat("\n")

  ## Calculate annotation p-values
  DataGraph_Stats <- read.csv(file = paste("TEMP/", NewFile, "_Stats.csv", sep=""),
                              header = TRUE, sep = ",", skipNul = TRUE, 
                              stringsAsFactors = TRUE)
  DataGraph_Stats <- 
    DataGraph_Stats[ , which(names(DataGraph_Stats) %in% 
                               c("ExptName", "Transcript", "pValue"))]
  DataGraph_Stats <- DataGraph_Stats[complete.cases (DataGraph_Stats$pValue),]
  DataGraph_Stats$pValue <- round(DataGraph_Stats$pValue, digits = 5)
  
  ## Order p-values by lowest to highest
  DataGraph_Stats <- DataGraph_Stats[order(DataGraph_Stats$pValue),]
  
  for (x in 1:length(TranscriptList)) {
    print(TranscriptList[x])
    
    for (k in 1:nrow(DataGraph_Stats)) {
      if (isTRUE(DataGraph_Stats$Transcript[k] == TranscriptList[x])) {
        print(DataGraph_Stats$pValue[k])
      }}}}

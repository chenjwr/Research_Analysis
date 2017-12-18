library("dplyr")
library("ggplot2")
library("tidyr")
library("reshape2")
library("stringr")

#################################################################################

qPCR_Analysis_Function <- function(ChemicalList, TranscriptList,
                                      Direction, TimePt) {
  
  for (k in 1:length(ChemicalList)) {
    Chemical = ChemicalList[k]
    
    InputFile = paste("~/Documents/", 
                      ChemicalList[k], ".csv", sep = "")
    
    Chemical_DF <- read.csv(file = InputFile,
                            header = TRUE, sep = ",", skipNul = TRUE, 
                            stringsAsFactors = FALSE)
    
    Bactin <- Chemical_DF[grep("Bactin", Chemical_DF$X),]
    
    for (x in 1:length(TranscriptList)) {
      Target = TranscriptList[x]
      
      TargetDF <- Chemical_DF[grep(Target, Chemical_DF$X),]
      Ct_Ref_Target <- rbind(Bactin, TargetDF)

      DF_long <- melt(Ct_Ref_Target[,2:5], id.vars=c("X.1"))
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
      
      FoldChange <- data.frame(matrix(ncol = 7, nrow = 0))
      colnames(FoldChange) <- c("Condition", "BRIndex", "Summary", 
                                "FoldChange", "TechRep", "Transcript", "TimePt")
      a=0
      
      BioRepCount = as.numeric(length(unique(DF_long$BRIndex)))
      
      for (j in 1:BioRepCount) {
        BioRep <- DF_long[grep(j, DF_long$BRIndex),]
        
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
        
        deltaCt_Control_Mean <- data.frame(matrix(ncol = ncol(deltaCt_Control), nrow = 0))
        deltaCt_Treatment_Mean <- data.frame(matrix(ncol = ncol(deltaCt_Treatment), nrow = 0))
        
        for (x in 1:ncol(deltaCt_Control)) {
          deltaCt_Control_Mean[1,x] <- mean(unlist(deltaCt_Control[,x]))
        }
        
        for (x in 1:ncol(deltaCt_Treatment)) {
          deltaCt_Treatment_Mean[1,x] <- mean(unlist(deltaCt_Treatment[,x]))
        }
        
        qPCR_Test <- t.test(unlist(deltaCt_Control_Mean), unlist(deltaCt_Treatment_Mean),
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
            deltadeltaCt_Control[i,x] <-
              (unlist(deltaCt_Control)[[x]] - unlist(deltaCt_Control)[[i]])
          }}
        
        deltadeltaCt_Treatment <- data.frame()
        for (x in 1:length(unlist(deltaCt_Treatment))) {
          for (i in 1:length(unlist(deltaCt_Control))) {
            deltadeltaCt_Treatment[i,x] <-
              (unlist(deltaCt_Treatment)[[x]] - unlist(deltaCt_Control)[[i]])
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
                                ExptName, TimePt) %>% summarise(Mean = mean(FoldChange),
                                                                Std = sd(FoldChange),
                                                                TechRep = mean(TechRep))
      
      FoldChange_Stats_Pre <- full_join(FoldChange_Stats_BioRep, PValueDF, 
                                        by = c('Condition' , 'BRIndex', 'Transcript'))
      
      FoldChange_Stats_All <- 
        FoldChange %>% group_by(Condition, Transcript,
                                Summary, TimePt) %>% summarise(Mean = mean(FoldChange), 
                                                               Std = sd(FoldChange), 
                                                               TechRep = mean(TechRep))
      
      FoldChange_Stats <- full_join(FoldChange_Stats_Pre, FoldChange_Stats_All,
                                    by = c("Transcript", "Condition", 
                                           "Mean", "Std", "TechRep", "TimePt"))
      
      FoldChange_Stats <- FoldChange_Stats %>% select(ExptName, TimePt, BRIndex,
                                                      Condition, Mean, Std, TechRep,
                                                      pValue, Direction,
                                                      Control_Mean, Control_StdDev,
                                                      Treatment_Mean, Treatment_StdDev,
                                                      DF, Method, Transcript, Percent95_CI)
      a <- nrow(ExptNames_Index) +1
      print(FoldChange_Stats[a:(2*nrow(ExptNames_Index)), 1:3])
      print(FoldChange_Stats [a:(2*nrow(ExptNames_Index)),7:13])
      
      FileName = paste("~/Documents/", Target, "_", Chemical, 
                       "_", TimePt, ".csv", 
                       sep = "")
      
      ## Append statistics to csv file
      write.table(FoldChange_Stats, file = FileName,
                  sep=",", row.names = FALSE, col.names = TRUE,
                  append = FALSE)
      
      DataGraph_Summary <- FoldChange_Stats_All
      
      for (x in 1:nrow(DataGraph_Summary)) {
        if (isTRUE(DataGraph_Summary$Condition[[x]] == "Control")) {
          DataGraph_Summary$Condition[[x]] <- "DMSO"
        } else {
          DataGraph_Summary$Condition[[x]] <- Chemical
        }}
      
      DataGraph_Summary$Condition = factor(DataGraph_Summary$Condition,
                                           levels = c("DMSO", Chemical))
      
      FileName = paste("~/Documents/", Chemical, "_Summary.csv", 
                       sep = "")
      
      ## Append statistics to csv file
      write.table(DataGraph_Summary, file = FileName,
                  sep=",", row.names = FALSE, col.names = FALSE,
                  append = TRUE)
      
      
    }}}

#################################################################################

for (k in 1:length(ChemicalList)) {
  DF <- data.frame()
  
  for (x in 1:length(TranscriptList)) {
    FileName = paste("~/Documents/", TranscriptList[x], "_",
                     ChemicalList[k], "_33-56h.csv", sep = "")
    Transcript_DF <- read.csv(file = FileName,
                              header = TRUE, sep = ",", skipNul = TRUE, 
                              stringsAsFactors = FALSE)
    
    DF <- rbind(DF, Transcript_DF)
    }
  
  FileName2 = paste("~/Documents/FoldChange_", ChemicalList[k], 
                    "_Stats.csv", sep = "")

  write.table(DF, file = FileName2,
              sep=",", row.names = FALSE, col.names = TRUE,
              append = FALSE)
  }

#################################################################################

for (k in 1:length(ChemicalList)) {
  
  FileName = paste("~/Documents/", ChemicalList[k], "_Summary.csv", sep = "")
  
  DataGraph_Summary <- read.csv(file = FileName,
                                header = FALSE, sep = ",", skipNul = TRUE, 
                                stringsAsFactors = FALSE)
  
  colnames(DataGraph_Summary) <- c("Condition", "Transcript", "Summary", "TimePt",
                                   "Mean", "Std", "TechRep")
  
  TreatmentName <- ChemicalList[k]
  TreatmentName <- TreatmentName %>% str_replace("-.*", "")
  
  for (j in 1:nrow(DataGraph_Summary)) {
    if (DataGraph_Summary$Condition[j] == ChemicalList[k]) {
      DataGraph_Summary$Condition[j] <- TreatmentName
    }}
  
  DataGraph_Summary$Condition = factor(DataGraph_Summary$Condition,
                                       levels = c("DMSO", TreatmentName))
  
  DataGraph_Summary$Name <- paste(DataGraph_Summary$Condition, 
                                  DataGraph_Summary$TimePt,
                                  sep = "_")
  
  TreatmentName2 = paste(TreatmentName, "_33-56h", sep = "")

  DataGraph_Summary$Name = factor(DataGraph_Summary$Name,
                                  levels = c("DMSO_33-56h", TreatmentName2))
  
  DataGraph_Summary$Transcript = factor(DataGraph_Summary$Transcript,
                                        levels = c("scxa", "sox9a", "col2a1", "myod1"))
  
  ## Input graph title, y-axis, save location 
  qPCR_Plot <- 
    ggplot(data=DataGraph_Summary,
           aes(x=DataGraph_Summary$Transcript,
               y=DataGraph_Summary$Mean,
               fill = Condition
           )) +
    geom_bar(stat="identity", 
             position=position_dodge(),
             colour="black",
             size = 0.3,
             show.legend = TRUE
    ) +
    geom_errorbar(aes(ymin = (DataGraph_Summary$Mean - DataGraph_Summary$Std), 
                      ymax = (DataGraph_Summary$Mean + DataGraph_Summary$Std)),
                  width = 0.1,
                  position = position_dodge(.9)) +
    ggtitle("(32-56 hpf) _ ZF") +
    xlab("") + 
    ylab("Relative Expression") + 
    theme_classic() +
    theme(legend.position="bottom",
          plot.title = element_text(hjust = 0.5, face="bold", size=12.5),
          axis.text.y = element_text(size=12.5, colour="black"),
          axis.text.x = element_text(size=12.5, colour="black"),
          axis.title.y = element_text(size=15)
    ) +
    coord_cartesian(ylim = c(0, 2))
  plot(qPCR_Plot)
  
  FileGraph = paste("~/Documents/FoldChange_", ChemicalList[k],
                    "_Stats.tiff", sep = "")
  
  tiff(file = FileGraph,
       width = 1600, height = 1800,
       units = "px",
       res = 350)
  plot(qPCR_Plot)
  dev.off()
}

#################################################################################
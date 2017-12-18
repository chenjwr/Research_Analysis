library("qqman")
library("gaston")
library("plyr")
library("ggplot2")
library("dplyr")
library("cowplot")

#################################################################

## Original GWAS of Gestational Age
Peru_Gest_MAF_LOCO_MLMA <- read.table(
  file="~/Documents/Peru.annot",
  header=TRUE, sep="")

manhattan(Peru_Gest_MAF_LOCO_MLMA, 
          chr = "Chr", bp = "bp", p = "p",
          suggestiveline = FALSE,
          col = c("#7094db", "#008080"),
          main = "Original Gestational Age GWAS", 
          cex = 0.5, cex.axis = 0.8)

## Resequencing analysis of Gestational Age
GestationalAge_LOCO_MLMA <- read.table(
  file="~/Documents/gest.mlma",
  header=TRUE, sep="")

manhattan((GestationalAge_LOCO_MLMA[complete.cases(GestationalAge_LOCO_MLMA$p),]), 
          chr = "Chr", bp = "bp", p = "p",
          suggestiveline = FALSE,
          col = c("#7094db", "#008080"),
          main = "Resequenced Gestational Age GWAS", 
          cex = 0.5, cex.axis = 0.8)

#################################################################

## Format as ANNOVAR input file
## CHR, START, END, REFerence allele, ALTernative allele, OPTIONAL

GestationalAge_ANNO_A1Ref <- data.frame(GestationalAge_LOCO_MLMA$Chr,
                                        GestationalAge_LOCO_MLMA$bp,
                                        GestationalAge_LOCO_MLMA$bp,
                                        GestationalAge_LOCO_MLMA$A1,
                                        GestationalAge_LOCO_MLMA$A2
                                        )

colnames(GestationalAge_ANNO_A1Ref) <- c("CHROM", "START", "END", "REF", "ALT")

Peru_Gest_ANNO_A1Ref <- data.frame(Peru_Gest_MAF_LOCO_MLMA$Chr,
                                   Peru_Gest_MAF_LOCO_MLMA$bp,
                                   Peru_Gest_MAF_LOCO_MLMA$bp,
                                   Peru_Gest_MAF_LOCO_MLMA$A1,
                                   Peru_Gest_MAF_LOCO_MLMA$A2
                                   )

colnames(Peru_Gest_ANNO_A1Ref) <- c("CHROM", "START", "END", "REF", "ALT")

GestationalAge_ANNO_A2Ref <- data.frame(GestationalAge_LOCO_MLMA$Chr,
                                        GestationalAge_LOCO_MLMA$bp,
                                        GestationalAge_LOCO_MLMA$bp,
                                        GestationalAge_LOCO_MLMA$A2,
                                        GestationalAge_LOCO_MLMA$A1
                                        )

colnames(GestationalAge_ANNO_A2Ref) <- c("CHROM", "START", "END", "REF", "ALT")

Peru_Gest_ANNO_A2Ref <- data.frame(Peru_Gest_MAF_LOCO_MLMA$Chr,
                                   Peru_Gest_MAF_LOCO_MLMA$bp,
                                   Peru_Gest_MAF_LOCO_MLMA$bp,
                                   Peru_Gest_MAF_LOCO_MLMA$A2,
                                   Peru_Gest_MAF_LOCO_MLMA$A1
                                   )

colnames(Peru_Gest_ANNO_A2Ref) <- c("CHROM", "START", "END", "REF", "ALT")

## Export as ANNOVAR input file
## Annotated versions saved August 10-11, 2017**

write.table(GestationalAge_ANNO_A1Ref,
            file = "~/Documents/Gest_A1Ref.avinput",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep="\t")

write.table(Peru_Gest_ANNO_A1Ref,
            file = "~/Documents/Peru_A1Ref.avinput",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep="\t")

write.table(GestationalAge_ANNO_A2Ref,
            file = "~/Documents/Gest_A2Ref.avinput",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep="\t")

write.table(Peru_Gest_ANNO_A2Ref,
            file = "~/Documents/Peru_A2Ref.avinput",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep="\t")

## Annotate with wANNOVAR [http://wannovar.wglab.org/index.php]
## Parameters: Reference Genome hg19 \ Input format ANNOVAR \ 
## Gene definition RefSeq Gene \ Individual Analysis All annotations

#################################################################

## Peru_Gest_ANNO_Exome has 11959 entries for both A1/A2 = REF allele
Peru_Gest_ANNO_A1Ref_Exome <- read.csv(
  file="~/Documents/Peru_Gest_ANNO_A1Ref_exome.csv",
  header=TRUE, sep=",", 
  na.strings = c("NA", "","."), stringsAsFactors = FALSE)

Peru_Gest_ANNO_A2Ref_Exome <- read.csv(
  file="~/Documents/Peru_Gest_ANNO_A2Ref_exome.csv",
  header=TRUE, sep=",", 
  na.strings = c("NA", "","."), stringsAsFactors = FALSE)

## Compare column contents for A1/A2 = REF allele
Peru_Gest_ANNO_Exome_Compare <- data.frame()

for (x in 1:length(colnames(Peru_Gest_ANNO_A1Ref_Exome))) {
  Peru_Gest_ANNO_Exome_Compare[x,1] <- 
    all(Peru_Gest_ANNO_A1Ref_Exome[[x]] == Peru_Gest_ANNO_A2Ref_Exome[[x]],
        na.rm = TRUE)
  }

Peru_Gest_ANNO_Exome_Diff <- which(Peru_Gest_ANNO_Exome_Compare == FALSE)

print(colnames(Peru_Gest_ANNO_A1Ref_Exome[Peru_Gest_ANNO_Exome_Diff]))

Peru_Gest_ANNO_A1Ref_Exome_P <- 
  Peru_Gest_ANNO_A1Ref_Exome %>% select (
    "Chr", "Start", "Func.refgene", "Gene.refgene", 
    "GeneDetail.refgene", "ExonicFunc.refgene", 
    "AAChange.refgene", "dbSNP")

colnames(Peru_Gest_ANNO_A1Ref_Exome_P) <- 
  c("Chr", "Start", "Func.refgene", "Gene.refgene", 
    "GeneDetail.refgene.A1REF_Ex", "ExonicFunc.refgene.A1REF_Ex", "AAChange.refgene.A1REF_Ex", "dbSNP")

Peru_Gest_ANNO_A2Ref_Exome_P <- 
  Peru_Gest_ANNO_A2Ref_Exome %>% select (
    "Chr", "Start", "GeneDetail.refgene", "ExonicFunc.refgene", "AAChange.refgene")

colnames(Peru_Gest_ANNO_A2Ref_Exome_P) <- 
  c("Chr", "Start", "GeneDetail.refgene.A2REF_Ex", "ExonicFunc.refgene.A2REF_Ex", 
    "AAChange.refgene.A2REF_Ex")

Peru_Gest_ANNO_Exome <- full_join(Peru_Gest_ANNO_A1Ref_Exome_P,
                                  Peru_Gest_ANNO_A2Ref_Exome_P, 
                                  by = c("Chr", "Start"))

## Peru_Gest_ANNO_Genome has 456942 entries for both A1/A2 = REF allele
Peru_Gest_ANNO_A1Ref_Genome <- read.csv(
  file="~/Documents/Peru_Gest_ANNO_A1Ref_genome.csv",
  header=TRUE, sep=",", 
  na.strings = c("NA", "","."), stringsAsFactors = FALSE)

Peru_Gest_ANNO_A2Ref_Genome <- read.csv(
  file="~/Documents/Peru_Gest_ANNO_A2Ref_genome.csv",
  header=TRUE, sep=",", 
  na.strings = c("NA", "","."), stringsAsFactors = FALSE)

## Compare column contents for A1/A2 = REF allele
Peru_Gest_ANNO_Genome_Compare <- data.frame()

for (x in 1:length(colnames(Peru_Gest_ANNO_A1Ref_Genome))) {
  Peru_Gest_ANNO_Genome_Compare[x,1] <- 
    all(Peru_Gest_ANNO_A1Ref_Genome[[x]] == Peru_Gest_ANNO_A2Ref_Genome[[x]],
        na.rm = TRUE)
  }

Peru_Gest_ANNO_Genome_Diff <- which(Peru_Gest_ANNO_Genome_Compare == FALSE)

print(colnames(Peru_Gest_ANNO_A1Ref_Genome[Peru_Gest_ANNO_Genome_Diff]))

Peru_Gest_ANNO_A1Ref_Genome_P <- 
  Peru_Gest_ANNO_A1Ref_Genome %>% select (
    "Chr", "Start", "Func.refgene", "Gene.refgene", 
    "GeneDetail.refgene","ExonicFunc.refgene", "AAChange.refgene", "dbSNP")

colnames(Peru_Gest_ANNO_A1Ref_Genome_P) <- 
  c("Chr", "Start", "Func.refgene", "Gene.refgene", 
    "GeneDetail.refgene.A1REF_Ge", "ExonicFunc.refgene.A1REF_Ge", 
    "AAChange.refgene.A1REF_Ge", "dbSNP")

Peru_Gest_ANNO_A2Ref_Genome_P <- 
  Peru_Gest_ANNO_A2Ref_Genome %>% select (
    "Chr", "Start", "GeneDetail.refgene", "ExonicFunc.refgene", "AAChange.refgene")

colnames(Peru_Gest_ANNO_A2Ref_Genome_P) <- 
  c("Chr", "Start", "GeneDetail.refgene.A2REF_Ge", 
    "ExonicFunc.refgene.A2REF_Ge", "AAChange.refgene.A2REF_Ge")

Peru_Gest_ANNO_Genome <- full_join(Peru_Gest_ANNO_A1Ref_Genome_P,
                                   Peru_Gest_ANNO_A2Ref_Genome_P, 
                                   by = c("Chr", "Start"))

## Merge ANNOVAR annotations for Peru_Gest_MAF_LOCO_MLMA 
Peru_Gest_ANNO_ExGe <- full_join(Peru_Gest_ANNO_Exome,
                                 Peru_Gest_ANNO_Genome,
                                 by = c("Chr", "Start", "Func.refgene", 
                                        "Gene.refgene", "dbSNP"))

save(Peru_Gest_ANNO_ExGe, 
     file="~/Documents/Peru_Gest_ANNO_ExGe.Rda")

#################################################################

## GestationalAge_ANNO_Genome has different REF alleles: A1 (342472) / A2 (341239)
GestationalAge_ANNO_A1Ref_Genome <- read.csv(
  file="~/Documents/GestationalAge_ANNO_A1Ref_genome.csv",
  header=TRUE, sep=",", 
  na.strings = c("NA", "","."), stringsAsFactors = FALSE)

GestationalAge_ANNO_A2Ref_Genome <- read.csv(
  file="~/Documents/GestationalAge_ANNO_A2Ref_genome.csv",
  header=TRUE, sep=",", 
  na.strings = c("NA", "","."), stringsAsFactors = FALSE)

GestationalAge_ANNO_Genome_OverLap <- inner_join(GestationalAge_ANNO_A1Ref_Genome, 
                                                 GestationalAge_ANNO_A2Ref_Genome,
                                                 by = c("Chr", "Start"))

GestationalAge_ANNO_Genome_OverLap_A1Ref <- GestationalAge_ANNO_Genome_OverLap[,3:127]

GestationalAge_ANNO_Genome_OverLap_A2Ref <- GestationalAge_ANNO_Genome_OverLap[,128:252]

GestationalAge_ANNO_Genome_OverLap_Compare <- data.frame()

for (x in 1:125) {
  GestationalAge_ANNO_Genome_OverLap_Compare[x,1] <- 
    all(GestationalAge_ANNO_Genome_OverLap_A1Ref[[x]] == 
          GestationalAge_ANNO_Genome_OverLap_A2Ref[[x]],
        na.rm = TRUE)
}

GestationalAge_ANNO_Genome_OverLap_Diff <- 
  which(GestationalAge_ANNO_Genome_OverLap_Compare == FALSE)

print(colnames(
  GestationalAge_ANNO_Genome_OverLap_A1Ref[GestationalAge_ANNO_Genome_OverLap_Diff]))

GestationalAge_ANNO_Genome_OverLap <- 
  GestationalAge_ANNO_Genome_OverLap %>% select(
    Chr, Start, Func.refgene.x, Gene.refgene.x, 
    GeneDetail.refgene.x, ExonicFunc.refgene.x, AAChange.refgene.x, 
    dbSNP.x, GeneDetail.refgene.y, ExonicFunc.refgene.y, AAChange.refgene.y)

colnames(GestationalAge_ANNO_Genome_OverLap) <- 
  c("Chr", "Start", "Func.refgene", "Gene.refgene", 
    "GeneDetail.refgene.A1REF_Ge", "ExonicFunc.refgene.A1REF_Ge", 
    "AAChange.refgene.A1REF_Ge", "dbSNP", 
    "GeneDetail.refgene.A2REF_Ge", "ExonicFunc.refgene.A2REF_Ge", 
    "AAChange.refgene.A2REF_Ge")

GestationalAge_ANNO_Genome_Single_A1Ref <- anti_join(GestationalAge_ANNO_A1Ref_Genome, 
                                                     GestationalAge_ANNO_A2Ref_Genome,
                                                     by = c("Chr", "Start"))

GestationalAge_ANNO_Genome_Single_A1Ref <- 
  GestationalAge_ANNO_Genome_Single_A1Ref %>% select(
    Chr, Start, Func.refgene, Gene.refgene, 
    GeneDetail.refgene, ExonicFunc.refgene, AAChange.refgene, dbSNP)

colnames(GestationalAge_ANNO_Genome_Single_A1Ref) <- 
  c("Chr", "Start", "Func.refgene","Gene.refgene", 
    "GeneDetail.refgene.A1REF_Ge", "ExonicFunc.refgene.A1REF_Ge", 
    "AAChange.refgene.A1REF_Ge","dbSNP")

GestationalAge_ANNO_Genome_Single_A2Ref <- anti_join(GestationalAge_ANNO_A2Ref_Genome, 
                                                     GestationalAge_ANNO_A1Ref_Genome,
                                                     by = c("Chr", "Start"))

GestationalAge_ANNO_Genome_Single_A2Ref <- 
  GestationalAge_ANNO_Genome_Single_A2Ref %>% select(
    Chr, Start, Func.refgene, Gene.refgene, GeneDetail.refgene, 
    ExonicFunc.refgene, AAChange.refgene, dbSNP)

colnames(GestationalAge_ANNO_Genome_Single_A2Ref) <- 
  c("Chr", "Start", "Func.refgene", "Gene.refgene", 
    "GeneDetail.refgene.A2REF_Ge", "ExonicFunc.refgene.A2REF_Ge", 
    "AAChange.refgene.A2REF_Ge", "dbSNP")

GestationalAge_ANNO_Genome_Single <- 
  full_join(GestationalAge_ANNO_Genome_Single_A1Ref, 
            GestationalAge_ANNO_Genome_Single_A2Ref,
            by = c("Chr", "Start", "Func.refgene", "Gene.refgene", "dbSNP"))

GestationalAge_ANNO_Genome <- 
  full_join(GestationalAge_ANNO_Genome_OverLap,
            GestationalAge_ANNO_Genome_Single,
            by = c("Chr", "Start", "Func.refgene", "Gene.refgene", 
                   "GeneDetail.refgene.A1REF_Ge", "ExonicFunc.refgene.A1REF_Ge", 
                   "AAChange.refgene.A1REF_Ge", "dbSNP", "GeneDetail.refgene.A2REF_Ge", 
                   "ExonicFunc.refgene.A2REF_Ge", "AAChange.refgene.A2REF_Ge"))

## GestationalAge_ANNO_Exome has different REF alleles: A1 (1983) / A2 (2067)
GestationalAge_ANNO_A1Ref_Exome <- read.csv(
  file="~/Documents/GestationalAge_ANNO_A1Ref_exome.csv",
  header=TRUE, sep=",", 
  na.strings = c("NA", "","."), stringsAsFactors = FALSE)

GestationalAge_ANNO_A2Ref_Exome <- read.csv(
  file="~/Documents/GestationalAge_ANNO_A2Ref_exome.csv",
  header=TRUE, sep=",", 
  na.strings = c("NA", "","."), stringsAsFactors = FALSE)

GestationalAge_ANNO_Exome_OverLap <- inner_join(GestationalAge_ANNO_A1Ref_Exome, 
                                                GestationalAge_ANNO_A2Ref_Exome,
                                                by = c("Chr", "Start"))

GestationalAge_ANNO_Exome_OverLap_A1Ref <- GestationalAge_ANNO_Exome_OverLap[,3:127]

GestationalAge_ANNO_Exome_OverLap_A2Ref <- GestationalAge_ANNO_Exome_OverLap[,128:252]

GestationalAge_ANNO_Exome_OverLap_Compare <- data.frame()

for (x in 1:125) {
  GestationalAge_ANNO_Exome_OverLap_Compare[x,1] <- 
    all(GestationalAge_ANNO_Exome_OverLap_A1Ref[[x]] == 
          GestationalAge_ANNO_Exome_OverLap_A2Ref[[x]],
        na.rm = TRUE)
}

GestationalAge_ANNO_Exome_OverLap_Diff <- 
  which(GestationalAge_ANNO_Exome_OverLap_Compare == FALSE)

print(colnames(
  GestationalAge_ANNO_Exome_OverLap_A1Ref[
    GestationalAge_ANNO_Exome_OverLap_Diff]))

GestationalAge_ANNO_Exome_OverLap <- 
  GestationalAge_ANNO_Exome_OverLap %>% select(
    Chr, Start, Func.refgene.x, Gene.refgene.x, 
    GeneDetail.refgene.x, ExonicFunc.refgene.x, AAChange.refgene.x, dbSNP.x, 
    GeneDetail.refgene.y, ExonicFunc.refgene.y, AAChange.refgene.y)

colnames(GestationalAge_ANNO_Exome_OverLap) <- 
  c("Chr", "Start", "Func.refgene", "Gene.refgene", 
    "GeneDetail.refgene.A1REF_Ex", "ExonicFunc.refgene.A1REF_Ex", 
    "AAChange.refgene.A1REF_Ex", "dbSNP", 
    "GeneDetail.refgene.A2REF_Ex", "ExonicFunc.refgene.A2REF_Ex", 
    "AAChange.refgene.A2REF_Ex")

GestationalAge_ANNO_Exome_Single_A1Ref <- anti_join(GestationalAge_ANNO_A1Ref_Exome, 
                                                    GestationalAge_ANNO_A2Ref_Exome,
                                                    by = c("Chr", "Start"))

GestationalAge_ANNO_Exome_Single_A1Ref <- 
  GestationalAge_ANNO_Exome_Single_A1Ref %>% select(
    Chr, Start, Func.refgene, Gene.refgene, 
    GeneDetail.refgene, ExonicFunc.refgene, AAChange.refgene, dbSNP)

colnames(GestationalAge_ANNO_Exome_Single_A1Ref) <- 
  c("Chr", "Start", "Func.refgene", "Gene.refgene", 
    "GeneDetail.refgene.A1REF_Ex", "ExonicFunc.refgene.A1REF_Ex", 
    "AAChange.refgene.A1REF_Ex", "dbSNP")

GestationalAge_ANNO_Exome_Single_A2Ref <- anti_join(GestationalAge_ANNO_A2Ref_Exome, 
                                                    GestationalAge_ANNO_A1Ref_Exome,
                                                    by = c("Chr", "Start"))

GestationalAge_ANNO_Exome_Single_A2Ref <- 
  GestationalAge_ANNO_Exome_Single_A2Ref %>% select(
    Chr, Start, Func.refgene, Gene.refgene, GeneDetail.refgene, 
    ExonicFunc.refgene, AAChange.refgene, dbSNP)

colnames(GestationalAge_ANNO_Exome_Single_A2Ref) <- 
  c("Chr", "Start", "Func.refgene", "Gene.refgene", 
    "GeneDetail.refgene.A2REF_Ex", "ExonicFunc.refgene.A2REF_Ex", 
    "AAChange.refgene.A2REF_Ex", "dbSNP")

GestationalAge_ANNO_Exome_Single <- 
  full_join(GestationalAge_ANNO_Exome_Single_A1Ref, 
            GestationalAge_ANNO_Exome_Single_A2Ref,
            by = c("Chr", "Start", "Func.refgene", "Gene.refgene", "dbSNP"))

GestationalAge_ANNO_Exome <- 
  full_join(GestationalAge_ANNO_Exome_OverLap,
            GestationalAge_ANNO_Exome_Single,
            by = c("Chr", "Start", "Func.refgene", "Gene.refgene", 
                   "GeneDetail.refgene.A1REF_Ex", "ExonicFunc.refgene.A1REF_Ex", 
                   "AAChange.refgene.A1REF_Ex", "dbSNP", "GeneDetail.refgene.A2REF_Ex", 
                   "ExonicFunc.refgene.A2REF_Ex", "AAChange.refgene.A2REF_Ex"))

## Merge ANNOVAR annotations for GestationalAge_LOCO_MLMA
GestationalAge_ANNO_ExGe <- 
  full_join(GestationalAge_ANNO_Exome,
            GestationalAge_ANNO_Genome,
            by = c("Chr", "Start", "Func.refgene", "Gene.refgene", "dbSNP"))

save(GestationalAge_ANNO_ExGe, 
     file="~/Documents/GestationalAge_ANNO_ExGe.Rda")

#################################################################

## Subset CHR 5 region
Peru_Gest_MAF_LOCO_MLMA_CHR5 <- subset(Peru_Gest_MAF_LOCO_MLMA, 
                                       Peru_Gest_MAF_LOCO_MLMA$Chr == 5)

GestationalAge_LOCO_MLMA_CHR5 <- subset(GestationalAge_LOCO_MLMA, 
                                        GestationalAge_LOCO_MLMA$Chr == 5)

## Merge ANNOVAR annotations for CHR 5 region
Peru_Gest_MAF_LOCO_MLMA_CHR5_ANNO <- left_join(Peru_Gest_MAF_LOCO_MLMA_CHR5, 
                                               Peru_Gest_ANNO_ExGe, 
                                               by= c("Chr", "bp"="Start"))

GestationalAge_LOCO_MLMA_CHR5_ANNO <- left_join(GestationalAge_LOCO_MLMA_CHR5, 
                                                GestationalAge_ANNO_ExGe, 
                                                by= c("Chr", "bp"="Start"))

#################################################################

## Remove NaN p-values (N/A for Peru_Gest)
Peru_Gest_MAF_LOCO_MLMA_CHR5_ANNO$Analysis <- "OriginalGWAS"

GestationalAge_LOCO_MLMA_CHR5_ANNO_NA <-  
  GestationalAge_LOCO_MLMA_CHR5_ANNO[complete.cases(GestationalAge_LOCO_MLMA_CHR5_ANNO$p),]

GestationalAge_LOCO_MLMA_CHR5_ANNO_NA$Analysis <- "Resequencing"
GestationalAge_LOCO_MLMA_CHR5_ANNO_NA$RSID <- NA

GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL <- rbind(Peru_Gest_MAF_LOCO_MLMA_CHR5_ANNO,
                                                 GestationalAge_LOCO_MLMA_CHR5_ANNO_NA)

save(GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL,
     file="~/Documents/GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL.Rda")

#################################################################

## Subset out FGF1 annotations 
GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1 <- 
  GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL[(
    grepl("FGF1", GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL$Gene.refgene) &
      !grepl ("FGF10", GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL$Gene.refgene) &
      !grepl ("FGF18", GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL$Gene.refgene)), ]

table(GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1$Gene.refgene)

## Remove empty columns 
GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1 <- 
  GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1[
    , colSums(is.na(GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1))
    < nrow(GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1)]

save(GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1,
     file="~/Documents/GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1.Rda")

table(GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1$Func.refgene)

## Plot FGF1 annotated region
FGF1_plot <-
  ggplot(data = GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1,
         aes(x=GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1$bp, 
             y=-log10(GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1$p))) +
  scale_color_manual(values=c("#db4dff",    ## downstream
                              "#00cc66",    ## intergenic
                              "#ff1a1a",    ## intronic
                              "#66b3ff",    ## upstream
                              "#204060",    ## UTR3
                              "#808000"     ## UTR5
                              )) +
  geom_point(size=2.5, 
             aes(colour = factor(Func.refgene),
                 shape = GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1$Analysis),
             alpha = 0.6) +
  ggtitle("FGF1 Annotated Region (chromosome 5)") +
  xlab("Physical Position") + 
  ylab("-log10 P-value") + 
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face="bold", size=20),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) 
plot(FGF1_plot)

tiff(file = 
       "~/Documents/GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1.tiff", 
     width = 9000, height = 5000, 
     units = "px", 
     res = 800) 
plot(FGF1_plot)
dev.off()

legend <- get_legend(FGF1_plot)
tiff(file = 
       "~/Documents/GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1_Legend.tiff", 
     width = 12000, height = 1000, 
     units = "px", 
     res = 800)
plot(legend)
dev.off()

#################################################################

## gunzip *.gz

## set working directory
## GTEx_Analysis_v6p_eQTL annotations
## Single-Tissue cis-eQTL Data
## eGene and significant snp-gene associations based on permutations
setwd("~/Documents/GTEx_Analysis_v6p_eQTL/")

## create a list from these files
list.filenames <- list.files(pattern=".v6p.egenes.txt$")

## create an empty list that will serve as a container to receive the incoming files
list.data <- list()

## create a loop to read in your data
## append filename as new column
load(file="~/Documents/list.data.Rda")

#################################################################

## Original GWAS of Gestational Age
Peru_Gest_MAF_LOCO_MLMA <- read.table(
  file="~/Documents/Peru.annot",
  header=TRUE, sep="")

Peru_Gest_MAF_LOCO_MLMA$Chr <- as.character(Peru_Gest_MAF_LOCO_MLMA$Chr)

## Filter out overlap btwn annotations and sequence data 
Peru_Gest_GTEx_Overlap <- list()

for (x in 1:length(list.data)) {
  Peru_Gest_GTEx_Overlap[[x]] <- inner_join(Peru_Gest_MAF_LOCO_MLMA,
                                            list.data[[x]],
                                            by= c("Chr"="chr", "bp"="snp_pos"))
}

Peru_Gest_GTEx_Overlap_DF <- 
  do.call(rbind.data.frame, Peru_Gest_GTEx_Overlap)

Peru_Gest_GTEx_Overlap_DF$Analysis <- c("OriginalGWAS")

## Resequencing analysis of Gestational Age
GestationalAge_LOCO_MLMA <- read.table(
  file="~/Documents/gestational-age-resequencing.loco.mlma",
  header=TRUE, sep="")

GestationalAge_LOCO_MLMA$Chr <- as.character(GestationalAge_LOCO_MLMA$Chr)

## Filter out overlap btwn annotations and sequence data 
GestationalAge_GTEx_Overlap <- list()

for (x in 1:length(list.data)) {
  GestationalAge_GTEx_Overlap[[x]] <- inner_join(GestationalAge_LOCO_MLMA,
                                                 list.data[[x]],
                                                 by= c("Chr"="chr", "bp"="snp_pos"))
}

GestationalAge_GTEx_Overlap_DF <- do.call(rbind.data.frame, GestationalAge_GTEx_Overlap)

GestationalAge_GTEx_Overlap_DF$Analysis <- c("Resequencing")
GestationalAge_GTEx_Overlap_DF$RSID <- NA

## Merge GTEx annotations
GestationalAge_LOCO_MLMA_GTEx_Overlap_FULL <- rbind(Peru_Gest_GTEx_Overlap_DF,
                                                    GestationalAge_GTEx_Overlap_DF)

save(GestationalAge_LOCO_MLMA_GTEx_Overlap_FULL,
     file="~/Documents/GestationalAge_LOCO_MLMA_GTEx_Overlap_FULL.Rda")

## Remove NaN p-values
GestationalAge_LOCO_MLMA_GTEx_Overlap_FULL_NA <-  
  GestationalAge_LOCO_MLMA_GTEx_Overlap_FULL[complete.cases(
    GestationalAge_LOCO_MLMA_GTEx_Overlap_FULL$p),]

## Subset CHR 5 region
GestationalAge_LOCO_MLMA_CHR5_GTEx_FULL <- 
  subset(GestationalAge_LOCO_MLMA_GTEx_Overlap_FULL_NA, 
         GestationalAge_LOCO_MLMA_GTEx_Overlap_FULL_NA$Chr == 5)

GestationalAge_LOCO_MLMA_CHR5_GTEx_FULL$`list.filenames[[i]]` <- 
  as.character(GestationalAge_LOCO_MLMA_CHR5_GTEx_FULL$`list.filenames[[i]]`)

## Subset out FGF1 annotations 
GestationalAge_LOCO_MLMA_CHR5_GTEx_FULL_FGF1 <- 
  GestationalAge_LOCO_MLMA_CHR5_GTEx_FULL[(
    grepl("FGF1", GestationalAge_LOCO_MLMA_CHR5_GTEx_FULL$gene_name) &
      !grepl ("FGF10", GestationalAge_LOCO_MLMA_CHR5_GTEx_FULL$gene_name) &
      !grepl ("FGF18", GestationalAge_LOCO_MLMA_CHR5_GTEx_FULL$gene_name)), ]

table(GestationalAge_LOCO_MLMA_CHR5_GTEx_FULL_FGF1$gene_name)

save(GestationalAge_LOCO_MLMA_CHR5_GTEx_FULL_FGF1,
     file="~/Documents/GestationalAge_LOCO_MLMA_CHR5_GTEx_FULL_FGF1.Rda")

#################################################################

## Merge with wANNOVAR annotations 

GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1$Chr <- 
  as.character(GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1$Chr)

GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1 <- 
  full_join(GestationalAge_LOCO_MLMA_CHR5_GTEx_FULL_FGF1,
            GestationalAge_LOCO_MLMA_CHR5_ANNO_FULL_FGF1,
            by = c("Chr", "SNP", "bp", "A1", "A2", "Freq", 
                   "b", "se", "p", "RSID", "Analysis"))

GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1 <-
  GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1[order(
    GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1$bp, decreasing = TRUE),]

## Merge entries from multiple tissues
## GTEx annotation parameters differ for tissues **

GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique <- 
  subset(GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1,
         (duplicated(GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1$bp) == FALSE))

GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Duplicates <- 
  subset(GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1,
         (duplicated(GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1$bp) == TRUE))

for (k in 1:nrow(GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Duplicates)) {
  
  x = which(GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$bp == 
              GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Duplicates$bp[[k]])
  
  DupName <- 
    GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Duplicates$`list.filenames[[i]]`[[k]]
  CurrName <- 
    GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$`list.filenames[[i]]`[[x]]
  
  GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$`list.filenames[[i]]`[[x]] <- 
    paste(CurrName, DupName, sep=";")
  
  DupAnalysis <- 
    GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Duplicates$Analysis[[k]]
  CurrAnalysis <- 
    GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$Analysis[[x]]
  
  GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$Analysis[[x]] <- 
    paste(CurrAnalysis, DupAnalysis, sep=";")
}

for (k in 1:nrow(GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique)) {
  
  if (is.na(GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$Func.refgene[[k]]) 
      == TRUE) {
    GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$Func.refgene[[k]] <- 
      c("GTEx_unknown")
  }}

## Indicate if present in OriginalGWAS, otherwise is Resequenced
GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$AnalysisV2 <- NA

Original <- which(
  grepl("OriginalGWAS", 
        GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$Analysis) == TRUE)

for (k in 1:length(Original)) {
  
  AnalysisNum <- Original[[k]]
  GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$AnalysisV2[[AnalysisNum]] <- 
    "OriginalGWAS"
}

Resequenced <- which(
  grepl("OriginalGWAS",
        GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$Analysis) == FALSE)

for (k in 1:length(Resequenced)) {
  
  AnalysisNum <- Resequenced[[k]]
  GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$AnalysisV2[[AnalysisNum]] <- 
    "Resequenced"
}

write.csv(GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique,
          file="~/Documents/GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique.csv")

save(GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique,
     file="~/Documents/GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique.Rda")

## Plot FGF1 annotated region

FGF1_plot <-
  ggplot(data = GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique,
         aes(x=GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$bp,
             y=-log10(GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$p))) +
  scale_color_manual(values=c("#db4dff",    ## downstream
                              "#ffb3b3",    ## GTEx_unknown
                              "#00cc66",    ## intergenic
                              "#ff1a1a",    ## intronic
                              "#66b3ff",    ## upstream
                              "#204060",    ## UTR3
                              "#808000"     ## UTR5
  )) +
  geom_point(size=2.5, 
             aes(colour = factor(Func.refgene),
                 shape = GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique$AnalysisV2),
             alpha = 0.6) +
  ggtitle("FGF1 Annotated Region (chromosome 5)") +
  xlab("Physical Position") + 
  ylab("-log10 P-value") + 
  theme(legend.position="bottom",
        plot.title = element_text(hjust = 0.5, face="bold", size=20),
        axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15)) 
plot(FGF1_plot)

tiff(file = 
       "~/Documents/GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique.tiff", 
     width = 9000, height = 5000, 
     units = "px", 
     res = 800) 
plot(FGF1_plot)
dev.off()

legend <- get_legend(FGF1_plot)
tiff(file = 
       "~/Documents/GestationalAge_LOCO_MLMA_CHR5_GTEx_ANNO_FULL_FGF1_Unique_Legend.tiff", 
     width = 15000, height = 1000, 
     units = "px", 
     res = 800)
plot(legend)
dev.off()

#################################################################

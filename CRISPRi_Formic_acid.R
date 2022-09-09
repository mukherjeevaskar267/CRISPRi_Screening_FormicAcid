#Import_metadata

Metadata_CRISPRi <- read.csv("COMPILED_DATA/library_keyfile1536.csv", na.strings = "#N/A", stringsAsFactors = FALSE)
str(Metadata_CRISPRi)


#Import of titration data
#Titration assay is normally performed only with CRISPRi library plate 8

## Extract Metadata for plate 8

Metadata_CRISPRi_p8 <- Metadata_CRISPRi[which(Metadata_CRISPRi$SOURCEPLATEID=="R2877.H.008"), ]
str(Metadata_CRISPRi_p8)

## Import dataset from the RAW_DATA/FORMIC_TITRATION/ folder

### GENERATE THE **ABSOLUTE** DATASET

file.names<- list.files("RAW_DATA/FORMIC_TITRATION/", pattern = ".phenotypes.Absolute", full.names = FALSE)
Path <- vector(mode = "character", length = 0)
FA_titration_df_Y <- data.frame(Metadata_CRISPRi_p8[, 2])
FA_titration_df_GT <- data.frame(Metadata_CRISPRi_p8[, 2])
temp_df<-data.frame()
for(i in 1:length(file.names)){
  Path <- dir("RAW_DATA/FORMIC_TITRATION", pattern = file.names[i], full.names = TRUE)
  temp_df <- read.csv(Path, na.strings = "NoGrowth")
  FA_titration_df_Y <- cbind(FA_titration_df_Y, temp_df[,14])
  FA_titration_df_GT <- cbind(FA_titration_df_GT, temp_df[,15])
}
colnames(FA_titration_df_Y) <- c("gRNA_name", paste0(sub("_phenotypes.*csv", "", file.names), "_Y"))
colnames(FA_titration_df_GT) <- c("gRNA_name", paste0(sub("_phenotypes.*csv", "", file.names), "_GT"))
str(FA_titration_df_Y)
str(FA_titration_df_GT)

#Install packages: ggplot2, reshape, ggridges
#Prepare the data in the format requisite for ggplot2 package using reshape
library(ggplot2)
library(reshape)
library(ggridges)
library(tidyverse)

FA_titration_df_Y_reshape <- reshape(data=FA_titration_df_Y, idvar="gRNA_name",
                                     varying = colnames(FA_titration_df_Y)[2:ncol(FA_titration_df_Y)],
                                     v.name=c("Yield"),
                                     new.row.names = 1:30000,
                                     direction="long",
                                     timevar = "Condition",
                                     times = colnames(FA_titration_df_Y)[2:ncol(FA_titration_df_Y)])

#NOT ESSENTIAL : This is just to reorder the lane in the desired format
FA_titration_df_Y_reshape <- FA_titration_df_Y_reshape %>% mutate(Condition = fct_relevel(Condition,
                                                                                             levels = "FA_0_UL_Y",
                                                                             "FA_50_UL_Y", 
                                                                             "FA_90_UL_Y", 
                                                                             "FA_100_UL_Y", 
                                                                             "FA_110_UL_Y",
                                                                             "FA_120_UL_Y",
                                                                             "FA_130_UL_Y",
                                                                             "FA_140_UL_Y",
                                                                             "FA_150_UL_Y",
                                                                             "FA_150_UL2_Y",
                                                                             "FA_160_UL_Y",
                                                                             "FA_170_UL_Y"))


FA_titration_df_GT_reshape <- reshape(data=FA_titration_df_GT, idvar="gRNA_name",
                                     varying = colnames(FA_titration_df_GT)[2:ncol(FA_titration_df_GT)],
                                     v.name=c("Generation_time"),
                                     new.row.names = 1:30000,
                                     direction="long",
                                     timevar = "Condition",
                                     times = colnames(FA_titration_df_GT)[2:ncol(FA_titration_df_GT)])

#NOT ESSENTIAL : This is just to reorder the lane in the desired format
FA_titration_df_GT_reshape <- FA_titration_df_GT_reshape %>% mutate(Condition = fct_relevel(Condition,
                                                                                          levels = "FA_0_UL_GT",
                                                                                          "FA_50_UL_GT", 
                                                                                          "FA_90_UL_GT", 
                                                                                          "FA_100_UL_GT", 
                                                                                          "FA_110_UL_GT",
                                                                                          "FA_120_UL_GT",
                                                                                          "FA_130_UL_GT",
                                                                                          "FA_140_UL_GT",
                                                                                          "FA_150_UL_GT",
                                                                                          "FA_150_UL2_GT",
                                                                                          "FA_160_UL_GT",
                                                                                          "FA_170_UL_GT"))

## Plot the Ridgeline plots: A nice way to compare the density trace of multiple dataset

pdf("FA_absolute_yield_20210906.pdf")
plt0 <- ggplot(FA_titration_df_Y_reshape, aes(x = Yield, y = Condition, height = stat(density))) + 
  geom_density_ridges2(stat = "binline", bins = 200, scale = 1, draw_baseline = FALSE)+
  theme_ridges()
suppressWarnings(print(plt0))
dev.off()

pdf("FA_absolute_GT_20210906.pdf")
plt1 <- ggplot(FA_titration_df_GT_reshape, aes(x = Generation_time, y = Condition, height = stat(density))) + 
  geom_density_ridges2(stat = "binline", bins = 200, scale = 1, draw_baseline = FALSE)+
  theme_ridges()
suppressWarnings(print(plt1))
dev.off()

### GENERATE THE **NORMALIZED** DATASET

file.names_Norm<- list.files("RAW_DATA/FORMIC_TITRATION/", pattern = ".phenotypes.Nor", full.names = FALSE)
Path <- vector(mode = "character", length = 0)
FA_titration_df_Norm_Y <- data.frame(Metadata_CRISPRi_p8[, 2])
FA_titration_df_Norm_GT <- data.frame(Metadata_CRISPRi_p8[, 2])
temp_df<-data.frame()
for(i in 1:length(file.names)){
  Path <- dir("RAW_DATA/FORMIC_TITRATION", pattern = file.names_Norm[i], full.names = TRUE)
  temp_df <- read.csv(Path, na.strings = "NoGrowth")
  FA_titration_df_Norm_Y <- cbind(FA_titration_df_Norm_Y, temp_df[,4])
  FA_titration_df_Norm_GT <- cbind(FA_titration_df_Norm_GT, temp_df[,5])
}
colnames(FA_titration_df_Norm_Y) <- c("gRNA_name", paste0(sub("_phenotypes.*csv", "", file.names), "_Y"))
colnames(FA_titration_df_Norm_GT) <- c("gRNA_name", paste0(sub("_phenotypes.*csv", "", file.names), "_GT"))
str(FA_titration_df_Norm_Y)
str(FA_titration_df_Norm_GT)

#Install packages: ggplot2, reshape, ggridges
#Prepare the data in the format requisite for ggplot2 package using reshape

FA_titration_df_Norm_Y_reshape <- reshape(data=FA_titration_df_Norm_Y, idvar="gRNA_name",
                                     varying = colnames(FA_titration_df_Norm_Y)[2:ncol(FA_titration_df_Norm_Y)],
                                     v.name=c("Yield"),
                                     new.row.names = 1:30000,
                                     direction="long",
                                     timevar = "Condition",
                                     times = colnames(FA_titration_df_Norm_Y)[2:ncol(FA_titration_df_Norm_Y)])

#NOT ESSENTIAL : This is just to reorder the lane in the desired format
FA_titration_df_Norm_Y_reshape <- FA_titration_df_Norm_Y_reshape %>% mutate(Condition = fct_relevel(Condition,
                                                                                          levels = "FA_0_UL_Y",
                                                                                          "FA_50_UL_Y", 
                                                                                          "FA_90_UL_Y", 
                                                                                          "FA_100_UL_Y", 
                                                                                          "FA_110_UL_Y",
                                                                                          "FA_120_UL_Y",
                                                                                          "FA_130_UL_Y",
                                                                                          "FA_140_UL_Y",
                                                                                          "FA_150_UL_Y",
                                                                                          "FA_150_UL2_Y",
                                                                                          "FA_160_UL_Y",
                                                                                          "FA_170_UL_Y"))

FA_titration_df_Norm_GT_reshape <- reshape(data=FA_titration_df_Norm_GT, idvar="gRNA_name",
                                      varying = colnames(FA_titration_df_Norm_GT)[2:ncol(FA_titration_df_Norm_GT)],
                                      v.name=c("Generation_time"),
                                      new.row.names = 1:30000,
                                      direction="long",
                                      timevar = "Condition",
                                      times = colnames(FA_titration_df_Norm_GT)[2:ncol(FA_titration_df_Norm_GT)])

#NOT ESSENTIAL : This is just to reorder the lane in the desired format
FA_titration_df_Norm_GT_reshape <- FA_titration_df_Norm_GT_reshape %>% mutate(Condition = fct_relevel(Condition,
                                                                                            levels = "FA_0_UL_GT",
                                                                                            "FA_50_UL_GT", 
                                                                                            "FA_90_UL_GT", 
                                                                                            "FA_100_UL_GT", 
                                                                                            "FA_110_UL_GT",
                                                                                            "FA_120_UL_GT",
                                                                                            "FA_130_UL_GT",
                                                                                            "FA_140_UL_GT",
                                                                                            "FA_150_UL_GT",
                                                                                            "FA_150_UL2_GT",
                                                                                            "FA_160_UL_GT",
                                                                                            "FA_170_UL_GT"))

## Plot the Ridgeline plots: A nice way to compare the density trace of multiple dataset
library(ggplot2)
library(reshape)
library(ggridges)

pdf("FA_Norm_yield_20210906.pdf")
plt0_N <- ggplot(FA_titration_df_Norm_Y_reshape, aes(x = Yield, y = Condition, height = stat(density))) + 
  geom_density_ridges2(stat = "binline", bins = 200, scale = 1, draw_baseline = FALSE)+
  theme_ridges()
suppressWarnings(print(plt0_N))
dev.off()

pdf("FA_Norm_GT_20210906.pdf")
plt1_N <- ggplot(FA_titration_df_Norm_GT_reshape, aes(x = Generation_time, y = Condition, height = stat(density))) + 
  geom_density_ridges2(stat = "binline", bins = 200, scale = 1, draw_baseline = FALSE)+
  theme_ridges()
suppressWarnings(print(plt1_N))
dev.off()

###############
#FORMIC ACID SCREENING ROUND1

### GENERATE THE DATASET (IMPORT DATA)

file.names_Abs_R1 <- list.files("RAW_DATA/20211105_FA_ROUND1/", recursive = TRUE, pattern = ".phenotypes.Absolute", full.names = FALSE)
file.names_Norm_R1 <- list.files("RAW_DATA/20211105_FA_ROUND1/", recursive = TRUE, pattern = ".phenotypes.Nor", full.names = FALSE)

Path_basal_Norm <- vector(mode = "character", length = 0)
Path_FA_Norm <- vector(mode = "character", length = 0)
Path_basal_ABS <- vector(mode = "character", length = 0)
Path_FA_ABS <- vector(mode = "character", length = 0)

temp_df_basal_Norm<-data.frame()
temp_df_FA_Norm<-data.frame()
temp_df_basal_ABS<-data.frame()
temp_df_FA_ABS<-data.frame()

temp_df_cbind <- data.frame()
whole_data_CRISPRi_FA_ROUND1 <- data.frame()

file.names_norm_basal <- read.csv("RAW_DATA/20211105_FA_ROUND1/File_names_basal.txt", header = TRUE, sep = "\t")
file.names_norm_FA <- read.csv("RAW_DATA/20211105_FA_ROUND1/File_names_FA.txt", header = TRUE, sep = "\t")
file.names_ABS_basal <- read.csv("RAW_DATA/20211105_FA_ROUND1/File_names_basal_ABS.txt", header = TRUE, sep = "\t")
file.names_ABS_FA <- read.csv("RAW_DATA/20211105_FA_ROUND1/File_names_FA_ABS.txt", header = TRUE, sep = "\t")


for(i in 1:nrow(file.names_norm_basal)){
  Path_basal_Norm <- dir("RAW_DATA/20211105_FA_ROUND1/", recursive = TRUE, pattern = as.character(file.names_norm_basal$FILE_NAME[i]), full.names = TRUE)
  Path_FA_Norm <- dir("RAW_DATA/20211105_FA_ROUND1/", recursive = TRUE, pattern = as.character(file.names_norm_FA$FILE_NAME[i]), full.names = TRUE)
  Path_basal_ABS <- dir("RAW_DATA/20211105_FA_ROUND1/", recursive = TRUE, pattern = as.character(file.names_ABS_basal$FILE_NAME[i]), full.names = TRUE)
  Path_FA_ABS <- dir("RAW_DATA/20211105_FA_ROUND1/", recursive = TRUE, pattern = as.character(file.names_ABS_FA$FILE_NAME[i]), full.names = TRUE)
  
  temp_df_basal_Norm <- read.csv(Path_basal_Norm, na.strings = "NoGrowth")
  temp_df_FA_Norm <- read.csv(Path_FA_Norm, na.strings = "NoGrowth")
  temp_df_basal_ABS <- read.csv(Path_basal_ABS, na.strings = "NoGrowth")
  temp_df_FA_ABS <- read.csv(Path_FA_ABS, na.strings = "NoGrowth")
  
  temp_df_cbind <- cbind(temp_df_basal_ABS[, 14:15], temp_df_FA_ABS[, 14:15], temp_df_basal_Norm[, 4:5], temp_df_FA_Norm[, 4:5])
  whole_data_CRISPRi_FA_ROUND1 <- rbind(whole_data_CRISPRi_FA_ROUND1, temp_df_cbind)
}
R <- rep("1", 36864)
Round_ID <- data.frame(R, stringsAsFactors = FALSE)
whole_data_CRISPRi_FA_ROUND1 <- cbind(Metadata_CRISPRi, Round_ID, whole_data_CRISPRi_FA_ROUND1)

colnames(whole_data_CRISPRi_FA_ROUND1) <- c(colnames(Metadata_CRISPRi), 
                                            "Round_ID", 
                                            "BASAL_Y_ABS",
                                            "BASAL_GT_ABS",
                                            "FA_Y_ABS",
                                            "FA_GT_ABS",
                                            "BASAL_Y_NORM",
                                            "BASAL_GT_NORM",
                                            "FA_Y_NORM",
                                            "FA_GT_NORM")

str(whole_data_CRISPRi_FA_ROUND1)

#### RELATIVE GT AND Y CALCULATION

whole_data_CRISPRi_FA_ROUND1[, 21] <- whole_data_CRISPRi_FA_ROUND1[, 19]-whole_data_CRISPRi_FA_ROUND1[, 17]
whole_data_CRISPRi_FA_ROUND1[, 22] <- whole_data_CRISPRi_FA_ROUND1[, 20]-whole_data_CRISPRi_FA_ROUND1[, 18]
colnames(whole_data_CRISPRi_FA_ROUND1)[21] <- "LPI_Y"
colnames(whole_data_CRISPRi_FA_ROUND1)[22] <- "LPI_GT"

### PLATE WISE BATCH CORRECTION

plate_ID <- as.character(unique(whole_data_CRISPRi_FA_ROUND1$SOURCEPLATEID))

whole_data_CRISPRi_FA_ROUND1_corrected <- whole_data_CRISPRi_FA_ROUND1


for(i in 1:24){
  whole_data_CRISPRi_FA_ROUND1_corrected$BASAL_GT_NORM[which(whole_data_CRISPRi_FA_ROUND1_corrected$SOURCEPLATEID==plate_ID[i])] <- 
    whole_data_CRISPRi_FA_ROUND1$BASAL_GT_NORM[which(whole_data_CRISPRi_FA_ROUND1$SOURCEPLATEID==plate_ID[i])] -
    median(na.omit(whole_data_CRISPRi_FA_ROUND1$BASAL_GT_NORM[which(whole_data_CRISPRi_FA_ROUND1$SOURCEPLATEID==plate_ID[i])]))
  
  whole_data_CRISPRi_FA_ROUND1_corrected$FA_GT_NORM[which(whole_data_CRISPRi_FA_ROUND1_corrected$SOURCEPLATEID==plate_ID[i])] <- 
    whole_data_CRISPRi_FA_ROUND1$FA_GT_NORM[which(whole_data_CRISPRi_FA_ROUND1$SOURCEPLATEID==plate_ID[i])] -
    median(na.omit(whole_data_CRISPRi_FA_ROUND1$FA_GT_NORM[which(whole_data_CRISPRi_FA_ROUND1$SOURCEPLATEID==plate_ID[i])]))
  
  whole_data_CRISPRi_FA_ROUND1_corrected$BASAL_Y_NORM[which(whole_data_CRISPRi_FA_ROUND1_corrected$SOURCEPLATEID==plate_ID[i])] <- 
    whole_data_CRISPRi_FA_ROUND1$BASAL_Y_NORM[which(whole_data_CRISPRi_FA_ROUND1$SOURCEPLATEID==plate_ID[i])] -
    median(na.omit(whole_data_CRISPRi_FA_ROUND1$BASAL_Y_NORM[which(whole_data_CRISPRi_FA_ROUND1$SOURCEPLATEID==plate_ID[i])]))
  
  whole_data_CRISPRi_FA_ROUND1_corrected$FA_Y_NORM[which(whole_data_CRISPRi_FA_ROUND1_corrected$SOURCEPLATEID==plate_ID[i])] <- 
    whole_data_CRISPRi_FA_ROUND1$FA_Y_NORM[which(whole_data_CRISPRi_FA_ROUND1$SOURCEPLATEID==plate_ID[i])] -
    median(na.omit(whole_data_CRISPRi_FA_ROUND1$FA_Y_NORM[which(whole_data_CRISPRi_FA_ROUND1$SOURCEPLATEID==plate_ID[i])]))
  
}

### ANALYSIS

#### RELATIVE GT AND Y CALCULATION BATCH CORRECTED DATASET

whole_data_CRISPRi_FA_ROUND1_corrected[, 21] <- whole_data_CRISPRi_FA_ROUND1_corrected[, 19]-whole_data_CRISPRi_FA_ROUND1_corrected[, 17]
whole_data_CRISPRi_FA_ROUND1_corrected[, 22] <- whole_data_CRISPRi_FA_ROUND1_corrected[, 20]-whole_data_CRISPRi_FA_ROUND1_corrected[, 18]
colnames(whole_data_CRISPRi_FA_ROUND1_corrected)[21] <- "LPI_Y"
colnames(whole_data_CRISPRi_FA_ROUND1_corrected)[22] <- "LPI_GT"

#### REMOVING SPATIAL CONTROL DATA AND MAKING A NEW DATAFRAME

Data_CRISPRi_FA <- whole_data_CRISPRi_FA_ROUND1_corrected[-which(whole_data_CRISPRi_FA_ROUND1_corrected$gRNA_name=="SP_Ctrl_CC23"), ]

#### CREATE A TABLE OF UNIQUE gRNA

unique_gRNA <- data.frame(table(Data_CRISPRi_FA$gRNA_name))

#### ARRANGE THE DATA IN THE DESIRED FORMAT

R1<-vector(mode = "integer", length = 0)
test2<-data.frame()
n<-nrow(unique_gRNA)
for(i in 1:n){
  R1 <- which(Data_CRISPRi_FA$gRNA_name==unique_gRNA$Var1[i])
  test1 <- Data_CRISPRi_FA[R1, ]
  test2[i, c(1:8)]<-test1[1, c(2:4, 6:7, 9:11)]
  test2[i, c(9:11)] <- test1$BASAL_GT_NORM
  test2[i, 12] <- mean(na.omit(test1$BASAL_GT_NORM[1:3]))
  test2[i, 13] <- sd(na.omit(test1$BASAL_GT_NORM[1:3]))
  test2[i, c(14:16)] <- test1$FA_GT_NORM
  test2[i, 17] <- mean(na.omit(test1$FA_GT_NORM[1:3]))
  test2[i, 18] <- sd(na.omit(test1$FA_GT_NORM[1:3]))
  test2[i, c(19:21)] <- test1$LPI_GT
  test2[i, 22] <- mean(na.omit(test1$LPI_GT[1:3]))
  test2[i, 23] <- sd(na.omit(test1$LPI_GT[1:3]))
  test2[i, 24] <- sum(!is.na(test1$BASAL_GT_NORM))
  test2[i, 25] <- sum(!is.na(test1$LPI_GT))
}

column_names <- read.table("COMPILED_DATA/column_names.txt", header = FALSE, sep = "\t", as.is = TRUE)
colnames(test2) <- column_names$V1
Analysis_CRISPRi_FA_R1 <- test2
str(Analysis_CRISPRi_FA_R1)

#### BOX PLOT - MEAN RELATIVE GENERATION TIME (LPI GT)


box_stat_LPI_GT_mean <- boxplot(Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN, cex=0.3)

### Null Hypothesis : µ~StrainX~(All_replicates_LPI_GT)- µ(InterquartileRange_LPI_GT) = 0

box_stat_LPI_GT_mean$stats 

# 25th Percentile = -0.02667832
# 75th Percentile = 0.05920068

#Therefore, extraction of the data points within IQR

Intermediate_50 <- Data_CRISPRi_FA$LPI_GT[which(Data_CRISPRi_FA$LPI_GT >=-0.02667832
                                                   &Data_CRISPRi_FA$LPI_GT<=0.05920068)]

### T-TEST
for(i in 1:nrow(Analysis_CRISPRi_FA_R1)){
  if(Analysis_CRISPRi_FA_R1$n_LPI[i] > 2){
    P.value <- t.test(Intermediate_50, as.numeric(Analysis_CRISPRi_FA_R1[i, 19:21]))
    Analysis_CRISPRi_FA_R1[i, 26] <- P.value$p.value
  } else {
    Analysis_CRISPRi_FA_R1[i, 26] <- NA
  }
}
colnames(Analysis_CRISPRi_FA_R1)[26] <- "P.value"

### BH CORRECTION

Analysis_CRISPRi_FA_R1[which(!is.na(Analysis_CRISPRi_FA_R1$P.value)), 27] <- p.adjust(Analysis_CRISPRi_FA_R1$P.value[which(!is.na(Analysis_CRISPRi_FA_R1$P.value))], 
                                                                              method = "BH", 
                                                                              n = length(Analysis_CRISPRi_FA_R1$P.value[which(!is.na(Analysis_CRISPRi_FA_R1$P.value))]))

colnames(Analysis_CRISPRi_FA_R1)[27] <- "P.adj"

### NUMBER OF SIGNIFICANT STRAINS
length(Analysis_CRISPRi_FA_R1$P.value[which(Analysis_CRISPRi_FA_R1$P.value<=0.1)])
# 2751

length(Analysis_CRISPRi_FA_R1$P.value[which(Analysis_CRISPRi_FA_R1$P.value<=0.01)])
# 568

length(Analysis_CRISPRi_FA_R1$P.adj[which(Analysis_CRISPRi_FA_R1$P.adj<=0.1)])
# 40

# Control strains

length(Analysis_CRISPRi_FA_R1$P.value[which(Analysis_CRISPRi_FA_R1$P.value<=0.1 & 
                                              Analysis_CRISPRi_FA_R1$Control.gRNA==1)])
# 5

length(Analysis_CRISPRi_FA_R1$P.value[which(Analysis_CRISPRi_FA_R1$P.value<=0.01 & 
                                              Analysis_CRISPRi_FA_R1$Control.gRNA==1)])
# 0

## BOXPLOT ONLY CONTROL STRAINS

box_stats_ctrl_gRNA <- boxplot(Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN[Analysis_CRISPRi_FA_R1$Control.gRNA==1])
box_stats_ctrl_gRNA$stats

#            [,1]
#[1,] -0.1069876160
#[2,] -0.0341112175
#[3,]  0.0006948636
#[4,]  0.0541541862
#[5,]  0.1007815673

min(na.omit(Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN[which(Analysis_CRISPRi_FA_R1$Control.gRNA==1)]))
# -0.106988

max(na.omit(Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN[which(Analysis_CRISPRi_FA_R1$Control.gRNA==1)]))
# 0.1007816

#TOLERANT STRAINS

candidate_p_val_0.1_FIT_M3 <- which(Analysis_CRISPRi_FA_R1$P.value< 0.1 & Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN < -0.106988)
length(candidate_p_val_0.1_FIT_M3)
# 131

candidate_p_val_0.01_FIT_M3 <- which(Analysis_CRISPRi_FA_R1$P.value< 0.01 & Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN < -0.106988)
length(candidate_p_val_0.01_FIT_M3)
# 22



#SENSITIVE STRAINS

super_sen_M3 <- which(!is.na(Analysis_CRISPRi_FA_R1$BASAL_GT_NORM_MEAN)
                      &(Analysis_CRISPRi_FA_R1$n_LPI<3)
                      &(
                        is.na(Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN)
                        |(Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN > 0.1007816)
                      ))
length(super_sen_M3)
# 116

candidate_p_val_0.1_SEN_M3 <- which((Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN > 0.1007816 & Analysis_CRISPRi_FA_R1$P.value<= 0.1))   
length(candidate_p_val_0.1_SEN_M3)
# 1053

candidate_p_val_0.01_SEN_M3 <- which((Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN > 0.1007816 & Analysis_CRISPRi_FA_R1$P.value<= 0.01))  
length(candidate_p_val_0.01_SEN_M3)
# 370

#### VIOLIN PLOT

Violin_LPI_Mean <- data.frame()
R <- length(which((Analysis_CRISPRi_FA_R1$Control.gRNA==0)
                  &(!is.na(Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN))
))
Violin_LPI_Mean[1:R, 1] <- Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN[which((Analysis_CRISPRi_FA_R1$Control.gRNA==0)
                                                                     &(!is.na(Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN)))]
Violin_LPI_Mean[1:R, 2] <- "ALL"
R2 <- length(which(Analysis_CRISPRi_FA_R1$Control.gRNA==1))
Violin_LPI_Mean[(R+1):(R+R2), 1] <- Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN[which(Analysis_CRISPRi_FA_R1$Control.gRNA==1)]
Violin_LPI_Mean[(R+1):(R+R2), 2] <- "CONTROL"
colnames(Violin_LPI_Mean)[1:2] <- c("Mean", "Label")

library(ggplot2)
pdf("violin_plot.pdf", height = 7, width = 5)
p_gg <- ggplot(Violin_LPI_Mean, aes(x=Label, y=Mean, fill=Label)) + 
  geom_violin(trim=FALSE) + 
  geom_boxplot(width=0.1, fill="white") +
  labs(title="Violin plot",x="Data Type", y = "LPI_GT") +
  scale_fill_manual(values=c("white", "green"))
p_gg + theme_classic()
dev.off()

#### SCATTER PLOT

pdf("Scatter_plot_formic_p_val0.1.pdf", width = 10, height = 10)
plot(Analysis_CRISPRi_FA_R1$BASAL_GT_NORM_MEAN, Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN, 
     pch = 16, 
     cex = 0.5, 
     col = "black", 
     main = "Selection of sensitive and tolerant strains (P-value<=0.1)", 
     xlab = "Normalized generation time (LSC GT) Basal.condition", 
     ylab = "Relative generation time (LPI GT) in 140mM Formic acid", 
     xlim = c(-0.5, 2.5), 
     ylim = c(-0.5, 2.5),
     yaxt="n",
     xaxt="n",
     cex.lab=1.5)
points(Analysis_CRISPRi_FA_R1$BASAL_GT_NORM_MEAN[candidate_p_val_0.1_FIT_M3], 
       Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN[candidate_p_val_0.1_FIT_M3], 
       pch = 16, 
       cex = 0.5, 
       col = "blue")
points(Analysis_CRISPRi_FA_R1$BASAL_GT_NORM_MEAN[c(super_sen_M3, candidate_p_val_0.1_SEN_M3)], 
       Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN[c(super_sen_M3, candidate_p_val_0.1_SEN_M3)], 
       pch = 16, 
       cex = 0.5, 
       col = "red")
points(Analysis_CRISPRi_FA_R1$BASAL_GT_NORM_MEAN[which(Analysis_CRISPRi_FA_R1$Control.gRNA==1)], 
       Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN[which(Analysis_CRISPRi_FA_R1$Control.gRNA==1)], 
       pch = 16, 
       cex = 0.6, 
       col = "green")
axis(side = 2, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
axis(side = 1, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
abline(h=c(-0.106988, 0.1007816), col="gray", lty=2, lwd=2)
dev.off()

pdf("Scatter_plot_formic_p_val0.01.pdf", width = 10, height = 10)
plot(Analysis_CRISPRi_FA_R1$BASAL_GT_NORM_MEAN, Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN, 
     pch = 16, 
     cex = 0.5, 
     col = "black", 
     main = "Selection of sensitive and tolerant strains (P-value<=0.01)", 
     xlab = "Normalized generation time (LSC GT) Basal.condition", 
     ylab = "Relative generation time (LPI GT) in 140mM Formic acid", 
     xlim = c(-0.5, 2.5), 
     ylim = c(-0.5, 2.5),
     yaxt="n",
     xaxt="n",
     cex.lab=1.5)
points(Analysis_CRISPRi_FA_R1$BASAL_GT_NORM_MEAN[candidate_p_val_0.01_FIT_M3], 
       Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN[candidate_p_val_0.01_FIT_M3], 
       pch = 16, 
       cex = 0.5, 
       col = "blue")
points(Analysis_CRISPRi_FA_R1$BASAL_GT_NORM_MEAN[c(super_sen_M3, candidate_p_val_0.01_SEN_M3)], 
       Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN[c(super_sen_M3, candidate_p_val_0.01_SEN_M3)], 
       pch = 16, 
       cex = 0.5, 
       col = "red")
points(Analysis_CRISPRi_FA_R1$BASAL_GT_NORM_MEAN[which(Analysis_CRISPRi_FA_R1$Control.gRNA==1)], 
       Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN[which(Analysis_CRISPRi_FA_R1$Control.gRNA==1)], 
       pch = 16, 
       cex = 0.6, 
       col = "green")
axis(side = 2, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
axis(side = 1, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
abline(h=c(-0.106988, 0.1007816), col="gray", lty=2, lwd=2)
dev.off()

#LIST FIT GENES pvalue 0.1

Fit_M3_complete_p0.1 <- Analysis_CRISPRi_FA_R1[candidate_p_val_0.1_FIT_M3, ]
Fit_M3_complete_p0.1 <- Fit_M3_complete_p0.1[order(Fit_M3_complete_p0.1$LPI_GT_NORM_MEAN, decreasing = FALSE), ]
str(Fit_M3_complete_p0.1)

whole_Gene_list_Final <- read.csv("COMPILED_DATA/Gene_List_CRISPRi_lib.csv", na.strings = "", stringsAsFactors = FALSE)
rownames(whole_Gene_list_Final) <- whole_Gene_list_Final$LIB_ID


length(candidate_p_val_0.1_FIT_M3)
# 131

Fit_all_M3_p0.1 <- data.frame(sort(table(Analysis_CRISPRi_FA_R1$GENE[candidate_p_val_0.1_FIT_M3]), decreasing = TRUE))
dim(Fit_all_M3_p0.1)
# 122

# THE INFORMATION OF MAX EFFECT (ALSO SIGNIFICANT) SIZE OBTAINED BY TARGETING A gRNA 
# IS ADDED. THAT IS THE EFFECT OBSERVED (MAX OR MIN LPI, MAX EFFECT INDUCING gRNA_name AND SEQUENCE)

test3 <- data.frame()
for(i in 1:nrow(Fit_all_M3_p0.1)){
  
  test3 <- Analysis_CRISPRi_FA_R1[which(Analysis_CRISPRi_FA_R1$P.value< 0.1
                                        & Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN < -0.106988 
                                        & Analysis_CRISPRi_FA_R1$GENE == Fit_all_M3_p0.1$Var1[i]), ]
  test3 <- test3[order(test3$LPI_GT_NORM_MEAN, decreasing = FALSE), ]
  
  Fit_all_M3_p0.1[i, 3] <- test3$LPI_GT_NORM_MEAN[1]
  Fit_all_M3_p0.1[i, 4] <- test3$gRNA_name[1]
  Fit_all_M3_p0.1[i, 5] <- test3$Seq[1]
}

colnames(Fit_all_M3_p0.1)[3:5] <- c("MAX_EFFECT_LPI", "MAX_EFFECT_gRNA", "MAX_EFFECT_gRNA_SEQ")

y <- as.character(Fit_all_M3_p0.1$Var1)
x <- whole_Gene_list_Final[y, ]

Fit_all_M3_description_pval_0.1 <- cbind(Fit_all_M3_p0.1, x[, -1])
str(Fit_all_M3_description_pval_0.1)
nrow(Fit_all_M3_description_pval_0.1)

write.csv(Fit_all_M3_description_pval_0.1, file = "COMPILED_DATA/FIT_STRAINS_R1_p0.1.csv")

# FOR GO ANALYSIS
GO_FIT_p0.1 <- Fit_all_M3_description_pval_0.1$SGD_DB_ID
write.table(GO_FIT_p0.1, 
            file = "GO_ANALYSIS/DATA_GO_FA_R1/GO_FIT_p0.1.txt", 
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

#LIST FIT GENES pvalue 0.01

Fit_M3_complete_p0.01 <- Analysis_CRISPRi_FA_R1[candidate_p_val_0.01_FIT_M3, ]
Fit_M3_complete_p0.01 <- Fit_M3_complete_p0.01[order(Fit_M3_complete_p0.01$LPI_GT_NORM_MEAN, decreasing = FALSE), ]
str(Fit_M3_complete_p0.01)

length(candidate_p_val_0.01_FIT_M3)
# [1] 22

Fit_all_M3_p0.01 <- data.frame(sort(table(Analysis_CRISPRi_FA_R1$GENE[candidate_p_val_0.01_FIT_M3]), decreasing = TRUE))
dim(Fit_all_M3_p0.01)
# [1] 22  5

test3 <- data.frame()
for(i in 1:nrow(Fit_all_M3_p0.01)){
  
  test3 <- Analysis_CRISPRi_FA_R1[which(Analysis_CRISPRi_FA_R1$P.value< 0.01
                                        & Analysis_CRISPRi_FA_R1$LPI_GT_NORM_MEAN < -0.106988 
                                        & Analysis_CRISPRi_FA_R1$GENE == Fit_all_M3_p0.01$Var1[i]), ]
  test3 <- test3[order(test3$LPI_GT_NORM_MEAN, decreasing = FALSE), ]
  
  Fit_all_M3_p0.01[i, 3] <- test3$LPI_GT_NORM_MEAN[1]
  Fit_all_M3_p0.01[i, 4] <- test3$gRNA_name[1]
  Fit_all_M3_p0.01[i, 5] <- test3$Seq[1]
}

colnames(Fit_all_M3_p0.01)[3:5] <- c("MAX_EFFECT_LPI", "MAX_EFFECT_gRNA", "MAX_EFFECT_gRNA_SEQ")

y <- as.character(Fit_all_M3_p0.01$Var1)
x <- whole_Gene_list_Final[y, ]

Fit_all_M3_description_pval_0.01 <- cbind(Fit_all_M3_p0.01, x[, -1])
str(Fit_all_M3_description_pval_0.01)
nrow(Fit_all_M3_description_pval_0.01)

write.csv(Fit_all_M3_description_pval_0.01, file = "COMPILED_DATA/FIT_STRAINS_R1_p0.01.csv")                                

# FOR GO ANALYSIS
GO_FIT_p0.01 <- Fit_all_M3_description_pval_0.01$SGD_DB_ID
write.table(GO_FIT_p0.01, 
            file = "GO_ANALYSIS/DATA_GO_FA_R1/GO_FIT_p0.01.txt", 
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

#LIST SENSITIVE GENES pvalue 0.1
Sen_M3_complete_p0.1 <- Analysis_CRISPRi_FA_R1[c(super_sen_M3, candidate_p_val_0.1_SEN_M3), ]
Sen_M3_complete_p0.1 <- Sen_M3_complete_p0.1[order(Sen_M3_complete_p0.1$LPI_GT_NORM_MEAN, decreasing = TRUE), ]
nrow(Sen_M3_complete_p0.1)

Sen_all_M3_p0.1 <- data.frame(sort(table(Analysis_CRISPRi_FA_R1$GENE[c(super_sen_M3, candidate_p_val_0.1_SEN_M3)]), decreasing = TRUE))
dim(Sen_all_M3_p0.1)
# [1] 718   2

test3 <- data.frame()
for(i in 1:nrow(Sen_all_M3_p0.1)){
  
  test3 <- Sen_M3_complete_p0.1[which(Sen_M3_complete_p0.1$GENE == Sen_all_M3_p0.1$Var1[i]), ]
  test3 <- test3[order(test3$LPI_GT_NORM_MEAN, decreasing = FALSE), ]
  Sen_all_M3_p0.1[i, 3] <- test3$LPI_GT_NORM_MEAN[nrow(test3)]
  Sen_all_M3_p0.1[i, 4] <- test3$gRNA_name[nrow(test3)]
  Sen_all_M3_p0.1[i, 5] <- test3$Seq[nrow(test3)]
}

colnames(Sen_all_M3_p0.1)[3:5] <- c("MAX_EFFECT_LPI", "MAX_EFFECT_gRNA", "MAX_EFFECT_gRNA_SEQ")

y <- as.character(Sen_all_M3_p0.1$Var1)
x <- whole_Gene_list_Final[y, ]

Sen_all_M3_description_pval_0.1 <- cbind(Sen_all_M3_p0.1, x[, -1])
str(Sen_all_M3_description_pval_0.1)
nrow(Sen_all_M3_description_pval_0.1)

write.csv(Sen_all_M3_description_pval_0.1, file = "COMPILED_DATA/Sen_STRAINS_R1_p0.1.csv")                                

# FOR GO ANALYSIS
GO_SEN_p0.1 <- Sen_all_M3_description_pval_0.1$SGD_DB_ID
write.table(GO_SEN_p0.1, 
            file = "GO_ANALYSIS/DATA_GO_FA_R1/GO_SEN_p0.1.txt", 
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

#LIST SENSITIVE GENES pvalue 0.01

Sen_M3_complete_p0.01 <- Analysis_CRISPRi_FA_R1[c(super_sen_M3, candidate_p_val_0.01_SEN_M3), ]
Sen_M3_complete_p0.01 <- Sen_M3_complete_p0.01[order(Sen_M3_complete_p0.01$LPI_GT_NORM_MEAN, decreasing = TRUE), ]
nrow(Sen_M3_complete_p0.01)
# [1] 486

Sen_all_M3_p0.01 <- data.frame(sort(table(Analysis_CRISPRi_FA_R1$GENE[c(super_sen_M3, candidate_p_val_0.01_SEN_M3)]), decreasing = TRUE))
dim(Sen_all_M3_p0.01)
# [1] 362   2

test3 <- data.frame()
for(i in 1:nrow(Sen_all_M3_p0.01)){
  
  test3 <- Sen_M3_complete_p0.01[which(Sen_M3_complete_p0.01$GENE == Sen_all_M3_p0.01$Var1[i]), ]
  test3 <- test3[order(test3$LPI_GT_NORM_MEAN, decreasing = FALSE), ]
  
  Sen_all_M3_p0.01[i, 3] <- test3$LPI_GT_NORM_MEAN[nrow(test3)]
  Sen_all_M3_p0.01[i, 4] <- test3$gRNA_name[nrow(test3)]
  Sen_all_M3_p0.01[i, 5] <- test3$Seq[nrow(test3)]
}

colnames(Sen_all_M3_p0.01)[3:5] <- c("MAX_EFFECT_LPI", "MAX_EFFECT_gRNA", "MAX_EFFECT_gRNA_SEQ")

y <- as.character(Sen_all_M3_p0.01$Var1)
x <- whole_Gene_list_Final[y, ]

Sen_all_M3_description_pval_0.01 <- cbind(Sen_all_M3_p0.01, x[, -1])
str(Sen_all_M3_description_pval_0.01)
nrow(Sen_all_M3_description_pval_0.01)

write.csv(Sen_all_M3_description_pval_0.01, file = "COMPILED_DATA/Sen_STRAINS_R1_p0.01.csv")     

# FOR GO ANALYSIS
GO_SEN_p0.01 <- Sen_all_M3_description_pval_0.01$SGD_DB_ID
write.table(GO_SEN_p0.01, 
            file = "GO_ANALYSIS/DATA_GO_FA_R1/GO_SEN_p0.01.txt", 
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

#BACKGROUND GENE SET
BACKGROUND_GENE_SET <- whole_Gene_list_Final$SGD_DB_ID
write.table(BACKGROUND_GENE_SET, 
            file = "GO_ANALYSIS/DATA_GO_FA_R1/BACKGROUND_GENE_SET.txt", 
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

###############
#FORMIC ACID SCREENING ROUND2

### GENERATE THE DATASET (IMPORT DATA)

file.names_Abs_R2 <- list.files("RAW_DATA/211202_08_220114_p8_3_211125/", recursive = TRUE, pattern = ".phenotypes.Absolute", full.names = FALSE)
file.names_Norm_R2 <- list.files("RAW_DATA/211202_08_220114_p8_3_211125/", recursive = TRUE, pattern = ".phenotypes.Nor", full.names = FALSE)

Path_basal_Norm <- vector(mode = "character", length = 0)
Path_FA_Norm <- vector(mode = "character", length = 0)
Path_basal_ABS <- vector(mode = "character", length = 0)
Path_FA_ABS <- vector(mode = "character", length = 0)

temp_df_basal_Norm<-data.frame()
temp_df_FA_Norm<-data.frame()
temp_df_basal_ABS<-data.frame()
temp_df_FA_ABS<-data.frame()

temp_df_cbind <- data.frame()
whole_data_CRISPRi_FA_ROUND2 <- data.frame()

file.names_norm_basal <- read.csv("RAW_DATA/211202_08_220114_p8_3_211125/File_names_basal.txt", header = TRUE, sep = "\t")
file.names_norm_FA <- read.csv("RAW_DATA/211202_08_220114_p8_3_211125/File_names_FA.txt", header = TRUE, sep = "\t")
file.names_ABS_basal <- read.csv("RAW_DATA/211202_08_220114_p8_3_211125/File_names_basal_ABS.txt", header = TRUE, sep = "\t")
file.names_ABS_FA <- read.csv("RAW_DATA/211202_08_220114_p8_3_211125/File_names_FA_ABS.txt", header = TRUE, sep = "\t")


for(i in 1:nrow(file.names_norm_basal)){
  Path_basal_Norm <- dir("RAW_DATA/211202_08_220114_p8_3_211125/", recursive = TRUE, pattern = as.character(file.names_norm_basal$FILE_NAME[i]), full.names = TRUE)
  Path_FA_Norm <- dir("RAW_DATA/211202_08_220114_p8_3_211125/", recursive = TRUE, pattern = as.character(file.names_norm_FA$FILE_NAME[i]), full.names = TRUE)
  Path_basal_ABS <- dir("RAW_DATA/211202_08_220114_p8_3_211125/", recursive = TRUE, pattern = as.character(file.names_ABS_basal$FILE_NAME[i]), full.names = TRUE)
  Path_FA_ABS <- dir("RAW_DATA/211202_08_220114_p8_3_211125/", recursive = TRUE, pattern = as.character(file.names_ABS_FA$FILE_NAME[i]), full.names = TRUE)
  
  temp_df_basal_Norm <- read.csv(Path_basal_Norm, na.strings = "NoGrowth")
  temp_df_FA_Norm <- read.csv(Path_FA_Norm, na.strings = "NoGrowth")
  temp_df_basal_ABS <- read.csv(Path_basal_ABS, na.strings = "NoGrowth")
  temp_df_FA_ABS <- read.csv(Path_FA_ABS, na.strings = "NoGrowth")
  
  temp_df_cbind <- cbind(temp_df_basal_ABS[, 14:15], temp_df_FA_ABS[, 14:15], temp_df_basal_Norm[, 4:5], temp_df_FA_Norm[, 4:5])
  whole_data_CRISPRi_FA_ROUND2 <- rbind(whole_data_CRISPRi_FA_ROUND2, temp_df_cbind)
}
R <- rep("2", 36864)
Round_ID <- data.frame(R, stringsAsFactors = FALSE)
whole_data_CRISPRi_FA_ROUND2 <- cbind(Metadata_CRISPRi, Round_ID, whole_data_CRISPRi_FA_ROUND2)

colnames(whole_data_CRISPRi_FA_ROUND2) <- c(colnames(Metadata_CRISPRi), 
                                            "Round_ID", 
                                            "BASAL_Y_ABS",
                                            "BASAL_GT_ABS",
                                            "FA_Y_ABS",
                                            "FA_GT_ABS",
                                            "BASAL_Y_NORM",
                                            "BASAL_GT_NORM",
                                            "FA_Y_NORM",
                                            "FA_GT_NORM")

str(whole_data_CRISPRi_FA_ROUND2)

#### RELATIVE GT AND Y CALCULATION

whole_data_CRISPRi_FA_ROUND2[, 21] <- whole_data_CRISPRi_FA_ROUND2[, 19]-whole_data_CRISPRi_FA_ROUND2[, 17]
whole_data_CRISPRi_FA_ROUND2[, 22] <- whole_data_CRISPRi_FA_ROUND2[, 20]-whole_data_CRISPRi_FA_ROUND2[, 18]
colnames(whole_data_CRISPRi_FA_ROUND2)[21] <- "LPI_Y"
colnames(whole_data_CRISPRi_FA_ROUND2)[22] <- "LPI_GT"

### PLATE WISE BATCH CORRECTION

plate_ID <- as.character(unique(whole_data_CRISPRi_FA_ROUND2$SOURCEPLATEID))

whole_data_CRISPRi_FA_ROUND2_corrected <- whole_data_CRISPRi_FA_ROUND2


for(i in 1:24){
  whole_data_CRISPRi_FA_ROUND2_corrected$BASAL_GT_NORM[which(whole_data_CRISPRi_FA_ROUND2_corrected$SOURCEPLATEID==plate_ID[i])] <- 
    whole_data_CRISPRi_FA_ROUND2$BASAL_GT_NORM[which(whole_data_CRISPRi_FA_ROUND2$SOURCEPLATEID==plate_ID[i])] -
    median(na.omit(whole_data_CRISPRi_FA_ROUND2$BASAL_GT_NORM[which(whole_data_CRISPRi_FA_ROUND2$SOURCEPLATEID==plate_ID[i])]))
  
  whole_data_CRISPRi_FA_ROUND2_corrected$FA_GT_NORM[which(whole_data_CRISPRi_FA_ROUND2_corrected$SOURCEPLATEID==plate_ID[i])] <- 
    whole_data_CRISPRi_FA_ROUND2$FA_GT_NORM[which(whole_data_CRISPRi_FA_ROUND2$SOURCEPLATEID==plate_ID[i])] -
    median(na.omit(whole_data_CRISPRi_FA_ROUND2$FA_GT_NORM[which(whole_data_CRISPRi_FA_ROUND2$SOURCEPLATEID==plate_ID[i])]))
  
  whole_data_CRISPRi_FA_ROUND2_corrected$BASAL_Y_NORM[which(whole_data_CRISPRi_FA_ROUND2_corrected$SOURCEPLATEID==plate_ID[i])] <- 
    whole_data_CRISPRi_FA_ROUND2$BASAL_Y_NORM[which(whole_data_CRISPRi_FA_ROUND2$SOURCEPLATEID==plate_ID[i])] -
    median(na.omit(whole_data_CRISPRi_FA_ROUND2$BASAL_Y_NORM[which(whole_data_CRISPRi_FA_ROUND2$SOURCEPLATEID==plate_ID[i])]))
  
  whole_data_CRISPRi_FA_ROUND2_corrected$FA_Y_NORM[which(whole_data_CRISPRi_FA_ROUND2_corrected$SOURCEPLATEID==plate_ID[i])] <- 
    whole_data_CRISPRi_FA_ROUND2$FA_Y_NORM[which(whole_data_CRISPRi_FA_ROUND2$SOURCEPLATEID==plate_ID[i])] -
    median(na.omit(whole_data_CRISPRi_FA_ROUND2$FA_Y_NORM[which(whole_data_CRISPRi_FA_ROUND2$SOURCEPLATEID==plate_ID[i])]))
  
}

### ANALYSIS

#### RELATIVE GT AND Y CALCULATION BATCH CORRECTED DATASET

whole_data_CRISPRi_FA_ROUND2_corrected[, 21] <- whole_data_CRISPRi_FA_ROUND2_corrected[, 19]-whole_data_CRISPRi_FA_ROUND2_corrected[, 17]
whole_data_CRISPRi_FA_ROUND2_corrected[, 22] <- whole_data_CRISPRi_FA_ROUND2_corrected[, 20]-whole_data_CRISPRi_FA_ROUND2_corrected[, 18]
colnames(whole_data_CRISPRi_FA_ROUND2_corrected)[21] <- "LPI_Y"
colnames(whole_data_CRISPRi_FA_ROUND2_corrected)[22] <- "LPI_GT"

#### REMOVING SPATIAL CONTROL DATA AND MAKING A NEW DATAFRAME

Data_CRISPRi_FA2 <- whole_data_CRISPRi_FA_ROUND2_corrected[-which(whole_data_CRISPRi_FA_ROUND2_corrected$gRNA_name=="SP_Ctrl_CC23"), ]


#### ARRANGE THE DATA IN THE DESIRED FORMAT

R1<-vector(mode = "integer", length = 0)
test2<-data.frame()
n<-nrow(unique_gRNA)
for(i in 1:n){
  R1 <- which(Data_CRISPRi_FA2$gRNA_name==unique_gRNA$Var1[i])
  test1 <- Data_CRISPRi_FA2[R1, ]
  test2[i, c(1:8)]<-test1[1, c(2:4, 6:7, 9:11)]
  test2[i, c(9:11)] <- test1$BASAL_GT_NORM
  test2[i, 12] <- mean(na.omit(test1$BASAL_GT_NORM[1:3]))
  test2[i, 13] <- sd(na.omit(test1$BASAL_GT_NORM[1:3]))
  test2[i, c(14:16)] <- test1$FA_GT_NORM
  test2[i, 17] <- mean(na.omit(test1$FA_GT_NORM[1:3]))
  test2[i, 18] <- sd(na.omit(test1$FA_GT_NORM[1:3]))
  test2[i, c(19:21)] <- test1$LPI_GT
  test2[i, 22] <- mean(na.omit(test1$LPI_GT[1:3]))
  test2[i, 23] <- sd(na.omit(test1$LPI_GT[1:3]))
  test2[i, 24] <- sum(!is.na(test1$BASAL_GT_NORM))
  test2[i, 25] <- sum(!is.na(test1$LPI_GT))
}

column_names <- read.table("COMPILED_DATA/column_names.txt", header = FALSE, sep = "\t", as.is = TRUE)
colnames(test2) <- column_names$V1
Analysis_CRISPRi_FA_R2 <- test2
str(Analysis_CRISPRi_FA_R2)

#### BOX PLOT - MEAN RELATIVE GENERATION TIME (LPI GT)


box_stat_LPI_GT_mean2 <- boxplot(Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN, cex=0.3)

### Null Hypothesis : µ~StrainX~(All_replicates_LPI_GT)- µ(InterquartileRange_LPI_GT) = 0

box_stat_LPI_GT_mean2$stats 

# 25th Percentile = -0.02382998
# 75th Percentile = 0.07384211

#Therefore, extraction of the data points within IQR

Intermediate_50_R2 <- Data_CRISPRi_FA2$LPI_GT[which(Data_CRISPRi_FA2$LPI_GT >=-0.02382998
                                                &Data_CRISPRi_FA2$LPI_GT<=0.07384211)]

### T-TEST
for(i in 1:nrow(Analysis_CRISPRi_FA_R2)){
  if(Analysis_CRISPRi_FA_R2$n_LPI[i] > 2){
    P.value <- t.test(Intermediate_50_R2, as.numeric(Analysis_CRISPRi_FA_R2[i, 19:21]))
    Analysis_CRISPRi_FA_R2[i, 26] <- P.value$p.value
  } else {
    Analysis_CRISPRi_FA_R2[i, 26] <- NA
  }
}
colnames(Analysis_CRISPRi_FA_R2)[26] <- "P.value"

### BH CORRECTION

Analysis_CRISPRi_FA_R2[which(!is.na(Analysis_CRISPRi_FA_R2$P.value)), 27] <- p.adjust(Analysis_CRISPRi_FA_R2$P.value[which(!is.na(Analysis_CRISPRi_FA_R2$P.value))], 
                                                                                      method = "BH", 
                                                                                      n = length(Analysis_CRISPRi_FA_R2$P.value[which(!is.na(Analysis_CRISPRi_FA_R2$P.value))]))

colnames(Analysis_CRISPRi_FA_R2)[27] <- "P.adj"

### NUMBER OF SIGNIFICANT STRAINS
length(Analysis_CRISPRi_FA_R2$P.value[which(Analysis_CRISPRi_FA_R2$P.value<=0.1)])
# 2881

length(Analysis_CRISPRi_FA_R2$P.value[which(Analysis_CRISPRi_FA_R2$P.value<=0.01)])
# 638

length(Analysis_CRISPRi_FA_R2$P.adj[which(Analysis_CRISPRi_FA_R2$P.adj<=0.1)])
# 283

# Control strains

length(Analysis_CRISPRi_FA_R2$P.value[which(Analysis_CRISPRi_FA_R2$P.value<=0.1 & 
                                              Analysis_CRISPRi_FA_R2$Control.gRNA==1)])
# 7

length(Analysis_CRISPRi_FA_R2$P.value[which(Analysis_CRISPRi_FA_R2$P.value<=0.01 & 
                                              Analysis_CRISPRi_FA_R2$Control.gRNA==1)])
# 1

## BOXPLOT ONLY CONTROL STRAINS

box_stats_ctrl_gRNA_R2 <- boxplot(Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN[Analysis_CRISPRi_FA_R2$Control.gRNA==1])
box_stats_ctrl_gRNA_R2$stats

#           [,1]
#[1,] -0.140971324
#[2,] -0.052167652
#[3,]  0.004285685
#[4,]  0.080207052
#[5,]  0.168919846

min(na.omit(Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN[which(Analysis_CRISPRi_FA_R2$Control.gRNA==1)]))
# -0.1409713

max(na.omit(Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN[which(Analysis_CRISPRi_FA_R2$Control.gRNA==1)]))
# 0.4300554
# One of the control strain had very poor had growth most likely due to faulty pinning or some other technical reasons. 
# Therefore, this extreme outlier data was ignored and the second last worst performer data (LPI GT mean = 0.168919846) 
# was set as the effect size threshold for sensitive strains. 

#TOLERANT STRAINS

candidate_p_val_0.1_FIT_M3_R2 <- which(Analysis_CRISPRi_FA_R2$P.value< 0.1 & Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN < -0.1409713)
length(candidate_p_val_0.1_FIT_M3_R2)
# 75

candidate_p_val_0.01_FIT_M3_R2 <- which(Analysis_CRISPRi_FA_R2$P.value< 0.01 & Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN < -0.1409713)
length(candidate_p_val_0.01_FIT_M3_R2)
# 11



#SENSITIVE STRAINS

super_sen_M3_R2 <- which(!is.na(Analysis_CRISPRi_FA_R2$BASAL_GT_NORM_MEAN)
                      &(Analysis_CRISPRi_FA_R2$n_LPI<3)
                      &(
                        is.na(Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN)
                        |(Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN > 0.168919846)
                      ))
length(super_sen_M3_R2)
# 91

candidate_p_val_0.1_SEN_M3_R2 <- which((Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN > 0.168919846 & Analysis_CRISPRi_FA_R2$P.value<= 0.1))   
length(candidate_p_val_0.1_SEN_M3_R2)
# 884

candidate_p_val_0.01_SEN_M3_R2 <- which((Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN > 0.168919846 & Analysis_CRISPRi_FA_R2$P.value<= 0.01))  
length(candidate_p_val_0.01_SEN_M3_R2)
# 388


#### VIOLIN PLOT R2

Violin_LPI_Mean_R2 <- data.frame()
R <- length(which((Analysis_CRISPRi_FA_R2$Control.gRNA==0)
                  &(!is.na(Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN))
))
Violin_LPI_Mean_R2[1:R, 1] <- Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN[which((Analysis_CRISPRi_FA_R2$Control.gRNA==0)
                                                                         &(!is.na(Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN)))]
Violin_LPI_Mean_R2[1:R, 2] <- "ALL"
R2 <- length(which(Analysis_CRISPRi_FA_R2$Control.gRNA==1))
Violin_LPI_Mean_R2[(R+1):(R+R2), 1] <- Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN[which(Analysis_CRISPRi_FA_R2$Control.gRNA==1)]
Violin_LPI_Mean_R2[(R+1):(R+R2), 2] <- "CONTROL"
colnames(Violin_LPI_Mean_R2)[1:2] <- c("Mean", "Label")

library(ggplot2)
pdf("violin_plot_R2.pdf", height = 7, width = 5)
p_gg_R2 <- ggplot(Violin_LPI_Mean_R2, aes(x=Label, y=Mean, fill=Label)) + 
  geom_violin(trim=FALSE) + 
  geom_boxplot(width=0.1, fill="white") +
  labs(title="Violin plot",x="Data Type", y = "LPI_GT") +
  scale_fill_manual(values=c("white", "green"))
p_gg_R2 + theme_classic()
dev.off()

#### SCATTER PLOT R2

pdf("Scatter_plot_formic_p_val0.1_R2.pdf", width = 10, height = 10)
plot(Analysis_CRISPRi_FA_R2$BASAL_GT_NORM_MEAN, Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN, 
     pch = 16, 
     cex = 0.5, 
     col = "black", 
     main = "Selection of sensitive and tolerant strains (P-value<=0.1)", 
     xlab = "Normalized generation time (LSC GT) Basal.condition", 
     ylab = "Relative generation time (LPI GT) in 140mM Formic acid", 
     xlim = c(-0.5, 2.5), 
     ylim = c(-0.5, 2.5),
     yaxt="n",
     xaxt="n",
     cex.lab=1.5)
points(Analysis_CRISPRi_FA_R2$BASAL_GT_NORM_MEAN[candidate_p_val_0.1_FIT_M3_R2], 
       Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN[candidate_p_val_0.1_FIT_M3_R2], 
       pch = 16, 
       cex = 0.5, 
       col = "blue")
points(Analysis_CRISPRi_FA_R2$BASAL_GT_NORM_MEAN[c(super_sen_M3_R2, candidate_p_val_0.1_SEN_M3_R2)], 
       Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN[c(super_sen_M3_R2, candidate_p_val_0.1_SEN_M3_R2)], 
       pch = 16, 
       cex = 0.5, 
       col = "red")
points(Analysis_CRISPRi_FA_R2$BASAL_GT_NORM_MEAN[which(Analysis_CRISPRi_FA_R2$Control.gRNA==1)], 
       Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN[which(Analysis_CRISPRi_FA_R2$Control.gRNA==1)], 
       pch = 16, 
       cex = 0.6, 
       col = "green")
axis(side = 2, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
axis(side = 1, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
abline(h=c(-0.1409713, 0.168919846), col="gray", lty=2, lwd=2)
dev.off()

pdf("Scatter_plot_formic_p_val0.01_R2.pdf", width = 10, height = 10)
plot(Analysis_CRISPRi_FA_R2$BASAL_GT_NORM_MEAN, Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN, 
     pch = 16, 
     cex = 0.5, 
     col = "black", 
     main = "Selection of sensitive and tolerant strains (P-value<=0.01)", 
     xlab = "Normalized generation time (LSC GT) Basal.condition", 
     ylab = "Relative generation time (LPI GT) in 140mM Formic acid", 
     xlim = c(-0.5, 2.5), 
     ylim = c(-0.5, 2.5),
     yaxt="n",
     xaxt="n",
     cex.lab=1.5)
points(Analysis_CRISPRi_FA_R2$BASAL_GT_NORM_MEAN[candidate_p_val_0.01_FIT_M3_R2], 
       Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN[candidate_p_val_0.01_FIT_M3_R2], 
       pch = 16, 
       cex = 0.5, 
       col = "blue")
points(Analysis_CRISPRi_FA_R2$BASAL_GT_NORM_MEAN[c(super_sen_M3_R2, candidate_p_val_0.01_SEN_M3_R2)], 
       Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN[c(super_sen_M3_R2, candidate_p_val_0.01_SEN_M3_R2)], 
       pch = 16, 
       cex = 0.5, 
       col = "red")
points(Analysis_CRISPRi_FA_R2$BASAL_GT_NORM_MEAN[which(Analysis_CRISPRi_FA_R2$Control.gRNA==1)], 
       Analysis_CRISPRi_FA_R2$LPI_GT_NORM_MEAN[which(Analysis_CRISPRi_FA_R2$Control.gRNA==1)], 
       pch = 16, 
       cex = 0.6, 
       col = "green")
axis(side = 2, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
axis(side = 1, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
abline(h=c(-0.1409713, 0.168919846), col="gray", lty=2, lwd=2)
dev.off()

#####################################
# COMBINED ANALYSIS FOR ROUND1 AND 2#
#####################################

# Data_CRISPRi_FA  is the round one batch corrected and spatial control eliminated data frame
# Data_CRISPRi_FA2 is the round two batch corrected and spatial control eliminated data frame

#For combined analysis, these two data frames are row bounded. 

Data_CRISPRi_FA_combined <- rbind(Data_CRISPRi_FA, Data_CRISPRi_FA2)

R1<-vector(mode = "integer", length = 0)
R2<-vector(mode = "integer", length = 0)
test2<-data.frame()
n<-nrow(unique_gRNA)
for(i in 1:n){
  R1 <- which(Data_CRISPRi_FA_combined$gRNA_name==unique_gRNA$Var1[i] & Data_CRISPRi_FA_combined$Round_ID=="1")
  R2 <- which(Data_CRISPRi_FA_combined$gRNA_name==unique_gRNA$Var1[i] & Data_CRISPRi_FA_combined$Round_ID=="2")
  test1 <- Data_CRISPRi_FA_combined[c(R1, R2), ]
  test2[i, c(1:8)]<-test1[1, c(2:4, 6:7, 9:11)]
  test2[i, c(9:14)] <- test1$BASAL_GT_NORM
  test2[i, 15] <- mean(na.omit(test1$BASAL_GT_NORM[1:3]))
  test2[i, 16] <- mean(na.omit(test1$BASAL_GT_NORM[4:6]))
  test2[i, 17] <- sd(na.omit(test1$BASAL_GT_NORM[1:3]))
  test2[i, 18] <- sd(na.omit(test1$BASAL_GT_NORM[4:6]))
  test2[i, 19] <- mean(na.omit(test1$BASAL_GT_NORM[1:6]))
  test2[i, 20] <- median(na.omit(test1$BASAL_GT_NORM[1:6]))
  test2[i, 21] <- sd(na.omit(test1$BASAL_GT_NORM[1:6]))
  test2[i, c(22:27)] <- test1$FA_GT_NORM
  test2[i, 28] <- mean(na.omit(test1$FA_GT_NORM[1:3]))
  test2[i, 29] <- mean(na.omit(test1$FA_GT_NORM[4:6]))
  test2[i, 30] <- sd(na.omit(test1$FA_GT_NORM[1:3]))
  test2[i, 31] <- sd(na.omit(test1$FA_GT_NORM[4:6]))
  test2[i, 32] <- mean(na.omit(test1$FA_GT_NORM[1:6]))
  test2[i, 33] <- median(na.omit(test1$FA_GT_NORM[1:6]))
  test2[i, 34] <- sd(na.omit(test1$FA_GT_NORM[1:6]))
  test2[i, c(35:40)] <- test1$LPI_GT
  test2[i, 41] <- mean(na.omit(test1$LPI_GT[1:3]))
  test2[i, 42] <- mean(na.omit(test1$LPI_GT[4:6]))
  test2[i, 43] <- sd(na.omit(test1$LPI_GT[1:3]))
  test2[i, 44] <- sd(na.omit(test1$LPI_GT[4:6]))
  test2[i, 45] <- mean(na.omit(test1$LPI_GT[1:6]))
  test2[i, 46] <- median(na.omit(test1$LPI_GT[1:6]))
  test2[i, 47] <- sd(na.omit(test1$LPI_GT[1:6]))
  test2[i, 48] <- sum(!is.na(test1$BASAL_GT_NORM))
  test2[i, 49] <- sum(!is.na(test1$LPI_GT))
}

column_names <- read.table("COMPILED_DATA/column_names_combined.txt", header = FALSE, sep = "\t", as.is = TRUE)
colnames(test2) <- column_names$V1
Analysis_CRISPRi_FA_Complete <- test2
str(Analysis_CRISPRi_FA_Complete)


#### BOX PLOT - MEAN RELATIVE GENERATION TIME (LPI GT) COMBINED


box_stat_LPI_GT_mean_comb <- boxplot(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN, cex=0.3)

### Null Hypothesis : µ~StrainX~(All_replicates_LPI_GT)- µ(InterquartileRange_LPI_GT) = 0

box_stat_LPI_GT_mean_comb$stats 

# 25th Percentile = -0.01831571
# 75th Percentile = 0.06050947

#Therefore, extraction of the data points within IQR

Intermediate_50_comb <- Data_CRISPRi_FA_combined$LPI_GT[which(Data_CRISPRi_FA_combined$LPI_GT >=-0.01831571
                                                    &Data_CRISPRi_FA_combined$LPI_GT<=0.06050947)]

### T-TEST
for(i in 1:nrow(Analysis_CRISPRi_FA_Complete)){
  if(Analysis_CRISPRi_FA_Complete$n_LPI[i] > 2){
    P.value <- t.test(Intermediate_50_comb, as.numeric(Analysis_CRISPRi_FA_Complete[i, 35:40]))
    Analysis_CRISPRi_FA_Complete[i, 50] <- P.value$p.value
  } else {
    Analysis_CRISPRi_FA_Complete[i, 50] <- NA
  }
}
colnames(Analysis_CRISPRi_FA_Complete)[50] <- "P.value"

### BH CORRECTION

Analysis_CRISPRi_FA_Complete[which(!is.na(Analysis_CRISPRi_FA_Complete$P.value)), 51] <- p.adjust(Analysis_CRISPRi_FA_Complete$P.value[which(!is.na(Analysis_CRISPRi_FA_Complete$P.value))], 
                                                                                      method = "BH", 
                                                                                      n = length(Analysis_CRISPRi_FA_Complete$P.value[which(!is.na(Analysis_CRISPRi_FA_Complete$P.value))]))

colnames(Analysis_CRISPRi_FA_Complete)[51] <- "P.adj"


### NUMBER OF SIGNIFICANT STRAINS
length(Analysis_CRISPRi_FA_Complete$P.value[which(Analysis_CRISPRi_FA_Complete$P.value<=0.1)])
# 3481

length(Analysis_CRISPRi_FA_Complete$P.value[which(Analysis_CRISPRi_FA_Complete$P.value<=0.01)])
# 1166

length(Analysis_CRISPRi_FA_Complete$P.adj[which(Analysis_CRISPRi_FA_Complete$P.adj<=0.1)])
# 1551

# Control strains

length(Analysis_CRISPRi_FA_Complete$P.value[which(Analysis_CRISPRi_FA_Complete$P.value<=0.1 & 
                                              Analysis_CRISPRi_FA_Complete$Control.gRNA==1)])
# 9

length(Analysis_CRISPRi_FA_Complete$P.value[which(Analysis_CRISPRi_FA_Complete$P.value<=0.01 & 
                                              Analysis_CRISPRi_FA_Complete$Control.gRNA==1)])
# 2

length(Analysis_CRISPRi_FA_Complete$P.value[which(Analysis_CRISPRi_FA_Complete$P.adj<=0.1 & 
                                                    Analysis_CRISPRi_FA_Complete$Control.gRNA==1)])

# 2

## P-Value Diagnostics

pdf("P_value_diagnostics_combined.pdf", width = 8, height = 8)
par(mfrow=c(2,2))
hist(Analysis_CRISPRi_FA_Complete$P.value,
     breaks = 20,
     xlab = "P-value", 
     ylab = "Frequency", 
     main = "P-value of all strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 3000),
     cex.lab= 1.5)
hist(Analysis_CRISPRi_FA_Complete$P.value[which(Analysis_CRISPRi_FA_Complete$Control.gRNA==1)],
     breaks = 20,
     xlab = "P-value", 
     ylab = "Frequency", 
     main = "P-value of control strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 10),
     cex.lab= 1.5)
hist(Analysis_CRISPRi_FA_Complete$P.adj,
     breaks = 20,
     xlab = "P.value adjusted", 
     ylab = "Frequency", 
     main = "P.adjusted values of all strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 3000),
     cex.lab= 1.5)
hist(Analysis_CRISPRi_FA_Complete$P.adj[which(Analysis_CRISPRi_FA_Complete$Control.gRNA==1)],
     breaks = 20,
     xlab = "P.value adjusted", 
     ylab = "Frequency", 
     main = "P.adjusted values of control strains", 
     col = "skyblue",
     xlim = c(0, 1),
     ylim = c(0, 10),
     cex.lab= 1.5)
dev.off()

## BOXPLOT ONLY CONTROL STRAINS COMBINED

box_stats_ctrl_gRNA_comb <- boxplot(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[Analysis_CRISPRi_FA_Complete$Control.gRNA==1])
box_stats_ctrl_gRNA_comb$stats

#             [,1]
#[1,] -0.06891051
#[2,] -0.04401365
#[3,]  0.01358809
#[4,]  0.04755488
#[5,]  0.11255980

min(na.omit(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[which(Analysis_CRISPRi_FA_Complete$Control.gRNA==1)]))
# -0.06891051

max(na.omit(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[which(Analysis_CRISPRi_FA_Complete$Control.gRNA==1)]))
# 0.2654185

#TOLERANT STRAINS

candidate_p_val_0.1_FIT_M3_comb <- which(Analysis_CRISPRi_FA_Complete$P.value< 0.1 & Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN < -0.06891051)
length(candidate_p_val_0.1_FIT_M3_comb)
# 352

candidate_p_val_0.01_FIT_M3_comb <- which(Analysis_CRISPRi_FA_Complete$P.value< 0.01 & Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN < -0.06891051)
length(candidate_p_val_0.01_FIT_M3_comb)
# 143

candidate_padj_val_0.1_FIT_M3_comb <- which(Analysis_CRISPRi_FA_Complete$P.adj< 0.1 & Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN < -0.06891051)
length(candidate_padj_val_0.1_FIT_M3_comb)
# 196

candidate_padj_val_0.01_FIT_M3_comb <- which(Analysis_CRISPRi_FA_Complete$P.adj< 0.01 & Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN < -0.06891051)
length(candidate_padj_val_0.01_FIT_M3_comb)
# 7

#SENSITIVE STRAINS

super_sen_M3_comb <- which(!is.na(Analysis_CRISPRi_FA_Complete$BASAL_GT_RND1_2_MEAN)
                         &(Analysis_CRISPRi_FA_Complete$n_LPI<3)
                         &(
                           is.na(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN)
                           |(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN > 0.2654185)
                         ))
length(super_sen_M3_comb)
# 43

candidate_p_val_0.1_SEN_M3_comb <- which((Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN > 0.2654185 & Analysis_CRISPRi_FA_Complete$P.value<= 0.1))   
length(candidate_p_val_0.1_SEN_M3_comb)
# 506

candidate_p_val_0.01_SEN_M3_comb <- which((Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN > 0.2654185 & Analysis_CRISPRi_FA_Complete$P.value<= 0.01))   
length(candidate_p_val_0.01_SEN_M3_comb)
# 348

candidate_padj_val_0.1_SEN_M3_comb <- which((Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN > 0.2654185 & Analysis_CRISPRi_FA_Complete$P.adj<= 0.1))   
length(candidate_padj_val_0.1_SEN_M3_comb)
# 402

candidate_padj_val_0.01_SEN_M3_comb <- which((Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN > 0.2654185 & Analysis_CRISPRi_FA_Complete$P.adj<= 0.01))   
length(candidate_padj_val_0.01_SEN_M3_comb)
# 70

#### VIOLIN PLOT R2

Violin_LPI_Mean_comb <- data.frame()
R <- length(which((Analysis_CRISPRi_FA_Complete$Control.gRNA==0)
                  &(!is.na(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN))
))
Violin_LPI_Mean_comb[1:R, 1] <- Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[which((Analysis_CRISPRi_FA_Complete$Control.gRNA==0)
                                                                            &(!is.na(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN)))]
Violin_LPI_Mean_comb[1:R, 2] <- "ALL"
R2 <- length(which(Analysis_CRISPRi_FA_Complete$Control.gRNA==1))
Violin_LPI_Mean_comb[(R+1):(R+R2), 1] <- Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[which(Analysis_CRISPRi_FA_Complete$Control.gRNA==1)]
Violin_LPI_Mean_comb[(R+1):(R+R2), 2] <- "CONTROL"
colnames(Violin_LPI_Mean_comb)[1:2] <- c("Mean", "Label")

library(ggplot2)
pdf("violin_plot_comb.pdf", height = 7, width = 5)
p_gg_comb <- ggplot(Violin_LPI_Mean_comb, aes(x=Label, y=Mean, fill=Label)) + 
  geom_violin(trim=FALSE) + 
  geom_boxplot(width=0.1, fill="white") +
  labs(title="Violin plot",x="Data Type", y = "LPI_GT") +
  scale_fill_manual(values=c("white", "green"))
p_gg_comb + theme_classic()
dev.off()


pdf("Scatter_plot_formic_padj_0.1_comb.pdf", width = 10, height = 10)
plot(Analysis_CRISPRi_FA_Complete$BASAL_GT_RND1_2_MEAN, Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN, 
     pch = 16, 
     cex = 0.5, 
     col = "black", 
     main = "Selection of sensitive and tolerant strains (P.adj<=0.1)", 
     xlab = "Normalized generation time (LSC GT) Basal.condition", 
     ylab = "Relative generation time (LPI GT) in 140mM Formic acid", 
     xlim = c(-0.5, 2.5), 
     ylim = c(-0.5, 2.5),
     yaxt="n",
     xaxt="n",
     cex.lab=1.5)
points(Analysis_CRISPRi_FA_Complete$BASAL_GT_RND1_2_MEAN[candidate_padj_val_0.1_FIT_M3_comb], 
       Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[candidate_padj_val_0.1_FIT_M3_comb], 
       pch = 16, 
       cex = 0.5, 
       col = "blue")
points(Analysis_CRISPRi_FA_Complete$BASAL_GT_RND1_2_MEAN[c(super_sen_M3_comb, candidate_padj_val_0.1_SEN_M3_comb)], 
       Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[c(super_sen_M3_comb, candidate_padj_val_0.1_SEN_M3_comb)], 
       pch = 16, 
       cex = 0.5, 
       col = "red")
points(Analysis_CRISPRi_FA_Complete$BASAL_GT_RND1_2_MEAN[which(Analysis_CRISPRi_FA_Complete$Control.gRNA==1)], 
       Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[which(Analysis_CRISPRi_FA_Complete$Control.gRNA==1)], 
       pch = 16, 
       cex = 0.6, 
       col = "green")
axis(side = 2, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
axis(side = 1, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
abline(h=c(-0.06891051, 0.2654185), col="gray", lty=2, lwd=2)
dev.off()

#####

#LIST FIT GENES p.adj value 0.1

Fit_M3_complete_padj0.1_comb <- Analysis_CRISPRi_FA_Complete[candidate_padj_val_0.1_FIT_M3_comb, ]
Fit_M3_complete_padj0.1_comb <- Fit_M3_complete_padj0.1_comb[order(Fit_M3_complete_padj0.1_comb$LPI_GT_RND1_2_MEAN, decreasing = FALSE), ]
str(Fit_M3_complete_padj0.1_comb)

dim(Fit_M3_complete_padj0.1_comb)
# 196

write.csv(Fit_M3_complete_padj0.1_comb, file = "COMPILED_DATA/FIT_STRAINS_COMB_padj0.1.csv")

Fit_all_M3_padj0.1_comb <- data.frame(sort(table(Analysis_CRISPRi_FA_Complete$GENE[candidate_padj_val_0.1_FIT_M3_comb]), decreasing = TRUE))
dim(Fit_all_M3_padj0.1_comb)
# 179

# THE INFORMATION OF MAX EFFECT (ALSO SIGNIFICANT) SIZE OBTAINED BY TARGETING A gRNA 
# IS ADDED. THAT IS THE EFFECT OBSERVED (MAX OR MIN LPI, MAX EFFECT INDUCING gRNA_name AND SEQUENCE)

test3 <- data.frame()
for(i in 1:nrow(Fit_all_M3_padj0.1_comb)){
  
  test3 <- Analysis_CRISPRi_FA_Complete[which(Analysis_CRISPRi_FA_Complete$P.adj< 0.1
                                        & Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN < -0.06891051 
                                        & Analysis_CRISPRi_FA_Complete$GENE == Fit_all_M3_padj0.1_comb$Var1[i]), ]
  test3 <- test3[order(test3$LPI_GT_RND1_2_MEAN, decreasing = FALSE), ]
  
  Fit_all_M3_padj0.1_comb[i, 3] <- test3$LPI_GT_RND1_2_MEAN[1]
  Fit_all_M3_padj0.1_comb[i, 4] <- test3$gRNA_name[1]
  Fit_all_M3_padj0.1_comb[i, 5] <- test3$Seq[1]
}

colnames(Fit_all_M3_padj0.1_comb)[3:5] <- c("MAX_EFFECT_LPI", "MAX_EFFECT_gRNA", "MAX_EFFECT_gRNA_SEQ")

y <- as.character(Fit_all_M3_padj0.1_comb$Var1)
x <- whole_Gene_list_Final[y, ]

Fit_all_M3_description_padj_val_0.1_comb <- cbind(Fit_all_M3_padj0.1_comb, x[, -1])
str(Fit_all_M3_description_padj_val_0.1_comb)
nrow(Fit_all_M3_description_padj_val_0.1_comb)

write.csv(Fit_all_M3_description_padj_val_0.1_comb, file = "COMPILED_DATA/FIT_GENES_COMB_padj0.1.csv")

# FOR GO ANALYSIS
GO_FIT_padj0.1_comb <- Fit_all_M3_description_padj_val_0.1_comb$SGD_DB_ID
write.table(GO_FIT_padj0.1_comb, 
            file = "GO_ANALYSIS/DATA_GO_FA_COMB/GO_FIT_padj0.1_com.txt", 
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)


#LIST SENSITIVE GENES padj 0.1
Sen_M3_complete_padj0.1_comb <- Analysis_CRISPRi_FA_Complete[c(super_sen_M3_comb, candidate_padj_val_0.1_SEN_M3_comb), ]
Sen_M3_complete_padj0.1_comb <- Sen_M3_complete_padj0.1_comb[order(Sen_M3_complete_padj0.1_comb$LPI_GT_RND1_2_MEAN, decreasing = TRUE), ]

dim(Sen_M3_complete_padj0.1_comb)

write.csv(Sen_M3_complete_padj0.1_comb, file = "COMPILED_DATA/SEN_STRAINS_COMB_padj0.1.csv")

Sen_all_M3_padj0.1_comb <- data.frame(sort(table(Analysis_CRISPRi_FA_Complete$GENE[c(super_sen_M3_comb, candidate_padj_val_0.1_SEN_M3_comb)]), decreasing = TRUE))
dim(Sen_all_M3_padj0.1_comb)
# [1] 319   2

test3 <- data.frame()
for(i in 1:nrow(Sen_all_M3_padj0.1_comb)){
  
  test3 <- Sen_M3_complete_padj0.1_comb[which(Sen_M3_complete_padj0.1_comb$GENE == Sen_all_M3_padj0.1_comb$Var1[i]), ]
  test3 <- test3[order(test3$LPI_GT_RND1_2_MEAN, decreasing = FALSE), ]
  Sen_all_M3_padj0.1_comb[i, 3] <- test3$LPI_GT_RND1_2_MEAN[nrow(test3)]
  Sen_all_M3_padj0.1_comb[i, 4] <- test3$gRNA_name[nrow(test3)]
  Sen_all_M3_padj0.1_comb[i, 5] <- test3$Seq[nrow(test3)]
}

colnames(Sen_all_M3_padj0.1_comb)[3:5] <- c("MAX_EFFECT_LPI", "MAX_EFFECT_gRNA", "MAX_EFFECT_gRNA_SEQ")

y <- as.character(Sen_all_M3_padj0.1_comb$Var1)
x <- whole_Gene_list_Final[y, ]

Sen_all_M3_description_padj_val_0.1_comb <- cbind(Sen_all_M3_padj0.1_comb, x[, -1])
str(Sen_all_M3_description_padj_val_0.1_comb)
nrow(Sen_all_M3_description_padj_val_0.1_comb)

write.csv(Sen_all_M3_description_padj_val_0.1_comb, file = "COMPILED_DATA/SEN_GENES_COMB_padj0.1.csv")                                

# FOR GO ANALYSIS
GO_SEN_padj0.1_comb <- Sen_all_M3_description_padj_val_0.1_comb$SGD_DB_ID
write.table(GO_SEN_padj0.1_comb, 
            file = "GO_ANALYSIS/DATA_GO_FA_COMB/GO_SEN_padj0.1_com.txt", 
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)


#BACKGROUND GENE SET
BACKGROUND_GENE_SET <- whole_Gene_list_Final$SGD_DB_ID
write.table(BACKGROUND_GENE_SET, 
            file = "GO_ANALYSIS/DATA_GO_FA_COMB/BACKGROUND_GENE_SET.txt", 
            quote = FALSE, 
            col.names = FALSE, 
            row.names = FALSE)

#############################
#CORREALTION BETWEEN TWO ROUNDS
#############################

#LPI data correlation between two rounds

pdf("Scatter_plot_CORREALTION_RND1vsRND2.pdf", width = 10, height = 10)
plot(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_MEAN, Analysis_CRISPRi_FA_Complete$LPI_GT_RND2_MEAN, 
     pch = 16, 
     cex = 0.5, 
     col = "gray", 
     main = "Correlation between mean relative generation time (LPI GT) of Round 1 and 2", 
     xlab = "Relative generation time (LPI GT) in ROUND1", 
     ylab = "Relative generation time (LPI GT) in ROUND2", 
     xlim = c(-0.5, 2.5), 
     ylim = c(-0.5, 2.5),
     yaxt="n",
     xaxt="n",
     cex.lab=1.5)
points(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_MEAN[candidate_padj_val_0.1_FIT_M3_comb], 
       Analysis_CRISPRi_FA_Complete$LPI_GT_RND2_MEAN[candidate_padj_val_0.1_FIT_M3_comb], 
       pch = 16, 
       cex = 0.5, 
       col = "blue")
points(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_MEAN[c(super_sen_M3_comb, candidate_padj_val_0.1_SEN_M3_comb)], 
       Analysis_CRISPRi_FA_Complete$LPI_GT_RND2_MEAN[c(super_sen_M3_comb, candidate_padj_val_0.1_SEN_M3_comb)], 
       pch = 16, 
       cex = 0.5, 
       col = "red")
points(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_MEAN[which(AA_Analysis_data$Control.gRNA==1)], 
       Analysis_CRISPRi_FA_Complete$LPI_GT_RND2_MEAN[which(AA_Analysis_data$Control.gRNA==1)], 
       pch = 16, 
       cex = 0.6, 
       col = "green")
axis(side = 2, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
axis(side = 1, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
stats_LPI_GT_Mean_RND1vsRND2_FA <- lm(LPI_GT_RND2_MEAN ~ LPI_GT_RND1_MEAN, 
                                      data = Analysis_CRISPRi_FA_Complete)
stats_LPI_GT_Mean_RND1vsRND2_FA_selected <- lm(LPI_GT_RND2_MEAN[c(candidate_padj_val_0.1_FIT_M3_comb, 
                                                                    super_sen_M3_comb, 
                                                                    candidate_padj_val_0.1_SEN_M3_comb)] ~ 
                                                 LPI_GT_RND1_MEAN[c(candidate_padj_val_0.1_FIT_M3_comb, 
                                                                    super_sen_M3_comb, 
                                                                    candidate_padj_val_0.1_SEN_M3_comb)], 
                                               data = Analysis_CRISPRi_FA_Complete)
abline(stats_LPI_GT_Mean_RND1vsRND2_FA, lty=2, lwd=2, col="darkgray")
abline(stats_LPI_GT_Mean_RND1vsRND2_FA_selected, lty=2, lwd=2, col="black")
dev.off()

summary(stats_LPI_GT_Mean_RND1vsRND2_FA)
summary(stats_LPI_GT_Mean_RND1vsRND2_FA_selected)


cor(Analysis_CRISPRi_FA_Complete$LPI_GT_RND2_MEAN,
    Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_MEAN, 
    method = "pearson", 
    use = "complete.obs")
#R = 0.5645077

cor(Analysis_CRISPRi_FA_Complete$LPI_GT_RND2_MEAN[c(candidate_padj_val_0.1_FIT_M3_comb, 
                                                    super_sen_M3_comb, 
                                                    candidate_padj_val_0.1_SEN_M3_comb)],
    Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_MEAN[c(candidate_padj_val_0.1_FIT_M3_comb, 
                                                    super_sen_M3_comb, 
                                                    candidate_padj_val_0.1_SEN_M3_comb)], 
    method = "pearson", 
    use = "complete.obs")
#R= 0.7235362

#############################
#CORRELATION WITH ACETIC ACID
#############################

AA_Analysis_data <- read.csv("COMPILED_DATA/ACETIC_ACID_DATA.csv", na.strings = "#N/A", stringsAsFactors = FALSE)

str(AA_Analysis_data)

#EXTRACT INDEX OF TOLERANT STRAINS AA
candidate_padj_0.1_FIT_M3_AA <- which((AA_Analysis_data$LPI_GT_Mean_all < -0.03680838 & AA_Analysis_data$P.adjusted_M3<= 0.1))
length(candidate_padj_0.1_FIT_M3_AA)

#EXTRACT INDEX OF SENSITIVE STRAINS AA
super_sen_AA <- which(!is.na(AA_Analysis_data$CTRL_GT_Mean_all)
                                       &(AA_Analysis_data$n_LPI<3)
                                       &(
                                         is.na(AA_Analysis_data$LPI_GT_Mean_all)
                                         |(AA_Analysis_data$LPI_GT_Mean_all> 0.165662)))
length(super_sen_AA)

candidate_padj_0.1_SEN_M3_AA <- which((AA_Analysis_data$LPI_GT_Mean_all > 0.165662 & AA_Analysis_data$P.adjusted_M3<= 0.1))
length(candidate_padj_0.1_SEN_M3_AA)

#MAKING SCATTER PLOT

# AA significant strains performance (LPI MEAN) comparison with FA vs AA
pdf("Scatter_plot_CORREALTION_AA_SIGN_STRAINS.pdf", width = 10, height = 10)
plot(AA_Analysis_data$LPI_GT_Mean_all, Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN, 
     pch = 16, 
     cex = 0.5, 
     col = "gray", 
     main = "Selection of AA sensitive and tolerant strains (P.adj<=0.1)", 
     xlab = "Relative generation time (LPI GT) in 150mM Acetic acid", 
     ylab = "Relative generation time (LPI GT) in 140mM Formic acid", 
     xlim = c(-0.5, 2.5), 
     ylim = c(-0.5, 2.5),
     yaxt="n",
     xaxt="n",
     cex.lab=1.5)
points(AA_Analysis_data$LPI_GT_Mean_all[candidate_padj_0.1_FIT_M3_AA], 
       Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[candidate_padj_0.1_FIT_M3_AA], 
       pch = 16, 
       cex = 0.5, 
       col = "blue")
points(AA_Analysis_data$LPI_GT_Mean_all[c(super_sen_AA, candidate_padj_0.1_SEN_M3_AA)], 
       Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[c(super_sen_AA, candidate_padj_0.1_SEN_M3_AA)], 
       pch = 16, 
       cex = 0.5, 
       col = "red")
points(AA_Analysis_data$LPI_GT_Mean_all[which(AA_Analysis_data$Control.gRNA==1)], 
       Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[which(AA_Analysis_data$Control.gRNA==1)], 
       pch = 16, 
       cex = 0.6, 
       col = "green")
axis(side = 2, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
axis(side = 1, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
abline(h=c(-0.06891051, 0.2654185), col="orange", lty=2, lwd=2)
abline(v=c(-0.03680838, 0.165662), col="magenta", lty=2, lwd=2)
stats_LPI_GT_Mean_FAvsAA_AA <- lm(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN ~ AA_Analysis_data$LPI_GT_Mean_all)
stats_LPI_GT_Mean_FAvsAA_AA_selected <- lm(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[c(candidate_padj_0.1_FIT_M3_AA, 
                                                                                             candidate_padj_0.1_SEN_M3_AA, 
                                                                                             super_sen_AA)] ~ 
                                             AA_Analysis_data$LPI_GT_Mean_all[c(candidate_padj_0.1_FIT_M3_AA, 
                                                                                candidate_padj_0.1_SEN_M3_AA, 
                                                                                super_sen_AA)])
abline(stats_LPI_GT_Mean_FAvsAA_AA, lty=2, lwd=2, col="darkgray")
abline(stats_LPI_GT_Mean_FAvsAA_AA_selected, lty=2, lwd=2, col="black")
dev.off()

summary(stats_LPI_GT_Mean_FAvsAA_AA)
summary(stats_LPI_GT_Mean_FAvsAA_AA_selected)


cor(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN,
    AA_Analysis_data$LPI_GT_Mean_all, 
    method = "pearson", 
    use = "complete.obs")
#R = 0.5917216

cor(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[c(candidate_padj_0.1_FIT_M3_AA, 
                                                      candidate_padj_0.1_SEN_M3_AA, 
                                                      super_sen_AA)], 
    AA_Analysis_data$LPI_GT_Mean_all[c(candidate_padj_0.1_FIT_M3_AA, 
                                       candidate_padj_0.1_SEN_M3_AA, 
                                       super_sen_AA)],  
    method = "pearson", 
    use = "complete.obs")
#R = 0.6947635

# FA significant strains performance (LPI MEAN) comparison with FA vs AA
pdf("Scatter_plot_CORREALTION_FA_SIGN_STRAINS.pdf", width = 10, height = 10)
plot(AA_Analysis_data$LPI_GT_Mean_all, Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN, 
     pch = 16, 
     cex = 0.5, 
     col = "gray", 
     main = "Selection of FA sensitive and tolerant strains (P.adj<=0.1)", 
     xlab = "Relative generation time (LPI GT) in 150mM Acetic acid", 
     ylab = "Relative generation time (LPI GT) in 140mM Formic acid", 
     xlim = c(-0.5, 2.5), 
     ylim = c(-0.5, 2.5),
     yaxt="n",
     xaxt="n",
     cex.lab=1.5)
points(AA_Analysis_data$LPI_GT_Mean_all[candidate_padj_val_0.1_FIT_M3_comb], 
       Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[candidate_padj_val_0.1_FIT_M3_comb], 
       pch = 16, 
       cex = 0.5, 
       col = "blue")
points(AA_Analysis_data$LPI_GT_Mean_all[c(super_sen_M3_comb, candidate_padj_val_0.1_SEN_M3_comb)], 
       Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[c(super_sen_M3_comb, candidate_padj_val_0.1_SEN_M3_comb)], 
       pch = 16, 
       cex = 0.5, 
       col = "red")
points(AA_Analysis_data$LPI_GT_Mean_all[which(Analysis_CRISPRi_FA_Complete$Control.gRNA==1)], 
       Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[which(Analysis_CRISPRi_FA_Complete$Control.gRNA==1)], 
       pch = 16, 
       cex = 0.6, 
       col = "green")
axis(side = 2, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
axis(side = 1, 
     at = c(-0.5, 0, 0.5, 1, 1.5, 2.5),
     cex.axis = 1.2,
     labels = c("-0.5", "0", "0.5", "1", "1.5", "2.5"), 
     tick = 0.05)
abline(h=c(-0.06891051, 0.2654185), col="orange", lty=2, lwd=2)
abline(v=c(-0.03680838, 0.165662), col="magenta", lty=2, lwd=2)
stats_LPI_GT_Mean_FAvsAA_AA <- lm(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN ~ AA_Analysis_data$LPI_GT_Mean_all)
stats_LPI_GT_Mean_FAvsAA_FA_selected <- lm(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[c(candidate_padj_val_0.1_FIT_M3_comb, 
                                                                                             super_sen_M3_comb, 
                                                                                             candidate_padj_val_0.1_SEN_M3_comb)] ~ 
                                             AA_Analysis_data$LPI_GT_Mean_all[c(candidate_padj_val_0.1_FIT_M3_comb, 
                                                                                super_sen_M3_comb, 
                                                                                candidate_padj_val_0.1_SEN_M3_comb)])
abline(stats_LPI_GT_Mean_FAvsAA_AA, lty=2, lwd=2, col="darkgray")
abline(stats_LPI_GT_Mean_FAvsAA_FA_selected, lty=2, lwd=2, col="black")
dev.off()

summary(stats_LPI_GT_Mean_FAvsAA_AA)
summary(stats_LPI_GT_Mean_FAvsAA_FA_selected)

cor(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN,
    AA_Analysis_data$LPI_GT_Mean_all, 
    method = "pearson", 
    use = "complete.obs")
#R = 0.5917216

cor(Analysis_CRISPRi_FA_Complete$LPI_GT_RND1_2_MEAN[c(candidate_padj_val_0.1_FIT_M3_comb, 
                                                      super_sen_M3_comb, 
                                                      candidate_padj_val_0.1_SEN_M3_comb)], 
    AA_Analysis_data$LPI_GT_Mean_all[c(candidate_padj_val_0.1_FIT_M3_comb, 
                                       super_sen_M3_comb, 
                                       candidate_padj_val_0.1_SEN_M3_comb)],  
    method = "pearson", 
    use = "complete.obs")
#R = 0.6341759

#SET ANALYSIS

intersect_FAvsAA <- Analysis_CRISPRi_FA_Complete[intersect(candidate_padj_0.1_FIT_M3_AA, candidate_padj_val_0.1_FIT_M3_comb), ]
write.csv(intersect_FAvsAA, file = "COMPILED_DATA/FIT_STRAINS_INTERSECT_FAvsAA.csv")


TOP_FA_str_wo_intersect_str <- Analysis_CRISPRi_FA_Complete[setdiff(candidate_padj_val_0.1_FIT_M3_comb, intersect(candidate_padj_0.1_FIT_M3_AA, candidate_padj_val_0.1_FIT_M3_comb)), ]

str(TOP_FA_str_wo_intersect_str)
write.csv(TOP_FA_str_wo_intersect_str, file = "COMPILED_DATA/FIT_STRAINS_SETDIFF_FAvsAA.csv")


intersect_genes_FAvsAA <- intersect(Analysis_CRISPRi_FA_Complete$GENE[candidate_padj_val_0.1_FIT_M3_comb], 
                                    AA_Analysis_data$GENE[candidate_padj_0.1_FIT_M3_AA])



intersect_SEN_FAvsAA <- Analysis_CRISPRi_FA_Complete[intersect(c(candidate_padj_0.1_SEN_M3_AA, super_sen_AA), c(super_sen_M3_comb, candidate_padj_val_0.1_SEN_M3_comb)), ]
write.csv(intersect_SEN_FAvsAA, file = "COMPILED_DATA/SEN_STRAINS_INTERSECT_FAvsAA.csv")

#CONTROL STRAINS
control_dataset <- Analysis_CRISPRi_FA_Complete[which(Analysis_CRISPRi_FA_Complete$Control.gRNA==1), ]
write.csv(control_dataset, file = "COMPILED_DATA/CONTROL_STRAINS.csv")

#PRINT DATA FOR IBAI
write.csv(Analysis_CRISPRi_FA_Complete, file = "COMPILED_DATA/Analysis_CRISPRi_FA_Complete.csv")





# Script to compute the proteomics data to be used, from raw files obtained from contact at Gothenburg University

# Load data
library(readr)
data_anaerobic <- read_delim("TMT_anaerobic_cultures.csv", locale = locale(decimal_mark = ","), 
                             ";", escape_double = FALSE, trim_ws = TRUE)
data_aerobic <- read_delim("TMT_aerobic_cultures.csv", locale = locale(decimal_mark = ","), 
                           ";", escape_double = FALSE, trim_ws = TRUE)
iBAQ <- read_delim("ref2824_iBAQ_unique_peptides.csv", locale = locale(decimal_mark = ","), 
                   ";", escape_double = FALSE, trim_ws = TRUE)

# load uniprot S288c (needed later to re-map the Protein column)
S288c <- read_csv("20180314 Uniprot S288c reviewed & fixed.csv")
S288c <- S288c[,c(1,7)]

colnames(data_anaerobic) <- c('Accession','Description','Protein_anaerobic','iBAQ_anaerobic_not_the_reference_pool',
                              'Min_1_glc_anaerobic','Min_2_glc_anaerobic','Min_3_glc_anaerobic',
                              'Rich_1_glc_anaerobic','Rich_2_glc_anaerobic','Rich_3_glc_anaerobic',
                              'Seq_coverage','Ident_peptides','PSMs','Ident_uniq_peptides','MW_kDa_anaerobic')
colnames(data_aerobic) <- c('Accession','Description','Seq_coverage','Ident_peptides','PSMs',
                            'Ident_uniq_peptides','MW_kDa_aerobic',
                            'Rich_1_glc_aerobic','Rich_2_glc_aerobic','Rich_3_glc_aerobic',
                            'Min_1_glc_aerobic','Min_2_glc_aerobic','Min_3_glc_aerobic')
colnames(iBAQ) <- c('Accession','fmol_ug_in_ref','Gene','Protein')
colnames(S288c) <- c("Accession", "Protein")

# Change one name that is imported as a date incorrectly
data_anaerobic$Protein_anaerobic <- as.character(data_anaerobic$Protein_anaerobic)
data_anaerobic[data_anaerobic=='Oct-01'] <- as.character('OCT1')

iBAQ$Protein <- as.character(iBAQ$Protein)
iBAQ[iBAQ=='01-Oct'] <- as.character('OCT1')

S288c[S288c=="01-Oct"] <- as.character('OCT1')



######## ANAEROBIC
### Merge the TMT-data with iBAQ data and filter based on TMT-cutoffs
########

data_anaerobic_2 <- merge(data_anaerobic,iBAQ,by='Accession',all.y = TRUE)

# Remove proteins with only 1 PSM (to reduce worst effect of variability between samples for a protein), based on consultance
data_anaerobic_2 <- data_anaerobic_2[data_anaerobic_2$PSMs > 1, ]

# We also use cut-offs for the TMT-value of a protein (level in sample compared to reference)
# As a lower cutoff we just make sure that the protein does not have 0.01, which is the lowest value per default
# these are seen as non-existant since they are very low and we cannot now the exact value
# This cutoff is done after separating data-set into min and rich, as we believe that this
# protein may be present in one of the samples although it is regarded as non-present in the other
# As a higher cutoff we use the limit of 10, after consulting with proteomics consultant
# This was done because the values get significantly less reliable as they get larger, since we will compare sums of proteins
# between rich and minimal, we do not want the same protein to take over in one of the samples
# with a quantified mass percentage that can be very much over-estimated in one of rich or min

# We do this cutoff strict, in the sense that we remove a protein from all samples if it does not lie
# within the high-cutoff. To always compute protein abundancies from triplicate values 

TMT_cutoff_high <- 10
#
data_anaerobic_2 <- data_anaerobic_2[data_anaerobic_2$Min_1_glc_anaerobic < TMT_cutoff_high, ]
data_anaerobic_2 <- data_anaerobic_2[data_anaerobic_2$Min_2_glc_anaerobic < TMT_cutoff_high, ]
data_anaerobic_2 <- data_anaerobic_2[data_anaerobic_2$Min_3_glc_anaerobic < TMT_cutoff_high, ]
data_anaerobic_2 <- data_anaerobic_2[data_anaerobic_2$Rich_1_glc_anaerobic < TMT_cutoff_high, ]
data_anaerobic_2 <- data_anaerobic_2[data_anaerobic_2$Rich_2_glc_anaerobic < TMT_cutoff_high, ]
data_anaerobic_2 <- data_anaerobic_2[data_anaerobic_2$Rich_3_glc_anaerobic < TMT_cutoff_high, ]

# ALSO, manually check so there are no occurences of a protein with high TMT-values that has a NaN for just a individual replicate
# Since proteins with NaN will be removed separately in a later stage. They are seen as not present in the sample then, and should
# Therefore not be removed across the board. (more precisely, if a protein has a high TMT-value for all samples except 1, we would
# potentially bias the analysis by removing the protein completely in either Minimal or Rich alone, and not in the other).

# Split the data set into Minimal and Rich samples
data_anaerobic_2_min <- data_anaerobic_2[,c(1,5,6,7,15,16,17,18)]
data_anaerobic_2_rich <- data_anaerobic_2[,c(1,8,9,10,15,16,17,18)]

# Remove protein that has a TMT-value that is the lower default value
TMT_cutoff_low <- 0.01
#
data_anaerobic_2_min <- data_anaerobic_2_min[data_anaerobic_2_min$Min_1_glc_anaerobic > TMT_cutoff_low, ]
data_anaerobic_2_min <- data_anaerobic_2_min[data_anaerobic_2_min$Min_2_glc_anaerobic > TMT_cutoff_low, ]
data_anaerobic_2_min <- data_anaerobic_2_min[data_anaerobic_2_min$Min_3_glc_anaerobic > TMT_cutoff_low, ]
data_anaerobic_2_rich <- data_anaerobic_2_rich[data_anaerobic_2_rich$Rich_1_glc_anaerobic > TMT_cutoff_low, ]
data_anaerobic_2_rich <- data_anaerobic_2_rich[data_anaerobic_2_rich$Rich_2_glc_anaerobic > TMT_cutoff_low, ]
data_anaerobic_2_rich <- data_anaerobic_2_rich[data_anaerobic_2_rich$Rich_3_glc_anaerobic > TMT_cutoff_low, ]

# Remove any undetected protein within Minimal and Rich separately. 
data_anaerobic_2_min$sum <- rowSums(data_anaerobic_2_min[,2:4], na.rm=T)
data_anaerobic_2_min <- data_anaerobic_2_min[data_anaerobic_2_min$sum>0,]
data_anaerobic_2_min$sum <- NULL

data_anaerobic_2_rich$sum <- rowSums(data_anaerobic_2_rich[,2:4], na.rm=T)
data_anaerobic_2_rich <- data_anaerobic_2_rich[data_anaerobic_2_rich$sum>0,]
data_anaerobic_2_rich$sum <- NULL


#####
## We now can calculate the allocation for both sets of min & rich, and do this individually, so that the sum of all
## mass percentages, within each replicates equals 1
#####

j <- length(data_anaerobic_2_min[,1])
jj <- length(data_anaerobic_2_rich[,1])

glu_proteomics_anaerobic_min <- as.data.frame(matrix(0,nrow = j,ncol = 6))
glu_proteomics_anaerobic_rich <- as.data.frame(matrix(0,nrow = jj,ncol = 6))

# We calculate TMT-value * iBAQ-value -> fmol/ug_sample_injected (in our actual sample)
# We then multiply with 10^-15 to get mol/ug
# We finally multiply with molecular weight*1000 (g/mol) to get protein abundancy in g/ug_sample_injected

# Minimal samples
glu_proteomics_anaerobic_min[,1] <- as.character(data_anaerobic_2_min[,1])
glu_proteomics_anaerobic_min[,2] <- data_anaerobic_2_min[,2]*data_anaerobic_2_min[,6]*(10^-15)*data_anaerobic_2_min[,5]*1000
glu_proteomics_anaerobic_min[,3] <- data_anaerobic_2_min[,3]*data_anaerobic_2_min[,6]*(10^-15)*data_anaerobic_2_min[,5]*1000
glu_proteomics_anaerobic_min[,4] <- data_anaerobic_2_min[,4]*data_anaerobic_2_min[,6]*(10^-15)*data_anaerobic_2_min[,5]*1000
glu_proteomics_anaerobic_min[,5] <- as.character(data_anaerobic_2_min[,7])
glu_proteomics_anaerobic_min[,6] <- as.character(data_anaerobic_2_min[,8])
#
# We extract the total sums of proteins in individual replicates, in order to be able to normalize into mass percentages
sum_anaerobic_min_1 <- sum(glu_proteomics_anaerobic_min[,2])
sum_anaerobic_min_2 <- sum(glu_proteomics_anaerobic_min[,3])
sum_anaerobic_min_3 <- sum(glu_proteomics_anaerobic_min[,4])
#
glu_proteomics_anaerobic_min[,2] <- glu_proteomics_anaerobic_min[,2]/sum_anaerobic_min_1
glu_proteomics_anaerobic_min[,3] <- glu_proteomics_anaerobic_min[,3]/sum_anaerobic_min_2
glu_proteomics_anaerobic_min[,4] <- glu_proteomics_anaerobic_min[,4]/sum_anaerobic_min_3
#
colnames(glu_proteomics_anaerobic_min) <- c('Accession','Min_1','Min_2','Min_3','Gene','Protein')

### We do the same for Rich samples
glu_proteomics_anaerobic_rich[,1] <- as.character(data_anaerobic_2_rich[,1])
glu_proteomics_anaerobic_rich[,2] <- data_anaerobic_2_rich[,2]*data_anaerobic_2_rich[,6]*(10^-15)*data_anaerobic_2_rich[,5]*1000
glu_proteomics_anaerobic_rich[,3] <- data_anaerobic_2_rich[,3]*data_anaerobic_2_rich[,6]*(10^-15)*data_anaerobic_2_rich[,5]*1000
glu_proteomics_anaerobic_rich[,4] <- data_anaerobic_2_rich[,4]*data_anaerobic_2_rich[,6]*(10^-15)*data_anaerobic_2_rich[,5]*1000
glu_proteomics_anaerobic_rich[,5] <- as.character(data_anaerobic_2_rich[,7])
glu_proteomics_anaerobic_rich[,6] <- as.character(data_anaerobic_2_rich[,8])
#
# We extract the total sums of proteins in individual replicates, in order to be able to normalize into mass percentages
sum_anaerobic_rich_1 <- sum(glu_proteomics_anaerobic_rich[,2])
sum_anaerobic_rich_2 <- sum(glu_proteomics_anaerobic_rich[,3])
sum_anaerobic_rich_3 <- sum(glu_proteomics_anaerobic_rich[,4])
#
glu_proteomics_anaerobic_rich[,2] <- glu_proteomics_anaerobic_rich[,2]/sum_anaerobic_rich_1
glu_proteomics_anaerobic_rich[,3] <- glu_proteomics_anaerobic_rich[,3]/sum_anaerobic_rich_2
glu_proteomics_anaerobic_rich[,4] <- glu_proteomics_anaerobic_rich[,4]/sum_anaerobic_rich_3
#
colnames(glu_proteomics_anaerobic_rich) <- c('Accession','Rich_1','Rich_2','Rich_3','Gene','Protein')

#remap the protein column 
glu_proteomics_anaerobic_min <- merge(glu_proteomics_anaerobic_min[,1:5], S288c, all.x=T)
glu_proteomics_anaerobic_rich <- merge(glu_proteomics_anaerobic_rich[,1:5], S288c, all.x=T)


# Finally, to prepare for analysis of fold change and its significance level, we combine rich and min again
# And in this data set we only use proteins that are present in both samples. Since we cannot get the fold change of proteins
# that are only present in one of Rich or Minimal
glu_proteomics_anaerobic <- merge(glu_proteomics_anaerobic_min[,c(1,5,6,2:4)],glu_proteomics_anaerobic_rich[,1:4],
                                  by='Accession',all = TRUE)
# Make sure that all proteins occur in both sets
glu_proteomics_anaerobic <- na.omit(glu_proteomics_anaerobic)

#### We now do the exact same analysis on proteomics from aerobic samples

######## AEROBIC
### Merge the TMT-data with iBAQ data and filter based on TMT-cutoffs
########

data_aerobic_2 <- merge(data_aerobic,iBAQ,by='Accession',all.y = TRUE)

# Remove proteins with only 1 PSM (to reduce worst effect of variability between samples for a protein), based on consultance
data_aerobic_2 <- data_aerobic_2[data_aerobic_2$PSMs > 1, ]

# We also use cut-offs for the TMT-value of a protein (level in sample compared to reference)
# As a lower cutoff we just make sure that the protein does not have 0.01, which is the lowest value per default
# these are seen as non-existant since they are very low and we cannot now the exact value
# This cutoff is done after separating data-set into min and rich, as we believe that this
# protein may be present in one of the samples although it is regarded as non-present in the other
# As a higher cutoff we use the limit of 10, after consulting with proteomics consultant
# This was done because the values get significantly less reliable as they get larger, since we will compare sums of proteins
# between rich and minimal, we do not want the same protein to take over in one of the samples
# with a quantified mass percentage that can be very much over-estimated in one of rich or min

# We do this cutoff strict, in the sense that we remove a protein from all samples if it does not lie
# within the high-cutoff. To always compute protein abundancies from triplicate values 

TMT_cutoff_high <- 10
#
data_aerobic_2 <- data_aerobic_2[data_aerobic_2$Min_1_glc_aerobic < TMT_cutoff_high, ]
data_aerobic_2 <- data_aerobic_2[data_aerobic_2$Min_2_glc_aerobic < TMT_cutoff_high, ]
data_aerobic_2 <- data_aerobic_2[data_aerobic_2$Rich_1_glc_aerobic < TMT_cutoff_high, ]
data_aerobic_2 <- data_aerobic_2[data_aerobic_2$Rich_2_glc_aerobic < TMT_cutoff_high, ]

# ALSO, manually check so there are no occurences of a protein with high TMT-values that has a NaN for just a individual replicate
# Since proteins with NaN will be removed separately in a later stage. They are seen as not present in the sample then, and should
# Therefore not be removed across the board. (more precisely, if a protein has a high TMT-value for all samples except 1, we would
# potentially bias the analysis by removing the protein completely in either Minimal or Rich alone, and not in the other).

# Split the data set into Minimal and Rich samples
data_aerobic_2_min <- data_aerobic_2[,c(1,11,12,7,14,15,16)]
data_aerobic_2_rich <- data_aerobic_2[,c(1,8,9,7,14,15,16)]


# Remove protein that has a TMT-value that is the lower default value
TMT_cutoff_low <- 0.01
#
data_aerobic_2_min <- data_aerobic_2_min[data_aerobic_2_min$Min_1_glc_aerobic > TMT_cutoff_low, ]
data_aerobic_2_min <- data_aerobic_2_min[data_aerobic_2_min$Min_2_glc_aerobic > TMT_cutoff_low, ]
data_aerobic_2_rich <- data_aerobic_2_rich[data_aerobic_2_rich$Rich_1_glc_aerobic > TMT_cutoff_low, ]
data_aerobic_2_rich <- data_aerobic_2_rich[data_aerobic_2_rich$Rich_2_glc_aerobic > TMT_cutoff_low, ]


# Remove any undetected protein within Minimal and Rich separately. 
data_aerobic_2_min$sum <- rowSums(data_aerobic_2_min[,2:3], na.rm=T)
data_aerobic_2_min <- data_aerobic_2_min[data_aerobic_2_min$sum>0,]
data_aerobic_2_min$sum <- NULL

data_aerobic_2_rich$sum <- rowSums(data_aerobic_2_rich[,2:3], na.rm=T)
data_aerobic_2_rich <- data_aerobic_2_rich[data_aerobic_2_rich$sum>0,]
data_aerobic_2_rich$sum <- NULL


#####
## We now can calculate the allocation for both sets of min & rich, and do this individually, so that the sum of all
## mass percentages, within each replicates equals 1
#####

j <- length(data_aerobic_2_min[,1])
jj <- length(data_aerobic_2_rich[,1])

glu_proteomics_aerobic_min <- as.data.frame(matrix(0,nrow = j,ncol = 5))
glu_proteomics_aerobic_rich <- as.data.frame(matrix(0,nrow = jj,ncol = 5))


# We calculate TMT-value * iBAQ-value -> fmol/ug_sample_injected (in our actual sample)
# We then multiply with 10^-15 to get mol/ug
# We finally multiply with molecular weight*1000 (g/mol) to get protein abundancy in g/ug_sample_injected

# Minimal samples
glu_proteomics_aerobic_min[,1] <- as.character(data_aerobic_2_min[,1])
glu_proteomics_aerobic_min[,2] <- data_aerobic_2_min[,2]*data_aerobic_2_min[,5]*(10^-15)*data_aerobic_2_min[,4]*1000
glu_proteomics_aerobic_min[,3] <- data_aerobic_2_min[,3]*data_aerobic_2_min[,5]*(10^-15)*data_aerobic_2_min[,4]*1000
glu_proteomics_aerobic_min[,4] <- as.character(data_aerobic_2_min[,6])
glu_proteomics_aerobic_min[,5] <- as.character(data_aerobic_2_min[,7])
#
# We extract the total sums of proteins in individual replicates, in order to be able to normalize into mass percentages
sum_aerobic_min_1 <- sum(glu_proteomics_aerobic_min[,2])
sum_aerobic_min_2 <- sum(glu_proteomics_aerobic_min[,3])
#
glu_proteomics_aerobic_min[,2] <- glu_proteomics_aerobic_min[,2]/sum_aerobic_min_1
glu_proteomics_aerobic_min[,3] <- glu_proteomics_aerobic_min[,3]/sum_aerobic_min_2
#
colnames(glu_proteomics_aerobic_min) <- c('Accession','Min_1','Min_2','Gene','Protein')

### We do the same for Rich samples
glu_proteomics_aerobic_rich[,1] <- as.character(data_aerobic_2_rich[,1])
glu_proteomics_aerobic_rich[,2] <- data_aerobic_2_rich[,2]*data_aerobic_2_rich[,5]*(10^-15)*data_aerobic_2_rich[,4]*1000
glu_proteomics_aerobic_rich[,3] <- data_aerobic_2_rich[,3]*data_aerobic_2_rich[,5]*(10^-15)*data_aerobic_2_rich[,4]*1000
glu_proteomics_aerobic_rich[,4] <- as.character(data_aerobic_2_rich[,6])
glu_proteomics_aerobic_rich[,5] <- as.character(data_aerobic_2_rich[,7])
#
# We extract the total sums of proteins in individual replicates, in order to be able to normalize into mass percentages
sum_aerobic_rich_1 <- sum(glu_proteomics_aerobic_rich[,2])
sum_aerobic_rich_2 <- sum(glu_proteomics_aerobic_rich[,3])
#
glu_proteomics_aerobic_rich[,2] <- glu_proteomics_aerobic_rich[,2]/sum_aerobic_rich_1
glu_proteomics_aerobic_rich[,3] <- glu_proteomics_aerobic_rich[,3]/sum_aerobic_rich_2
#
colnames(glu_proteomics_aerobic_rich) <- c('Accession','Rich_1','Rich_2','Gene','Protein')

#remap the protein column 
glu_proteomics_aerobic_min <- merge(glu_proteomics_aerobic_min[,1:4], S288c, all.x=T)
glu_proteomics_aerobic_rich <- merge(glu_proteomics_aerobic_rich[,1:4], S288c, all.x=T)


# Finally, to prepare for analysis of fold change and its significance level, we combine rich and min again
# And in this data set we only use proteins that are present in both samples. Since we cannot get the fold change of proteins
# that are only present in one of Rich or Minimal

glu_proteomics_aerobic <- merge(glu_proteomics_aerobic_min[,c(1,4,5,2:3)],glu_proteomics_aerobic_rich[,1:3],
                                by='Accession',all = TRUE)
# Make sure that all proteins occur in both sets
glu_proteomics_aerobic <- na.omit(glu_proteomics_aerobic)


print('done')







###write tables
Dataset_Proteomics <- merge(glu_proteomics_aerobic_min[,c(1,4,5,2:3)], glu_proteomics_aerobic_rich[,c(1,4,5,2:3)], all=T)
colnames(Dataset_Proteomics)[4:7] <- c("Min_aerobic_1", "Min_aerobic_2", "Rich_aerobic_1", "Rich_aerobic_2")
Dataset_Proteomics <- merge(Dataset_Proteomics, glu_proteomics_anaerobic_min[,c(1,5,6,2:4)], all=T)
colnames(Dataset_Proteomics)[8:10] <- c("Min_anaerobic_1", "Min_anaerobic_2", "Min_anaerobic_3")
Dataset_Proteomics <- merge(Dataset_Proteomics, glu_proteomics_anaerobic_rich[,c(1,5,6,2:4)], all=T)
colnames(Dataset_Proteomics)[11:13] <- c("Rich_anaerobic_1", "Rich_anaerobic_2", "Rich_anaerobic_3")

write.table(glu_proteomics_aerobic_min, file = "glu_proteomics_aerobic_min.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
write.table(glu_proteomics_aerobic_rich, file = "glu_proteomics_aerobic_rich.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
write.table(glu_proteomics_anaerobic_min, file = "glu_proteomics_anaerobic_min.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
write.table(glu_proteomics_anaerobic_rich, file = "glu_proteomics_anaerobic_rich.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
write.table(glu_proteomics_aerobic, file = "glu_proteomics_aerobic.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)
write.table(glu_proteomics_anaerobic, file = "glu_proteomics_anaerobic.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)

write.table(Dataset_Proteomics, file = "20200505_Dataset_Proteomics.csv", 
            sep = ",", col.names = TRUE, row.names = FALSE)











# p-value and FC

anaerobic_DE <- glu_proteomics_anaerobic
anaerobic_DE$p.value <- NA
anaerobic_DE$log2.FC.rich.over.min.ana <- NA

for (i in 1:nrow(anaerobic_DE)){
  anaerobic_DE[i,"p.value"] <- t.test(log2(anaerobic_DE[i,4:6]), log2(anaerobic_DE[i,7:9]))$p.value
  anaerobic_DE[i,"log2.FC.rich.over.min.ana"] <- log2(rowMeans(anaerobic_DE[i,7:9]) / rowMeans(anaerobic_DE[i,4:6]))
}
anaerobic_DE$FDR.ana <- p.adjust(anaerobic_DE$p.value, method ="fdr")

aerobic_DE <- glu_proteomics_aerobic
aerobic_DE$p.value.aer <- NA
aerobic_DE$log2.FC.rich.over.min.aer <- NA

for(i in 1:nrow(aerobic_DE)){
  aerobic_DE[i,"p.value.aer"] <- t.test(log2(aerobic_DE[i,4:5]), log2(aerobic_DE[i,6:7]))$p.value
  aerobic_DE[i,"log2.FC.rich.over.min.aer"] <- log2( rowMeans(aerobic_DE[i,6:7]) / rowMeans(aerobic_DE[i,4:5]) )
}





upreg_aerobic <- aerobic_DE[aerobic_DE$log2.FC.rich.over.min>1 & aerobic_DE$p.value<0.05,] #18 proteins
upreg_anaerobic <- anaerobic_DE[anaerobic_DE$log2.FC.rich.over.min>1 & anaerobic_DE$FDR<0.05,] #12 proteins
upreg_both <- merge(upreg_aerobic[,c(1:3,9,8)], upreg_anaerobic[,c(1:3,11,12)]) #9 proteins
upreg_both <- upreg_both[,c(1:3,4,6,5,7)]



downreg_aerobic <- aerobic_DE[aerobic_DE$log2.FC.rich.over.min < (-1) & aerobic_DE$p.value<0.05,] #59
downreg_anaerobic <- anaerobic_DE[anaerobic_DE$log2.FC.rich.over.min < (-1) & anaerobic_DE$FDR<0.05,] #69
downreg_both <- merge(downreg_aerobic[,c(1:3,9,8)], downreg_anaerobic[,c(1:3,11,12)]) #46
downreg_both <- downreg_both[,c(1:3,4,6,5,7)]



DE_Supp_Figs_no_FC_cutoff <- aerobic_DE[aerobic_DE$p.value<0.05,]
DE_Supp_Figs_no_FC_cutoff <- merge(DE_Supp_Figs_no_FC_cutoff[,c(1:3,8,9)], 
                                   anaerobic_DE[anaerobic_DE$FDR<0.05,], all=T)
DE_Supp_Figs_no_FC_cutoff <- DE_Supp_Figs_no_FC_cutoff[,c(1:5,13,14)]



###Fig 3C-D

barplot(rev(c(rowMeans(Dataset_Proteomics[42,6:7]), #rich aerobic
              rowMeans(Dataset_Proteomics[42,4:5]),  # min aerobic
          rowMeans(Dataset_Proteomics[42,11:13]),  #rich anaerobic
          rowMeans(Dataset_Proteomics[42,8:10]))*100), #min anaerobic
        xlim=c(0,0.03), main=Dataset_Proteomics[42,3],
        horiz = T, names.arg = "")
#abline(v=c(0.01,0.02,0.03), col="grey")

barplot(rev(c(rowMeans(Dataset_Proteomics[127,6:7]), #rich aerobic
              rowMeans(Dataset_Proteomics[127,4:5]),  # min aerobic
              rowMeans(Dataset_Proteomics[127,11:13]),  #rich anaerobic
              rowMeans(Dataset_Proteomics[127,8:10]))*100), #min anaerobic
        xlim=c(0,2.5), main=Dataset_Proteomics[127,3],
        horiz = T, names.arg = "", 
        col=rev(c("orange4", "orange1", "blue3", "cadetblue2")))
#abline(v=c(.5,1,1.5,2,2.5), col="grey")



#Fig 3E-F

# Match all individual samples against YeastGenome (GO)-slim-mapper terms (obtained in March 2019).
# This is a way to analyze differences in protein groups, rather only individual proteins
# 


file <- "Slim_mapper_proteins_iBAQ.csv"
test_slim_mapper_proteins <- read.csv(file,header=TRUE,sep=";",dec=',',na.strings=c("", 'NA'))

i <- nrow(test_slim_mapper_proteins)
j <- ncol(test_slim_mapper_proteins)

column_names_matrices <- as.data.frame(matrix(0,nrow = 1,ncol = j*2))


for(i in seq(1,j,by=1)){
  
  column_names_matrices[1,(i*2-1)] <- as.character(paste0(as.character(colnames(test_slim_mapper_proteins)[i]),'_Rich'))
  column_names_matrices[1,(i*2)] <- as.character(paste0(as.character(colnames(test_slim_mapper_proteins)[i]),'_Min'))
  
}
##################################################################

######### PERFORM CALCULATIONS BASED ON THE GROUPS ABOVE #########

##################################################################

aerobic_slim <- merge(glu_proteomics_aerobic_min,glu_proteomics_aerobic_rich,by='Accession',all=TRUE)
anaerobic_slim <- merge(glu_proteomics_anaerobic_min,glu_proteomics_anaerobic_rich,by='Accession',all=TRUE)
aerobic_slim[is.na(aerobic_slim)] <- 0
anaerobic_slim[is.na(anaerobic_slim)] <- 0

Input_aer_groups <- aerobic_slim[,c(1,4,5,6,7,2,3)]
colnames(Input_aer_groups) <- c('Accession','Gene','Protein','Rich_1','Rich_2','Min_1','Min_2')
Input_anaer_groups <- anaerobic_slim[,c(1,5,6,7,8,9,2,3,4)]
colnames(Input_anaer_groups) <- c('Accession','Gene','Protein','Rich_1','Rich_2','Rich_3','Min_1','Min_2','Min_3')


Pvalue_groups_aerobic <- as.data.frame(matrix(0,nrow = 1,ncol = j))
colnames(Pvalue_groups_aerobic) <- colnames(test_slim_mapper_proteins)
Stdev_groups_aerobic <- as.data.frame(matrix(0,nrow = 1,ncol = j*2))
colnames(Stdev_groups_aerobic) <- column_names_matrices
Mass_allocation_groups_aerobic <- as.data.frame(matrix(0,nrow = 1,ncol = j*2))
colnames(Mass_allocation_groups_aerobic) <- column_names_matrices
Pvalue_groups_anaerobic <- as.data.frame(matrix(0,nrow = 1,ncol = j))
colnames(Pvalue_groups_anaerobic) <- colnames(test_slim_mapper_proteins)
Stdev_groups_anaerobic <- as.data.frame(matrix(0,nrow = 1,ncol = j*2))
colnames(Stdev_groups_anaerobic) <- column_names_matrices
Mass_allocation_groups_anaerobic <- as.data.frame(matrix(0,nrow = 1,ncol = j*2))
colnames(Mass_allocation_groups_anaerobic) <- column_names_matrices

# For later purposes we also need to save the number of identified proteins in each group that is summarized
n_and_p <- as.data.frame(matrix(0,nrow = 1,ncol = 3))
colnames(n_and_p) <- c("Process", "aerobic.red.n", "anaerobic.blue.n")

for(r in seq(1,j,by=1)){
  
  Proteins_to_match <- as.data.frame(test_slim_mapper_proteins[,r])
  Proteins_to_match <- na.omit(Proteins_to_match)
  
  # We make sure that the protein OCT1 gets its correct name, instead of the automatic 01-Oct, at import
  Proteins_to_match[,1] <- as.character(Proteins_to_match[,1])
  Proteins_to_match[Proteins_to_match=='01-Oct'] <- as.character('OCT1')
  
  # WHEN MERGING WE SEARCH FIRST AGAINST PROTEIN AND THEN AGAINST GENE, SINCE BOTH NAME TYPES OCCUR IN Proteins_to_match
  
  colnames(Proteins_to_match) <- c('Protein')
  aerobic_1 <- merge(Input_aer_groups,Proteins_to_match,by='Protein',all.x = FALSE)
  anaerobic_1 <- merge(Input_anaer_groups,Proteins_to_match,by='Protein',all.x = FALSE)
  colnames(Proteins_to_match) <- c('Gene')
  aerobic_2 <- merge(Input_aer_groups,Proteins_to_match,by='Gene',all.x = FALSE)
  anaerobic_2 <- merge(Input_anaer_groups,Proteins_to_match,by='Gene',all.x = FALSE)
  
  # Take out the actual data and rbind _1 and _2 together for both cases
  
  aerobic_1 <- aerobic_1[,c(4,5,6,7)]
  aerobic_2 <- aerobic_2[,c(4,5,6,7)]
  anaerobic_1 <- anaerobic_1[,c(4,5,6,7,8,9)]
  anaerobic_2 <- anaerobic_2[,c(4,5,6,7,8,9)]
  
  aerobic_data <- rbind(aerobic_1,aerobic_2)
  anaerobic_data <- rbind(anaerobic_1,anaerobic_2)
  
  # For later purposes we save the number of identified proteins in each group that is summarized
  n_and_p[r,1] <- colnames(test_slim_mapper_proteins)[r]
  n_and_p[r,2] <- nrow(aerobic_data)
  n_and_p[r,3] <- nrow(anaerobic_data)
  
  
  
  #### T_TEST
  # T-test on log2-transformed values
  ###
  ttesting_aerobic <- t.test(c(as.numeric(log2(sum(aerobic_data[,1]))),as.numeric(log2(sum(aerobic_data[,2])))),c(as.numeric(log2(sum(aerobic_data[,3]))),as.numeric(log2(sum(aerobic_data[,4])))),alternative='t',paired=FALSE)
  pvalue_aerobic <- ttesting_aerobic$p.value
  ttesting_anaerobic <- t.test(c(as.numeric(log2(sum(anaerobic_data[,1]))),as.numeric(log2(sum(anaerobic_data[,2]))),as.numeric(log2(sum(anaerobic_data[,3])))),c(as.numeric(log2(sum(anaerobic_data[,4]))),as.numeric(log2(sum(anaerobic_data[,5]))),as.numeric(log2(sum(anaerobic_data[,6])))),alternative='t',paired=FALSE)
  pvalue_anaerobic <- ttesting_anaerobic$p.value
  
  Pvalue_groups_aerobic[1,r] <- pvalue_aerobic
  Pvalue_groups_anaerobic[1,r] <- pvalue_anaerobic
  
  #### Summarized average mass-% for each group
  #
  ####
  Mass_allocation_groups_aerobic[1,r*2-1] <- mean(c(as.numeric(sum(aerobic_data[,1])*100),as.numeric(sum(aerobic_data[,2])*100)))
  Mass_allocation_groups_aerobic[1,r*2] <- mean(c(as.numeric(sum(aerobic_data[,3])*100),as.numeric(sum(aerobic_data[,4])*100)))
  Mass_allocation_groups_anaerobic[1,r*2-1] <- mean(c(as.numeric(sum(anaerobic_data[,1])*100),as.numeric(sum(anaerobic_data[,2])*100),as.numeric(sum(anaerobic_data[,3])*100)))
  Mass_allocation_groups_anaerobic[1,r*2] <- mean(c(as.numeric(sum(anaerobic_data[,4])*100),as.numeric(sum(anaerobic_data[,5])*100),as.numeric(sum(anaerobic_data[,6])*100)))
  
  ####Standard deviation of mass-% summation for each group
  #
  ####
  Stdev_groups_aerobic[1,r*2-1] <- sd(c(as.numeric(sum(aerobic_data[,1])*100),as.numeric(sum(aerobic_data[,2])*100)))
  Stdev_groups_aerobic[1,r*2] <- sd(c(as.numeric(sum(aerobic_data[,3])*100),as.numeric(sum(aerobic_data[,4])*100)))
  Stdev_groups_anaerobic[1,r*2-1] <- sd(c(as.numeric(sum(anaerobic_data[,1])*100),as.numeric(sum(anaerobic_data[,2])*100),as.numeric(sum(anaerobic_data[,3])*100)))
  Stdev_groups_anaerobic[1,r*2] <- sd(c(as.numeric(sum(anaerobic_data[,4])*100),as.numeric(sum(anaerobic_data[,5])*100),as.numeric(sum(anaerobic_data[,6])*100)))
  
}


check_collection_aerobic <- as.data.frame(matrix(nrow=202,ncol=5,0))
check_collection_anaerobic <- as.data.frame(matrix(nrow=202,ncol=5,0))

for(i in seq(1,202,by=1)){
  
  # Collect names of processes
  check_collection_aerobic[i,1] <- as.character(colnames(Mass_allocation_groups_aerobic)[i])
  check_collection_anaerobic[i,1] <- as.character(colnames(Mass_allocation_groups_anaerobic)[i])
  
  # Collect average_mass, stdev_mass and p_value corresponding to each group
  # For aerobic samples
  check_collection_aerobic[i,2] <- as.numeric(Mass_allocation_groups_aerobic[1,i])
  check_collection_aerobic[i,3] <- as.numeric(Stdev_groups_aerobic[1,i])
  check_collection_aerobic[i,5] <- as.numeric(i)
  if(i<102){
    check_collection_aerobic[((i*2)-1),4] <- as.numeric(Pvalue_groups_aerobic[1,i])
    check_collection_aerobic[(i*2),4] <- as.numeric(Pvalue_groups_aerobic[1,i])
  }
  # Collect average_mass, stdev_mass and p_value corresponding to each group
  # For anaerobic samples
  check_collection_anaerobic[i,2] <- as.numeric(Mass_allocation_groups_anaerobic[1,i])
  check_collection_anaerobic[i,3] <- as.numeric(Stdev_groups_anaerobic[1,i])
  check_collection_anaerobic[i,5] <- as.numeric(i)
  if(i<102){
    check_collection_anaerobic[((i*2)-1),4] <- as.numeric(Pvalue_groups_anaerobic[1,i])
    check_collection_anaerobic[(i*2),4] <- as.numeric(Pvalue_groups_anaerobic[1,i])
  }
  
}

colnames(check_collection_aerobic) <- c('Process_and_sample','Average_mass','Stdev_mass','P_value','Identifier')
colnames(check_collection_anaerobic) <- c('Process_and_sample','Average_mass','Stdev_mass','P_value','Identifier')


# Subset the data and collect processes that have p<0.05

check_collection_aerobic_for_filter <- check_collection_aerobic[check_collection_aerobic$P_value < 0.05,]
check_collection_anaerobic_for_filter <- check_collection_anaerobic[check_collection_anaerobic$P_value < 0.05,]
#

# Filter on average_mass. If either rich sample or minimal sample shows >2% average mass, the process is saved

i1 <- nrow(check_collection_aerobic_for_filter)/2
i2 <- nrow(check_collection_anaerobic_for_filter)/2

for(i in seq(1,i1,by=1)){
  
  max_average_mass_aerobic <- max(c(as.numeric(check_collection_aerobic_for_filter[((i*2)-1),2]),as.numeric(check_collection_aerobic_for_filter[i*2,2])))
  
  if(max_average_mass_aerobic<2){
    check_collection_aerobic_for_filter[((i*2)-1),] <- NA
    check_collection_aerobic_for_filter[(i*2),] <- NA
  }
}

for(i in seq(1,i2,by=1)){
  
  max_average_mass_anaerobic <- max(c(as.numeric(check_collection_anaerobic_for_filter[((i*2)-1),2]),as.numeric(check_collection_anaerobic_for_filter[i*2,2])))
  
  if(max_average_mass_anaerobic<2){
    check_collection_anaerobic_for_filter[((i*2)-1),] <- NA
    check_collection_anaerobic_for_filter[(i*2),] <- NA
  }
}

check_collection_aerobic_for_filter <- na.omit(check_collection_aerobic_for_filter) #24
check_collection_anaerobic_for_filter <- na.omit(check_collection_anaerobic_for_filter) #44
#




# Collect and unite the identifiers of all positive matches, and remove duplicates
identifier_aerobic <- as.data.frame(check_collection_aerobic_for_filter$Identifier)
identifier_anaerobic <- as.data.frame(check_collection_anaerobic_for_filter$Identifier)
colnames(identifier_aerobic) <- c('Identifier')
colnames(identifier_anaerobic) <- c('Identifier')
identifier_vector <- rbind(identifier_aerobic,identifier_anaerobic)
identifier_vector <- unique(identifier_vector)

Plot_these_from_aerobic <- merge(identifier_vector,check_collection_aerobic,by='Identifier',all.x=FALSE)
Plot_these_from_anaerobic <- merge(identifier_vector,check_collection_anaerobic,by='Identifier',all.x=FALSE)


# Set-up the data in the correct way and Check the descending order of mass for minimal samples in anaerobic condition and add this information

Plot_these_from_aerobic_final <- as.data.frame(matrix(nrow=nrow(Plot_these_from_aerobic),ncol=4,0))
Plot_these_from_anaerobic_final <- as.data.frame(matrix(nrow=nrow(Plot_these_from_anaerobic),ncol=4,0))

# Aerobic and anaerobic have the same number of rows, since we have matched to the same identifiers
# We run the collection in the same loop
for(i in seq(1,nrow(Plot_these_from_aerobic)/2,by=1)){
  
  # This will be used for grouping
  Plot_these_from_aerobic_final[(i*2-1),1] <- as.character('Rich')
  Plot_these_from_aerobic_final[i*2,1] <- as.character('Min')
  Plot_these_from_anaerobic_final[(i*2-1),1] <- as.character('Rich')
  Plot_these_from_anaerobic_final[i*2,1] <- as.character('Min')
  
  # Subset the name of the process so they show the same for minimal and rich
  name_aerobic_rich <- Plot_these_from_aerobic[(i*2-1),2]
  name_aerobic_min <- Plot_these_from_aerobic[(i*2),2]
  name_anaerobic_rich <- Plot_these_from_aerobic[(i*2-1),2]
  name_anaerobic_min <- Plot_these_from_aerobic[(i*2),2]
  #
  Plot_these_from_aerobic_final[(i*2-1),2] <- as.character(substr(name_aerobic_rich,1,(nchar(name_aerobic_rich)-5)))
  Plot_these_from_aerobic_final[(i*2),2] <- as.character(substr(name_aerobic_min,1,(nchar(name_aerobic_min)-4)))
  Plot_these_from_anaerobic_final[(i*2-1),2] <- as.character(substr(name_anaerobic_rich,1,(nchar(name_anaerobic_rich)-5)))
  Plot_these_from_anaerobic_final[(i*2),2] <- as.character(substr(name_anaerobic_min,1,(nchar(name_anaerobic_min)-4)))
  
  # Collect the average masses again
  Plot_these_from_aerobic_final[(i*2-1),3] <- Plot_these_from_aerobic[(i*2-1),3]
  Plot_these_from_aerobic_final[(i*2),3] <- Plot_these_from_aerobic[(i*2),3]
  Plot_these_from_anaerobic_final[(i*2-1),3] <- Plot_these_from_anaerobic[(i*2-1),3]
  Plot_these_from_anaerobic_final[(i*2),3] <- Plot_these_from_anaerobic[(i*2),3]
  
  # Collect average Stdev of masses again
  Plot_these_from_aerobic_final[(i*2-1),4] <- Plot_these_from_aerobic[(i*2-1),4]
  Plot_these_from_aerobic_final[(i*2),4] <- Plot_these_from_aerobic[(i*2),4]
  Plot_these_from_anaerobic_final[(i*2-1),4] <- Plot_these_from_anaerobic[(i*2-1),4]
  Plot_these_from_anaerobic_final[(i*2),4] <- Plot_these_from_anaerobic[(i*2),4]
}


# Input order factor (described above)

Plot_these_from_aerobic_final[,5] <- as.numeric(c(12,12,18,18,22,22,3,3,7,7,10,10,27,27,23,23,21,21,8,8,6,6,4,4,25,25,26,26,5,5,2,2,17,17,24,24,14,14,20,20,16,16,1,1,9,9,13,13,15,15,19,19,11,11))
Plot_these_from_anaerobic_final[,5] <- as.numeric(c(12,12,18,18,22,22,3,3,7,7,10,10,27,27,23,23,21,21,8,8,6,6,4,4,25,25,26,26,5,5,2,2,17,17,24,24,14,14,20,20,16,16,1,1,9,9,13,13,15,15,19,19,11,11))

colnames(Plot_these_from_aerobic_final) <- c('Sample','Process','Average_mass','Stdev_average_mass','Order')
colnames(Plot_these_from_anaerobic_final) <- c('Sample','Process','Average_mass','Stdev_average_mass','Order')


# PLOT FOR AEROBIC CONDITIONS
tiff("Fig3e.tif", width=4, height=6, units = 'in', res = 300)
ggplot(data=Plot_these_from_aerobic_final, aes(x=reorder(Process,Order), y=Average_mass, fill=Sample)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),size=0.1) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  scale_fill_manual(values = c('orange','brown4')) +
  scale_y_continuous(limits=c(0,35),breaks=c(0,10,20,30)) +
  labs(x=" ") +
  labs(y="Mass ratio of proteome") +
  labs(fill=' ') +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(face = c('bold'),size=7)) +
  theme(axis.text.y = element_text(face = c('bold'),size=7)) +
  geom_errorbar(aes(ymin=`Average_mass`-`Stdev_average_mass`,
                    ymax=`Average_mass`+`Stdev_average_mass`), 
                width=.4,position=position_dodge(.9),size=0.3)


dev.off()


# PLOT FOR ANAEROBIC CONDITIONS
tiff("Fig3f.tif", width=4, height=6, units = 'in', res = 300)
ggplot(data=Plot_these_from_anaerobic_final, aes(x=reorder(Process,Order), y=Average_mass, fill=Sample)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),size=0.1) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  scale_fill_manual(values = c('deepskyblue2','purple4')) +
  scale_y_continuous(limits=c(0,35),breaks=c(0,10,20,30)) +
  labs(x=" ") +
  labs(y="Mass ratio of proteome") +
  labs(fill=' ') +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(face = c('bold'),size=7)) +
  theme(axis.text.y = element_text(face = c('bold'),size=7)) +
  geom_errorbar(aes(ymin=`Average_mass`-`Stdev_average_mass`, 
                    ymax=`Average_mass`+`Stdev_average_mass`), 
                width=.4,position=position_dodge(.9),size=0.3)

dev.off()




#get numbers and p-values

Fig3E_n_and_p <- Plot_these_from_aerobic_final
Fig3E_n_and_p <- unique(Fig3E_n_and_p[,c(2,5)])

Fig3E_n_and_p <- merge(Fig3E_n_and_p, n_and_p)
Fig3E_n_and_p <- Fig3E_n_and_p[,c(2,1,4,3)]


Fig3E_n_and_p$anaerobic.blue.p <- NA
Fig3E_n_and_p$aerobic.red.p <- NA


#* p-value < 5e-2, ** p-value < 1e-2, *** p-value < 1e-3, **** p-value < 1e-4.

Fig3E_n_and_p[Fig3E_n_and_p$anaerobic.blue.p<1e-4, "anaerobic.blue.p"] <- 4
Fig3E_n_and_p[Fig3E_n_and_p$anaerobic.blue.p<1e-3, "anaerobic.blue.p"] <- 3
Fig3E_n_and_p[Fig3E_n_and_p$anaerobic.blue.p<1e-2, "anaerobic.blue.p"] <- 2
Fig3E_n_and_p[Fig3E_n_and_p$anaerobic.blue.p<5e-2, "anaerobic.blue.p"] <- 1
Fig3E_n_and_p[Fig3E_n_and_p$anaerobic.blue.p<1, "anaerobic.blue.p"] <- NA

Fig3E_n_and_p[Fig3E_n_and_p$aerobic.red.p<1e-4, "aerobic.red.p"] <- 4
Fig3E_n_and_p[Fig3E_n_and_p$aerobic.red.p<1e-3, "aerobic.red.p"] <- 3
Fig3E_n_and_p[Fig3E_n_and_p$aerobic.red.p<1e-2, "aerobic.red.p"] <- 2
Fig3E_n_and_p[Fig3E_n_and_p$aerobic.red.p<5e-2, "aerobic.red.p"] <- 1
Fig3E_n_and_p[Fig3E_n_and_p$aerobic.red.p<1, "aerobic.red.p"] <- NA

Fig3E_n_and_p <- Fig3E_n_and_p[rev(rank(Fig3E_n_and_p$Order)),]







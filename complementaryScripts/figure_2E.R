# To map significantly changing proteins (both on FC and p-value) to metabolism of amino acids,
# we pull out the proteins listed in pathways of synthesis degradation of each amino acid, from yeastgenome database (SGD).

library(ggplot2)
library(dplyr)

# Define proteins involved in amino acid metabolism
# biosynthesis
Arginine_biosynthesis <- data.frame(c('ARG2','ARG7','ARG5,6','ARG8','ARG3','ARG1','ARG4','CPA1','CPA2'))
colnames(Arginine_biosynthesis) <- c('Protein')
Aspartate_biosynthesis <- data.frame(c('AAT1','AAT2'))
colnames(Aspartate_biosynthesis) <- c('Protein')
Glutamate_biosynthesis <- data.frame(c('GLN1','GLT1','GDH1','GDH3'))
colnames(Glutamate_biosynthesis) <- c('Protein')
Glycine_biosynthesis <- data.frame(c('AGX1','GLY1','SHM1','SHM2'))
colnames(Glycine_biosynthesis) <- c('Protein')
Histidine_biosynthesis <- data.frame(c('HIS1','HIS2','HIS3','HIS4','HIS5','HIS6','HIS7'))
colnames(Histidine_biosynthesis) <- c('Protein')
Isoleucine_biosynthesis <- data.frame(c('ILV1','ILV2','ILV3','ILV5','ILV6','BAT1','BAT2'))
colnames(Isoleucine_biosynthesis) <- c('Protein')
Leucine_biosynthesis <- data.frame(c('LEU1','LEU2','LEU4','LEU9','BAT1','BAT2'))
colnames(Leucine_biosynthesis) <- c('Protein')
Lysine_biosynthesis <- data.frame(c('LYS1','LYS2','LYS4','LYS9','LYS12','LYS20','LYS21'))
colnames(Lysine_biosynthesis) <- c('Protein')
Methionine_biosynthesis <- data.frame(c('MET1','MET2','MET3','MET5','MET6','MET7','MET8','MET10','MET12','MET13','MET14','MET16','MET17','MET18',
                                        'MET22','MET28','HOM2','HOM3','HOM6'))
colnames(Methionine_biosynthesis) <- c('Protein')
Phenylalanine_biosynthesis <- data.frame(c('ARO7','ARO8','ARO9','PHA2'))
colnames(Phenylalanine_biosynthesis) <- c('Protein')
Threonine_biosynthesis <- data.frame(c('HOM2','HOM3','HOM6','THR1','THR4'))
colnames(Threonine_biosynthesis) <- c('Protein')
Tryptophan_biosynthesis <- data.frame(c('TRP1','TRP2','TRP3','TRP4','TRP5'))
colnames(Tryptophan_biosynthesis) <- c('Protein')
Tyrosine_biosynthesis <- data.frame(c('ARO7','ARO8','ARO9','TYR1'))
colnames(Tyrosine_biosynthesis) <- c('Protein')
Valine_biosynthesis <- data.frame(c('ILV2','ILV3','ILV5','ILV6','BAT1','BAT2'))
colnames(Valine_biosynthesis) <- c('Protein')

# degradation
Arginine_degradation <- data.frame(c('CAR1','CAR2','PRO3','DUR1,2'))
colnames(Arginine_degradation) <- c('Protein')
Aspartate_degradation <- data.frame(c('NO'))
colnames(Aspartate_degradation) <- c('Protein')
Glutamate_degradation <- data.frame(c('GDH2','GAD1','UGA1','UGA2')) 
colnames(Glutamate_degradation) <- c('Protein')
Glycine_degradation <- data.frame(c('LPD1','GCV1','GCV2','GCV3'))
colnames(Glycine_degradation) <- c('Protein')
Histidine_degradation <- data.frame(c('NO'))
colnames(Histidine_degradation) <- c('Protein')
Isoleucine_degradation <- data.frame(c('THI3','ARO10','PDC1','PDC5','PDC6','SFA1','ADH5','ADH4'))
colnames(Isoleucine_degradation) <- c('Protein')
Leucine_degradation <- data.frame(c('THI3','ARO10','SFA1','ADH5','ADH4'))
colnames(Leucine_degradation) <- c('Protein')
Lysine_degradation <- data.frame(c('NO'))
colnames(Lysine_degradation) <- c('Protein')
Methionine_degradation <- data.frame(c('ADI1','UTR4','MDE1','MRI1','MEU1','SPE2','SPE3','SPE4','SAM1','SAM2','ARO8','ARO9','ARO10')) 
colnames(Methionine_degradation) <- c('Protein')
Phenylalanine_degradation <- data.frame(c('ARO9','ARO8','HIS5','ARO10','PDC1','PDC5','PDC6','SFA1','ADH4','ADH5'))
colnames(Phenylalanine_degradation) <- c('Protein')
Threonine_degradation <- data.frame(c('CHA1'))
colnames(Threonine_degradation) <- c('Protein')
Tryptophan_degradation <- data.frame(c('ARO8','ARO9','ARO10','PDC1','PDC5','PDC6','ADH1','ADH2','ADH3','BNA1','BNA2','BNA3','BNA4','BNA5','BNA6','BNA7'))
colnames(Tryptophan_degradation) <- c('Protein')
Tyrosine_degradation <- data.frame(c('ARO8','ARO9','PDC1','PDC5','PDC6','ADH1','ADH2','ADH3','ADH4','ADH5'))
colnames(Tyrosine_degradation) <- c('Protein')
Valine_degradation <- data.frame(c('BAT1','BAT2','PDC1','PDC5','PDC6','SFA1','ADH5','ADH4'))
colnames(Valine_degradation) <- c('Protein')

Arginine_metabolism <- unique(rbind(Arginine_biosynthesis,Arginine_degradation))
Aspartate_metabolism <- unique(rbind(Aspartate_biosynthesis,Aspartate_degradation))
Glutamate_metabolism <- unique(rbind(Glutamate_biosynthesis,Glutamate_degradation))
Glycine_metabolism <- unique(rbind(Glycine_biosynthesis,Glycine_degradation))
Histidine_metabolism <- unique(rbind(Histidine_biosynthesis,Histidine_degradation))
Isoleucine_metabolism <- unique(rbind(Isoleucine_biosynthesis,Isoleucine_degradation))
Leucine_metabolism <- unique(rbind(Leucine_biosynthesis,Leucine_degradation))
Lysine_metabolism <- unique(rbind(Lysine_biosynthesis,Lysine_degradation))
Methionine_metabolism <- unique(rbind(Methionine_biosynthesis,Methionine_degradation))
Phenylalanine_metabolism <- unique(rbind(Phenylalanine_biosynthesis,Phenylalanine_degradation))
Threonine_metabolism <- unique(rbind(Threonine_biosynthesis,Threonine_degradation))
Tryptophan_metabolism <- unique(rbind(Tryptophan_biosynthesis,Tryptophan_degradation))
Tyrosine_metabolism <- unique(rbind(Tyrosine_biosynthesis,Tyrosine_degradation))
Valine_metabolism <- unique(rbind(Valine_biosynthesis,Valine_degradation))

# Perform statistical testing - aerobic data
# Load data
input_data_aerobic <- read.csv('glu_proteomics_aerobic.csv',header = TRUE,sep = ',',dec = '.')[,c(1,6,7,4,5,2,3)]

nr_rows_aerobic <- nrow(input_data_aerobic)
all_aerobic <- as.data.frame(matrix(0,nrow=nr_rows_aerobic,ncol=6))

# Define from which column in the input-data set that different names related to the examined protein should be taken
accession_variable <- 1
prot_name_variable <- 7
gene_name_variable <- 6

for (i in c(1:nr_rows_aerobic)){
  
  # T-test the difference in abundancy of a protein between minimal and rich, in aerobic
  ttesting_row_aerobic <- t.test(c(as.numeric(log2(input_data_aerobic[i,2:3]))),c(as.numeric(log2(input_data_aerobic[i,4:5]))),alternative='t',paired=FALSE)
  
  # Calculate the log2(FC)
  rich_mean_aerobic <- mean(c(as.numeric(input_data_aerobic[i,2]),as.numeric(input_data_aerobic[i,3])))
  min_mean_aerobic <- mean(c(as.numeric(input_data_aerobic[i,4]),as.numeric(input_data_aerobic[i,5])))
  ratio_rich_min_aerobic <- rich_mean_aerobic/min_mean_aerobic
  log2_fold_change_aerobic <- log2(ratio_rich_min_aerobic)
  
  # Collect the average mass of the protein, in minimal medium samples
  mass_average_min_aerobic <- mean(c(as.numeric(input_data_aerobic[i,4])*100,as.numeric(input_data_aerobic[i,5])*100))
  
  # Collect the information computed above
  all_aerobic[i,1] <- as.character(input_data_aerobic[i,accession_variable])
  all_aerobic[i,2] <- as.character(input_data_aerobic[i,prot_name_variable])
  all_aerobic[i,3] <- as.character(input_data_aerobic[i,gene_name_variable])
  all_aerobic[i,4] <- log2_fold_change_aerobic
  all_aerobic[i,5] <- ttesting_row_aerobic$p.value
  all_aerobic[i,6] <- mass_average_min_aerobic
  
}

colnames(all_aerobic) <- c('Accession','Protein','Gene','log2FC','P_value','Average_mass_minimal')

upreg_aerobic <- all_aerobic[all_aerobic$log2FC > 1 & all_aerobic$P_value < 0.05,]
downreg_aerobic <- all_aerobic[all_aerobic$log2FC < (-1) & all_aerobic$P_value < 0.05,]

# Get proteins with significantly regulated expression
Arginine_signUp <- merge(upreg_aerobic,Arginine_metabolism)
Arginine_signDown <- merge(downreg_aerobic,Arginine_metabolism)
Aspartate_signUp <- merge(upreg_aerobic,Aspartate_metabolism)
Aspartate_signDown <- merge(downreg_aerobic,Aspartate_metabolism)
Glutamate_signUp <- merge(upreg_aerobic,Glutamate_metabolism)
Glutamate_signDown <- merge(downreg_aerobic,Glutamate_metabolism)
Glycine_signUp <- merge(upreg_aerobic,Glycine_metabolism)
Glycine_signDown <- merge(downreg_aerobic,Glycine_metabolism)
Histidine_signUp <- merge(upreg_aerobic,Histidine_metabolism)
Histidine_signDown <- merge(downreg_aerobic,Histidine_metabolism)
Isoleucine_signUp <- merge(upreg_aerobic,Isoleucine_metabolism)
Isoleucine_signDown <- merge(downreg_aerobic,Isoleucine_metabolism)
Leucine_signUp <- merge(upreg_aerobic,Leucine_metabolism)
Leucine_signDown <- merge(downreg_aerobic,Leucine_metabolism)
Lysine_signUp <- merge(upreg_aerobic,Lysine_metabolism)
Lysine_signDown <- merge(downreg_aerobic,Lysine_metabolism)
Methionine_signUp <- merge(upreg_aerobic,Methionine_metabolism)
Methionine_signDown <- merge(downreg_aerobic,Methionine_metabolism)
Phenylalanine_signUp <- merge(upreg_aerobic,Phenylalanine_metabolism)
Phenylalanine_signDown <- merge(downreg_aerobic,Phenylalanine_metabolism)
Threonine_signUp <- merge(upreg_aerobic,Threonine_metabolism)
Threonine_signDown <- merge(downreg_aerobic,Threonine_metabolism)
Tryptophan_signUp <- merge(upreg_aerobic,Tryptophan_metabolism)
Tryptophan_signDown <- merge(downreg_aerobic,Tryptophan_metabolism)
Tyrosine_signUp <- merge(upreg_aerobic,Tyrosine_metabolism)
Tyrosine_signDown <- merge(downreg_aerobic,Tyrosine_metabolism)
Valine_signUp <- merge(upreg_aerobic,Valine_metabolism)
Valine_signDown <- merge(downreg_aerobic,Valine_metabolism)

# Create dataset for plotting
aminoAcid <-  c(rep('Arg',2), rep('Asp',2), rep('Glu',2), rep('Gly',2), rep('His',2),
              rep('Ile',2), rep('Leu',2), rep('Lys',2), rep('Met',2), rep('Phe',2),
              rep('Thr',2), rep('Trp',2), rep('Tyr',2), rep('Val',2))
condition <- rep(c('UP','DOWN'),14)
values <- c(nrow(Arginine_signUp),nrow(Arginine_signDown),nrow(Aspartate_signUp),nrow(Aspartate_signDown),
            nrow(Glutamate_signUp),nrow(Glutamate_signDown),nrow(Glycine_signUp),nrow(Glycine_signDown),
            nrow(Histidine_signUp),nrow(Histidine_signDown),nrow(Isoleucine_signUp),nrow(Isoleucine_signDown),
            nrow(Leucine_signUp),nrow(Leucine_signDown),nrow(Lysine_signUp),nrow(Lysine_signDown),
            nrow(Methionine_signUp),nrow(Methionine_signDown),nrow(Phenylalanine_signUp),nrow(Phenylalanine_signDown),
            nrow(Threonine_signUp),nrow(Threonine_signDown),nrow(Tryptophan_signUp),nrow(Tryptophan_signDown),
            nrow(Tyrosine_signUp),nrow(Tyrosine_signDown),nrow(Valine_signUp),nrow(Valine_signDown))
AAmetabolism_aerobic <- data.frame(aminoAcid,condition,values)
AAmetabolism_aerobic[,4] <- as.numeric(c(9,9,5,5,6,6,13,13,7,7,12,12,10,10,8,8,14,14,1,1,4,4,3,3,2,2,11,11))
colnames(AAmetabolism_aerobic)[4] <- c('Order')

ggplot(AAmetabolism_aerobic, aes(fill=condition, y=values, x=reorder(aminoAcid, Order))) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  scale_y_continuous(breaks=c(0,5,10,15,20),limits=c(0,20)) +
  scale_fill_manual(values = c(rgb(26,35,104, max=255),rgb(123,19,21,max=255)))

# Perform statistical testing - anaerobic data
# Load data
input_data_anaerobic <- read.csv('glu_proteomics_anaerobic.csv',header = TRUE,sep = ',',dec = '.')[,c(1,7,8,9,4,5,6,2,3)]

nr_rows_anaerobic <- nrow(input_data_anaerobic)
all_anaerobic <- as.data.frame(matrix(0,nrow=nr_rows_anaerobic,ncol=6))


# Define from which column in the input-data set that different names related to the examined protein should be taken
accession_variable <- 1
prot_name_variable <- 9
gene_name_variable <- 8

for (i in c(1:nr_rows_anaerobic)){
  
  # T-test the difference in abundancy of a protein between minimal and rich, in anaerobic
  ttesting_row_anaerobic <- t.test(c(as.numeric(log2(input_data_anaerobic[i,2:4]))),c(as.numeric(log2(input_data_anaerobic[i,5:7]))),alternative='t',paired=FALSE)
  
  # Calculate the log2(FC)
  rich_mean_anaerobic <- mean(c(as.numeric(input_data_anaerobic[i,2]),as.numeric(input_data_anaerobic[i,3]),as.numeric(input_data_anaerobic[i,4])))
  min_mean_anaerobic <- mean(c(as.numeric(input_data_anaerobic[i,5]),as.numeric(input_data_anaerobic[i,6]),as.numeric(input_data_anaerobic[i,7])))
  ratio_rich_min_anaerobic <- rich_mean_anaerobic/min_mean_anaerobic
  log2_fold_change_anaerobic <- log2(ratio_rich_min_anaerobic)
  
  # Collect the average mass of the protein, in minimal medium samples
  mass_average_min_anaerobic <- mean(c(as.numeric(input_data_anaerobic[i,5])*100,as.numeric(input_data_anaerobic[i,6])*100,as.numeric(input_data_anaerobic[i,7])*100))
  
  # Collect the information computed above
  all_anaerobic[i,1] <- as.character(input_data_anaerobic[i,accession_variable])
  all_anaerobic[i,2] <- as.character(input_data_anaerobic[i,prot_name_variable])
  all_anaerobic[i,3] <- as.character(input_data_anaerobic[i,gene_name_variable])
  all_anaerobic[i,4] <- log2_fold_change_anaerobic
  all_anaerobic[i,5] <- ttesting_row_anaerobic$p.value
  all_anaerobic[i,6] <- mass_average_min_anaerobic
  
}

colnames(all_anaerobic) <- c('Accession','Protein','Gene','log2FC','P_value','Average_mass_minimal')

# Perform FDR correction using the Benjamini-Hochberg method
adj_pvalues_anaerobic <- p.adjust(all_anaerobic$P_value, method = 'fdr', n = nrow(all_anaerobic))
all_anaerobic <- all_anaerobic %>%
  mutate(adj_p_values = adj_pvalues_anaerobic)

upreg_anaerobic <- all_anaerobic[all_anaerobic$log2FC > 1 & all_anaerobic$adj_p_values < 0.05,]
downreg_anaerobic <- all_anaerobic[all_anaerobic$log2FC < (-1) & all_anaerobic$adj_p_values < 0.05,]

# Get proteins with significantly regulated expression
Arginine_signUp <- merge(upreg_anaerobic,Arginine_metabolism)
Arginine_signDown <- merge(downreg_anaerobic,Arginine_metabolism)
Aspartate_signUp <- merge(upreg_anaerobic,Aspartate_metabolism)
Aspartate_signDown <- merge(downreg_anaerobic,Aspartate_metabolism)
Glutamate_signUp <- merge(upreg_anaerobic,Glutamate_metabolism)
Glutamate_signDown <- merge(downreg_anaerobic,Glutamate_metabolism)
Glycine_signUp <- merge(upreg_anaerobic,Glycine_metabolism)
Glycine_signDown <- merge(downreg_anaerobic,Glycine_metabolism)
Histidine_signUp <- merge(upreg_anaerobic,Histidine_metabolism)
Histidine_signDown <- merge(downreg_anaerobic,Histidine_metabolism)
Isoleucine_signUp <- merge(upreg_anaerobic,Isoleucine_metabolism)
Isoleucine_signDown <- merge(downreg_anaerobic,Isoleucine_metabolism)
Leucine_signUp <- merge(upreg_anaerobic,Leucine_metabolism)
Leucine_signDown <- merge(downreg_anaerobic,Leucine_metabolism)
Lysine_signUp <- merge(upreg_anaerobic,Lysine_metabolism)
Lysine_signDown <- merge(downreg_anaerobic,Lysine_metabolism)
Methionine_signUp <- merge(upreg_anaerobic,Methionine_metabolism)
Methionine_signDown <- merge(downreg_anaerobic,Methionine_metabolism)
Phenylalanine_signUp <- merge(upreg_anaerobic,Phenylalanine_metabolism)
Phenylalanine_signDown <- merge(downreg_anaerobic,Phenylalanine_metabolism)
Threonine_signUp <- merge(upreg_anaerobic,Threonine_metabolism)
Threonine_signDown <- merge(downreg_anaerobic,Threonine_metabolism)
Tryptophan_signUp <- merge(upreg_anaerobic,Tryptophan_metabolism)
Tryptophan_signDown <- merge(downreg_anaerobic,Tryptophan_metabolism)
Tyrosine_signUp <- merge(upreg_aerobic,Tyrosine_metabolism)
Tyrosine_signDown <- merge(downreg_anaerobic,Tyrosine_metabolism)
Valine_signUp <- merge(upreg_anaerobic,Valine_metabolism)
Valine_signDown <- merge(downreg_anaerobic,Valine_metabolism)

# Create dataset for plotting
aminoAcid <-  c(rep('Arg',2), rep('Asp',2), rep('Glu',2), rep('Gly',2), rep('His',2),
                rep('Ile',2), rep('Leu',2), rep('Lys',2), rep('Met',2), rep('Phe',2),
                rep('Thr',2), rep('Trp',2), rep('Tyr',2), rep('Val',2))
condition <- rep(c('UP','DOWN'),14)
values <- c(nrow(Arginine_signUp),nrow(Arginine_signDown),nrow(Aspartate_signUp),nrow(Aspartate_signDown),
            nrow(Glutamate_signUp),nrow(Glutamate_signDown),nrow(Glycine_signUp),nrow(Glycine_signDown),
            nrow(Histidine_signUp),nrow(Histidine_signDown),nrow(Isoleucine_signUp),nrow(Isoleucine_signDown),
            nrow(Leucine_signUp),nrow(Leucine_signDown),nrow(Lysine_signUp),nrow(Lysine_signDown),
            nrow(Methionine_signUp),nrow(Methionine_signDown),nrow(Phenylalanine_signUp),nrow(Phenylalanine_signDown),
            nrow(Threonine_signUp),nrow(Threonine_signDown),nrow(Tryptophan_signUp),nrow(Tryptophan_signDown),
            nrow(Tyrosine_signUp),nrow(Tyrosine_signDown),nrow(Valine_signUp),nrow(Valine_signDown))
AAmetabolism_anaerobic <- data.frame(aminoAcid,condition,values)
AAmetabolism_anaerobic[,4] <- as.numeric(c(9,9,5,5,6,6,13,13,7,7,12,12,10,10,8,8,14,14,1,1,4,4,3,3,2,2,11,11))
colnames(AAmetabolism_anaerobic)[4] <- c('Order')

ggplot(AAmetabolism_anaerobic, aes(fill=condition, y=values, x=reorder(aminoAcid, Order))) + 
  geom_bar(position="stack", stat="identity") +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  scale_y_continuous(breaks=c(0,5,10,15,20),limits=c(0,20)) +
  scale_fill_manual(values = c(rgb(26,35,104, max=255),rgb(123,19,21,max=255)))







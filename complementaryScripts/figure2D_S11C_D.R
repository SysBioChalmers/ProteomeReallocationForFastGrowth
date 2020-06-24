# This script will filter individual proteomic data samples through matching to specified proteins,that have been deduced
# from SGD pathways, to be part in the biosynthesis of supplemented amino acids

# Create data-sets for the individual sample conditions, from the supporting data-set.
# Input the data according to the colnames specified sirectly below, one protein per row in the data.frame, for aerobic samples only use the 2 replicates for rich and min respectively
# Example; Accession Rich_repl_1_mass_ratio  Rich_repl_2_mass_ratio Rich_repl_3_mass_ratio Min_repl_1_mass_ratio Min_repl_2_mass_ratio Min_repl_3_mass_ratio Gene_Name Protein_Name
# Example; PXXXXX          0.0056567               0.005343                0.004943           0.005143                  0.005243              0.005003        YXXXXXX    ABC1    

# Load data
input_data_aerobic_min <-  read.csv('glu_proteomics_aerobic_min.csv',header=TRUE,sep = ',',dec = '.')[,c(1,2,3,4,5)]
input_data_aerobic_rich <- read.csv('glu_proteomics_aerobic_rich.csv',header=TRUE,sep = ',',dec = '.')[,c(1,2,3,4,5)]
input_data_anaerobic_min <- read.csv('glu_proteomics_anaerobic_min.csv',header=TRUE,sep = ',',dec = '.')[,c(1,2,3,4,5,6)]
input_data_anaerobic_rich <- read.csv('glu_proteomics_anaerobic_rich.csv',header=TRUE,sep = ',',dec = '.')[,c(1,2,3,4,5,6)]

# Set-up data frames of proteins relating to each amino acid biosynthesis
# Can be found in supporting data-set
#
Arginine_proteins <- data.frame(c('ARG2','ARG7','ARG5,6','ARG8','ARG3','ARG1','ARG4','CPA1','CPA2'))
colnames(Arginine_proteins) <- c('Protein')
Aspartate_proteins <- data.frame(c('AAT1','AAT2'))
colnames(Aspartate_proteins) <- c('Protein')
Glutamate_proteins <- data.frame(c('GLN1','GLT1','GDH1','GDH3'))
colnames(Glutamate_proteins) <- c('Protein')
Glycine_proteins <- data.frame(c('AGX1','GLY1','SHM1','SHM2'))
colnames(Glycine_proteins) <- c('Protein')
Histidine_proteins <- data.frame(c('HIS1','HIS2','HIS3','HIS4','HIS5','HIS6','HIS7'))
colnames(Histidine_proteins) <- c('Protein')
Isoleucine_proteins <- data.frame(c('ILV1','ILV2','ILV3','ILV5','ILV6','BAT1','BAT2'))
colnames(Isoleucine_proteins) <- c('Protein')
Leucine_proteins <- data.frame(c('LEU1','LEU2','LEU4','LEU9','BAT1','BAT2'))
colnames(Leucine_proteins) <- c('Protein')
Lysine_proteins <- data.frame(c('LYS1','LYS2','LYS4','LYS9','LYS12','LYS20','LYS21'))
colnames(Lysine_proteins) <- c('Protein')
Methionine_proteins <- data.frame(c('MET2','MET3','MET5','MET6','MET7','MET10','MET12','MET13','MET14','MET16','MET17','MET22'))
colnames(Methionine_proteins) <- c('Protein')
Phenylalanine_proteins <- data.frame(c('ARO7','ARO8','ARO9','PHA2'))
colnames(Phenylalanine_proteins) <- c('Protein')
Threonine_proteins <- data.frame(c('HOM2','HOM3','HOM6','THR1','THR4'))
colnames(Threonine_proteins) <- c('Protein')
Tryptophan_proteins <- data.frame(c('TRP1','TRP2','TRP3','TRP4','TRP5'))
colnames(Tryptophan_proteins) <- c('Protein')
Tyrosine_proteins <- data.frame(c('ARO7','ARO8','ARO9','TYR1'))
colnames(Tyrosine_proteins) <- c('Protein')
Valine_proteins <- data.frame(c('ILV2','ILV3','ILV5','ILV6','BAT1','BAT2'))
colnames(Valine_proteins) <- c('Protein')

# Sub-set data-sets of proteins to the identified enzymes in amino acid biosynthesis (above)

# Arginine
Arginine_matches_min_anaerobic <- merge(Arginine_proteins,input_data_anaerobic_min,by='Protein',all.x=FALSE)
Arginine_matches_rich_anaerobic <- merge(Arginine_proteins,input_data_anaerobic_rich,by='Protein',all.x=FALSE)
Arginine_matches_min_aerobic <- merge(Arginine_proteins,input_data_aerobic_min,by='Protein',all.x=FALSE)
Arginine_matches_rich_aerobic <- merge(Arginine_proteins,input_data_aerobic_rich,by='Protein',all.x=FALSE)
# Aspartate
Aspartate_matches_min_anaerobic <- merge(Aspartate_proteins,input_data_anaerobic_min,by='Protein',all.x=FALSE)
Aspartate_matches_rich_anaerobic <- merge(Aspartate_proteins,input_data_anaerobic_rich,by='Protein',all.x=FALSE)
Aspartate_matches_min_aerobic <- merge(Aspartate_proteins,input_data_aerobic_min,by='Protein',all.x=FALSE)
Aspartate_matches_rich_aerobic <- merge(Aspartate_proteins,input_data_aerobic_rich,by='Protein',all.x=FALSE)
# Glutamate
Glutamate_matches_min_anaerobic <- merge(Glutamate_proteins,input_data_anaerobic_min,by='Protein',all.x=FALSE)
Glutamate_matches_rich_anaerobic <- merge(Glutamate_proteins,input_data_anaerobic_rich,by='Protein',all.x=FALSE)
Glutamate_matches_min_aerobic <- merge(Glutamate_proteins,input_data_aerobic_min,by='Protein',all.x=FALSE)
Glutamate_matches_rich_aerobic <- merge(Glutamate_proteins,input_data_aerobic_rich,by='Protein',all.x=FALSE)
# Glycine
Glycine_matches_min_anaerobic <- merge(Glycine_proteins,input_data_anaerobic_min,by='Protein',all.x=FALSE)
Glycine_matches_rich_anaerobic <- merge(Glycine_proteins,input_data_anaerobic_rich,by='Protein',all.x=FALSE)
Glycine_matches_min_aerobic <- merge(Glycine_proteins,input_data_aerobic_min,by='Protein',all.x=FALSE)
Glycine_matches_rich_aerobic <- merge(Glycine_proteins,input_data_aerobic_rich,by='Protein',all.x=FALSE)
# Histidine
Histidine_matches_min_anaerobic <- merge(Histidine_proteins,input_data_anaerobic_min,by='Protein',all.x=FALSE)
Histidine_matches_rich_anaerobic <- merge(Histidine_proteins,input_data_anaerobic_rich,by='Protein',all.x=FALSE)
Histidine_matches_min_aerobic <- merge(Histidine_proteins,input_data_aerobic_min,by='Protein',all.x=FALSE)
Histidine_matches_rich_aerobic <- merge(Histidine_proteins,input_data_aerobic_rich,by='Protein',all.x=FALSE)
# Isoleucine
Isoleucine_matches_min_anaerobic <- merge(Isoleucine_proteins,input_data_anaerobic_min,by='Protein',all.x=FALSE)
Isoleucine_matches_rich_anaerobic <- merge(Isoleucine_proteins,input_data_anaerobic_rich,by='Protein',all.x=FALSE)
Isoleucine_matches_min_aerobic <- merge(Isoleucine_proteins,input_data_aerobic_min,by='Protein',all.x=FALSE)
Isoleucine_matches_rich_aerobic <- merge(Isoleucine_proteins,input_data_aerobic_rich,by='Protein',all.x=FALSE)
# Leucine
Leucine_matches_min_anaerobic <- merge(Leucine_proteins,input_data_anaerobic_min,by='Protein',all.x=FALSE)
Leucine_matches_rich_anaerobic <- merge(Leucine_proteins,input_data_anaerobic_rich,by='Protein',all.x=FALSE)
Leucine_matches_min_aerobic <- merge(Leucine_proteins,input_data_aerobic_min,by='Protein',all.x=FALSE)
Leucine_matches_rich_aerobic <- merge(Leucine_proteins,input_data_aerobic_rich,by='Protein',all.x=FALSE)
# Lysine
Lysine_matches_min_anaerobic <- merge(Lysine_proteins,input_data_anaerobic_min,by='Protein',all.x=FALSE)
Lysine_matches_rich_anaerobic <- merge(Lysine_proteins,input_data_anaerobic_rich,by='Protein',all.x=FALSE)
Lysine_matches_min_aerobic <- merge(Lysine_proteins,input_data_aerobic_min,by='Protein',all.x=FALSE)
Lysine_matches_rich_aerobic <- merge(Lysine_proteins,input_data_aerobic_rich,by='Protein',all.x=FALSE)
# Methionine
Methionine_matches_min_anaerobic <- merge(Methionine_proteins,input_data_anaerobic_min,by='Protein',all.x=FALSE)
Methionine_matches_rich_anaerobic <- merge(Methionine_proteins,input_data_anaerobic_rich,by='Protein',all.x=FALSE)
Methionine_matches_min_aerobic <- merge(Methionine_proteins,input_data_aerobic_min,by='Protein',all.x=FALSE)
Methionine_matches_rich_aerobic <- merge(Methionine_proteins,input_data_aerobic_rich,by='Protein',all.x=FALSE)
# Phenylalanine
Phenylalanine_matches_min_anaerobic <- merge(Phenylalanine_proteins,input_data_anaerobic_min,by='Protein',all.x=FALSE)
Phenylalanine_matches_rich_anaerobic <- merge(Phenylalanine_proteins,input_data_anaerobic_rich,by='Protein',all.x=FALSE)
Phenylalanine_matches_min_aerobic <- merge(Phenylalanine_proteins,input_data_aerobic_min,by='Protein',all.x=FALSE)
Phenylalanine_matches_rich_aerobic <- merge(Phenylalanine_proteins,input_data_aerobic_rich,by='Protein',all.x=FALSE)
# Threonine
Threonine_matches_min_anaerobic <- merge(Threonine_proteins,input_data_anaerobic_min,by='Protein',all.x=FALSE)
Threonine_matches_rich_anaerobic <- merge(Threonine_proteins,input_data_anaerobic_rich,by='Protein',all.x=FALSE)
Threonine_matches_min_aerobic <- merge(Threonine_proteins,input_data_aerobic_min,by='Protein',all.x=FALSE)
Threonine_matches_rich_aerobic <- merge(Threonine_proteins,input_data_aerobic_rich,by='Protein',all.x=FALSE)
# Tryptophan
Tryptophan_matches_min_anaerobic <- merge(Tryptophan_proteins,input_data_anaerobic_min,by='Protein',all.x=FALSE)
Tryptophan_matches_rich_anaerobic <- merge(Tryptophan_proteins,input_data_anaerobic_rich,by='Protein',all.x=FALSE)
Tryptophan_matches_min_aerobic <- merge(Tryptophan_proteins,input_data_aerobic_min,by='Protein',all.x=FALSE)
Tryptophan_matches_rich_aerobic <- merge(Tryptophan_proteins,input_data_aerobic_rich,by='Protein',all.x=FALSE)
# Tyrosine
Tyrosine_matches_min_anaerobic <- merge(Tyrosine_proteins,input_data_anaerobic_min,by='Protein',all.x=FALSE)
Tyrosine_matches_rich_anaerobic <- merge(Tyrosine_proteins,input_data_anaerobic_rich,by='Protein',all.x=FALSE)
Tyrosine_matches_min_aerobic <- merge(Tyrosine_proteins,input_data_aerobic_min,by='Protein',all.x=FALSE)
Tyrosine_matches_rich_aerobic <- merge(Tyrosine_proteins,input_data_aerobic_rich,by='Protein',all.x=FALSE)
# Valine
Valine_matches_min_anaerobic <- merge(Valine_proteins,input_data_anaerobic_min,by='Protein',all.x=FALSE)
Valine_matches_rich_anaerobic <- merge(Valine_proteins,input_data_anaerobic_rich,by='Protein',all.x=FALSE)
Valine_matches_min_aerobic <- merge(Valine_proteins,input_data_aerobic_min,by='Protein',all.x=FALSE)
Valine_matches_rich_aerobic <- merge(Valine_proteins,input_data_aerobic_rich,by='Protein',all.x=FALSE)

# Start to set up the data frames that will collect all information, for plotting results
# Mean sums of protein and corresponding standard deviation should be collected, together with
# name of amino acid, the sample condition name and how many proteins that was identified

Collected_results_from_aerobic_min <- data.frame(c(rep('Minimal_Aerobic',14)),c('Arg','Asp','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Thr','Trp','Tyr','Val'),c(rep(0,14)),c(rep(0,14)),c(rep(0,14)))
Collected_results_from_aerobic_rich <- data.frame(c(rep('Rich_Aerobic',14)),c('Arg','Asp','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Thr','Trp','Tyr','Val'),c(rep(0,14)),c(rep(0,14)),c(rep(0,14)))
Collected_results_from_anaerobic_min <- data.frame(c(rep('Minimal_Anaerobic',14)),c('Arg','Asp','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Thr','Trp','Tyr','Val'),c(rep(0,14)),c(rep(0,14)),c(rep(0,14)))
Collected_results_from_anaerobic_rich <- data.frame(c(rep('Rich_Anaerobic',14)),c('Arg','Asp','Glu','Gly','His','Ile','Leu','Lys','Met','Phe','Thr','Trp','Tyr','Val'),c(rep(0,14)),c(rep(0,14)),c(rep(0,14)))
colnames(Collected_results_from_aerobic_min) <- c('Sample','Amino_acid','Mean_summed_mass','Stdev_summed_mass','Number_of_proteins')
colnames(Collected_results_from_aerobic_rich) <- c('Sample','Amino_acid','Mean_summed_mass','Stdev_summed_mass','Number_of_proteins')
colnames(Collected_results_from_anaerobic_min) <- c('Sample','Amino_acid','Mean_summed_mass','Stdev_summed_mass','Number_of_proteins')
colnames(Collected_results_from_anaerobic_rich) <- c('Sample','Amino_acid','Mean_summed_mass','Stdev_summed_mass','Number_of_proteins')

# Start collecting the data

# Mean, stdev and number of identified proteins for minimal aerobic
# Mean
Collected_results_from_aerobic_min[1,3] <- mean(c(sum(Arginine_matches_min_aerobic[,3]*100),sum(Arginine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[2,3] <- mean(c(sum(Aspartate_matches_min_aerobic[,3]*100),sum(Aspartate_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[3,3] <- mean(c(sum(Glutamate_matches_min_aerobic[,3]*100),sum(Glutamate_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[4,3] <- mean(c(sum(Glycine_matches_min_aerobic[,3]*100),sum(Glycine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[5,3] <- mean(c(sum(Histidine_matches_min_aerobic[,3]*100),sum(Histidine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[6,3] <- mean(c(sum(Isoleucine_matches_min_aerobic[,3]*100),sum(Isoleucine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[7,3] <- mean(c(sum(Leucine_matches_min_aerobic[,3]*100),sum(Leucine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[8,3] <- mean(c(sum(Lysine_matches_min_aerobic[,3]*100),sum(Lysine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[9,3] <- mean(c(sum(Methionine_matches_min_aerobic[,3]*100),sum(Methionine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[10,3] <- mean(c(sum(Phenylalanine_matches_min_aerobic[,3]*100),sum(Phenylalanine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[11,3] <- mean(c(sum(Threonine_matches_min_aerobic[,3]*100),sum(Threonine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[12,3] <- mean(c(sum(Tryptophan_matches_min_aerobic[,3]*100),sum(Tryptophan_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[13,3] <- mean(c(sum(Tyrosine_matches_min_aerobic[,3]*100),sum(Tyrosine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[14,3] <- mean(c(sum(Valine_matches_min_aerobic[,3]*100),sum(Valine_matches_min_aerobic[,4]*100)))
# Stdev
Collected_results_from_aerobic_min[1,4] <- sd(c(sum(Arginine_matches_min_aerobic[,3]*100),sum(Arginine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[2,4] <- sd(c(sum(Aspartate_matches_min_aerobic[,3]*100),sum(Aspartate_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[3,4] <- sd(c(sum(Glutamate_matches_min_aerobic[,3]*100),sum(Glutamate_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[4,4] <- sd(c(sum(Glycine_matches_min_aerobic[,3]*100),sum(Glycine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[5,4] <- sd(c(sum(Histidine_matches_min_aerobic[,3]*100),sum(Histidine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[6,4] <- sd(c(sum(Isoleucine_matches_min_aerobic[,3]*100),sum(Isoleucine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[7,4] <- sd(c(sum(Leucine_matches_min_aerobic[,3]*100),sum(Leucine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[8,4] <- sd(c(sum(Lysine_matches_min_aerobic[,3]*100),sum(Lysine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[9,4] <- sd(c(sum(Methionine_matches_min_aerobic[,3]*100),sum(Methionine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[10,4] <- sd(c(sum(Phenylalanine_matches_min_aerobic[,3]*100),sum(Phenylalanine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[11,4] <- sd(c(sum(Threonine_matches_min_aerobic[,3]*100),sum(Threonine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[12,4] <- sd(c(sum(Tryptophan_matches_min_aerobic[,3]*100),sum(Tryptophan_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[13,4] <- sd(c(sum(Tyrosine_matches_min_aerobic[,3]*100),sum(Tyrosine_matches_min_aerobic[,4]*100)))
Collected_results_from_aerobic_min[14,4] <- sd(c(sum(Valine_matches_min_aerobic[,3]*100),sum(Valine_matches_min_aerobic[,4]*100)))
# Number of identified proteins
Collected_results_from_aerobic_min[1,5] <- nrow(Arginine_matches_min_aerobic)
Collected_results_from_aerobic_min[2,5] <- nrow(Aspartate_matches_min_aerobic)
Collected_results_from_aerobic_min[3,5] <- nrow(Glutamate_matches_min_aerobic)
Collected_results_from_aerobic_min[4,5] <- nrow(Glycine_matches_min_aerobic)
Collected_results_from_aerobic_min[5,5] <- nrow(Histidine_matches_min_aerobic)
Collected_results_from_aerobic_min[6,5] <- nrow(Isoleucine_matches_min_aerobic)
Collected_results_from_aerobic_min[7,5] <- nrow(Leucine_matches_min_aerobic)
Collected_results_from_aerobic_min[8,5] <- nrow(Lysine_matches_min_aerobic)
Collected_results_from_aerobic_min[9,5] <- nrow(Methionine_matches_min_aerobic)
Collected_results_from_aerobic_min[10,5] <- nrow(Phenylalanine_matches_min_aerobic)
Collected_results_from_aerobic_min[11,5] <- nrow(Threonine_matches_min_aerobic)
Collected_results_from_aerobic_min[12,5] <- nrow(Tryptophan_matches_min_aerobic)
Collected_results_from_aerobic_min[13,5] <- nrow(Tyrosine_matches_min_aerobic)
Collected_results_from_aerobic_min[14,5] <- nrow(Valine_matches_min_aerobic)

# Mean, stdev and number of identified proteins for rich aerobic
# Mean
Collected_results_from_aerobic_rich[1,3] <- mean(c(sum(Arginine_matches_rich_aerobic[,3]*100),sum(Arginine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[2,3] <- mean(c(sum(Aspartate_matches_rich_aerobic[,3]*100),sum(Aspartate_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[3,3] <- mean(c(sum(Glutamate_matches_rich_aerobic[,3]*100),sum(Glutamate_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[4,3] <- mean(c(sum(Glycine_matches_rich_aerobic[,3]*100),sum(Glycine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[5,3] <- mean(c(sum(Histidine_matches_rich_aerobic[,3]*100),sum(Histidine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[6,3] <- mean(c(sum(Isoleucine_matches_rich_aerobic[,3]*100),sum(Isoleucine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[7,3] <- mean(c(sum(Leucine_matches_rich_aerobic[,3]*100),sum(Leucine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[8,3] <- mean(c(sum(Lysine_matches_rich_aerobic[,3]*100),sum(Lysine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[9,3] <- mean(c(sum(Methionine_matches_rich_aerobic[,3]*100),sum(Methionine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[10,3] <- mean(c(sum(Phenylalanine_matches_rich_aerobic[,3]*100),sum(Phenylalanine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[11,3] <- mean(c(sum(Threonine_matches_rich_aerobic[,3]*100),sum(Threonine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[12,3] <- mean(c(sum(Tryptophan_matches_rich_aerobic[,3]*100),sum(Tryptophan_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[13,3] <- mean(c(sum(Tyrosine_matches_rich_aerobic[,3]*100),sum(Tyrosine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[14,3] <- mean(c(sum(Valine_matches_rich_aerobic[,3]*100),sum(Valine_matches_rich_aerobic[,4]*100)))
# Stdev
Collected_results_from_aerobic_rich[1,4] <- sd(c(sum(Arginine_matches_rich_aerobic[,3]*100),sum(Arginine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[2,4] <- sd(c(sum(Aspartate_matches_rich_aerobic[,3]*100),sum(Aspartate_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[3,4] <- sd(c(sum(Glutamate_matches_rich_aerobic[,3]*100),sum(Glutamate_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[4,4] <- sd(c(sum(Glycine_matches_rich_aerobic[,3]*100),sum(Glycine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[5,4] <- sd(c(sum(Histidine_matches_rich_aerobic[,3]*100),sum(Histidine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[6,4] <- sd(c(sum(Isoleucine_matches_rich_aerobic[,3]*100),sum(Isoleucine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[7,4] <- sd(c(sum(Leucine_matches_rich_aerobic[,3]*100),sum(Leucine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[8,4] <- sd(c(sum(Lysine_matches_rich_aerobic[,3]*100),sum(Lysine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[9,4] <- sd(c(sum(Methionine_matches_rich_aerobic[,3]*100),sum(Methionine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[10,4] <- sd(c(sum(Phenylalanine_matches_rich_aerobic[,3]*100),sum(Phenylalanine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[11,4] <- sd(c(sum(Threonine_matches_rich_aerobic[,3]*100),sum(Threonine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[12,4] <- sd(c(sum(Tryptophan_matches_rich_aerobic[,3]*100),sum(Tryptophan_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[13,4] <- sd(c(sum(Tyrosine_matches_rich_aerobic[,3]*100),sum(Tyrosine_matches_rich_aerobic[,4]*100)))
Collected_results_from_aerobic_rich[14,4] <- sd(c(sum(Valine_matches_rich_aerobic[,3]*100),sum(Valine_matches_rich_aerobic[,4]*100)))
# Number of identified proteins
Collected_results_from_aerobic_rich[1,5] <- nrow(Arginine_matches_rich_aerobic)
Collected_results_from_aerobic_rich[2,5] <- nrow(Aspartate_matches_rich_aerobic)
Collected_results_from_aerobic_rich[3,5] <- nrow(Glutamate_matches_rich_aerobic)
Collected_results_from_aerobic_rich[4,5] <- nrow(Glycine_matches_rich_aerobic)
Collected_results_from_aerobic_rich[5,5] <- nrow(Histidine_matches_rich_aerobic)
Collected_results_from_aerobic_rich[6,5] <- nrow(Isoleucine_matches_rich_aerobic)
Collected_results_from_aerobic_rich[7,5] <- nrow(Leucine_matches_rich_aerobic)
Collected_results_from_aerobic_rich[8,5] <- nrow(Lysine_matches_rich_aerobic)
Collected_results_from_aerobic_rich[9,5] <- nrow(Methionine_matches_rich_aerobic)
Collected_results_from_aerobic_rich[10,5] <- nrow(Phenylalanine_matches_rich_aerobic)
Collected_results_from_aerobic_rich[11,5] <- nrow(Threonine_matches_rich_aerobic)
Collected_results_from_aerobic_rich[12,5] <- nrow(Tryptophan_matches_rich_aerobic)
Collected_results_from_aerobic_rich[13,5] <- nrow(Tyrosine_matches_rich_aerobic)
Collected_results_from_aerobic_rich[14,5] <- nrow(Valine_matches_rich_aerobic)

# Mean, stdev and number of identified proteins for minimal anaerobic
# Mean
Collected_results_from_anaerobic_min[1,3] <- mean(c(sum(Arginine_matches_min_anaerobic[,3]*100),sum(Arginine_matches_min_anaerobic[,4]*100),sum(Arginine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[2,3] <- mean(c(sum(Aspartate_matches_min_anaerobic[,3]*100),sum(Aspartate_matches_min_anaerobic[,4]*100),sum(Aspartate_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[3,3] <- mean(c(sum(Glutamate_matches_min_anaerobic[,3]*100),sum(Glutamate_matches_min_anaerobic[,4]*100),sum(Glutamate_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[4,3] <- mean(c(sum(Glycine_matches_min_anaerobic[,3]*100),sum(Glycine_matches_min_anaerobic[,4]*100),sum(Glycine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[5,3] <- mean(c(sum(Histidine_matches_min_anaerobic[,3]*100),sum(Histidine_matches_min_anaerobic[,4]*100),sum(Histidine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[6,3] <- mean(c(sum(Isoleucine_matches_min_anaerobic[,3]*100),sum(Isoleucine_matches_min_anaerobic[,4]*100),sum(Isoleucine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[7,3] <- mean(c(sum(Leucine_matches_min_anaerobic[,3]*100),sum(Leucine_matches_min_anaerobic[,4]*100),sum(Leucine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[8,3] <- mean(c(sum(Lysine_matches_min_anaerobic[,3]*100),sum(Lysine_matches_min_anaerobic[,4]*100),sum(Lysine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[9,3] <- mean(c(sum(Methionine_matches_min_anaerobic[,3]*100),sum(Methionine_matches_min_anaerobic[,4]*100),sum(Methionine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[10,3] <- mean(c(sum(Phenylalanine_matches_min_anaerobic[,3]*100),sum(Phenylalanine_matches_min_anaerobic[,4]*100),sum(Phenylalanine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[11,3] <- mean(c(sum(Threonine_matches_min_anaerobic[,3]*100),sum(Threonine_matches_min_anaerobic[,4]*100),sum(Threonine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[12,3] <- mean(c(sum(Tryptophan_matches_min_anaerobic[,3]*100),sum(Tryptophan_matches_min_anaerobic[,4]*100),sum(Tryptophan_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[13,3] <- mean(c(sum(Tyrosine_matches_min_anaerobic[,3]*100),sum(Tyrosine_matches_min_anaerobic[,4]*100),sum(Tyrosine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[14,3] <- mean(c(sum(Valine_matches_min_anaerobic[,3]*100),sum(Valine_matches_min_anaerobic[,4]*100),sum(Valine_matches_min_anaerobic[,5]*100)))
# Stdev
Collected_results_from_anaerobic_min[1,4] <- sd(c(sum(Arginine_matches_min_anaerobic[,3]*100),sum(Arginine_matches_min_anaerobic[,4]*100),sum(Arginine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[2,4] <- sd(c(sum(Aspartate_matches_min_anaerobic[,3]*100),sum(Aspartate_matches_min_anaerobic[,4]*100),sum(Aspartate_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[3,4] <- sd(c(sum(Glutamate_matches_min_anaerobic[,3]*100),sum(Glutamate_matches_min_anaerobic[,4]*100),sum(Glutamate_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[4,4] <- sd(c(sum(Glycine_matches_min_anaerobic[,3]*100),sum(Glycine_matches_min_anaerobic[,4]*100),sum(Glycine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[5,4] <- sd(c(sum(Histidine_matches_min_anaerobic[,3]*100),sum(Histidine_matches_min_anaerobic[,4]*100),sum(Histidine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[6,4] <- sd(c(sum(Isoleucine_matches_min_anaerobic[,3]*100),sum(Isoleucine_matches_min_anaerobic[,4]*100),sum(Isoleucine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[7,4] <- sd(c(sum(Leucine_matches_min_anaerobic[,3]*100),sum(Leucine_matches_min_anaerobic[,4]*100),sum(Leucine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[8,4] <- sd(c(sum(Lysine_matches_min_anaerobic[,3]*100),sum(Lysine_matches_min_anaerobic[,4]*100),sum(Lysine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[9,4] <- sd(c(sum(Methionine_matches_min_anaerobic[,3]*100),sum(Methionine_matches_min_anaerobic[,4]*100),sum(Methionine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[10,4] <- sd(c(sum(Phenylalanine_matches_min_anaerobic[,3]*100),sum(Phenylalanine_matches_min_anaerobic[,4]*100),sum(Phenylalanine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[11,4] <- sd(c(sum(Threonine_matches_min_anaerobic[,3]*100),sum(Threonine_matches_min_anaerobic[,4]*100),sum(Threonine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[12,4] <- sd(c(sum(Tryptophan_matches_min_anaerobic[,3]*100),sum(Tryptophan_matches_min_anaerobic[,4]*100),sum(Tryptophan_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[13,4] <- sd(c(sum(Tyrosine_matches_min_anaerobic[,3]*100),sum(Tyrosine_matches_min_anaerobic[,4]*100),sum(Tyrosine_matches_min_anaerobic[,5]*100)))
Collected_results_from_anaerobic_min[14,4] <- sd(c(sum(Valine_matches_min_anaerobic[,3]*100),sum(Valine_matches_min_anaerobic[,4]*100),sum(Valine_matches_min_anaerobic[,5]*100)))
# Number of identified proteins
Collected_results_from_anaerobic_min[1,5] <- nrow(Arginine_matches_min_anaerobic)
Collected_results_from_anaerobic_min[2,5] <- nrow(Aspartate_matches_min_anaerobic)
Collected_results_from_anaerobic_min[3,5] <- nrow(Glutamate_matches_min_anaerobic)
Collected_results_from_anaerobic_min[4,5] <- nrow(Glycine_matches_min_anaerobic)
Collected_results_from_anaerobic_min[5,5] <- nrow(Histidine_matches_min_anaerobic)
Collected_results_from_anaerobic_min[6,5] <- nrow(Isoleucine_matches_min_anaerobic)
Collected_results_from_anaerobic_min[7,5] <- nrow(Leucine_matches_min_anaerobic)
Collected_results_from_anaerobic_min[8,5] <- nrow(Lysine_matches_min_anaerobic)
Collected_results_from_anaerobic_min[9,5] <- nrow(Methionine_matches_min_anaerobic)
Collected_results_from_anaerobic_min[10,5] <- nrow(Phenylalanine_matches_min_anaerobic)
Collected_results_from_anaerobic_min[11,5] <- nrow(Threonine_matches_min_anaerobic)
Collected_results_from_anaerobic_min[12,5] <- nrow(Tryptophan_matches_min_anaerobic)
Collected_results_from_anaerobic_min[13,5] <- nrow(Tyrosine_matches_min_anaerobic)
Collected_results_from_anaerobic_min[14,5] <- nrow(Valine_matches_min_anaerobic)

# Mean, stdev and number of identified proteins for rich anaerobic
# Mean
Collected_results_from_anaerobic_rich[1,3] <- mean(c(sum(Arginine_matches_rich_anaerobic[,3]*100),sum(Arginine_matches_rich_anaerobic[,4]*100),sum(Arginine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[2,3] <- mean(c(sum(Aspartate_matches_rich_anaerobic[,3]*100),sum(Aspartate_matches_rich_anaerobic[,4]*100),sum(Aspartate_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[3,3] <- mean(c(sum(Glutamate_matches_rich_anaerobic[,3]*100),sum(Glutamate_matches_rich_anaerobic[,4]*100),sum(Glutamate_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[4,3] <- mean(c(sum(Glycine_matches_rich_anaerobic[,3]*100),sum(Glycine_matches_rich_anaerobic[,4]*100),sum(Glycine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[5,3] <- mean(c(sum(Histidine_matches_rich_anaerobic[,3]*100),sum(Histidine_matches_rich_anaerobic[,4]*100),sum(Histidine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[6,3] <- mean(c(sum(Isoleucine_matches_rich_anaerobic[,3]*100),sum(Isoleucine_matches_rich_anaerobic[,4]*100),sum(Isoleucine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[7,3] <- mean(c(sum(Leucine_matches_rich_anaerobic[,3]*100),sum(Leucine_matches_rich_anaerobic[,4]*100),sum(Leucine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[8,3] <- mean(c(sum(Lysine_matches_rich_anaerobic[,3]*100),sum(Lysine_matches_rich_anaerobic[,4]*100),sum(Lysine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[9,3] <- mean(c(sum(Methionine_matches_rich_anaerobic[,3]*100),sum(Methionine_matches_rich_anaerobic[,4]*100),sum(Methionine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[10,3] <- mean(c(sum(Phenylalanine_matches_rich_anaerobic[,3]*100),sum(Phenylalanine_matches_rich_anaerobic[,4]*100),sum(Phenylalanine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[11,3] <- mean(c(sum(Threonine_matches_rich_anaerobic[,3]*100),sum(Threonine_matches_rich_anaerobic[,4]*100),sum(Threonine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[12,3] <- mean(c(sum(Tryptophan_matches_rich_anaerobic[,3]*100),sum(Tryptophan_matches_rich_anaerobic[,4]*100),sum(Tryptophan_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[13,3] <- mean(c(sum(Tyrosine_matches_rich_anaerobic[,3]*100),sum(Tyrosine_matches_rich_anaerobic[,4]*100),sum(Tyrosine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[14,3] <- mean(c(sum(Valine_matches_rich_anaerobic[,3]*100),sum(Valine_matches_rich_anaerobic[,4]*100),sum(Valine_matches_rich_anaerobic[,5]*100)))
# Stdev
Collected_results_from_anaerobic_rich[1,4] <- sd(c(sum(Arginine_matches_rich_anaerobic[,3]*100),sum(Arginine_matches_rich_anaerobic[,4]*100),sum(Arginine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[2,4] <- sd(c(sum(Aspartate_matches_rich_anaerobic[,3]*100),sum(Aspartate_matches_rich_anaerobic[,4]*100),sum(Aspartate_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[3,4] <- sd(c(sum(Glutamate_matches_rich_anaerobic[,3]*100),sum(Glutamate_matches_rich_anaerobic[,4]*100),sum(Glutamate_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[4,4] <- sd(c(sum(Glycine_matches_rich_anaerobic[,3]*100),sum(Glycine_matches_rich_anaerobic[,4]*100),sum(Glycine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[5,4] <- sd(c(sum(Histidine_matches_rich_anaerobic[,3]*100),sum(Histidine_matches_rich_anaerobic[,4]*100),sum(Histidine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[6,4] <- sd(c(sum(Isoleucine_matches_rich_anaerobic[,3]*100),sum(Isoleucine_matches_rich_anaerobic[,4]*100),sum(Isoleucine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[7,4] <- sd(c(sum(Leucine_matches_rich_anaerobic[,3]*100),sum(Leucine_matches_rich_anaerobic[,4]*100),sum(Leucine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[8,4] <- sd(c(sum(Lysine_matches_rich_anaerobic[,3]*100),sum(Lysine_matches_rich_anaerobic[,4]*100),sum(Lysine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[9,4] <- sd(c(sum(Methionine_matches_rich_anaerobic[,3]*100),sum(Methionine_matches_rich_anaerobic[,4]*100),sum(Methionine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[10,4] <- sd(c(sum(Phenylalanine_matches_rich_anaerobic[,3]*100),sum(Phenylalanine_matches_rich_anaerobic[,4]*100),sum(Phenylalanine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[11,4] <- sd(c(sum(Threonine_matches_rich_anaerobic[,3]*100),sum(Threonine_matches_rich_anaerobic[,4]*100),sum(Threonine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[12,4] <- sd(c(sum(Tryptophan_matches_rich_anaerobic[,3]*100),sum(Tryptophan_matches_rich_anaerobic[,4]*100),sum(Tryptophan_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[13,4] <- sd(c(sum(Tyrosine_matches_rich_anaerobic[,3]*100),sum(Tyrosine_matches_rich_anaerobic[,4]*100),sum(Tyrosine_matches_rich_anaerobic[,5]*100)))
Collected_results_from_anaerobic_rich[14,4] <- sd(c(sum(Valine_matches_rich_anaerobic[,3]*100),sum(Valine_matches_rich_anaerobic[,4]*100),sum(Valine_matches_rich_anaerobic[,5]*100)))
# Number of identified proteins
Collected_results_from_anaerobic_rich[1,5] <- nrow(Arginine_matches_rich_anaerobic)
Collected_results_from_anaerobic_rich[2,5] <- nrow(Aspartate_matches_rich_anaerobic)
Collected_results_from_anaerobic_rich[3,5] <- nrow(Glutamate_matches_rich_anaerobic)
Collected_results_from_anaerobic_rich[4,5] <- nrow(Glycine_matches_rich_anaerobic)
Collected_results_from_anaerobic_rich[5,5] <- nrow(Histidine_matches_rich_anaerobic)
Collected_results_from_anaerobic_rich[6,5] <- nrow(Isoleucine_matches_rich_anaerobic)
Collected_results_from_anaerobic_rich[7,5] <- nrow(Leucine_matches_rich_anaerobic)
Collected_results_from_anaerobic_rich[8,5] <- nrow(Lysine_matches_rich_anaerobic)
Collected_results_from_anaerobic_rich[9,5] <- nrow(Methionine_matches_rich_anaerobic)
Collected_results_from_anaerobic_rich[10,5] <- nrow(Phenylalanine_matches_rich_anaerobic)
Collected_results_from_anaerobic_rich[11,5] <- nrow(Threonine_matches_rich_anaerobic)
Collected_results_from_anaerobic_rich[12,5] <- nrow(Tryptophan_matches_rich_anaerobic)
Collected_results_from_anaerobic_rich[13,5] <- nrow(Tyrosine_matches_rich_anaerobic)
Collected_results_from_anaerobic_rich[14,5] <- nrow(Valine_matches_rich_anaerobic)

# Combine the data-sets, for aerobic and anaerobic respectively, before plotting the data

Collected_results_from_aerobic <- rbind(Collected_results_from_aerobic_min,Collected_results_from_aerobic_rich)
Collected_results_from_anaerobic <- rbind(Collected_results_from_anaerobic_min,Collected_results_from_anaerobic_rich)

# Before plotting we add a column that will make us decide the order of amino acids in the plot
# We organize them based on groups of amino acids, with the knowledge of sum of allocation
# between conditions, for that amino acid
Collected_results_from_aerobic[,6] <- as.numeric(c(9,5,6,13,7,12,10,8,14,1,4,3,2,11))
Collected_results_from_anaerobic[,6] <- as.numeric(c(9,5,6,13,7,12,10,8,14,1,4,3,2,11))
colnames(Collected_results_from_aerobic)[6] <- c('Order')
colnames(Collected_results_from_anaerobic)[6] <- c('Order')


# Plot results

# Important!
# The number of identified proteins for each individual sample is added manually to the y-axis title, afterwards,
# corresponding to the value of matched proteins for each individual replicate

library(ggplot2)

# Plot summed abundance for aerobic cultures
# Explanation of colors and what samples they represent
# Aerobic: brown = rich samples,  orange = minimal samples

ggplot(data=Collected_results_from_aerobic, aes(x=reorder(Amino_acid, Order), y=Mean_summed_mass, fill=Sample)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),size=0.1) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  scale_fill_manual(values = c('darkorange','darkorange4')) +
  labs(x=" ") +
  labs(y="Mass ratio of proteome") +
  labs(fill=' ') +
  scale_y_continuous(breaks=c(0,1,2,3,4,5),limits=c(0,5)) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(face = c('bold'),size=10)) +
  theme(axis.title.x = element_text(face = c('bold'),size=10)) +
  theme(axis.text.y = element_text(colour='black',face = c('bold'),size=10)) +
  theme(axis.title.y=element_text(color = "black",size=10,face=c('bold'))) +
  geom_errorbar(aes(ymin=Collected_results_from_aerobic$Mean_summed_mass-Collected_results_from_aerobic$Stdev_summed_mass, ymax=Collected_results_from_aerobic$Mean_summed_mass+Collected_results_from_aerobic$Stdev_summed_mass), width=0.4,position=position_dodge(.9),size=0.3)


# Plot summed abundance for anaerobic cultures
# Explanation of colors and what samples they represent
# Aerobic: purple = rich samples,  light blue = minimal samples
ggplot(data=Collected_results_from_anaerobic, aes(x=reorder(Amino_acid, Order), y=Mean_summed_mass, fill=Sample)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),size=0.1) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor.x = element_blank()) +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  scale_fill_manual(values = c('deepskyblue1','darkorchid4')) +
  labs(x=" ") +
  labs(y="Mass ratio of proteome") +
  labs(fill=' ') +
  scale_y_continuous(breaks=c(0,1,2,3,4,5),limits=c(0,5)) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(face = c('bold'),size=10)) +
  theme(axis.title.x = element_text(face = c('bold'),size=10)) +
  theme(axis.text.y = element_text(colour='black',face = c('bold'),size=10)) +
  theme(axis.title.y=element_text(color = "black",size=10,face=c('bold'))) +
  geom_errorbar(aes(ymin=Collected_results_from_anaerobic$Mean_summed_mass-Collected_results_from_anaerobic$Stdev_summed_mass, ymax=Collected_results_from_anaerobic$Mean_summed_mass+Collected_results_from_anaerobic$Stdev_summed_mass), width=0.4,position=position_dodge(.9),size=0.3)

# Figure S11 C-D
# We also plot the average enzyme mass-allocation within each individual amino acid
#
# After investigation we see that we can actually divide the mean and standard deviations already obtained
# by the number of identified proteins, for each amino acid and sample, and get the same results as if
# we calculate the mean and standard deviation, per protein, based on sums of mass normalized per protein

Collected_results_from_aerobic_normalized_per_protein <- Collected_results_from_aerobic
Collected_results_from_anaerobic_normalized_per_protein <- Collected_results_from_anaerobic
Collected_results_from_aerobic_normalized_per_protein[,3] <- Collected_results_from_aerobic_normalized_per_protein[,3]/Collected_results_from_aerobic_normalized_per_protein[,5]
Collected_results_from_aerobic_normalized_per_protein[,4] <- Collected_results_from_aerobic_normalized_per_protein[,4]/Collected_results_from_aerobic_normalized_per_protein[,5]
Collected_results_from_anaerobic_normalized_per_protein[,3] <- Collected_results_from_anaerobic_normalized_per_protein[,3]/Collected_results_from_anaerobic_normalized_per_protein[,5]
Collected_results_from_anaerobic_normalized_per_protein[,4] <- Collected_results_from_anaerobic_normalized_per_protein[,4]/Collected_results_from_anaerobic_normalized_per_protein[,5]


# Plot average abundance for aerobic cultures
# Explanation of colors and what samples they represent
# Aerobic: brown = rich samples,  orange = minimal samples
ggplot(data=Collected_results_from_aerobic_normalized_per_protein, aes(x=reorder(Amino_acid, Mean_summed_mass), y=Mean_summed_mass, fill=Sample)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),size=0.1) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  scale_fill_manual(values = c('darkorange','darkorange4')) +
  labs(x=" ") +
  labs(y="Mass ratio of proteome") +
  labs(fill=' ') +
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6),limits=c(0,0.6)) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(face = c('bold'),size=10)) +
  theme(axis.title.x = element_text(face = c('bold'),size=10)) +
  theme(axis.text.y = element_text(colour='black',face = c('bold'),size=10)) +
  theme(axis.title.y=element_text(color = "black",size=10,face=c('bold'))) +
  geom_errorbar(aes(ymin=Collected_results_from_aerobic_normalized_per_protein$Mean_summed_mass-Collected_results_from_aerobic_normalized_per_protein$Stdev_summed_mass, ymax=Collected_results_from_aerobic_normalized_per_protein$Mean_summed_mass+Collected_results_from_aerobic_normalized_per_protein$Stdev_summed_mass), width=0.4,position=position_dodge(.9),size=0.3)

# Plot average abundance for anaerobic cultures
# Explanation of colors and what samples they represent
# Aerobic: purple = rich samples,  light blue = minimal samples
ggplot(data=Collected_results_from_anaerobic_normalized_per_protein, aes(x=reorder(Amino_acid, Mean_summed_mass), y=Mean_summed_mass, fill=Sample)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),size=0.1) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  scale_fill_manual(values = c('deepskyblue1','darkorchid4')) +
  labs(x=" ") +
  labs(y="Mass ratio of proteome") +
  labs(fill=' ') +
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6),limits=c(0,0.6)) +
  theme(legend.position = 'none') +
  theme(axis.text.x = element_text(face = c('bold'),size=10)) +
  theme(axis.title.x = element_text(face = c('bold'),size=10)) +
  theme(axis.text.y = element_text(colour='black',face = c('bold'),size=10)) +
  theme(axis.title.y=element_text(color = "black",size=10,face=c('bold'))) +
  geom_errorbar(aes(ymin=Collected_results_from_anaerobic_normalized_per_protein$Mean_summed_mass-Collected_results_from_anaerobic_normalized_per_protein$Stdev_summed_mass, ymax=Collected_results_from_anaerobic_normalized_per_protein$Mean_summed_mass+Collected_results_from_anaerobic_normalized_per_protein$Stdev_summed_mass), width=0.4,position=position_dodge(.9),size=0.3)



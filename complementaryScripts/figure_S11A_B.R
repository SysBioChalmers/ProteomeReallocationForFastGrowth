library(ggplot2)

# load the file containing kcats, from enzymes in specific reactions, and those enzymes molecular weight.
# Note that some enzymes catalyze more than 1 reaction in the specific amino acid biosynthetic pathway, and occur more than 1 time.
# The focus of this analysis is to visualize reaction properties, through the responsible enzymes. Therefore, if an enzyme occurs two times
# in the pathway its molecular weight will be plotted twice.
# For identified isoenzymes, through the model, kcats are not plotted twice when the isoenzyes are reported with the same kcat-value.
# This does not add new information about the reaction. However, both molecular weights of the isoenzymes are plotted, since this adds information.

# Load values of kcat and MWs for investigated enzymes related to reactions
file <- 'kcats_AA_biosynthetic_enzymes.csv'
kcats_and_MWs <- read.csv(file,header=TRUE,sep=';',dec=',')
colnames(kcats_and_MWs) <- c('Protein','kcat','MW')

#
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


# Merge the data sets and divide into specific amino acid biosynthetic pathway
ARGININE_kcat_mw <- merge(kcats_and_MWs,Arginine_proteins,by='Protein',all.x=FALSE)
ASPARTATE_kcat_mw <- merge(kcats_and_MWs,Aspartate_proteins,by='Protein',all.x=FALSE)
GLUTAMATE_kcat_mw <- merge(kcats_and_MWs,Glutamate_proteins,by='Protein',all.x=FALSE)
GLYCINE_kcat_mw <- merge(kcats_and_MWs,Glycine_proteins,by='Protein',all.x=FALSE)
HISTIDINE_kcat_mw <- merge(kcats_and_MWs,Histidine_proteins,by='Protein',all.x=FALSE)
ISOLEUCINE_kcat_mw <- merge(kcats_and_MWs,Isoleucine_proteins,by='Protein',all.x=FALSE)
LEUCINE_kcat_mw <- merge(kcats_and_MWs,Leucine_proteins,by='Protein',all.x=FALSE)
LYSINE_kcat_mw <- merge(kcats_and_MWs,Lysine_proteins,by='Protein',all.x=FALSE)
METHIONINE_kcat_mw <- merge(kcats_and_MWs,Methionine_proteins,by='Protein',all.x=FALSE)
PHENYLALANINE_kcat_mw <- merge(kcats_and_MWs,Phenylalanine_proteins,by='Protein',all.x=FALSE)
THREONINE_kcat_mw <- merge(kcats_and_MWs,Threonine_proteins,by='Protein',all.x=FALSE)
TRYPTOPHAN_kcat_mw <- merge(kcats_and_MWs,Tryptophan_proteins,by='Protein',all.x=FALSE)
TYROSINE_kcat_mw <- merge(kcats_and_MWs,Tyrosine_proteins,by='Protein',all.x=FALSE)
VALINE_kcat_mw <- merge(kcats_and_MWs,Valine_proteins,by='Protein',all.x=FALSE)

# Get lengths of each data frame
a1 <- length(ARGININE_kcat_mw[,1])
a2 <- length(ASPARTATE_kcat_mw[,1])
a3 <- length(GLUTAMATE_kcat_mw[,1])
a4 <- length(GLYCINE_kcat_mw[,1])
a5 <- length(HISTIDINE_kcat_mw[,1])
a6 <- length(ISOLEUCINE_kcat_mw[,1])
a7 <- length(LEUCINE_kcat_mw[,1])
a8 <- length(LYSINE_kcat_mw[,1])
a9 <- length(METHIONINE_kcat_mw[,1])
a10 <- length(PHENYLALANINE_kcat_mw[,1])
a11 <- length(THREONINE_kcat_mw[,1])
a12 <- length(TRYPTOPHAN_kcat_mw[,1])
a13 <- length(TYROSINE_kcat_mw[,1])
a14 <- length(VALINE_kcat_mw[,1])

# Use lengths to add amino acid name in each data frame (to be used in grouping while plotting)
ARGININE_kcat_mw[,4] <- c(rep('Arg',a1))
ASPARTATE_kcat_mw[,4] <- c(rep('Asp',a2))
GLUTAMATE_kcat_mw[,4] <- c(rep('Glu',a3))
GLYCINE_kcat_mw[,4] <- c(rep('Gly',a4))
HISTIDINE_kcat_mw[,4] <- c(rep('His',a5))
ISOLEUCINE_kcat_mw[,4] <- c(rep('Ile',a6))
LEUCINE_kcat_mw[,4] <- c(rep('Leu',a7))
LYSINE_kcat_mw[,4] <- c(rep('Lys',a8))
METHIONINE_kcat_mw[,4] <- c(rep('Met',a9))
PHENYLALANINE_kcat_mw[,4] <- c(rep('Phe',a10))
THREONINE_kcat_mw[,4] <- c(rep('Thr',a11))
TRYPTOPHAN_kcat_mw[,4] <- c(rep('Trp',a12))
TYROSINE_kcat_mw[,4] <- c(rep('Tyr',a13))
VALINE_kcat_mw[,4] <- c(rep('Val',a14))

# Name columns
colnames(ARGININE_kcat_mw) <- c('Protein','kcat','MW','amino_acid')
colnames(ASPARTATE_kcat_mw) <- c('Protein','kcat','MW','amino_acid')
colnames(GLUTAMATE_kcat_mw) <- c('Protein','kcat','MW','amino_acid')
colnames(GLYCINE_kcat_mw) <- c('Protein','kcat','MW','amino_acid')
colnames(HISTIDINE_kcat_mw) <- c('Protein','kcat','MW','amino_acid')
colnames(ISOLEUCINE_kcat_mw) <- c('Protein','kcat','MW','amino_acid')
colnames(LEUCINE_kcat_mw) <- c('Protein','kcat','MW','amino_acid')
colnames(LYSINE_kcat_mw) <- c('Protein','kcat','MW','amino_acid')
colnames(METHIONINE_kcat_mw) <- c('Protein','kcat','MW','amino_acid')
colnames(PHENYLALANINE_kcat_mw) <- c('Protein','kcat','MW','amino_acid')
colnames(THREONINE_kcat_mw) <- c('Protein','kcat','MW','amino_acid')
colnames(TRYPTOPHAN_kcat_mw) <- c('Protein','kcat','MW','amino_acid')
colnames(TYROSINE_kcat_mw) <- c('Protein','kcat','MW','amino_acid')
colnames(VALINE_kcat_mw) <- c('Protein','kcat','MW','amino_acid')

# Merge all the information
collected_kcat_mw <- rbind(ARGININE_kcat_mw,ASPARTATE_kcat_mw,GLUTAMATE_kcat_mw,GLYCINE_kcat_mw,HISTIDINE_kcat_mw,ISOLEUCINE_kcat_mw,LEUCINE_kcat_mw,LYSINE_kcat_mw,METHIONINE_kcat_mw,PHENYLALANINE_kcat_mw,THREONINE_kcat_mw,TRYPTOPHAN_kcat_mw,TYROSINE_kcat_mw,VALINE_kcat_mw)

# For plotting of kcats, we want to remove duplicate values, so we get the kcat of individual reactions (see information in beginning)
kcat_plot <- collected_kcat_mw[-c(8,10,11,17,27,30,34,37,39,43,48,69,78,80),]

# Get the median values of both data-types
MW_median <- median(collected_kcat_mw$MW)
kcat_median <- median(kcat_plot$kcat)

# Plot kcat values
ggplot(data=kcat_plot,aes(x=amino_acid,y=log2(kcat_plot$kcat))) +
  geom_hline(yintercept = log2(kcat_median),linetype='dashed',colour='black') +
  scale_y_continuous(breaks=c(-6,-4,-2,0,2,4,6,8,10,12),limits=c(-6,12)) +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(),legend.position = 'none') +
  geom_jitter(data=kcat_plot,aes(x=amino_acid,y=log2(kcat_plot$kcat),color=kcat_plot$amino_acid),width=0.15) +
  scale_color_manual(values=c('gray71','gray68','gray65','gray62','gray59','gray56','gray53','gray50','gray47','gray44','gray41','gray38','gray35','gray32'))

# Plot MW values
ggplot(data=collected_kcat_mw,aes(x=amino_acid,y=log2(collected_kcat_mw$MW))) +
  geom_hline(yintercept = log2(MW_median),linetype='dashed',colour='black') +
  scale_y_continuous(breaks=c(4,5,6,7,8),limits=c(4,8)) +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank(),legend.position = 'none') +
  geom_jitter(data=collected_kcat_mw,aes(x=amino_acid,y=log2(collected_kcat_mw$MW),color=collected_kcat_mw$amino_acid),width=0.15) +
  scale_color_manual(values=c('gray71','gray68','gray65','gray62','gray59','gray56','gray53','gray50','gray47','gray44','gray41','gray38','gray35','gray32'))


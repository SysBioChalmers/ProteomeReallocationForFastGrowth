# 17 May 2019
# Figure S10 of "Proteome re-allocation from amino acid biosynthesis to ribosomes enables yeast to grow faster in rich media"

#setwd
setwd(dir="/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Kate")

#load libraries
library(ggplot2)
library(data.table)
library(gridExtra)


# load data
Rich_anaer<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/glu_proteomics_anaerobic_rich.csv")
Min_anaer<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/glu_proteomics_anaerobic_min.csv")
Rich_aer<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/glu_proteomics_aerobic_rich.csv")
Min_aer<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/glu_proteomics_aerobic_min.csv")

Translation<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Translation_gene_names_from_elife.csv",
                      sep=";")
Transcription<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Transcription_gene_names_from_elife.csv",
                      sep=";")
Stress<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Stress_gene_names_from_elife.csv",
                        sep=";")
Ribosomes<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Ribosomes_gene_names_from_elife.csv",
                    sep=";")
Nucleotides<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Nucleotides_gene_names_from_elife.csv",
                      sep=";")
Mitochondria<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Mitochondria_gene_names_from_elife.csv",
                       sep=";")
Lipids<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Lipids_gene_names_from_elife.csv",
                 sep=";")
Glycolysis<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Glycolysis_gene_names_from_elife.csv",
                     sep=";")
Energy<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Energy_gene_names_from_elife.csv",
                 sep=";")
ER<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/ER_gene_names_from_elife.csv",
                 sep=";")
Chaperones<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Chaperons_gene_names_from_elife.csv",
             sep=";")
AA_metabolism<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/AA_gene_names_from_elife.csv",
                        sep=";")

########################################################################################
##########################         FIGURE 4A              ##############################
########################################################################################

Rich_anaer$Type<-"Rich_anaer"
Min_anaer$Type<-"Min_anaer"
Rich_aer$Type<-"Rich_aer"
Min_aer$Type<-"Min_aer"

# Merge processes
Translation2<-Translation
Transcription2<-Transcription
Stress2<-Stress
Ribosomes2<-Ribosomes
Nucleotides2<-Nucleotides
Mitochondria2<-Mitochondria
Lipids2<-Lipids
Glycolysis2<-Glycolysis
Energy2<-Energy
ER2<-ER
Chaperones2<-Chaperones
AA_metabolism2<-AA_metabolism


Translation2$Process<-"Translation"
Transcription2$Process<-"Transcription"
Stress2$Process<-"Stress"
Ribosomes2$Process<-"Ribosomes"
Nucleotides2$Process<-"Nucleotides"
Mitochondria2$Process<-"Mitochondria"
Lipids2$Process<-"Lipids"
Glycolysis2$Process<-"Glycolysis"
Energy2$Process<-"Energy"
ER2$Process<-"ER"
Chaperones2$Process<-"Chaperones"
AA_metabolism2$Process<-"AA_metabolism"

Chaperones2<-Chaperones2[,c(1,2,4)]
AA_metabolism2<-AA_metabolism2[,c(1,2,4)]

names(Ribosomes2)[1]<-"ORF"
names(Glycolysis2)[1]<-"ORF"

eLife_all<-rbind(Translation2, Transcription2, Stress2, Ribosomes2,
                 Nucleotides2, Mitochondria2, Lipids2, Glycolysis2, Energy2,
                 ER2, Chaperones2, AA_metabolism2)
unique(eLife_all$Process)

# Rename columns to unify names used
names(eLife_all)[1:2]<-c("Gene","Accession")


######## ANAEROBIC ############################################################## 

# Combine Rich and Min for Anaerobic samples
names(Rich_anaer)[2:4]<-c("rep1","rep2","rep3")
names(Min_anaer)[2:4]<-c("rep1","rep2","rep3")
Anaer<-rbind(Rich_anaer, Min_anaer)

# Determine mean for reps
Anaer$Mean_anaer<-apply(Anaer[,2:4], 1, mean)

# Determine stdev for reps
Anaer$Stdev_anaer<-apply(Anaer[,2:4], 1, sd)

# Extract protein abundance for genese related to AA metabolisma and Translation
Anaer_subset<-merge(Anaer, eLife_all, by=c("Gene","Accession"))

# Summarise
Anaer_subset_sum_means<- as.data.frame(as.data.table(Anaer_subset)[, sum(Mean_anaer), by = .(Type, Process)])
Anaer_subset_sum_stdev<- as.data.frame(as.data.table(Anaer_subset)[, sd(Mean_anaer), by = .(Type, Process)])

names(Anaer_subset_sum_means)[3]<-"Sum"
names(Anaer_subset_sum_stdev)[3]<-"Stdev"

Anaer_subset2<-merge(Anaer_subset_sum_means, Anaer_subset_sum_stdev)

# Convert to %
Anaer_subset2$Sum<-Anaer_subset2$Sum*100
Anaer_subset2$Stdev<-Anaer_subset2$Stdev*100


# Plot for anaerobic samples 

Anaer_S10A<-ggplot(Anaer_subset2, aes(Process, Sum, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge(0.9), colour="black") +
  geom_errorbar(aes(ymin=Sum-Stdev, ymax=Sum+Stdev),
                position=position_dodge(0.9), width=0.2)+
  theme_bw(base_family= "serif") +
  theme(text = element_text(colour = "black",face="bold"),
        axis.text= element_text(colour = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title="Fig S10A",
       x="",
       y="Mass fraction of proteome (%)",
       fill="")  +
  scale_y_continuous(limits=c(0,50)) +
  scale_fill_manual(values=c("skyblue","purple"))  +
  theme(legend.position="bottom") +
  coord_flip() 


Anaer_S10A 

######## AEROBIC ############################################################## 

# Combine Rich and Min for Aerobic samples
names(Rich_aer)[2:3]<-c("rep1","rep2")
names(Min_aer)[2:3]<-c("rep1","rep2")
Aer<-rbind(Rich_aer, Min_aer)

# Dtermine mean for reps
Aer$Mean_Aer<-apply(Aer[,2:3], 1, mean)

# Determine stdev for reps
Aer$Stdev_Aer<-apply(Aer[,2:3], 1, sd)

# Extract protein abundance for genese related to AA metabolisma and Translation
Aer_subset<-merge(Aer, eLife_all, by=c("Gene","Accession"))

# Summarise
Aer_subset_sum_means<- as.data.frame(as.data.table(Aer_subset)[, sum(Mean_Aer), by = .(Type, Process)])
Aer_subset_sum_stdev<- as.data.frame(as.data.table(Aer_subset)[, sd(Mean_Aer), by = .(Type, Process)])

names(Aer_subset_sum_means)[3]<-"Sum"
names(Aer_subset_sum_stdev)[3]<-"Stdev"

Aer_subset2<-merge(Aer_subset_sum_means, Aer_subset_sum_stdev)

# Convert to %
Aer_subset2$Sum<-Aer_subset2$Sum*100
Aer_subset2$Stdev<-Aer_subset2$Stdev*100


# Plot for Aerobic samples 

Aer_S10B<-ggplot(Aer_subset2, aes(Process, Sum, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge(0.9), colour="black") +
  geom_errorbar(aes(ymin=Sum-Stdev, ymax=Sum+Stdev),
                position=position_dodge(0.9), width=0.2)+
  theme_bw(base_family= "serif") +
  theme(text = element_text(colour = "black",face="bold"),
        axis.text= element_text(colour = "black"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank()) +
  labs(title="Fig S10B",
         x="",
       y="Mass fraction of proteome (%)",
       fill="")  +
  scale_y_continuous(limits=c(0,50)) +
  scale_fill_manual(values=c("orange", "tan4")) +
  theme(legend.position="bottom") +
  coord_flip() 

Aer_S10B 

# Combine plots
grid.arrange(Anaer_S10A, Aer_S10B, nrow = 1)

# t-tests

################################ Aerobic ################################ 

# Translation
Rich_aer_Translation<-Rich_aer[which(Rich_aer$Gene %in% Translation$ORF),]
Min_aer_Translation<-Min_aer[which(Min_aer$Gene %in% Translation$ORF),]
t.test(Rich_aer_Translation[,2:3],Min_aer_Translation[,2:3]) # not significant

# Transcription
Rich_aer_Transcription<-Rich_aer[which(Rich_aer$Gene %in% Transcription$ORF),]
Min_aer_Transcription<-Min_aer[which(Min_aer$Gene %in% Transcription$ORF),]
t.test(Rich_aer_Transcription[,2:3],Min_aer_Transcription[,2:3]) # not significant

# Stress
Rich_aer_Stress<-Rich_aer[which(Rich_aer$Gene %in% Stress$ORF),]
Min_aer_Stress<-Min_aer[which(Min_aer$Gene %in% Stress$ORF),]
t.test(Rich_aer_Stress[,2:3],Min_aer_Stress[,2:3]) # not significant

# Ribosomes
Rich_aer_Ribosomes<-Rich_aer[which(Rich_aer$Gene %in% Ribosomes2$ORF),]
Min_aer_Ribosomes<-Min_aer[which(Min_aer$Gene %in% Ribosomes2$ORF),]
t.test(Rich_aer_Ribosomes[,2:3],Min_aer_Ribosomes[,2:3]) # not significant

# Nucleotides
Rich_aer_Nucleotides<-Rich_aer[which(Rich_aer$Gene %in% Nucleotides$ORF),]
Min_aer_Nucleotides<-Min_aer[which(Min_aer$Gene %in% Nucleotides$ORF),]
t.test(Rich_aer_Nucleotides[,2:3],Min_aer_Nucleotides[,2:3]) # not significant

# Mitochondria
Rich_aer_Mitochondria<-Rich_aer[which(Rich_aer$Gene %in% Mitochondria$ORF),]
Min_aer_Mitochondria<-Min_aer[which(Min_aer$Gene %in% Mitochondria$ORF),]
t.test(Rich_aer_Mitochondria[,2:3],Min_aer_Mitochondria[,2:3]) # not significant

# Lipids
Rich_aer_Lipids<-Rich_aer[which(Rich_aer$Gene %in% Lipids$ORF),]
Min_aer_Lipids<-Min_aer[which(Min_aer$Gene %in% Lipids$ORF),]
t.test(Rich_aer_Lipids[,2:3],Min_aer_Lipids[,2:3]) # not significant

# Glycolysis
Rich_aer_Glycolysis<-Rich_aer[which(Rich_aer$Gene %in% Glycolysis2$ORF),]
Min_aer_Glycolysis<-Min_aer[which(Min_aer$Gene %in% Glycolysis2$ORF),]
t.test(Rich_aer_Glycolysis[,2:3],Min_aer_Glycolysis[,2:3]) # not significant

# ER
Rich_aer_ER<-Rich_aer[which(Rich_aer$Gene %in% ER$ORF),]
Min_aer_ER<-Min_aer[which(Min_aer$Gene %in% ER$ORF),]
t.test(Rich_aer_ER[,2:3],Min_aer_ER[,2:3]) # not significant

# Energy
Rich_aer_Energy<-Rich_aer[which(Rich_aer$Gene %in% Energy$ORF),]
Min_aer_Energy<-Min_aer[which(Min_aer$Gene %in% Energy$ORF),]
t.test(Rich_aer_Energy[,2:3],Min_aer_Energy[,2:3]) # not significant

# Chaperones
Rich_aer_Chaperones<-Rich_aer[which(Rich_aer$Gene %in% Chaperones2$ORF),]
Min_aer_Chaperones<-Min_aer[which(Min_aer$Gene %in% Chaperones2$ORF),]
t.test(Rich_aer_Chaperones[,2:3],Min_aer_Chaperones[,2:3]) # not significant

# AA_metabolism
Rich_aer_AA_metabolism<-Rich_aer[which(Rich_aer$Gene %in% AA_metabolism2$ORF),]
Min_aer_AA_metabolism<-Min_aer[which(Min_aer$Gene %in% AA_metabolism2$ORF),]
t.test(Rich_aer_AA_metabolism[,2:3],Min_aer_AA_metabolism[,2:3]) # < 0.05 *


################################ Ananaerobic ################################ 

# Translation
Rich_anaer_Translation<-Rich_anaer[which(Rich_anaer$Gene %in% Translation$ORF),]
Min_anaer_Translation<-Min_anaer[which(Min_anaer$Gene %in% Translation$ORF),]
t.test(Rich_anaer_Translation[,2:3],Min_anaer_Translation[,2:3]) # not significant

# Transcription
Rich_anaer_Transcription<-Rich_anaer[which(Rich_anaer$Gene %in% Transcription$ORF),]
Min_anaer_Transcription<-Min_anaer[which(Min_anaer$Gene %in% Transcription$ORF),]
t.test(Rich_anaer_Transcription[,2:3],Min_anaer_Transcription[,2:3]) # not significant

# Stress
Rich_anaer_Stress<-Rich_anaer[which(Rich_anaer$Gene %in% Stress$ORF),]
Min_anaer_Stress<-Min_anaer[which(Min_anaer$Gene %in% Stress$ORF),]
t.test(Rich_anaer_Stress[,2:3],Min_anaer_Stress[,2:3]) # not significant

# Ribosomes
Rich_anaer_Ribosomes<-Rich_anaer[which(Rich_anaer$Gene %in% Ribosomes2$ORF),]
Min_anaer_Ribosomes<-Min_anaer[which(Min_anaer$Gene %in% Ribosomes2$ORF),]
t.test(Rich_anaer_Ribosomes[,2:3],Min_anaer_Ribosomes[,2:3]) # not significant

# Nucleotides
Rich_anaer_Nucleotides<-Rich_anaer[which(Rich_anaer$Gene %in% Nucleotides$ORF),]
Min_anaer_Nucleotides<-Min_anaer[which(Min_anaer$Gene %in% Nucleotides$ORF),]
t.test(Rich_anaer_Nucleotides[,2:3],Min_anaer_Nucleotides[,2:3]) # not significant

# Mitochondria
Rich_anaer_Mitochondria<-Rich_anaer[which(Rich_anaer$Gene %in% Mitochondria$ORF),]
Min_anaer_Mitochondria<-Min_anaer[which(Min_anaer$Gene %in% Mitochondria$ORF),]
t.test(Rich_anaer_Mitochondria[,2:3],Min_anaer_Mitochondria[,2:3]) # not significant

# Lipids
Rich_anaer_Lipids<-Rich_anaer[which(Rich_anaer$Gene %in% Lipids$ORF),]
Min_anaer_Lipids<-Min_anaer[which(Min_anaer$Gene %in% Lipids$ORF),]
t.test(Rich_anaer_Lipids[,2:3],Min_anaer_Lipids[,2:3]) # not significant

# Glycolysis
Rich_anaer_Glycolysis<-Rich_anaer[which(Rich_anaer$Gene %in% Glycolysis2$ORF),]
Min_anaer_Glycolysis<-Min_anaer[which(Min_anaer$Gene %in% Glycolysis2$ORF),]
t.test(Rich_anaer_Glycolysis[,2:3],Min_anaer_Glycolysis[,2:3]) # not significant

# ER
Rich_anaer_ER<-Rich_anaer[which(Rich_anaer$Gene %in% ER$ORF),]
Min_anaer_ER<-Min_anaer[which(Min_anaer$Gene %in% ER$ORF),]
t.test(Rich_anaer_ER[,2:3],Min_anaer_ER[,2:3]) # not significant

# Energy
Rich_anaer_Energy<-Rich_anaer[which(Rich_anaer$Gene %in% Energy$ORF),]
Min_anaer_Energy<-Min_anaer[which(Min_anaer$Gene %in% Energy$ORF),]
t.test(Rich_anaer_Energy[,2:3],Min_anaer_Energy[,2:3]) # not significant

# Chaperones
Rich_anaer_Chaperones<-Rich_anaer[which(Rich_anaer$Gene %in% Chaperones2$ORF),]
Min_anaer_Chaperones<-Min_anaer[which(Min_anaer$Gene %in% Chaperones2$ORF),]
t.test(Rich_anaer_Chaperones[,2:3],Min_anaer_Chaperones[,2:3]) # not significant

# AA_metabolism
Rich_anaer_AA_metabolism<-Rich_anaer[which(Rich_anaer$Gene %in% AA_metabolism2$ORF),]
Min_anaer_AA_metabolism<-Min_anaer[which(Min_anaer$Gene %in% AA_metabolism2$ORF),]
t.test(Rich_anaer_AA_metabolism[,2:3],Min_anaer_AA_metabolism[,2:3]) # not significant


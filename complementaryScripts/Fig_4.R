# 10 May 2019
# Figure 4 of "Proteome re-allocation from amino acid biosynthesis to ribosomes enables yeast to grow faster in rich media"

#setwd
setwd(dir="/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Kate")

#load libraries
library(ggplot2)
library(dplyr)
library(data.table)
library(grid)
library(gridExtra)
library(reshape2)

# load data
Rich_anaer<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/glu_proteomics_anaerobic_rich.csv")
Min_anaer<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/glu_proteomics_anaerobic_min.csv")
Rich_aer<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/glu_proteomics_aerobic_rich.csv")
Min_aer<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/glu_proteomics_aerobic_min.csv")
AA_metabolism<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/AA_gene_names_from_elife.csv",
                        sep=";")
Translation<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Translation_gene_names_from_elife.csv",
                      sep=";")
growth_rates <- read.csv("/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Kate/growth_rates.csv",
                         row.names = 1)
Ribosomes<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Ribosomes_gene_names_from_elife.csv",
                    sep=";")
Lipids<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Lipids_gene_names_from_elife.csv",
                 sep=";")
Glycolysis<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Glycolysis_gene_names_from_elife.csv",
                     sep=";")
Energy<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Energy_gene_names_from_elife.csv",
                 sep=";")
Nucleotides<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Nucleotides_gene_names_from_elife.csv",
                      sep=";")
Mitochondria<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Supporting_documents/Mitochondria_gene_names_from_elife.csv",
                       sep=";")

########################################################################################
##########################         FIGURE 4A              ##############################
########################################################################################

Rich_anaer$Type<-"Rich_anaer"
Min_anaer$Type<-"Min_anaer"
Rich_aer$Type<-"Rich_aer"
Min_aer$Type<-"Min_aer"

# Merge processes
AA_metabolism2<-AA_metabolism
AA_metabolism2$Process<-"AA_metabolism"
AA_metabolism3<-AA_metabolism2[,c(1,2,4)]
Translation$Process<-"Translation"

AA_metabolism_Translation<-rbind(AA_metabolism3, Translation)

# Rename columns to unify names used
names(AA_metabolism_Translation)[1:2]<-c("Gene","Accession")


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
Anaer_subset<-merge(Anaer, AA_metabolism_Translation, by=c("Gene","Accession"))

# Summarise
Anaer_subset_sum_means<- as.data.frame(as.data.table(Anaer_subset)[, sum(Mean_anaer), by = .(Type, Process)])
Anaer_subset_sum_stdev<- as.data.frame(as.data.table(Anaer_subset)[, sd(Mean_anaer), by = .(Type, Process)])

names(Anaer_subset_sum_means)[3]<-"Sum"
names(Anaer_subset_sum_stdev)[3]<-"Stdev"

Anaer_subset2<-merge(Anaer_subset_sum_means, Anaer_subset_sum_stdev)

# Convert to %
Anaer_subset2$Sum<-Anaer_subset2$Sum*100
Anaer_subset2$Stdev<-Anaer_subset2$Stdev*100

# Find differences between Min and Rich for each process

# AA metabolism
subset (Anaer_subset2$Sum, Anaer_subset2$Type=="Min_anaer" & Anaer_subset2$Process=="AA_metabolism")
subset (Anaer_subset2$Sum, Anaer_subset2$Type=="Rich_anaer" & Anaer_subset2$Process=="AA_metabolism")

(subset (Anaer_subset2$Sum, Anaer_subset2$Type=="Min_anaer" & Anaer_subset2$Process=="AA_metabolism"))-
  (subset (Anaer_subset2$Sum, Anaer_subset2$Type=="Rich_anaer" & Anaer_subset2$Process=="AA_metabolism"))  

# Translation
subset (Anaer_subset2$Sum, Anaer_subset2$Type=="Rich_anaer" & Anaer_subset2$Process=="Translation")
subset (Anaer_subset2$Sum, Anaer_subset2$Type=="Min_anaer" & Anaer_subset2$Process=="Translation")

(subset (Anaer_subset2$Sum, Anaer_subset2$Type=="Rich_anaer" & Anaer_subset2$Process=="Translation"))-
  (subset (Anaer_subset2$Sum, Anaer_subset2$Type=="Min_anaer" & Anaer_subset2$Process=="Translation"))  

# Plot for anaerobic samples 

Anaer_4A<-ggplot(Anaer_subset2, aes(Process, Sum, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge(0.9), colour="black") +
  geom_errorbar(aes(ymin=Sum-Stdev, ymax=Sum+Stdev),
                position=position_dodge(0.9), width=0.2)+
  theme_bw(base_family= "serif") +
  theme(text = element_text(colour = "black",face="bold"),
          axis.text= element_text(colour = "black"),
        panel.grid.major.y = element_blank()) +
  labs(x="",
       y="Mass fraction of proteome (%)",
       fill="")  +
  scale_x_discrete(limits=c("Translation", "AA_metabolism")) +
  scale_y_continuous(limits=c(0,50)) +
  scale_fill_manual(values=c("deepskyblue2", "deepskyblue4")) +
  coord_flip() + 
  annotate("text", x = 2.2, y = 14, label = "D 5.59%", family= "serif") +
  annotate("text", x = 0.7, y = 35, label = "D 5.42%",family= "serif") +
  theme(legend.position = c(0.8, 0.7))

Anaer_4A 

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
Aer_subset<-merge(Aer, AA_metabolism_Translation, by=c("Gene","Accession"))

# Summarise
Aer_subset_sum_means<- as.data.frame(as.data.table(Aer_subset)[, sum(Mean_Aer), by = .(Type, Process)])
Aer_subset_sum_stdev<- as.data.frame(as.data.table(Aer_subset)[, sd(Mean_Aer), by = .(Type, Process)])

names(Aer_subset_sum_means)[3]<-"Sum"
names(Aer_subset_sum_stdev)[3]<-"Stdev"

Aer_subset2<-merge(Aer_subset_sum_means, Aer_subset_sum_stdev)

# Convert to %
Aer_subset2$Sum<-Aer_subset2$Sum*100
Aer_subset2$Stdev<-Aer_subset2$Stdev*100

# Find differences between Min and Rich for each process

# AA metabolism
subset (Aer_subset2$Sum, Aer_subset2$Type=="Min_aer" & Aer_subset2$Process=="AA_metabolism")
subset (Aer_subset2$Sum, Aer_subset2$Type=="Rich_aer" & Aer_subset2$Process=="AA_metabolism")

(subset (Aer_subset2$Sum, Aer_subset2$Type=="Min_aer" & Aer_subset2$Process=="AA_metabolism"))-
  (subset (Aer_subset2$Sum, Aer_subset2$Type=="Rich_aer" & Aer_subset2$Process=="AA_metabolism"))  

# Translation
subset (Aer_subset2$Sum, Aer_subset2$Type=="Rich_aer" & Aer_subset2$Process=="Translation")
subset (Aer_subset2$Sum, Aer_subset2$Type=="Min_aer" & Aer_subset2$Process=="Translation")

(subset (Aer_subset2$Sum, Aer_subset2$Type=="Rich_aer" & Aer_subset2$Process=="Translation"))-
  (subset (Aer_subset2$Sum, Aer_subset2$Type=="Min_aer" & Aer_subset2$Process=="Translation"))  

# Plot for Aerobic samples 

Aer_4A<-ggplot(Aer_subset2, aes(Process, Sum, fill=Type)) +
  geom_bar(stat="identity", position=position_dodge(0.9), colour="black") +
  geom_errorbar(aes(ymin=Sum-Stdev, ymax=Sum+Stdev),
                position=position_dodge(0.9), width=0.2)+
  theme_bw(base_family= "serif") +
  theme(text = element_text(colour = "black",face="bold"),
        axis.text= element_text(colour = "black"),
        panel.grid.major.y = element_blank()) +
  labs(x="",
       y="Mass fraction of proteome (%)",
       fill="")  +
  scale_x_discrete(limits=c("Translation", "AA_metabolism")) +
  scale_y_continuous(limits=c(0,50)) +
  scale_fill_manual(values=c("orange", "tan4")) +
  coord_flip() + 
  annotate("text", x = 2.2, y = 14, label = "D 6.51%", family= "serif") +
  annotate("text", x = 0.7, y = 45, label = "D 4.55%",family= "serif") +
  theme(legend.position = c(0.8, 0.7))

Aer_4A 

# Combine plots
grid.arrange(Aer_4A, Anaer_4A, nrow = 2, left = textGrob("Figure 4A"))


########################################################################################
##########################         FIGURE 4B              ##############################
########################################################################################

# To determine metabolic enzymes group > genes were taken from eLife csv files:
# "Energy", "AA metabolism","Lipids","Nucleotides","Mitochondria" and "Glycolysis

head(AA_metabolism,1)
# remove "Dalton"
AA_metabolism4<-AA_metabolism[,-3]

# merge all metabolic genes together and then remove duplicates
head(AA_metabolism4,1)
head(Energy,1)
head(Lipids,1)
head(Nucleotides,1)
head(Mitochondria,1)
head(Glycolysis,1)
names(Glycolysis)[1]<-"ORF"

Metabolic_enzymes<-rbind(AA_metabolism4, Energy,Lipids,Nucleotides, Mitochondria, Glycolysis)
# remove duplicates
Metabolic_enzymes2<-Metabolic_enzymes[!duplicated(Metabolic_enzymes$ORF), ]

Metabolic_enzymes2$Process<-"Metabolic_enzymes"

head(Translation,1)
head(Metabolic_enzymes2,1)

Metabolic_enzymes_translation<-rbind(Metabolic_enzymes2, Translation)

names(Metabolic_enzymes_translation)[1:2]<-c("Gene","Accession")

######## ANAEROBIC ############################################################## 

# Extract protein abundance for genese related to AA metabolisma and Translation
Anaer_subset3<-merge(Anaer, Metabolic_enzymes_translation, by=c("Gene","Accession"))

# Summarise
Anaer_subset_sum_means2<- as.data.frame(as.data.table(Anaer_subset3)[, sum(Mean_anaer), by = .(Type, Process)])
Anaer_subset_sum_stdev2<- as.data.frame(as.data.table(Anaer_subset3)[, sd(Mean_anaer), by = .(Type, Process)])

names(Anaer_subset_sum_means2)[3]<-"Sum"
names(Anaer_subset_sum_stdev2)[3]<-"Stdev"

Anaer_subset4<-merge(Anaer_subset_sum_means2, Anaer_subset_sum_stdev2)

# Convert to %
Anaer_subset4$Sum<-Anaer_subset4$Sum*100
Anaer_subset4$Stdev<-Anaer_subset4$Stdev*100

######## AEROBIC ############################################################## 


# Extract protein abundance for genese related to AA metabolisma and Translation
Aer_subset3<-merge(Aer, Metabolic_enzymes_translation, by=c("Gene","Accession"))

# Summarise
Aer_subset_sum_means2<- as.data.frame(as.data.table(Aer_subset3)[, sum(Mean_Aer), by = .(Type, Process)])
Aer_subset_sum_stdev2<- as.data.frame(as.data.table(Aer_subset3)[, sd(Mean_Aer), by = .(Type, Process)])

names(Aer_subset_sum_means2)[3]<-"Sum"
names(Aer_subset_sum_stdev2)[3]<-"Stdev"

Aer_subset4<-merge(Aer_subset_sum_means2, Aer_subset_sum_stdev2)

# Convert to %
Aer_subset4$Sum<-Aer_subset4$Sum*100
Aer_subset4$Stdev<-Aer_subset4$Stdev*100

# Combine results

Combined<-rbind(Aer_subset4, Anaer_subset4)
split<-as.data.table(split(Combined, Combined$Process))

# Plot 

Fig_4B<-ggplot(split, aes(split$Metabolic_enzymes.Sum, split$Translation.Sum, colour=split$Metabolic_enzymes.Type)) +
  geom_point(size=4)+
  theme_bw(base_family= "serif") +
  theme(text = element_text(colour = "black",face="bold"),
        axis.text= element_text(colour = "black"),
        panel.grid.major.y = element_blank()) +
  scale_x_continuous(limits=c(0,60)) +
  scale_y_continuous(limits=c(0,50)) +
  scale_colour_manual(values=c("orange", "deepskyblue2", "tan4", "deepskyblue4")) +
  labs(x="Metabolic enzyme fraction \nof proteome (%)",
       y="Translational fraction \nof proteome (%)",
       colour="",
       title="Figure 4B") +
  geom_errorbar(aes(ymin=split$Translation.Sum-split$Translation.Stdev, 
                    ymax=split$Translation.Sum+split$Translation.Stdev),
                width=.2,colour="black") +
  geom_errorbarh(aes(xmin=split$Metabolic_enzymes.Sum-split$Metabolic_enzymes.Stdev, 
                    xmax=split$Metabolic_enzymes.Sum+split$Metabolic_enzymes.Stdev),
                 height=.2, colour="black")
  
Fig_4B

########################################################################################
##########################         FIGURE 4C left         ##############################
########################################################################################

# growth rates taken from '20200505 datasets_updated.xlsx' under 'Physiological parameters'/sheet 2
# for row starting with 'µ_cDW (1/h)' and converted to .csv file

head(Translation,1)
head(Ribosomes,1)

Ribosomes$Process<-"Ribosome"

names(Translation)[1:2]<-c("Gene","Accession")
names(Ribosomes)[1:2]<-c("Gene","Accession")

Ribosomes_translation<-rbind(Ribosomes, Translation)

######## ANAEROBIC ############################################################## 

# Extract protein abundance for genese related to AA metabolisma and Translation
Anaer_subset5<-merge(Anaer, Ribosomes_translation, by=c("Gene","Accession"))

# Summarise
Anaer_subset_sum_means3<- as.data.frame(as.data.table(Anaer_subset5)[, sum(Mean_anaer), by = .(Type, Process)])
Anaer_subset_sum_stdev3<- as.data.frame(as.data.table(Anaer_subset5)[, sd(Mean_anaer), by = .(Type, Process)])

names(Anaer_subset_sum_means3)[3]<-"Sum"
names(Anaer_subset_sum_stdev3)[3]<-"Stdev"

Anaer_subset6<-merge(Anaer_subset_sum_means3, Anaer_subset_sum_stdev3)

# Convert to %
Anaer_subset6$Sum<-Anaer_subset6$Sum*100
Anaer_subset6$Stdev<-Anaer_subset6$Stdev*100

######## AEROBIC ############################################################## 


# Extract protein abundance for genese related to AA metabolisma and Translation
Aer_subset5<-merge(Aer, Ribosomes_translation, by=c("Gene","Accession"))

# Summarise
Aer_subset_sum_means3<- as.data.frame(as.data.table(Aer_subset5)[, sum(Mean_Aer), by = .(Type, Process)])
Aer_subset_sum_stdev3<- as.data.frame(as.data.table(Aer_subset5)[, sd(Mean_Aer), by = .(Type, Process)])

names(Aer_subset_sum_means3)[3]<-"Sum"
names(Aer_subset_sum_stdev3)[3]<-"Stdev"

Aer_subset6<-merge(Aer_subset_sum_means3, Aer_subset_sum_stdev3)

# Convert to %
Aer_subset6$Sum<-Aer_subset6$Sum*100
Aer_subset6$Stdev<-Aer_subset6$Stdev*100

# Combine results

Combined2<-rbind(Aer_subset6, Anaer_subset6)
split2<-as.data.table(split(Combined2, Combined2$Process))

# merge with growth rates
split2$Ribosome.Type
growth_rates$Ribosome.Type<-c("Min_aer", "Rich_aer", "Min_anaer", "Rich_anaer")
split3<-merge(split2, growth_rates, by="Ribosome.Type")



# Plot 

cor(split3$Translation.Sum, split3$growth_rates)

Fig_4C_left<-ggplot(split3, aes(split3$growth_rates, split3$Translation.Sum, colour=split3$Translation.Type, group=1)) +
  geom_point(size=4)+
  theme_bw(base_family= "serif") +
  theme(text = element_text(colour = "black",face="bold"),
        axis.text= element_text(colour = "black"),
        panel.grid.major.y = element_blank(),
        legend.position = "none") +
 scale_x_continuous(limits=c(0,0.5)) +
  scale_y_continuous(limits=c(0,50)) +
  labs(x="Growth rate (1/h)",
       y="Translational fraction \nof proteome (%)",
      colour="",
      title="Figure 4C left \nTranslation") +
  geom_errorbar(aes(ymin=split3$Translation.Sum-split3$Translation.Stdev, 
                    ymax=split3$Translation.Sum+split3$Translation.Stdev),
                width=0,colour="black") +
  geom_smooth(method=lm,   
              se=FALSE, colour="black", size=0.5, fullrange=TRUE) +
  scale_colour_manual(values=c("orange", "deepskyblue2", "tan4", "deepskyblue4")) +
  annotate("text", x = 0.1, y = 40, label = expression(R^"2":0.931), family= "serif")

Fig_4C_left


########################################################################################
##########################         FIGURE 4C right              ##############################
########################################################################################


cor(split3$Ribosome.Sum, split3$growth_rates)

Fig_4C_right<-ggplot(split3, aes(split3$growth_rates, split3$Ribosome.Sum, colour=split3$Ribosome.Type, group=1)) +
  geom_point(size=4)+
  theme_bw(base_family= "serif") +
  theme(text = element_text(colour = "black",face="bold"),
        axis.text= element_text(colour = "black"),
        panel.grid.major.y = element_blank(),
        legend.position = "none") +
  scale_x_continuous(limits=c(0,0.5)) +
  scale_y_continuous(limits=c(0,40)) +
  labs(x="Growth rate (1/h)",
       y="Ribosomal fraction \nof proteome (%)",
       colour="",
       title="Figure 4C right \nRibosomes") +
  geom_errorbar(aes(ymin=split3$Ribosome.Sum-split3$Ribosome.Stdev, 
                    ymax=split3$Ribosome.Sum+split3$Ribosome.Stdev),
                width=0,colour="black") +
  geom_smooth(method=lm,   
              se=FALSE, colour="black", size=0.5, fullrange=TRUE) +
  scale_colour_manual(values=c("orange", "deepskyblue2", "tan4", "deepskyblue4")) +
  annotate("text", x = 0.1, y = 30, label = expression(R^"2":0.971), family= "serif")

Fig_4C_right

grid.arrange(Fig_4C_left, Fig_4C_right, nrow=1)


########################################################################################
##########################         FIGURE 4D left         ##############################
########################################################################################

head(Anaer_subset5, 1)
Anaer_subset7<-Anaer_subset5[which(Anaer_subset5$Process =="Ribosome"),]
Anaer_subset8<-Anaer_subset7
Anaer_subset8$Mean_anaer<-Anaer_subset8$Mean_anaer*100
Anaer_subset7$Stdev_anaer<-Anaer_subset7$Stdev_anaer*100

split4<-as.data.table(split(Anaer_subset8, Anaer_subset8$Type))

cor(split4$Rich_anaer.Mean_anaer, split4$Min_anaer.Mean_anaer)

Fig_4D_left<-ggplot(split4, aes(split4$Min_anaer.Mean_anaer, split4$Rich_anaer.Mean_anaer, colour=split4$Rich_anaer.Type, group=1)) +
  geom_point(size=4)+
  theme_bw(base_family= "serif") +
  theme(text = element_text(colour = "black",face="bold"),
        axis.text= element_text(colour = "black"),
        panel.grid.major.y = element_blank(),
        legend.position = "none") +
  scale_x_continuous(limits=c(0,2)) +
  scale_y_continuous(limits=c(0,2)) +
  labs(x="Ribosomal fraction \nof proteome (%) (minimal)",
       y="Ribosomal fraction \nof proteome (%) (rich)",
       colour="",
       title="Figure 4D left \nAnaerobic") +
  geom_errorbar(aes(ymin=split4$Rich_anaer.Mean_anaer-split4$Rich_anaer.Stdev_anaer, 
                    ymax=split4$Rich_anaer.Mean_anaer+split4$Rich_anaer.Stdev_anaer),
                width=0.1,colour="black", size=0.3) +
  geom_errorbarh(aes(xmin=split4$Min_anaer.Mean_anaer-split4$Min_anaer.Stdev_anaer, 
                    xmax=split4$Min_anaer.Mean_anaer+split4$Min_anaer.Stdev_anaer),
                height=0.1,colour="black",size=0.3) +
  geom_smooth(method=lm,   
              se=FALSE, colour="black", size=0.5, fullrange=TRUE) +
  scale_colour_manual(values=c("deepskyblue2")) +
  annotate("text", x = 1.5, y = 0.2, label = expression(R^"2":0.999), family= "serif") 

Fig_4D_left

########################################################################################
##########################         FIGURE 4D right        ##############################
########################################################################################

head(Aer_subset5, 1)
Aer_subset7<-Aer_subset5[which(Aer_subset5$Process =="Ribosome"),]
Aer_subset8<-Aer_subset7
Aer_subset8$Mean_Aer<-Aer_subset8$Mean_Aer*100
Aer_subset8$Stdev_Aer<-Aer_subset8$Stdev_Aer*100

split5<-as.data.table(split(Aer_subset8, Aer_subset8$Type))

cor(split5$Rich_aer.Mean_Aer, split5$Min_aer.Mean_Aer)

Fig_4D_right<-ggplot(split5, aes(split5$Min_aer.Mean_Aer, split5$Rich_aer.Mean_Aer, colour=split5$Rich_aer.Type, group=1)) +
  geom_point(size=4)+
  theme_bw(base_family= "serif") +
  theme(text = element_text(colour = "black",face="bold"),
        axis.text= element_text(colour = "black"),
        panel.grid.major.y = element_blank(),
        legend.position = "none") +
  scale_x_continuous(limits=c(0,2)) +
  scale_y_continuous(limits=c(0,2)) +
  labs(x="Ribosomal fraction \nof proteome (%) (minimal)",
       y="Ribosomal fraction \nof proteome (%) (rich)",
       colour="",
       title="Figure 4D right \nAerobic") +
  geom_errorbar(aes(ymin=split5$Rich_aer.Mean_Aer-split5$Rich_aer.Stdev_Aer, 
                    ymax=split5$Rich_aer.Mean_Aer+split5$Rich_aer.Stdev_Aer),
                width=0.1,colour="black", size=0.3) +
  geom_errorbarh(aes(xmin=split5$Min_aer.Mean_Aer-split5$Min_aer.Stdev_Aer, 
                     xmax=split5$Min_aer.Mean_Aer+split5$Min_aer.Stdev_Aer),
                 height=0.1,colour="black",size=0.3) +
  geom_smooth(method=lm,   
              se=FALSE, colour="black", size=0.5, fullrange=TRUE) +
  scale_colour_manual(values=c("tan3")) +
  annotate("text", x = 1.5, y = 0.2, label = expression(R^"2":0.996), family= "serif") 

Fig_4D_right

grid.arrange(Fig_4D_left, Fig_4D_right, nrow=1)


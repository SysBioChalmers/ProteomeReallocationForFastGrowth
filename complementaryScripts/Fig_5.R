

library(ggplot2)

Glycolysis_enzymes <- data.frame(c('PGI1','PFK2','PFK1','FBA1','TPI1','TDH1','TDH2','TDH3','PGK1','GPM1','ENO1','ENO2','PYK2','CDC19','HXK2','HXK1','GLK1'))
TCA_enzymes <- data.frame(c('PYC1','PYC2','CIT1','CIT3','ACO1','ACO2','IDH1','IDH2','KGD1','KGD2','LPD1','LSC1','LSC2','SDH1','SDH2','SDH3','SDH4','YJL045W','FUM1','MDH1'))
PPP_enzymes <- data.frame(c('ZWF1','SOL4','SOL3','GND2','GND1','RPE1','RKI1','TKL1','TKL2','TAL1'))
colnames(Glycolysis_enzymes) <- c('Protein')
colnames(TCA_enzymes) <- c('Protein')
colnames(PPP_enzymes) <- c('Protein')


Input_aerobic_glytcappp <- merge(glu_proteomics_aerobic_min,glu_proteomics_aerobic_rich,by='Accession',all=TRUE)
Input_anaerobic_glytcappp <- merge(glu_proteomics_anaerobic_min,glu_proteomics_anaerobic_rich,by='Accession',all=TRUE)

colnames(Input_aerobic_glytcappp) <- c('Accession','Min_1','Min_2','Gene','Protein','Rich_1','Rich_2','Gene_2','Protein_2')
colnames(Input_anaerobic_glytcappp) <- c('Accession','Min_1','Min_2','Min_3','Gene','Protein','Rich_1','Rich_2','Rich_3','Gene_2','Protein_2')

Glycolysis_abundancies_aer <- merge(Input_aerobic_glytcappp,Glycolysis_enzymes,by='Protein',all.x=FALSE)
Glycolysis_abundancies_anaer <- merge(Input_anaerobic_glytcappp,Glycolysis_enzymes,by='Protein',all.x=FALSE)
TCA_abundancies_aer <- merge(Input_aerobic_glytcappp,TCA_enzymes,by='Protein',all.x=FALSE)
TCA_abundancies_anaer <- merge(Input_anaerobic_glytcappp,TCA_enzymes,by='Protein',all.x=FALSE)

TCA_extra_enzyme <- data.frame(c('YJL045W'))
colnames(TCA_extra_enzyme) <- c('Gene')
TCA_extra_aer <- merge(Input_aerobic_glytcappp,TCA_extra_enzyme,by='Gene',all.x=FALSE)
TCA_extra_anaer <- merge(Input_anaerobic_glytcappp,TCA_extra_enzyme,by='Gene',all.x=FALSE)

j_aer <- nrow(TCA_abundancies_aer)+1
j_anaer <- nrow(TCA_abundancies_anaer)+1

TCA_abundancies_aer[j_aer,] <- TCA_extra_aer
TCA_abundancies_anaer[j_anaer,] <- TCA_extra_anaer


glycolysis_name_aer_min <- as.data.frame(rep('Minimal Aerobic',13))
colnames(glycolysis_name_aer_min) <- 'Name'
glycolysis_name_aer_rich <- as.data.frame(rep('Rich Aerobic',13))
colnames(glycolysis_name_aer_rich) <- 'Name'
glycolysis_sample_names_aer <- rbind(glycolysis_name_aer_min,glycolysis_name_aer_rich)
#
glycolysis_name_anaer_min <- as.data.frame(rep('Minimal Anaerobic',13))
colnames(glycolysis_name_anaer_min) <- 'Name'
glycolysis_name_anaer_rich <- as.data.frame(rep('Rich Anaerobic',13))
colnames(glycolysis_name_anaer_rich) <- 'Name'
glycolysis_sample_names_anaer <- rbind(glycolysis_name_anaer_min,glycolysis_name_anaer_rich)
#
tca_name_aer_min <- as.data.frame(rep('Minimal Aerobic',20))
colnames(tca_name_aer_min) <- 'Name'
tca_name_aer_rich <- as.data.frame(rep('Rich Aerobic',20))
colnames(tca_name_aer_rich) <- 'Name'
tca_sample_names_aer <- rbind(tca_name_aer_min,tca_name_aer_rich)
#
tca_name_anaer_min <- as.data.frame(rep('Minimal Anaerobic',20))
colnames(tca_name_anaer_min) <- 'Name'
tca_name_anaer_rich <- as.data.frame(rep('Rich Anaerobic',20))
colnames(tca_name_anaer_rich) <- 'Name'
tca_sample_names_anaer <- rbind(tca_name_anaer_min,tca_name_anaer_rich)
#
#
#
Collected_glycolysis_aer <- glycolysis_sample_names_aer
Collected_glycolysis_anaer <- glycolysis_sample_names_anaer
Collected_tca_aer <- tca_sample_names_aer
Collected_tca_anaer <- tca_sample_names_anaer
#
#
#
glycolysis_enzyme_names <- rbind(as.data.frame(Glycolysis_abundancies_aer[,1]),as.data.frame(Glycolysis_abundancies_aer[,1]))
tca_enzyme_names <- rbind(as.data.frame(TCA_abundancies_aer[,1]),as.data.frame(TCA_abundancies_aer[,1]))
#
Collected_glycolysis_aer <- as.data.frame(matrix(0,nrow=nrow(Glycolysis_abundancies_aer),ncol=5))
Collected_glycolysis_anaer <- as.data.frame(matrix(0,nrow=nrow(Glycolysis_abundancies_anaer),ncol=5))
Collected_tca_aer <- as.data.frame(matrix(0,nrow=nrow(TCA_abundancies_aer),ncol=5))
Collected_tca_anaer <- as.data.frame(matrix(0,nrow=nrow(TCA_abundancies_anaer),ncol=5))
#
Collected_glycolysis_aer[,1] <- as.data.frame(Glycolysis_abundancies_aer[,1])
Collected_glycolysis_anaer[,1] <- as.data.frame(Glycolysis_abundancies_anaer[,1])
Collected_tca_aer[,1] <- as.data.frame(TCA_abundancies_aer[,1])
Collected_tca_anaer[,1] <- as.data.frame(TCA_abundancies_anaer[,1])
#
Collected_glycolysis_aer[,3] <- 0
Collected_glycolysis_aer[,4] <- 0
Collected_glycolysis_anaer[,3] <- 0
Collected_glycolysis_anaer[,4] <- 0
Collected_tca_aer[,3] <- 0
Collected_tca_aer[,4] <- 0
Collected_tca_anaer[,3] <- 0
Collected_tca_anaer[,4] <- 0
#
j <- nrow(Collected_glycolysis_aer)
for(i in seq(1,j,by=1)){
  
  Collected_glycolysis_aer[i,2] <- mean(c(as.numeric(Glycolysis_abundancies_aer[i,3])*100,as.numeric(Glycolysis_abundancies_aer[i,4])*100))
  Collected_glycolysis_aer[i,4] <- mean(c(as.numeric(Glycolysis_abundancies_aer[i,6])*100,as.numeric(Glycolysis_abundancies_aer[i,7])*100))
  Collected_glycolysis_aer[i,3] <- sd(c(as.numeric(Glycolysis_abundancies_aer[i,3])*100,as.numeric(Glycolysis_abundancies_aer[i,4])*100))
  Collected_glycolysis_aer[i,5] <- sd(c(as.numeric(Glycolysis_abundancies_aer[i,6])*100,as.numeric(Glycolysis_abundancies_aer[i,7])*100))
}

j <- nrow(Collected_glycolysis_anaer)
for(i in seq(1,j,by=1)){
  
  Collected_glycolysis_anaer[i,2] <- mean(c(as.numeric(Glycolysis_abundancies_anaer[i,3])*100,as.numeric(Glycolysis_abundancies_anaer[i,4])*100,as.numeric(Glycolysis_abundancies_anaer[i,5])*100))
  Collected_glycolysis_anaer[i,4] <- mean(c(as.numeric(Glycolysis_abundancies_anaer[i,7])*100,as.numeric(Glycolysis_abundancies_anaer[i,8])*100,as.numeric(Glycolysis_abundancies_anaer[i,9])*100))
  Collected_glycolysis_anaer[i,3] <- sd(c(as.numeric(Glycolysis_abundancies_anaer[i,3])*100,as.numeric(Glycolysis_abundancies_anaer[i,4])*100,as.numeric(Glycolysis_abundancies_anaer[i,5])*100))
  Collected_glycolysis_anaer[i,5] <- sd(c(as.numeric(Glycolysis_abundancies_anaer[i,7])*100,as.numeric(Glycolysis_abundancies_anaer[i,8])*100,as.numeric(Glycolysis_abundancies_anaer[i,9])*100))
}

colnames(Collected_glycolysis_aer) <- c('Protein','Min_mass','Min_mass_stdev','Rich_mass','Rich_mass_stdev')
colnames(Collected_glycolysis_anaer) <- c('Protein','Min_mass','Min_mass_stdev','Rich_mass','Rich_mass_stdev')
#
#
#
j <- nrow(Collected_tca_aer)
for(i in seq(1,j,by=1)){
  
  Collected_tca_aer[i,2] <- mean(c(as.numeric(TCA_abundancies_aer[i,3])*100,as.numeric(TCA_abundancies_aer[i,4])*100))
  Collected_tca_aer[i,4] <- mean(c(as.numeric(TCA_abundancies_aer[i,6])*100,as.numeric(TCA_abundancies_aer[i,7])*100))
  Collected_tca_aer[i,3] <- sd(c(as.numeric(TCA_abundancies_aer[i,3])*100,as.numeric(TCA_abundancies_aer[i,4])*100))
  Collected_tca_aer[i,5] <- sd(c(as.numeric(TCA_abundancies_aer[i,6])*100,as.numeric(TCA_abundancies_aer[i,7])*100))
}

j <- nrow(Collected_tca_anaer)
for(i in seq(1,j,by=1)){
  Collected_tca_anaer[i,2] <- mean(c(as.numeric(TCA_abundancies_anaer[i,3])*100,as.numeric(TCA_abundancies_anaer[i,4])*100,as.numeric(TCA_abundancies_anaer[i,5])*100))
  Collected_tca_anaer[i,4] <- mean(c(as.numeric(TCA_abundancies_anaer[i,7])*100,as.numeric(TCA_abundancies_anaer[i,8])*100,as.numeric(TCA_abundancies_anaer[i,9])*100))
  Collected_tca_anaer[i,3] <- sd(c(as.numeric(TCA_abundancies_anaer[i,3])*100,as.numeric(TCA_abundancies_anaer[i,4])*100,as.numeric(TCA_abundancies_anaer[i,5])*100))
  Collected_tca_anaer[i,5] <- sd(c(as.numeric(TCA_abundancies_anaer[i,7])*100,as.numeric(TCA_abundancies_anaer[i,8])*100,as.numeric(TCA_abundancies_anaer[i,9])*100))
  
}

colnames(Collected_tca_aer) <- c('Protein','Min_mass','Min_mass_stdev','Rich_mass','Rich_mass_stdev')
colnames(Collected_tca_anaer) <- c('Protein','Min_mass','Min_mass_stdev','Rich_mass','Rich_mass_stdev')
#


### Collect mean and sd of sums of glycolysis and tca and ppp for all samples

Sums_glycolysis <- data.frame(c('Minimal Anaerobic','Rich Anaerobic','Minimal Aerobic','Rich Aerobic'),c('2 Glycolysis Anaerobic','2 Glycolysis Anaerobic','1 Glycolysis Aerobic','1 Glycolysis Aerobic'))
Sums_tca <- data.frame(c('Minimal Anaerobic','Rich Anaerobic','Minimal Aerobic','Rich Aerobic'),c('2 TCA-cycle Anaerobic','2 TCA-cycle Anaerobic','1 TCA-cycle Aerobic','1 TCA-cycle Aerobic'))
Sums_glycolysis[,3] <- 0
Sums_glycolysis[,4] <- 0
Sums_tca[,3] <- 0
Sums_tca[,4] <- 0
colnames(Sums_glycolysis) <- c('Sample','Process','Mass','Mass_Stdev')
colnames(Sums_tca) <- c('Sample','Process','Mass','Mass_Stdev')
#
Sums_glycolysis[1,3] <- mean(c(sum(Glycolysis_abundancies_anaer[,3])*100,sum(Glycolysis_abundancies_anaer[,4])*100,sum(Glycolysis_abundancies_anaer[,5])*100))
Sums_glycolysis[2,3] <- mean(c(sum(Glycolysis_abundancies_anaer[,7])*100,sum(Glycolysis_abundancies_anaer[,8])*100,sum(Glycolysis_abundancies_anaer[,9])*100))
Sums_glycolysis[3,3] <- mean(c(sum(Glycolysis_abundancies_aer[,3])*100,sum(Glycolysis_abundancies_aer[,4])*100))
Sums_glycolysis[4,3] <- mean(c(sum(Glycolysis_abundancies_aer[,6])*100,sum(Glycolysis_abundancies_aer[,7])*100))
Sums_glycolysis[1,4] <- sd(c(sum(Glycolysis_abundancies_anaer[,3])*100,sum(Glycolysis_abundancies_anaer[,4])*100,sum(Glycolysis_abundancies_anaer[,5])*100))
Sums_glycolysis[2,4] <- sd(c(sum(Glycolysis_abundancies_anaer[,7])*100,sum(Glycolysis_abundancies_anaer[,8])*100,sum(Glycolysis_abundancies_anaer[,9])*100))
Sums_glycolysis[3,4] <- sd(c(sum(Glycolysis_abundancies_aer[,3])*100,sum(Glycolysis_abundancies_aer[,4])*100))
Sums_glycolysis[4,4] <- sd(c(sum(Glycolysis_abundancies_aer[,6])*100,sum(Glycolysis_abundancies_aer[,7])*100))
#
t_test_glycolysis_anaer <- t.test(c(as.numeric(log2(c(sum(Glycolysis_abundancies_anaer[,3])*100,sum(Glycolysis_abundancies_anaer[,4])*100,sum(Glycolysis_abundancies_anaer[,5])*100)))),c(as.numeric(log2(c(sum(Glycolysis_abundancies_anaer[,7])*100,sum(Glycolysis_abundancies_anaer[,8])*100,sum(Glycolysis_abundancies_anaer[,9])*100)))),alternative='t',paired=FALSE)
t_test_glycolysis_aer <- t.test(c(as.numeric(log2(c(sum(Glycolysis_abundancies_aer[,3])*100,sum(Glycolysis_abundancies_aer[,4])*100)))),c(as.numeric(log2(c(sum(Glycolysis_abundancies_aer[,6])*100,sum(Glycolysis_abundancies_aer[,7])*100)))),alternative='t',paired=FALSE)
#
##
#
Sums_tca[1,3] <- mean(c(sum(TCA_abundancies_anaer[,3])*100,sum(TCA_abundancies_anaer[,4])*100,sum(TCA_abundancies_anaer[,5])*100))
Sums_tca[2,3] <- mean(c(sum(TCA_abundancies_anaer[,7])*100,sum(TCA_abundancies_anaer[,8])*100,sum(TCA_abundancies_anaer[,9])*100))
Sums_tca[3,3] <- mean(c(sum(TCA_abundancies_aer[,3])*100,sum(TCA_abundancies_aer[,4])*100))
Sums_tca[4,3] <- mean(c(sum(TCA_abundancies_aer[,6])*100,sum(TCA_abundancies_aer[,7])*100))
Sums_tca[1,4] <- sd(c(sum(TCA_abundancies_anaer[,3])*100,sum(TCA_abundancies_anaer[,4])*100,sum(TCA_abundancies_anaer[,5])*100))
Sums_tca[2,4] <- sd(c(sum(TCA_abundancies_anaer[,7])*100,sum(TCA_abundancies_anaer[,8])*100,sum(TCA_abundancies_anaer[,9])*100))
Sums_tca[3,4] <- sd(c(sum(TCA_abundancies_aer[,3])*100,sum(TCA_abundancies_aer[,4])*100))
Sums_tca[4,4] <- sd(c(sum(TCA_abundancies_aer[,6])*100,sum(TCA_abundancies_aer[,7])*100))
#
t_test_tca_anaer <- t.test(c(as.numeric(log2(c(sum(TCA_abundancies_anaer[,3])*100,sum(TCA_abundancies_anaer[,4])*100,sum(TCA_abundancies_anaer[,5])*100)))),c(as.numeric(log2(c(sum(TCA_abundancies_anaer[,7])*100,sum(TCA_abundancies_anaer[,8])*100,sum(TCA_abundancies_anaer[,9])*100)))),alternative='t',paired=FALSE)
t_test_tca_aer <- t.test(c(as.numeric(log2(c(sum(TCA_abundancies_aer[,3])*100,sum(TCA_abundancies_aer[,4])*100)))),c(as.numeric(log2(c(sum(TCA_abundancies_aer[,6])*100,sum(TCA_abundancies_aer[,7])*100)))),alternative='t',paired=FALSE)
#
##

### PLOT MEAN AND SD OF SUMS

# PLOT GLYCOLYSIS SUMS

ggplot(data=Sums_glycolysis, aes(x=Process , y=Mass, fill=Sample)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5,size=0.3) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c('darkorange1','deepskyblue1','brown','darkorchid4')) +
  labs(x=" ") +
  labs(y="Mass ratio of proteome") +
  scale_y_continuous(breaks=c(0,5,10,15,20,25),limits=c(0,25)) +
  theme(plot.title = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y=element_text(color = "black",size=12,face=c('bold'))) +
  geom_errorbar(aes(ymin=Sums_glycolysis$Mass-Sums_glycolysis$Mass_Stdev, ymax=Sums_glycolysis$Mass+Sums_glycolysis$Mass_Stdev), width=.2,position=position_dodge(0.50),size=0.3) +
  theme(legend.position = 'none') +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  ggsave("glycolysis_sums.tiff", height=2, width=3, units='in', dpi=600)


# PLOT TCA SUMS

ggplot(data=Sums_tca, aes(x=Process , y=Mass, fill=Sample)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5,size=0.3) +
  coord_flip() +
  theme_bw() +
  scale_fill_manual(values = c('darkorange1','deepskyblue1','brown','darkorchid4')) +
  labs(x=" ") +
  labs(y="Mass ratio of proteome") +
  scale_y_continuous(breaks=c(0,0.2,0.4,0.6,0.8,1),limits=c(0,1)) +
  theme(plot.title = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.title.x = element_blank()) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y=element_text(color = "black",size=12,face=c('bold'))) +
  geom_errorbar(aes(ymin=Sums_tca$Mass-Sums_tca$Mass_Stdev, ymax=Sums_tca$Mass+Sums_tca$Mass_Stdev), width=.2,position=position_dodge(0.50),size=0.3) +
  theme(legend.position = 'none') +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  ggsave("tca_sums.tiff", height=2, width=3, units='in', dpi=600)


#### PLOT TCA ENZYMES
# Data to dra lines in the plot
line_to_plot <- data.frame(c(0,0.25),c(0,0.25),c(0,0.125))
colnames(line_to_plot) <- c('x','y1','y2')
# Add ribbon data to data frame
x1 <- seq(0,0.25,by=0.25/(nrow(Collected_tca_anaer)-1))
x2 <- seq(0,0.125,by=0.125/(nrow(Collected_tca_anaer)-1))
Collected_tca_anaer[,6] <- x1
Collected_tca_anaer[,7] <- x2


#
ggplot(data=Collected_tca_anaer,aes(x=Min_mass,y=Rich_mass)) +
  geom_point(size=2,color='deepskyblue3',alpha=0.7) +
  geom_line(data=line_to_plot,aes(x=x,y=y1)) +
  geom_line(data=line_to_plot,aes(x=x,y=y2),color='blue',alpha=0.5) +
  #geom_line(data=line_to_plot,aes(x=x,y=y3),color='orange',alpha=0.8) +
  geom_ribbon(aes(ymin=V7,ymax=V6,x=V6),fill='blue',alpha=0.2) +
  scale_x_continuous(breaks=c(0,0.05,0.10,0.15,0.20,0.25),limits=c(-0,0.25)) +
  scale_y_continuous(breaks=c(0,0.05,0.10,0.15,0.20,0.25),limits=c(-0,0.25)) +
  geom_errorbar(data=Collected_tca_anaer,aes(ymin=Collected_tca_anaer$Rich_mass-Collected_tca_anaer$Rich_mass_stdev,ymax=Collected_tca_anaer$Rich_mass+Collected_tca_anaer$Rich_mass_stdev),color='black',size=0.5,width=0) +
  geom_errorbarh(data=Collected_tca_anaer,aes(xmin=Collected_tca_anaer$Min_mass-Collected_tca_anaer$Min_mass_stdev,xmax=Collected_tca_anaer$Min_mass+Collected_tca_anaer$Min_mass_stdev),color='black',size=0.5,height=0) +
  theme_bw() +
  theme(panel.grid.minor=element_blank()) +
  ggsave("tca_enzymes.tiff", height=2.5, width=3, units='in', dpi=600)


#### PLOT TCA ENZYMES
# Data to dra lines in the plot
line_to_plot <- data.frame(c(0,0.25),c(0,0.25),c(0,0.125))
colnames(line_to_plot) <- c('x','y1','y2')
# Add ribbon data to data frame
x1 <- seq(0,0.25,by=0.25/(nrow(Collected_tca_aer)-1))
x2 <- seq(0,0.125,by=0.125/(nrow(Collected_tca_aer)-1))
Collected_tca_aer[,6] <- x1
Collected_tca_aer[,7] <- x2

#
ggplot(data=Collected_tca_aer,aes(x=Min_mass,y=Rich_mass)) +
  geom_point(size=2,color='darkorange3',alpha=0.7) +
  geom_line(data=line_to_plot,aes(x=x,y=y1)) +
  geom_line(data=line_to_plot,aes(x=x,y=y2),color='orange',alpha=0.5) +
  #geom_line(data=line_to_plot,aes(x=x,y=y3),color='orange',alpha=0.8) +
  geom_ribbon(aes(ymin=V7,ymax=V6,x=V6),fill='orange',alpha=0.2) +
  scale_x_continuous(breaks=c(0,0.05,0.10,0.15,0.20,0.25),limits=c(-0,0.25)) +
  scale_y_continuous(breaks=c(0,0.05,0.10,0.15,0.20,0.25),limits=c(-0,0.25)) +
  geom_errorbar(data=Collected_tca_aer,aes(ymin=Collected_tca_aer$Rich_mass-Collected_tca_aer$Rich_mass_stdev,ymax=Collected_tca_aer$Rich_mass+Collected_tca_aer$Rich_mass_stdev),color='black',size=0.5,width=0) +
  geom_errorbarh(data=Collected_tca_aer,aes(xmin=Collected_tca_aer$Min_mass-Collected_tca_aer$Min_mass_stdev,xmax=Collected_tca_aer$Min_mass+Collected_tca_aer$Min_mass_stdev),color='black',size=0.5,height=0) +
  theme_bw() +
  theme(panel.grid.minor=element_blank()) +
  ggsave("tca_enzymes_aer.tiff", height=2.5, width=3, units='in', dpi=600)

##### PLOT GLYCOLYSIS ENZYMES
# Data to dra lines in the plot
line_to_plot <- data.frame(c(0,8),c(0,8),c(0,12),c(0,4))
colnames(line_to_plot) <- c('x','y1','y2','y3')
# Add ribbon data to data frame
x1 <- seq(0,12,by=12/(nrow(Collected_glycolysis_aer)-1))
x2 <- seq(0,4,by=4/(nrow(Collected_glycolysis_aer)-1))
x3 <- seq(0,8,by=8/(nrow(Collected_glycolysis_aer)-1))
Collected_glycolysis_aer[,6] <- x1
Collected_glycolysis_aer[,7] <- x2
Collected_glycolysis_aer[,8] <- x3
#
ggplot(data=Collected_glycolysis_aer,aes(x=Min_mass,y=Rich_mass)) +
  geom_point(size=2,color='orange',alpha=0.7) +
  geom_line(data=line_to_plot,aes(x=x,y=y1)) +
  #geom_line(data=line_to_plot,aes(x=x,y=y2),color='orange',alpha=0.8) +
  #geom_line(data=line_to_plot,aes(x=x,y=y3),color='orange',alpha=0.8) +
  #geom_ribbon(aes(ymin=V7,ymax=V6,x=V8),fill='orange3',alpha=0.2) +
  scale_x_continuous(breaks=c(0,2,4,6,8),limits=c(-0,8)) +
  scale_y_continuous(breaks=c(0,2,4,6,8),limits=c(-0,8)) +
  geom_errorbar(data=Collected_glycolysis_aer,aes(ymin=Collected_glycolysis_aer$Rich_mass-Collected_glycolysis_aer$Rich_mass_stdev,ymax=Collected_glycolysis_aer$Rich_mass+Collected_glycolysis_aer$Rich_mass_stdev),color='black',size=0.3,width=0) +
  geom_errorbarh(data=Collected_glycolysis_aer,aes(xmin=Collected_glycolysis_aer$Min_mass-Collected_glycolysis_aer$Min_mass_stdev,xmax=Collected_glycolysis_aer$Min_mass+Collected_glycolysis_aer$Min_mass_stdev),color='black',size=0.3,height=0) +
  theme_bw() +
  theme(panel.grid.minor=element_blank()) +
  ggsave("glycolysis_enzymes.tiff", height=1.8, width=2, units='in', dpi=600)

------
  ##### PLOT GLYCOLYSIS ENZYMES
  # Data to dra lines in the plot
  line_to_plot <- data.frame(c(0,8),c(0,8),c(0,12),c(0,4))
colnames(line_to_plot) <- c('x','y1','y2','y3')
# Add ribbon data to data frame
x1 <- seq(0,12,by=12/(nrow(Collected_glycolysis_anaer)-1))
x2 <- seq(0,4,by=4/(nrow(Collected_glycolysis_anaer)-1))
x3 <- seq(0,8,by=8/(nrow(Collected_glycolysis_anaer)-1))
Collected_glycolysis_anaer[,6] <- x1
Collected_glycolysis_anaer[,7] <- x2
Collected_glycolysis_anaer[,8] <- x3
#
ggplot(data=Collected_glycolysis_anaer,aes(x=Min_mass,y=Rich_mass)) +
  geom_point(size=2,color='deepskyblue3',alpha=0.7) +
  geom_line(data=line_to_plot,aes(x=x,y=y1)) +
  #geom_line(data=line_to_plot,aes(x=x,y=y2),color='orange',alpha=0.8) +
  #geom_line(data=line_to_plot,aes(x=x,y=y3),color='orange',alpha=0.8) +
  #geom_ribbon(aes(ymin=V7,ymax=V6,x=V8),fill='orange3',alpha=0.2) +
  scale_x_continuous(breaks=c(0,2,4,6,8),limits=c(-0,8)) +
  scale_y_continuous(breaks=c(0,2,4,6,8),limits=c(-0,8)) +
  geom_errorbar(data=Collected_glycolysis_anaer,aes(ymin=Collected_glycolysis_anaer$Rich_mass-Collected_glycolysis_anaer$Rich_mass_stdev,ymax=Collected_glycolysis_anaer$Rich_mass+Collected_glycolysis_anaer$Rich_mass_stdev),color='black',size=0.3,width=0) +
  geom_errorbarh(data=Collected_glycolysis_anaer,aes(xmin=Collected_glycolysis_anaer$Min_mass-Collected_glycolysis_anaer$Min_mass_stdev,xmax=Collected_glycolysis_anaer$Min_mass+Collected_glycolysis_anaer$Min_mass_stdev),color='black',size=0.3,height=0) +
  theme_bw() +
  theme(panel.grid.minor=element_blank()) +
  ggsave("glycolysis_enzymes_anaer.tiff", height=1.8, width=2, units='in', dpi=600)



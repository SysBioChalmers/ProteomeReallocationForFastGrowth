
# Anaerobic
glu_proteomics_aerobic <- glu_proteomics_aerobic[,c(1,2,3,4,5,6,7)]
colnames(glu_proteomics_aerobic) <- c('Accession','Rich_1','Rich_2','Min_1','Min_2','Gene','Protein')

nr_rows_aerobic <- nrow(glu_proteomics_aerobic)

accession_variable <- 1
prot_name_variable <- 7
gene_name_variable <- 6

for (i in c(1:nr_rows_aerobic)){
  
  ttesting_row_aerobic <- t.test(c(as.numeric(log2(input_data_aerobic[i,2:3]))),c(as.numeric(log2(input_data_aerobic[i,4:5]))),alternative='t',paired=FALSE)
  
  rich_mean_aerobic <- mean(c(as.numeric(input_data_aerobic[i,2]),as.numeric(input_data_aerobic[i,3])))
  min_mean_aerobic <- mean(c(as.numeric(input_data_aerobic[i,5]),as.numeric(input_data_aerobic[i,6])))
  ratio_rich_min_aerobic <- rich_mean_anaerobic/min_mean_aerobic
  log2_fold_change_aerobic <- log2(ratio_rich_min_aerobic)
  
  mass_average_min_aerobic <- mean(c(as.numeric(input_data_anaerobic[i,5])*100,as.numeric(input_data_anaerobic[i,6])*100,as.numeric(input_data_anaerobic[i,7])*100))
  
  # Collect the information computed above
  all_aerobic[i,1] <- as.character(input_data_aerobic[i,accession_variable])
  all_aerobic[i,2] <- as.character(input_data_aerobic[i,prot_name_variable])
  all_aerobic[i,3] <- as.character(input_data_aerobic[i,gene_name_variable])
  all_aerobic[i,4] <- log2_fold_change_aerobic
  all_aerobic[i,5] <- ttesting_row_aerobic$p.value
  all_aerobic[i,6] <- mass_average_min_aerobic
  
}

colnames(all_aerobic) <- c('Accession','Protein','Gene','Fold_change','P_value','Average_mass_minimal')

Sign_up_regulated_aerobic <- all_anaerobic[all_aerobic$P_value < 0.05, ]
Sign_up_regulated_aerobic <- Sign_up_regulated_aerobic[Sign_up_regulated_aerobic$Fold_change > 1 , ]
Sign_down_regulated_aerobic <- all_aerobic[all_aerobic$P_value < 0.05, ]
Sign_down_regulated_aerobic <- Sign_down_regulated_aerobic[Sign_down_regulated_aerobic$Fold_change < (-1) , ]

Not_sign_aerobic_1 <- all_aerobic[all_aerobic$Fold_change >= (-1), ]
Not_sign_aerobic_1 <- Not_sign_aerobic_1[Not_sign_aerobic_1$Fold_change <= 1 , ]
#
Not_sign_aerobic_2 <- all_aerobic[all_aerobic$Fold_change < (-1), ]
Not_sign_aerobic_2 <- Not_sign_aerobic_2[Not_sign_aerobic_2$P_value >= 0.05 , ]
#
Not_sign_anaerobic_3 <- all_anaerobic[all_anaerobic$Fold_change > 1, ]
Not_sign_anaerobic_3 <- Not_sign_anaerobic_3[Not_sign_anaerobic_3$P_value >= 0.05 , ]
#
if(nrow(Not_sign_anaerobic_2)>0 && nrow(Not_sign_anaerobic_3)>0){
  Not_sign_anaerobic <- rbind(Not_sign_anaerobic_1,Not_sign_anaerobic_2)
  Not_sign_anaerobic <- rbind(Not_sign_anaerobic,Not_sign_anaerobic_3)
} else if (nrow(Not_sign_anaerobic_2)>0){
  Not_sign_anaerobic <- rbind(Not_sign_anaerobic_1,Not_sign_anaerobic_2)
} else if (nrow(Not_sign_anaerobic_3)>0){
  Not_sign_anaerobic <- rbind(Not_sign_anaerobic_1,Not_sign_anaerobic_3)
} else{
  Not_sign_anaerobic <- Not_sign_anaerobic_1
}

check_collection <- nrow(Sign_up_regulated_anaerobic) + nrow(Sign_down_regulated_anaerobic) + nrow(Not_sign_anaerobic)
check_collection <- check_collection - nrow(input_data_anaerobic)
if(check_collection==0){
  print('COLLECTION IS OK!')
} else {
  print('COLLECTION WENT WRONG')
}


print('done')

######## VISUALIZE THE VULCANO PLOT
# Anaerobic
size_dots <- 3

ggplot(Sign_up_regulated_anaerobic,aes(x=Fold_change,y=abs(log10(P_value)))) +
  geom_point(colour='red4',size=size_dots,alpha=0.6) +
  geom_point(data=Not_sign_anaerobic,aes(x=Fold_change,y=abs(log10(P_value))),colour='grey10',size=size_dots,alpha=0.4) +
  geom_point(data=Sign_down_regulated_anaerobic,aes(x=Fold_change,y=abs(log10(P_value))),colour='blue4',size=size_dots,alpha=0.6) +
  theme_bw() +
  theme(panel.grid.minor=element_blank()) +
  scale_x_continuous(breaks=c(-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6),limits=c(-6,6)) +
  scale_y_continuous(breaks=seq(0,8,1),limits=c(0,8)) +
  geom_hline(yintercept = 1.30103,colour='black',linetype='dashed',size=1) +
  geom_vline(xintercept = c(1),colour='red4',linetype='dashed',size=1) +
  geom_vline(xintercept = c(-1),colour='blue4',linetype='dashed',size=1)+
  ggsave("anaerobic_volcano.tiff", height=5, width=5, units='in', dpi=600)

##
##Aerobic
input_data_aerobic <- glu_proteomics_aerobic[,c(1,2,3,4,5,6,7)]
colnames(input_data_aerobic) <- c('Accession','Rich_1','Rich_2','Min_1','Min_2','Gene','Protein')
nr_rows_aerobic <- nrow(input_data_aerobic)
all_aerobic <- as.data.frame(matrix(0,nrow=nr_rows_aerobic,ncol=6))
#
#
accession_variable <- 1
prot_name_variable <- 7
gene_name_variable <- 6

for (i in c(1:nr_rows_aerobic)){
  
  ttesting_row_aerobic <- t.test(c(as.numeric(log2(input_data_aerobic[i,2:3]))),c(as.numeric(log2(input_data_aerobic[i,4:5]))),alternative='t',paired=FALSE)
  
  rich_mean_aerobic <- mean(c(as.numeric(input_data_aerobic[i,2]),as.numeric(input_data_aerobic[i,3])))
  min_mean_aerobic <- mean(c(as.numeric(input_data_aerobic[i,4]),as.numeric(input_data_aerobic[i,5])))
  ratio_rich_min_aerobic <- rich_mean_aerobic/min_mean_aerobic
  log2_fold_change_aerobic <- log2(ratio_rich_min_aerobic)
  
  mass_average_min_aerobic <- mean(c(as.numeric(input_data_aerobic[i,4])*100,as.numeric(input_data_aerobic[i,5])*100))
  
  all_aerobic[i,1] <- as.character(input_data_aerobic[i,accession_variable])
  all_aerobic[i,2] <- as.character(input_data_aerobic[i,prot_name_variable])
  all_aerobic[i,3] <- as.character(input_data_aerobic[i,gene_name_variable])
  all_aerobic[i,4] <- log2_fold_change_aerobic
  all_aerobic[i,5] <- ttesting_row_aerobic$p.value
  all_aerobic[i,6] <- mass_average_min_aerobic
  
}

colnames(all_aerobic) <- c('Accession','Protein','Gene','Fold_change','P_value','Average_mass_minimal')

Sign_up_regulated_aerobic <- all_aerobic[all_aerobic$P_value < 0.05, ]
Sign_up_regulated_aerobic <- Sign_up_regulated_aerobic[Sign_up_regulated_aerobic$Fold_change > 1 , ]
Sign_down_regulated_aerobic <- all_aerobic[all_aerobic$P_value < 0.05, ]
Sign_down_regulated_aerobic <- Sign_down_regulated_aerobic[Sign_down_regulated_aerobic$Fold_change < (-1) , ]
Not_sign_aerobic_1 <- all_aerobic[all_aerobic$Fold_change >= (-1), ]
Not_sign_aerobic_1 <- Not_sign_aerobic_1[Not_sign_aerobic_1$Fold_change <= 1 , ]
#
Not_sign_aerobic_2 <- all_aerobic[all_aerobic$Fold_change < (-1), ]
Not_sign_aerobic_2 <- Not_sign_aerobic_2[Not_sign_aerobic_2$P_value >= 0.05 , ]
#
Not_sign_aerobic_3 <- all_aerobic[all_aerobic$Fold_change > 1, ]
Not_sign_aerobic_3 <- Not_sign_aerobic_3[Not_sign_aerobic_3$P_value >= 0.05 , ]
#
if(nrow(Not_sign_aerobic_2)>0 && nrow(Not_sign_aerobic_3)>0){
  Not_sign_aerobic <- rbind(Not_sign_aerobic_1,Not_sign_aerobic_2)
  Not_sign_aerobic <- rbind(Not_sign_aerobic,Not_sign_aerobic_3)
} else if (nrow(Not_sign_aerobic_2)>0){
  Not_sign_aerobic <- rbind(Not_sign_aerobic_1,Not_sign_aerobic_2)
} else if (nrow(Not_sign_aerobic_3)>0){
  Not_sign_aerobic <- rbind(Not_sign_aerobic_1,Not_sign_aerobic_3)
} else{
  Not_sign_aerobic <- Not_sign_aerobic_1
}
check_collection <- nrow(Sign_up_regulated_aerobic) + nrow(Sign_down_regulated_aerobic) + nrow(Not_sign_aerobic)
check_collection <- check_collection - nrow(input_data_aerobic)
if(check_collection==0){
  print('COLLECTION IS OK!')
} else {
  print('COLLECTION WENT WRONG')
}

print('done')

######## VISUALIZE THE VULCANO PLOT
# Aerobic
size_dots <- 3

ggplot(Sign_up_regulated_aerobic,aes(x=Fold_change,y=abs(log10(P_value)))) +
  geom_point(colour='red4',size=size_dots,alpha=0.6) +
  geom_point(data=Not_sign_aerobic,aes(x=Fold_change,y=abs(log10(P_value))),colour='grey10',size=size_dots,alpha=0.4) +
  geom_point(data=Sign_down_regulated_aerobic,aes(x=Fold_change,y=abs(log10(P_value))),colour='blue4',size=size_dots,alpha=0.6) +
  theme_bw() +
  theme(panel.grid.minor=element_blank()) +
  scale_x_continuous(breaks=c(-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6),limits=c(-6,6)) +
  scale_y_continuous(breaks=seq(0,8,1),limits=c(0,8)) +
  geom_hline(yintercept = 1.30103,colour='black',linetype='dashed',size=1) +
  geom_vline(xintercept = c(1),colour='red4',linetype='dashed',size=1) +
  geom_vline(xintercept = c(-1),colour='blue4',linetype='dashed',size=1) +
  ggsave("aerobic_volcano.tiff", height=5, width=5, units='in', dpi=600)


# Analyze the overlap of significantly up- and down-regulated proteins between rich and min

up_both_conditions <- merge(Sign_up_regulated_aerobic,Sign_up_regulated_anaerobic,by='Accession',all=TRUE)
down_both_conditions <- merge(Sign_down_regulated_aerobic,Sign_down_regulated_anaerobic,by='Accession',all=TRUE)
up_both_conditions <- na.omit(up_both_conditions)
down_both_conditions <- na.omit(down_both_conditions)



#########################

Adh1p_aer_min <- glu_proteomics_aerobic_min[glu_proteomics_aerobic_min$Protein_Name == 'ADH1', ]
Adh1p_aer_rich <- glu_proteomics_aerobic_rich[glu_proteomics_aerobic_rich$Protein_Name == 'ADH1', ]
Adh1p_anaer_min <- glu_proteomics_anaerobic_min[glu_proteomics_anaerobic_min$Protein_Name == 'ADH1', ]
Adh1p_anaer_rich <- glu_proteomics_anaerobic_rich[glu_proteomics_anaerobic_rich$Protein_Name == 'ADH1', ]

Adh3p_aer_min <- glu_proteomics_aerobic_min[glu_proteomics_aerobic_min$Protein_Name == 'ADH3', ]
Adh3p_aer_rich <- glu_proteomics_aerobic_rich[glu_proteomics_aerobic_rich$Protein_Name == 'ADH3', ]
Adh3p_anaer_min <- glu_proteomics_anaerobic_min[glu_proteomics_anaerobic_min$Protein_Name == 'ADH3', ]
Adh3p_anaer_rich <- glu_proteomics_anaerobic_rich[glu_proteomics_anaerobic_rich$Protein_Name == 'ADH3', ]

Zrt1p_aer_min <- glu_proteomics_aerobic_min[glu_proteomics_aerobic_min$Protein_Name == 'ZRT1', ]
Zrt1p_aer_rich <- glu_proteomics_aerobic_rich[glu_proteomics_aerobic_rich$Protein_Name == 'ZRT1', ]
Zrt1p_anaer_min <- glu_proteomics_anaerobic_min[glu_proteomics_anaerobic_min$Protein_Name == 'ZRT1', ]
Zrt1p_anaer_rich <- glu_proteomics_anaerobic_rich[glu_proteomics_anaerobic_rich$Protein_Name == 'ZRT1', ]

Zps1p_aer_min <- glu_proteomics_aerobic_min[glu_proteomics_aerobic_min$Protein_Name == 'ZPS1', ]
Zps1p_aer_rich <- glu_proteomics_aerobic_rich[glu_proteomics_aerobic_rich$Protein_Name == 'ZPS1', ]
Zps1p_anaer_min <- glu_proteomics_anaerobic_min[glu_proteomics_anaerobic_min$Protein_Name == 'ZPS1', ]
Zps1p_anaer_rich <- glu_proteomics_anaerobic_rich[glu_proteomics_anaerobic_rich$Protein_Name == 'ZPS1', ]

# Produce data frames that will contain average mass and stdev of average mass
ADH1p <- as.data.frame(matrix(0,nrow=4,ncol=3))
ADH3p <- as.data.frame(matrix(0,nrow=4,ncol=3))
ZRT1p <- as.data.frame(matrix(0,nrow=4,ncol=3))
ZPS1p <- as.data.frame(matrix(0,nrow=4,ncol=3))

ADH1p[,1] <- c('4 Rich aerobic','3 Min aerobic','2 Rich anaerobic','1 Min anaerobic')
ADH3p[,1] <- c('4 Rich aerobic','3 Min aerobic','2 Rich anaerobic','1 Min anaerobic')
ZRT1p[,1] <- c('4 Rich aerobic','3 Min aerobic','2 Rich anaerobic','1 Min anaerobic')
ZPS1p[,1] <- c('4 Rich aerobic','3 Min aerobic','2 Rich anaerobic','1 Min anaerobic')

# Calculate and save average mass and stdev of average mass#
ADH1p[1,2] <- mean(c(as.numeric(Adh1p_aer_rich[2]*100),as.numeric(Adh1p_aer_rich[3]*100)))
ADH1p[2,2] <- mean(c(as.numeric(Adh1p_aer_min[2]*100),as.numeric(Adh1p_aer_min[3]*100)))
ADH1p[3,2] <- mean(c(as.numeric(Adh1p_anaer_rich[2]*100),as.numeric(Adh1p_anaer_rich[3]*100),as.numeric(Adh1p_anaer_rich[4]*100)))
ADH1p[4,2] <- mean(c(as.numeric(Adh1p_anaer_min[2]*100),as.numeric(Adh1p_anaer_min[3]*100),as.numeric(Adh1p_anaer_min[4]*100)))
ADH1p[1,3] <- sd(c(as.numeric(Adh1p_aer_rich[2]*100),as.numeric(Adh1p_aer_rich[3]*100)))
ADH1p[2,3] <- sd(c(as.numeric(Adh1p_aer_min[2]*100),as.numeric(Adh1p_aer_min[3]*100)))
ADH1p[3,3] <- sd(c(as.numeric(Adh1p_anaer_rich[2]*100),as.numeric(Adh1p_anaer_rich[3]*100),as.numeric(Adh1p_anaer_rich[4]*100)))
ADH1p[4,3] <- sd(c(as.numeric(Adh1p_anaer_min[2]*100),as.numeric(Adh1p_anaer_min[3]*100),as.numeric(Adh1p_anaer_min[4]*100)))
#
ADH3p[1,2] <- mean(c(as.numeric(Adh3p_aer_rich[2]*100),as.numeric(Adh3p_aer_rich[3]*100)))
ADH3p[2,2] <- mean(c(as.numeric(Adh3p_aer_min[2]*100),as.numeric(Adh3p_aer_min[3]*100)))
ADH3p[3,2] <- mean(c(as.numeric(Adh3p_anaer_rich[2]*100),as.numeric(Adh3p_anaer_rich[3]*100),as.numeric(Adh3p_anaer_rich[4]*100)))
ADH3p[4,2] <- mean(c(as.numeric(Adh3p_anaer_min[2]*100),as.numeric(Adh3p_anaer_min[3]*100),as.numeric(Adh3p_anaer_min[4]*100)))
ADH3p[1,3] <- sd(c(as.numeric(Adh3p_aer_rich[2]*100),as.numeric(Adh3p_aer_rich[3]*100)))
ADH3p[2,3] <- sd(c(as.numeric(Adh3p_aer_min[2]*100),as.numeric(Adh3p_aer_min[3]*100)))
ADH3p[3,3] <- sd(c(as.numeric(Adh3p_anaer_rich[2]*100),as.numeric(Adh3p_anaer_rich[3]*100),as.numeric(Adh3p_anaer_rich[4]*100)))
ADH3p[4,3] <- sd(c(as.numeric(Adh3p_anaer_min[2]*100),as.numeric(Adh3p_anaer_min[3]*100),as.numeric(Adh3p_anaer_min[4]*100)))
#
ZRT1p[1,2] <- mean(c(as.numeric(Zrt1p_aer_rich[2]*100),as.numeric(Zrt1p_aer_rich[3]*100)))
ZRT1p[2,2] <- mean(c(as.numeric(Zrt1p_aer_min[2]*100),as.numeric(Zrt1p_aer_min[3]*100)))
ZRT1p[3,2] <- mean(c(as.numeric(Zrt1p_anaer_rich[2]*100),as.numeric(Zrt1p_anaer_rich[3]*100),as.numeric(Zrt1p_anaer_rich[4]*100)))
ZRT1p[4,2] <- mean(c(as.numeric(Zrt1p_anaer_min[2]*100),as.numeric(Zrt1p_anaer_min[3]*100),as.numeric(Zrt1p_anaer_min[4]*100)))
ZRT1p[1,3] <- sd(c(as.numeric(Zrt1p_aer_rich[2]*100),as.numeric(Zrt1p_aer_rich[3]*100)))
ZRT1p[2,3] <- sd(c(as.numeric(Zrt1p_aer_min[2]*100),as.numeric(Zrt1p_aer_min[3]*100)))
ZRT1p[3,3] <- sd(c(as.numeric(Zrt1p_anaer_rich[2]*100),as.numeric(Zrt1p_anaer_rich[3]*100),as.numeric(Zrt1p_anaer_rich[4]*100)))
ZRT1p[4,3] <- sd(c(as.numeric(Zrt1p_anaer_min[2]*100),as.numeric(Zrt1p_anaer_min[3]*100),as.numeric(Zrt1p_anaer_min[4]*100)))
#
ZPS1p[1,2] <- mean(c(as.numeric(Zps1p_aer_rich[2]*100),as.numeric(Zps1p_aer_rich[3]*100)))
ZPS1p[2,2] <- mean(c(as.numeric(Zps1p_aer_min[2]*100),as.numeric(Zps1p_aer_min[3]*100)))
ZPS1p[3,2] <- mean(c(as.numeric(Zps1p_anaer_rich[2]*100),as.numeric(Zps1p_anaer_rich[3]*100),as.numeric(Zps1p_anaer_rich[4]*100)))
ZPS1p[4,2] <- mean(c(as.numeric(Zps1p_anaer_min[2]*100),as.numeric(Zps1p_anaer_min[3]*100),as.numeric(Zps1p_anaer_min[4]*100)))
ZPS1p[1,3] <- sd(c(as.numeric(Zps1p_aer_rich[2]*100),as.numeric(Zps1p_aer_rich[3]*100)))
ZPS1p[2,3] <- sd(c(as.numeric(Zps1p_aer_min[2]*100),as.numeric(Zps1p_aer_min[3]*100)))
ZPS1p[3,3] <- sd(c(as.numeric(Zps1p_anaer_rich[2]*100),as.numeric(Zps1p_anaer_rich[3]*100),as.numeric(Zps1p_anaer_rich[4]*100)))
ZPS1p[4,3] <- sd(c(as.numeric(Zps1p_anaer_min[2]*100),as.numeric(Zps1p_anaer_min[3]*100),as.numeric(Zps1p_anaer_min[4]*100)))
#
colnames(MET17p) <- c('Condition','Mass','Stdev_Mass')
colnames(CAR1p) <- c('Condition','Mass','Stdev_Mass')
colnames(ADH1p) <- c('Condition','Mass','Stdev_Mass')
colnames(ADH3p) <- c('Condition','Mass','Stdev_Mass')
colnames(ZRT1p) <- c('Condition','Mass','Stdev_Mass')
colnames(ZPS1p) <- c('Condition','Mass','Stdev_Mass')
#


### PLOT GATHERED DATA


# Plot Adh3p
ggplot(data=ADH3p, aes(x=Condition, y=Mass, fill=Condition)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5,size=0.2) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  scale_fill_manual(values = c('deepskyblue1','purple','darkorange1','darkorange4')) +
  labs(x=" ") +
  labs(y="Mass ratio of proteome (%)") +
  labs(fill=' ') +
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3),limits=c(0,0.3)) +
  theme(axis.text.x = element_text(face = c('bold'),size=7)) +
  theme(axis.title.x = element_text(face = c('bold'),size=7)) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y= element_blank()) +
  theme(legend.position = "none") +
  geom_errorbar(aes(ymin=ADH3p$Mass-ADH3p$Stdev_Mass, ymax=ADH3p$Mass+ADH3p$Stdev_Mass), width=.2,position=position_dodge(.9),size=0.3)+
  ggsave("plot_Adh3p.tiff", height=2, width=3, units='in', dpi=600)


# Plot Adh1p
ggplot(data=ADH1p, aes(x=Condition, y=Mass, fill=Condition)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5,size=0.2) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  scale_fill_manual(values = c('deepskyblue1','purple','darkorange1','darkorange4')) +
  labs(x=" ") +
  labs(y="Mass ratio of proteome (%)") +
  labs(fill=' ') +
  scale_y_continuous(breaks=c(0,1,2,3,4),limits=c(0,4)) +
  theme(axis.text.x = element_text(face = c('bold'),size=7)) +
  theme(axis.title.x = element_text(face = c('bold'),size=7)) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y= element_blank()) +
  theme(legend.position = "none") +
  geom_errorbar(aes(ymin=ADH1p$Mass-ADH1p$Stdev_Mass, ymax=ADH1p$Mass+ADH1p$Stdev_Mass), width=.2,position=position_dodge(.9),size=0.3)+
  ggsave("plot_Adh1p.tiff", height=2, width=3, units='in', dpi=600)


# Plot Zrt1p
ggplot(data=ZRT1p, aes(x=Condition, y=Mass, fill=Condition)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5,size=0.2) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  scale_fill_manual(values = c('deepskyblue1','purple','darkorange1','darkorange4')) +
  labs(x=" ") +
  labs(y="Mass ratio of proteome (%)") +
  labs(fill=' ') +
  scale_y_continuous(breaks=c(0,0.1,0.2,0.3),limits=c(0,0.3)) +
  theme(axis.text.x = element_text(face = c('bold'),size=7)) +
  theme(axis.title.x = element_text(face = c('bold'),size=7)) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y= element_blank()) +
  theme(legend.position = "none") +
  geom_errorbar(aes(ymin=ZRT1p$Mass-ZRT1p$Stdev_Mass, ymax=ZRT1p$Mass+ZRT1p$Stdev_Mass), width=.2,position=position_dodge(.9),size=0.3)+
  ggsave("plot_Zrt1p.tiff", height=2, width=3, units='in', dpi=600)


# Plot Zps1p
ggplot(data=ZPS1p, aes(x=Condition, y=Mass, fill=Condition)) +
  geom_bar(stat="identity", color="black", position=position_dodge(),width=0.5,size=0.2) +
  coord_flip() +
  theme_bw() +
  theme(panel.grid.minor.y = element_blank()) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.minor.x = element_blank()) +
  scale_fill_manual(values = c('deepskyblue1','purple','darkorange1','darkorange4')) +
  labs(x=" ") +
  labs(y="Mass ratio of proteome (%)") +
  labs(fill=' ') +
  scale_y_continuous(breaks=c(0,0.05,0.1,0.15),limits=c(0,0.15)) +
  theme(axis.text.x = element_text(face = c('bold'),size=7)) +
  theme(axis.title.x = element_text(face = c('bold'),size=7)) +
  theme(axis.text.y = element_blank()) +
  theme(axis.title.y= element_blank()) +
  theme(legend.position = "none") +
  geom_errorbar(aes(ymin=ZPS1p$Mass-ZPS1p$Stdev_Mass, ymax=ZPS1p$Mass+ZPS1p$Stdev_Mass), width=.2,position=position_dodge(.9),size=0.3)+
  ggsave("plot_Zps1p.tiff", height=2, width=3, units='in', dpi=600)






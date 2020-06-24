# 17 May 2019
# Figure S4 of "Proteome re-allocation from amino acid biosynthesis to ribosomes enables yeast to grow faster in rich media"

#setwd
setwd(dir="/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Kate")

#load libraries
library(pheatmap)
library(matrixTests)

# load data
Rich_anaer<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/glu_proteomics_anaerobic_rich.csv")
Min_anaer<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/glu_proteomics_anaerobic_min.csv")
Rich_aer<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/glu_proteomics_aerobic_rich.csv")
Min_aer<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/glu_proteomics_aerobic_min.csv")
GO_cell_AA_met_process<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Kate/cellular-amino-acid-metabolic-process.csv")
GO_AA_transport<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Kate/amino-acid-transport.csv")

# after calculating p-values
AA_properties<-read.csv(file = "/Users/katec/Box Sync/Manuscripts/in_progress/Johan/Rosemary/Kate/AA_properties.csv")

Rich_anaer$Type<-"Rich_anaer"
Min_anaer$Type<-"Min_anaer"
Rich_aer$Type<-"Rich_aer"
Min_aer$Type<-"Min_aer"

###### log10 transform each replicate

# Rich_anaer
Rich_anaer$log10_rich_anaer_rep1<-log10(Rich_anaer$Rich_1)
Rich_anaer$log10_rich_anaer_rep2<-log10(Rich_anaer$Rich_2)
Rich_anaer$log10_rich_anaer_rep3<-log10(Rich_anaer$Rich_3)

# Min_anaer
Min_anaer$log10_min_anaer_rep1<-log10(Min_anaer$Min_1)
Min_anaer$log10_min_anaer_rep2<-log10(Min_anaer$Min_2)
Min_anaer$log10_min_anaer_rep3<-log10(Min_anaer$Min_3)

# Rich_aer
Rich_aer$log10_rich_aer_rep1<-log10(Rich_aer$Rich_1)
Rich_aer$log10_rich_aer_rep2<-log10(Rich_aer$Rich_2)

# Min_aer
Min_aer$log10_min_aer_rep1<-log10(Min_aer$Min_1)
Min_aer$log10_min_aer_rep2<-log10(Min_aer$Min_2)

########################################################################################
##########################         FIGURE S4A              ##############################
########################################################################################


# Subset for cellular amino acid metabolic process genes

names(GO_cell_AA_met_process)<-"Protein"

Rich_anaer_AA_met_process<-Rich_anaer[which(Rich_anaer$Protein %in% GO_cell_AA_met_process$Protein),]
Min_anaer_AA_met_process<-Min_anaer[which(Min_anaer$Protein %in% GO_cell_AA_met_process$Protein),]
Rich_aer_AA_met_process<-Rich_aer[which(Rich_aer$Protein %in% GO_cell_AA_met_process$Protein),]
Min_aer_AA_met_process<-Min_aer[which(Min_aer$Protein %in% GO_cell_AA_met_process$Protein),]


Rich_aer_AA_met_process2<-Rich_aer_AA_met_process[,c(5,7:8)]
Min_aer_AA_met_process2<-Min_aer_AA_met_process[,c(5,7:8)]

Rich_min_aer_AA_met_process<-merge(Rich_aer_AA_met_process2, Min_aer_AA_met_process2, by="Protein")
row.names(Rich_min_aer_AA_met_process)<-Rich_min_aer_AA_met_process$Protein
Rich_min_aer_AA_met_process<-Rich_min_aer_AA_met_process[,-1]

##### determine p value for each row and subset based on this

row_t_welch<-row_t_welch(Rich_min_aer_AA_met_process[,1:2], Rich_min_aer_AA_met_process[,3:4])

Rich_min_aer_AA_met_process$pvalue<-row_t_welch$pvalue

Rich_min_aer_AA_met_process_SIG<-Rich_min_aer_AA_met_process[
  which(Rich_min_aer_AA_met_process$pvalue <0.01),]

Rich_min_aer_AA_met_process_SIG2<-Rich_min_aer_AA_met_process_SIG[,-5]

row.names(Rich_min_aer_AA_met_process_SIG2)

# format annotation for rows and columns 

AA_properties2<-AA_properties
AA_properties2<-as.data.frame(AA_properties2[,-1])
names(AA_properties2)<-"Properties"
row.names(AA_properties2)<-AA_properties$Gene

annotation_aer_col<-data.frame("Rep" = c(1:2,1:2),
                               "Medium" = c("Rich_aer","Rich_aer","Min_aer","Min_aer"),
                               "Name" = c("log10_rich_aer_rep1","log10_rich_aer_rep2",
                                          "log10_min_aer_rep1","log10_min_aer_rep2"))

row.names(annotation_aer_col)<-annotation_aer_col$Name
annotation_aer_col<-annotation_aer_col[,-3]
str(annotation_aer_col)
annotation_aer_col$Rep<-as.character(annotation_aer_col$Rep)

# plot
pheatmap(Rich_min_aer_AA_met_process_SIG2,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50),
         annotation_col=annotation_aer_col,
         annotation_row = AA_properties2,
         cutree_cols = 2,
         show_colnames=F,
         fontsize_row =5)

row.names(Rich_min_aer_AA_met_process_SIG2)



########################################################################################
##########################         FIGURE S4B              ##############################
########################################################################################


Rich_anaer_AA_met_process2<-Rich_anaer_AA_met_process[,c(6,8:10)]
Min_anaer_AA_met_process2<-Min_anaer_AA_met_process[,c(6,8:10)]

Rich_min_anaer_AA_met_process<-merge(Rich_anaer_AA_met_process2, Min_anaer_AA_met_process2, by="Protein")
row.names(Rich_min_anaer_AA_met_process)<-Rich_min_anaer_AA_met_process$Protein
Rich_min_anaer_AA_met_process<-Rich_min_anaer_AA_met_process[,-1]

##### determine p value for each row and subset based on this

row_t_welch<-row_t_welch(Rich_min_anaer_AA_met_process[,1:3], Rich_min_anaer_AA_met_process[,4:6])

Rich_min_anaer_AA_met_process$pvalue<-row_t_welch$pvalue

Rich_min_anaer_AA_met_process_SIG<-Rich_min_anaer_AA_met_process[
  which(Rich_min_anaer_AA_met_process$pvalue <0.001),]

Rich_min_anaer_AA_met_process_SIG2<-Rich_min_anaer_AA_met_process_SIG[,-7]

# format annotation for rows and columns 

AA_properties2<-AA_properties
AA_properties2<-as.data.frame(AA_properties2[,-1])
names(AA_properties2)<-"Properties"
row.names(AA_properties2)<-AA_properties$Gene

annotation_anaer_col<-data.frame("Rep" = c(1:3,1:3),
                             "Medium" = c("Rich_anaer","Rich_anaer","Rich_anaer","Min_anaer","Min_anaer","Min_anaer"),
                             "Name" = c("log10_rich_anaer_rep1","log10_rich_anaer_rep2","log10_rich_anaer_rep3",
                                        "log10_min_anaer_rep1","log10_min_anaer_rep2","log10_min_anaer_rep3"))

row.names(annotation_anaer_col)<-annotation_anaer_col$Name
annotation_anaer_col<-annotation_anaer_col[,-3]
str(annotation_anaer_col)
annotation_anaer_col$Rep<-as.character(annotation_anaer_col$Rep)

# plot
pheatmap(Rich_min_anaer_AA_met_process_SIG2,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50),
         annotation_col=annotation_anaer_col,
         annotation_row = AA_properties2,
         cutree_cols = 2,
         show_colnames=F,
         fontsize_row =5)

row.names(Rich_min_anaer_AA_met_process_SIG2)

########################################################################################
##########################         FIGURE S4C              ##############################
########################################################################################


# Subset for amino acid transport genes

names(GO_AA_transport)<-"Protein"

Rich_anaer_AA_transport<-Rich_anaer[which(Rich_anaer$Protein %in% GO_AA_transport$Protein),]
Min_anaer_AA_transport<-Min_anaer[which(Min_anaer$Protein %in% GO_AA_transport$Protein),]
Rich_aer_AA_transport<-Rich_aer[which(Rich_aer$Protein %in% GO_AA_transport$Protein),]
Min_aer_AA_transport<-Min_aer[which(Min_aer$Protein %in% GO_AA_transport$Protein),]

Rich_aer_AA_transport2<-Rich_aer_AA_transport[,c(5,7:8)]
Min_aer_AA_transport2<-Min_aer_AA_transport[,c(5,7:8)]

Rich_min_aer_AA_transport<-merge(Rich_aer_AA_transport2, Min_aer_AA_transport2, by="Protein")
row.names(Rich_min_aer_AA_transport)<-Rich_min_aer_AA_transport$Protein
Rich_min_aer_AA_transport<-Rich_min_aer_AA_transport[,-1]

##### determine p value for each row and subset based on this

row_t_welch<-row_t_welch(Rich_min_aer_AA_transport[,1:2], Rich_min_aer_AA_transport[,3:4])

Rich_min_aer_AA_transport$pvalue<-row_t_welch$pvalue

Rich_min_aer_AA_transport_SIG<-Rich_min_aer_AA_transport[
  which(Rich_min_aer_AA_transport$pvalue <0.05),]

Rich_min_aer_AA_transport_SIG2<-Rich_min_aer_AA_transport_SIG[,-5]

row.names(Rich_min_aer_AA_transport_SIG2)

annotation_aer_col<-data.frame("Rep" = c(1:2,1:2),
                               "Medium" = c("Rich_aer","Rich_aer","Min_aer","Min_aer"),
                               "Name" = c("log10_rich_aer_rep1","log10_rich_aer_rep2",
                                          "log10_min_aer_rep1","log10_min_aer_rep2"))

row.names(annotation_aer_col)<-annotation_aer_col$Name
annotation_aer_col<-annotation_aer_col[,-3]
str(annotation_aer_col)
annotation_aer_col$Rep<-as.character(annotation_aer_col$Rep)

# plot
pheatmap(Rich_min_aer_AA_transport_SIG2,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50),
         annotation_col=annotation_aer_col,
         cutree_cols = 2,
         show_colnames=F,
         fontsize_row =7)

row.names(Rich_min_aer_AA_transport_SIG2)

########################################################################################
##########################         FIGURE S4D              ##############################
########################################################################################


Rich_anaer_AA_transport2<-Rich_anaer_AA_transport[,c(6,8:10)]
Min_anaer_AA_transport2<-Min_anaer_AA_transport[,c(6,8:10)]

Rich_min_anaer_AA_transport<-merge(Rich_anaer_AA_transport2, Min_anaer_AA_transport2, by="Protein")
row.names(Rich_min_anaer_AA_transport)<-Rich_min_anaer_AA_transport$Protein
Rich_min_anaer_AA_transport<-Rich_min_anaer_AA_transport[,-1]

##### determine p value for each row and subset based on this

row_t_welch<-row_t_welch(Rich_min_anaer_AA_transport[,1:3], Rich_min_anaer_AA_transport[,4:6])

Rich_min_anaer_AA_transport$pvalue<-row_t_welch$pvalue

Rich_min_anaer_AA_transport_SIG<-Rich_min_anaer_AA_transport[
  which(Rich_min_anaer_AA_transport$pvalue <0.05),]

Rich_min_anaer_AA_transport_SIG2<-Rich_min_anaer_AA_transport_SIG[,-7]

# format annotation for columns 

annotation_anaer_col<-data.frame("Rep" = c(1:3,1:3),
                                 "Medium" = c("Rich_anaer","Rich_anaer","Rich_anaer","Min_anaer","Min_anaer","Min_anaer"),
                                 "Name" = c("log10_rich_anaer_rep1","log10_rich_anaer_rep2","log10_rich_anaer_rep3",
                                            "log10_min_anaer_rep1","log10_min_anaer_rep2","log10_min_anaer_rep3"))

row.names(annotation_anaer_col)<-annotation_anaer_col$Name
annotation_anaer_col<-annotation_anaer_col[,-3]
str(annotation_anaer_col)
annotation_anaer_col$Rep<-as.character(annotation_anaer_col$Rep)

# plot
pheatmap(Rich_min_anaer_AA_transport_SIG2,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         color = colorRampPalette(c("dodgerblue3", "white", "firebrick3"))(50),
         annotation_col=annotation_anaer_col,
         cutree_cols = 2,
         show_colnames=F,
         fontsize_row =7)

row.names(Rich_min_anaer_AA_transport_SIG2)




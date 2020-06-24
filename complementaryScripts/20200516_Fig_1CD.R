

#Fig 1C

#
# Creating the VennDiagram for overlap of quantified proteins in different samples
#

# The difference in numbers of proteins between different intersections has not been taken into account
# A tiff-file will be saved to a choosen directory, and the numbers between overlaps are correct
# ,however, manual correction of the image is required afterwards. This seemed alright, since we neglect
# to scale the areas to the actual numbers reported

# Notice that the position of samples is important, but the placement of samples in the actual picture does not
# Have the same order. Therefore the labels; R_aer (Rich aerobic), M_aer (Min aerobic), R_anaer (Rich anaerobic) and M_anaer (Min anaerobic)
# has been added to the output figure, to aid and not mix up the areas.

library(VennDiagram)

venn.diagram(
  x                  = list(a=as.character(glu_proteomics_aerobic_rich$Gene),c=as.character(glu_proteomics_anaerobic_min$Gene),b=as.character(glu_proteomics_aerobic_min$Gene),d=as.character(glu_proteomics_anaerobic_rich$Gene)),
  fill               = c('darkorange3',"deepskyblue1",'darkorange','darkorchid4'),
  main.fontface      = "plain",
  main.fontfamily    = "serif",
  main.col           = "black",
  alpha              = c(0.75, 0.65, 0.65, 0.75),
  cex                = 0.5,
  category.names     = c('R_aer','M_anaer','Min_aer','R_anaer'),
  cat.cex            = 0.5,
  imagetype          = "tiff",
  units              = "px",
  compression        = "lzw",
  height             = 1000,
  width              = 1000,
  resolution         = 500,
  lwd                = 1,
  
  filename           = "VennDiagram_of_quantified_proteins.tif"
)


#

#####PCA

#the below is the base plot, aesthetics adjustment is done elsewhere

print(prcomp(na.omit(Dataset_Proteomics[,4:13]), center = TRUE, scale. = TRUE)) #for standard deviations and rotations
summary(prcomp(na.omit(Dataset_Proteomics[,4:13]), center = TRUE, scale. = TRUE)) #for importance of components


plot(prcomp(na.omit(Dataset_Proteomics[,4:13]), center = TRUE, scale. = TRUE)$rotation[,1],
     -prcomp(na.omit(Dataset_Proteomics[,4:13]), center = TRUE, scale. = TRUE)$rotation[,2],
     xlab="PC1", ylab="PC2", xlim=c(0.31,0.321), pch=19)
text(prcomp(na.omit(Dataset_Proteomics[,4:13]), center = TRUE, scale. = TRUE)$rotation[,1],
     -prcomp(na.omit(Dataset_Proteomics[,4:13]), center = TRUE, scale. = TRUE)$rotation[,2],
     c(1,2,1,2,1:3,1:3), pos=2, cex=0.8)






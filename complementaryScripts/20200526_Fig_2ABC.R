


library(readxl)
AA_ana_std <- read_excel("20200525_anaerobic_AA_external_known.xlsx")
AA_ana_data_raw <- read_excel("20200525_anaerobic_AA_data.xlsx")

AA_aer_std <- read_excel("20200525_aerobic_AA_external_known.xlsx")
AA_aer_data_raw <- read_excel("20200525_aerobic_AA_data.xlsx")


### analyze anaerobic data

AA_ana_data <- AA_ana_data_raw[,1:2]
AA_ana_data[,3:16] <- NA
for(i in 3:16){
  #calculate sample ratio to IS
  AA_ana_data[,i] <- AA_ana_data_raw[,(i-1)*2]/AA_ana_data_raw[,((i-1)*2)-1]
  
  #build standard curve
  stdcurve_y <- na.omit(AA_ana_data[16:20,i])
  stdcurve_x <- na.omit(AA_ana_std[1:5,i-1])
  stdcurve_b <- lm(unlist(stdcurve_y) ~ unlist(stdcurve_x))$coefficients[1]
  stdcurve_m <- lm(unlist(stdcurve_y) ~ unlist(stdcurve_x))$coefficients[2]
  
  #calculate sample conc
  AA_ana_data[,i] <- (AA_ana_data[,i]-stdcurve_b)/stdcurve_m*2
  colnames(AA_ana_data)[i] <- paste(colnames(AA_ana_data_raw)[(i-1)*2], "uM", sep=".")
}

AA_ana_data <- AA_ana_data[1:15,]

#take mean of reps
AA_ana_data_mean <- AA_ana_data[1:5,]
AA_ana_data_mean[,2:16] <- NA
AA_ana_data_mean[,1] <- c("t1.mean", "t2.mean", "t3.mean", "t4.mean", "t5.mean")
for(i in 2:16){
  AA_ana_data_mean[,i] <- rowMeans(matrix(unlist(AA_ana_data[,i]), ncol=3,byrow=F))
}

#calculate yield
AA_ana_yield <- rep(0,14) 
names(AA_ana_yield) <- colnames(AA_ana_std)[2:15]
for(i in 1:14){
  AA_ana_yield[i] <- lm(unlist(AA_ana_data_mean[,i+2]) ~ unlist(AA_ana_data_mean[,2]))$coefficients[2]
}


###analyze aerobic data



AA_aer_data <- AA_aer_data_raw[,1:2]
AA_aer_data[,3:16] <- NA
for(i in 3:16){
  #calculate sample ratio to IS
  AA_aer_data[,i] <- AA_aer_data_raw[,(i-1)*2]/AA_aer_data_raw[,((i-1)*2)-1]
  
  #build standard curve
  stdcurve_y <- na.omit(AA_aer_data[16:19,i])
  stdcurve_x <- na.omit(AA_aer_std[1:5,i-1])
  stdcurve_b <- lm(unlist(stdcurve_y) ~ unlist(stdcurve_x))$coefficients[1]
  stdcurve_m <- lm(unlist(stdcurve_y) ~ unlist(stdcurve_x))$coefficients[2]
  
  #calculate sample conc
  AA_aer_data[,i] <- (AA_aer_data[,i]-stdcurve_b)/stdcurve_m*2
  colnames(AA_aer_data)[i] <- paste(colnames(AA_aer_data_raw)[(i-1)*2], "uM", sep=".")
}

AA_aer_data <- AA_aer_data[1:15,]

#take mean of reps
AA_aer_data_mean <- AA_aer_data[1:5,]
AA_aer_data_mean[,2:16] <- NA
AA_aer_data_mean[,1] <- c("t1.mean", "t2.mean", "t3.mean", "t4.mean", "t5.mean")
for(i in 2:16){
  AA_aer_data_mean[,i] <- rowMeans(matrix(unlist(AA_aer_data[,i]), ncol=3,byrow=F))
}

#calculate yield
AA_aer_yield <- rep(0,14) 
names(AA_aer_yield) <- colnames(AA_aer_std)[2:15]
for(i in 1:14){
  AA_aer_yield[i] <- lm(unlist(AA_aer_data_mean[,i+2]) ~ unlist(AA_aer_data_mean[,2]))$coefficients[2]
}



### Fig 2A

plot(AA_aer_yield, AA_ana_yield, xlim=c(-400,20), ylim=c(-900,50),type="n", axes=F, xlab="", ylab="")
abline(a=0,b=1)
abline(v=c(-400,-350,-300,-250, -200,-150,-100,-50,0), col="grey")
abline(h=c(-800,-700,-600,-500,-400,-300,-200,-100,0), col="grey")
abline(v=0, lty=5)
abline(h=0, lty=5)
par(new=T)
plot(AA_aer_yield, AA_ana_yield, xlim=c(-400,20), ylim=c(-900,50), las=1,
     pch=c("Q", "T", "K", "L", "R", "I", "M", "F", "V", "H", "W", "Y", "D", "G"))




### Fig 2B

AA_ecYeast_8 <- c(346.67, 219.86, 328.75, 340.47, 184.59, 221.35, 58.238, 153.81,
                  303.94, 76.157, 32.622, 117.16, 341.73, 333.58)
names(AA_ecYeast_8) <- colnames(AA_aer_std)[2:15]

plot(AA_ecYeast_8, AA_aer_yield, xlim=c(0,400), ylim=c(-800,50),type="n", axes=F, xlab="", ylab="")
abline(v=c(50,100,150,200,250,300,350,400), col="grey")
abline(h=c(-800,-700,-600,-500,-400,-300,-200,-100,0), col="grey")
abline(v=0, lty=5)
abline(h=0, lty=5)
abline(a=0,b=-1)
abline(a=0,b=-2, col="orange")
abline(a=0,b=-0.5, col="orange")
par(new=T)
plot(AA_ecYeast_8, AA_aer_yield, xlim=c(0,400), ylim=c(-800,50), las=1,
     pch=c("Q", "T", "K", "L", "R", "I", "M", "F", "V", "H", "W", "Y", "D", "G"))



### Fig 3C

plot(AA_ecYeast_8, AA_ana_yield, xlim=c(0,400), ylim=c(-800,50),type="n", axes=F, xlab="", ylab="")
abline(v=c(50,100,150,200,250,300,350,400), col="grey")
abline(h=c(-800,-700,-600,-500,-400,-300,-200,-100,0), col="grey")
abline(v=0, lty=5)
abline(h=0, lty=5)
abline(a=0,b=-1)
abline(a=0,b=-2, col="royalblue")
abline(a=0,b=-0.5, col="royalblue")
par(new=T)
plot(AA_ecYeast_8, AA_ana_yield, xlim=c(0,400), ylim=c(-800,50), las=1,
     pch=c("Q", "T", "K", "L", "R", "I", "M", "F", "V", "H", "W", "Y", "D", "G"))






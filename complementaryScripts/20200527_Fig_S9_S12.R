


par(oma=c(0,0,0,2)) #set margins
#hist(log10(iBAQ$fmol_ug_in_ref*1e6), xlim=c(2,10))
plot(density(log10(iBAQ$fmol_ug_in_ref*1e6)), xlim=c(2,10), ylim=c(-.032,.45), las=1, main="", 
     xlab="log10(iBAQ intensity)", ylab="density of detected proteins", col="red", lwd=2)
par(new=TRUE)
plot(UPS2$log10.mean.IBAQ, UPS2$log10.fmol, xlab="", ylab="",
     xlim=c(2,10), axes=F)

abline(a=summary(lm(UPS2$log10.fmol ~ UPS2$log10.mean.IBAQ))$coefficients[1,1],
       b=summary(lm(UPS2$log10.fmol ~ UPS2$log10.mean.IBAQ))$coefficients[2,1])
axis(4, las=1) #4th axis is the right y-axis
mtext("log10(fmol)", side=4, outer=TRUE)

par(oma=c(0,0,0,0)) #reset margins
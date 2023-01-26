Vm <- read.table(file="/home/drew/projects/kernik_2019/IPSC-model/iPSC_Baseline_Model_C_code/ys/Vm_3000ms.txt",header=TRUE)
t.ms <- Vm$t
Vm <- Vm$Vm
Xs.kernik <- read.table(file="/home/drew/projects/kernik_2019/IPSC-model/iPSC_Baseline_Model_C_code/ys/Xs_3000ms.txt",header=TRUE)
Xs.kernik <- Xs.kernik$Xs


Xs.test <- read.table("./test_Xs.txt",header=TRUE)
jgpplotter(t.ms,Xs.kernik,type="l",xlab="t (ms)",ylab="Xs")
lines(Xs.test$t,Xs.test$Xs_dclamp,col="red",lty=2)
legend(x="topleft",legend=c("Kernik","dclamp"),lty=c(1,2),col=c("black","red"),bty="n",text.col=c("black","red"))

tmp <- read.table("./test_data/test_IK1.txt",header=TRUE)
IK1.kernik <- tmp$IK1
IK1.predicted <- tmp$IK1_dclamp
rm(tmp)

ys <- read.table(file="~/projects/kernik_2019/IPSC-model/iPSC_Baseline_Model_C_code/ys/ys_labeled.txt",header=TRUE)

jgpplotter(t.ms,IK1.predicted,type="l",xlab="t (ms)",ylab=expression(I[K1]~nS/pF),col="red")
lines(t.ms,IK1.kernik,lwd=2)
lines(t.ms,ik1.func.EK,col="red",lty=2,lwd=1.5)
legend(x="topleft",legend=IK1.legend,lty=c(1,1,2),bty="n",col=c("black","red","red"),text.col=c("black","red","red"))

Xrx <- read.table("./test_data/test_Xrx_010720.txt",header=T)
jgpplotter(x=t.ms,y=Xrx$Xr1,type="l",xlab="t (ms)",ylab="Xr1 , Xr2",lwd=1.5)
lines(x=t.ms,y=Xrx$Xr1_dclamp,col="red",lty=2)
lines(x=t.ms,y=Xrx$Xr2,col="navy",lwd=1.5)
lines(x=t.ms,y=Xrx$Xr2_dclamp,col="orange",lty=2,lwd=1.5)

legend.text <- c("Kernik Xr1","Xr1-dclamp","Kernik Xr2","Xr2-dclamp")
legend.col <- c("black","red","navy","orange")
legend.lty <- c(1,2,1,2)

legend("topleft",legend=legend.text,col=legend.col,text.col=legend.col,bty="n",lty=legend.lty)


Ko <- 5.4
R <- 8.314472  # joule_per_mole_kelvin (in model_parameters)
T  <-  310.0   # kelvin (in model <- parameters)
F  <-  96.4853415


# multi axis plot #
jgpplotter(t.ms,Vm,type="l",axes=FALSE,xlab="t (ms)",ylab="",mar=c(5.1,5.1,5.1,5.1))
axis(side=1)
axis(side=2)
mtext(text="mV",side=2,line=2.5)
par(new=TRUE)
plot(x=t.ms,y=E.K,axes=FALSE,col="red",bty="n",type="l",ylab="",xlab="")
axis(side=4,col.axis="red",col="red")
mtext(text="E_K",side=4,line=2.5,col="red")


SR <- read.table("./test_data/test_SR.txt",header=T)
jgpplotter(x=t.ms,y=SR$s,type="l",xlab="t (ms)",ylab="S , R",lwd=1.5)
lines(x=t.ms,y=SR$S_dclamp,col="red",lty=2)
lines(x=t.ms,y=SR$r,col="navy",lwd=1.5)
lines(x=t.ms,y=SR$R_dclamp,col="orange",lty=2,lwd=1.5)

legend.text <- c("Kernik S","S-dclamp","Kernik R","R-dclamp")
legend.col <- c("black","red","navy","orange")
legend.lty <- c(1,2,1,2)

legend("topleft",legend=legend.text,col=legend.col,text.col=legend.col,bty="n",lty=legend.lty)
dev.copy(pdf,file="./plots/test_SR_022120.pdf")
dev.off()
dev.off()

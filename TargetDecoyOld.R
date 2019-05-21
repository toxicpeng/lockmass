###################
#this script is used for Target-Decoy Search and untargeted chemical analysis
###Hui, 20170921

####################for drinking water, cpp funtion, the oxygen number could be 1.2*carbon+3, cutoff e5, S/N>5, change is.BrCl
library(xcms)
library(MassSpecWavelet)
library(Rcpp)
library(RcppArmadillo)
library(isopat)
data(iso_list)
setwd("C:/Users/Steven/Dropbox/MSc/Mass Spec/Target Decoy")
source("TargetDecoyFun.r")
polarity<--1##if neg -1, if pos 1
LockMass<-c(131.9614,165.0188,172.956991,227.201104,255.232405,312.972263,411.12912)##Lock Mass, C8H5O4, C14H27O2, C16H31O2, humic acid
LockMass1<-206.8524
LockMass2<-216.8328
LockMass3<-216.8812
LockMass4<-222.8803
LockMass5<-232.9278
LockMass6<-248.8960
LockMass7<-250.8422
LockMass8<-256.8413
LockMass9<-258.9248



####################Mass Calibration###############
setwd("C:/Users/Steven/Dropbox/MSc/Mass Spec/Target Decoy/data")
msfiles<-list.files()
for (i in 1:length(msfiles)){
    setwd("C:/Users/Steven/Dropbox/MSc/Mass Spec/Target Decoy/data")
    xrawdata<-xcmsRaw(msfiles[i])
    xrawdata<-MassCal(xrawdata,LockMass)
    setwd("C:/Users/Steven/Dropbox/MSc/Mass Spec/Target Decoy/CalData")
    write.mzdata(xrawdata,msfiles[i])####save the calibrated data to new files
}

############read the library############
setwd("C:/procedure/R work/Target_Decoy/library")
Library<-read.table("library.csv",header=TRUE,sep=',')
TargetDatabase<-TargetConstruct(Library)
ImAdducts<-c(4.002603,7.016005,9.012183,11.009305,19.992439,26.981541,27.976928,39.962383,39.962591,44.955914,47.947947)##Only those ions with lesser than 50 m/z,implausible adducts
DecoyDatabase<-DecoyConstruct(TargetDatabase,ImAdducts)

############peaks detection##############################
setwd("C:/Users/skutarna/Dropbox/MSc/Mass Spec/Target Decoy/data/results")
xset.all<-xcmsSet(msfiles,method='centWave',ppm=2.5,peakwidth=c(5,10),snthresh=10,nSlaves=min(6,length(msfiles)))
xset_group<-group(xset.all,bw=30,minsamp=1,minfrac=0.1,mzwid=0.001)
len<-length(xset_group@groupidx)#group number
len2<-length(msfiles)#data files
Peak.init<-array(rep(0,len*(len2+3)),dim=c(len,(len2+3)))##columns are m/z, rt,sampleID 
for (i in 1:len){
    temp<-unlist(xset_group@groupidx[i])
    len3<-length(temp)
    for (j in 1:len3){
         index1<-xset_group@peaks[temp[j],11]
         Peak.init[i,1]<-xset_group@peaks[temp[j],1]##mz
         Peak.init[i,2]<-xset_group@peaks[temp[j],4]/60##rt
         Peak.init[i,3]<-xset_group@peaks[temp[j],11]###window id
         Peak.init[i,index1+3]<-max(xset_group@peaks[temp[j],9],Peak.init[i,index1+3])##intensity,select the maximum one if there are multiple peaks for one sample
         }}

for (i in 1:nrow(Peak.init)){####replace the sample id as the one with maximal abundance
    sample.id<-Peak.init[i,4:ncol(Peak.init)]
    Peak.init[i,3]<-which.max(sample.id)
}         

######################retention time alignment#################
setwd("C:/procedure/R work/Target_Decoy/library")
Reference<-read.table("Reference.csv",header=TRUE,sep=',')
ReferDatabase<-TargetConstruct(Reference)
setwd("C:/procedure/R work/Target_Decoy/CalData")
refermatch<-DatabaseMatch(Peak.init,ReferDatabase,2*10^(-6),msfiles)
index.refer<-which(refermatch[,10]>0)
refermatch<-refermatch[index.refer,]
for (i in 1:nrow(refermatch)){
     refermatch[i,3]<-ReferDatabase[refermatch[i,6],1] 
}

for (i in 1:nrow(refermatch)){
    reg<-lm(refermatch[,2]~refermatch[,3])
    coeff<-reg$coefficients
    res<-refermatch[,3]*coeff[2]+coeff[1]-refermatch[,2]
    index<-which.max(abs(res))##find the maximal residual, and delete it
    refermatch<-refermatch[-index,]
    if (coeff[2]>0.99||nrow(refermatch)<2){break}}
plot(refermatch[,3]*coeff[2]+coeff[1],refermatch[,2])
TargetDatabase[,1]<-TargetDatabase[,1]*coeff[2]+coeff[1]###correct the retention time of the database
DecoyDatabase[,1]<-DecoyDatabase[,1]*coeff[2]+coeff[1]    


############Database Match#########################
Target.score<-DatabaseMatch(Peak.init,TargetDatabase,2*10^(-6),msfiles)
Decoy.score<-DatabaseMatch(Peak.init,DecoyDatabase,2*10^(-6),msfiles)
index<-which(Target.score[,10]>1.5)
index1<-which(Decoy.score[,10]>1.5)
plot(Target.score[index,10])
points(Decoy.score[index1,10],col="red")

##########FDR calculation######################
score.save<-NULL
for (k in 1:100){##repeat the decoy database for 100 times, for random data
    print(paste('iteration...',k,sep=''))
    DecoyDatabase<-DecoyConstruct(TargetDatabase,ImAdducts)
    Decoy.score<-DatabaseMatch(Peak.init,DecoyDatabase,2*10^(-6),msfiles)
    for (j in 1:99){
        score.fdr<-0.03*j
        index<-which(Target.score[,10]>score.fdr)
        index1<-which(Decoy.score[,10]>score.fdr)
        if (length(index1)/length(index)<0.05){
        score.save<-c(score.save,score.fdr)
        break}
    }}
index<-which(Target.score[,10]>mean(score.save))
FinalID<-Target.score[index,]

########delete the duplicate information, the isotopic peaks##############
score.save<-NULL
ratio<-NULL
for (j in 1:99){
    score.fdr<-0.03*j
    score.save<-c(score.save,score.fdr)
    index<-which(Target.score[,10]>score.fdr)
    index1<-which(Decoy.score[,10]>score.fdr)
    ratio<-c(ratio,1-length(index1)/length(index))
}
plot(score.save,ratio)

  
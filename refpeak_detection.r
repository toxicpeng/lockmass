findref<-function(xset.input,ppm,btw){

ppm<-ppm/10^6

for (i in 1:length(xset.input@groupidx)){
  index<-unlist(xset.input@groupidx[i])
  sampleid<-xset.input@peaks[index,11]
  mz.value<-mean(xset.input@peaks[index,1])
  rt.value<-mean(xset.input@peaks[index,4])
  for (j in 1:length(unlist(phenoData(xset_pos)))){
    index<-which(sampleid==j)
    if (length(index)<1){
      idsave<-rbind(idsave,c(i,j,mz.value,rt.value))##save groupidx,sampleid
    }}}
for (n in 1:length(index)){
  mz.value<-idsave[index[n],3]
  mz.min<-mz.value-mz.value*ppm
  mz.min<-max(mz.min, xraw@mzrange[1])
  mz.max<-mz.value+mz.value*ppm
  mz.max<-min(mz.max,xraw@mzrange[2])
  rt.min<-min(idsave[,4])
  rt.max<-max(idsave[,4])

  peak<-getEIC(xraw,mzrange=cbind(mz.min,mz.max),rtrange=cbind(rt.min, rt.max),step=0.001)
  peak<-unlist(peak@eic[[1]])
  peak<-peak[((length(peak)/2)+1):length(peak)]
  for (k in 1:length(peak))
    
oil<-scan("oil.dat")
#oil data -remove outliers
est<-rep(0,4050)
for(i in 1:(4050)){
 s<-max(1,i-40)
 t<-min(4050,i+40)
est[i]<-median(oil[s:t])
residual<-oil-est
index<-(1:4050)[abs(residual)> 8000]
}
y<-oil
y<-y[-index] #remove outliers for simplicity
##y[index]<- -1 ##this would be better - treat outliers as missing data; initially just removing outliers is simplest

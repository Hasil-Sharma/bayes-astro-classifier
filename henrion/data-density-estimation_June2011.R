#####################
## QUICK DATA LOAD ##
#####################


## the model is fitted in 4 steps
## i) fit GM model for ClassStat noise distribution by ML
## ii) use UKIDSS pipeline posterior star proba to roughly identify galaxy locus (idx.gal = prob_star<0.1) and fit galaxy locus parameters by LS
## iii) use UKIDSS pipeline posterior star proba to get number counts (idx.gal = prob_star<0.5) and fit Number Count model by LS
## iv) fit magnitude distribution by ML

  temp.seed<-unlist(strsplit(date(),split=" "))
  temp.seed2<-unlist(strsplit(temp.seed[4],split=":"))
  temp.seed2<-sum(as.numeric(temp.seed2))
  temp.seed<-as.numeric(temp.seed[3])+as.numeric(temp.seed[5])
  temp.seed<-temp.seed*temp.seed2
  set.seed(temp.seed)
  rm(temp.seed,temp.seed2)

  source("mixfit.R")

  lasdat<-read.csv(file="sample_with-colours_no-miss_std.csv")
  datcl<-read.csv(file="SDSS Stripe 82/sources_classified_sdss.dat",sep=" ")
  #datsdss<-read.csv(file="SDSS Stripe 82/sources_sdss.dat",sep=",",nrows=75000)
  n<-nrow(lasdat)
  m<-ncol(lasdat)
  
  idx.miss.y<-(1:n)[lasdat$miss.y.bin==1]
  idx.miss.j<-(1:n)[lasdat$miss.j.bin==1]
  idx.miss.h<-(1:n)[lasdat$miss.h.bin==1]
  idx.miss.k<-(1:n)[lasdat$miss.k.bin==1]

##############################

## taking random sample
idx.rv<-sample(1:nrow(lasdat),size=20000)
lasdatrv<-lasdat[idx.rv,]

##############################





#########################################################
## mixture model for star locus and galaxy locus noise ##
#########################################################


  starpars<-vector("list",4)
  names(starpars)<-c("y","j","h","k")

## define dstar
  dstar<-function(x,starpars){
    res<-(starpars[1]*dnorm(x+starpars[2],sd=1) + (1-starpars[1])*dnorm(x,mean=starpars[3],sd=starpars[4]))
    res
  }

  # Yband
    dtemp<-lasdat[lasdat$mapp_y>13 & lasdat$mapp_y<17 & abs(lasdat$cs_y)<10 & lasdat$miss.y.bin==0 & lasdat$cs_y!=0,]
    starloglik<-function(para){
      res<-para[1]*dnorm(dtemp$cs_y+para[2],sd=1) + (1-para[1])*dnorm(dtemp$cs_y,mean=para[3],sd=para[4])
      res[res==0]<-10^(-200)
      res<-(-sum(log(res)))
      res
    }
    starpars$y<-optim(fn=starloglik,par=c(0.95,0.1,4,4),lower=c(0,-50,-50,0.01),upper=c(1,50,50,50),method="L-BFGS-B")
    #hist(lasdat$cs_y[lasdat$mapp_y>13 & lasdat$mapp_y<15 & lasdat$miss.y.bin==0 & abs(lasdat$cs_y)<20],breaks=100,freq=F,xlim=c(-10,20))
    #lines(mgvect,yvvect,col="red")
    #lines(mgvect,dnorm(mgvect),col="blue")

         ## ALTERNATIVE WAY -- GETTING STANDARD ERRORS
           library(MASS)
           dstaralt<-function(x,par1,par2,par3,par4){dstar(x,starpars=c(par1,par2,par3,par4))}
           startls<-vector("list",4); names(startls)<-c("par1","par2","par3","par4")
           startls$par1<-0.95; startls$par2<-0.1; startls$par3<-4; startls$par4<-4
           idx.t<-sample(x=(1:nrow(dtemp)),size=3000,replace=F)
           fitdistr(x=dtemp$cs_y[idx.t],densfun=dstaralt,start=startls)

  # Jband
    dtemp<-lasdat[lasdat$mapp_j>13 & lasdat$mapp_j<15 & lasdat$miss.j.bin==0 & abs(lasdat$cs_j)<12 & lasdat$cs_j!=0,]
    starloglik<-function(para){
      res<-para[1]*dnorm(dtemp$cs_j+para[2],sd=1) + (1-para[1])*dnorm(dtemp$cs_j,mean=para[3],sd=para[4])
      res[res==0]<-10^(-200)
      res<-(-sum(log(res)))
      res
    }
    starpars$j<-optim(fn=starloglik,par=c(0.95,0.1,4,4),lower=c(0,-50,-50,0.01),upper=c(1,50,50,50),method="L-BFGS-B")

         ## ALTERNATIVE WAY -- GETTING STANDARD ERRORS
           library(MASS)
           dstaralt<-function(x,par1,par2,par3,par4){dstar(x,starpars=c(par1,par2,par3,par4))}
           startls<-vector("list",4); names(startls)<-c("par1","par2","par3","par4")
           startls$par1<-0.95; startls$par2<-0.1; startls$par3<-4; startls$par4<-4
           idx.t<-sample(x=(1:nrow(dtemp)),size=3000,replace=F)
           fitdistr(x=dtemp$cs_j[idx.t],densfun=dstaralt,start=startls)

  # Hband
    dtemp<-lasdat[lasdat$mapp_h>12.5 & lasdat$mapp_h<15 & lasdat$miss.h.bin==0 & abs(lasdat$cs_h)<12 & lasdat$cs_h!=0,]
    starloglik<-function(para){
      res<-para[1]*dnorm(dtemp$cs_h+para[2],sd=1) + (1-para[1])*dnorm(dtemp$cs_h,mean=para[3],sd=para[4])
      res[res==0]<-10^(-200)
      res<-(-sum(log(res)))
      res
    }
    starpars$h<-optim(fn=starloglik,par=c(0.95,0.1,4,4),lower=c(0,-50,-50,0.01),upper=c(1,50,50,50),method="L-BFGS-B")

         ## ALTERNATIVE WAY -- GETTING STANDARD ERRORS
           library(MASS)
           dstaralt<-function(x,par1,par2,par3,par4){dstar(x,starpars=c(par1,par2,par3,par4))}
           startls<-vector("list",4); names(startls)<-c("par1","par2","par3","par4")
           startls$par1<-0.95; startls$par2<-0.1; startls$par3<-4; startls$par4<-4
           idx.t<-sample(x=(1:nrow(dtemp)),size=3000,replace=F)
           fitdistr(x=dtemp$cs_h[idx.t],densfun=dstaralt,start=startls)

  # Kband
    dtemp<-lasdat[lasdat$mapp_k>12.5 & lasdat$mapp_k<15 & lasdat$miss.k.bin==0 & abs(lasdat$cs_k)<12 & lasdat$cs_k!=0,]
    starloglik<-function(para){
      res<-para[1]*dnorm(dtemp$cs_k+para[2],sd=1) + (1-para[1])*dnorm(dtemp$cs_k,mean=para[3],sd=para[4])
      res[res==0]<-10^(-200)
      res<-(-sum(log(res)))
      res
    }
    starpars$k<-optim(fn=starloglik,par=c(0.95,0.1,4,4),lower=c(0,-50,-50,0.01),upper=c(1,50,50,50),method="L-BFGS-B")

         ## ALTERNATIVE WAY -- GETTING STANDARD ERRORS
           library(MASS)
           dstaralt<-function(x,par1,par2,par3,par4){dstar(x,starpars=c(par1,par2,par3,par4))}
           startls<-vector("list",4); names(startls)<-c("par1","par2","par3","par4")
           startls$par1<-0.95; startls$par2<-0.05; startls$par3<-4; startls$par4<-4
           idx.t<-sample(x=(1:nrow(dtemp)),size=3000,replace=F)
           fitdistr(x=dtemp$cs_k[idx.t],densfun=dstaralt,start=startls)
    
  # check
    hist(lasdat$cs_h[lasdat$mapp_h>12.5 & lasdat$mapp_h<15 & lasdat$miss.h.bin==0 & abs(lasdat$cs_h)<35 & lasdat$cs_h!=0],freq=F,breaks=100)
    cs<-seq(-10,40,length=500)
    funval<-(starpars$h$par[1]*dnorm(cs+starpars$h$par[2],sd=1) + (1-starpars$h$par[1])*dnorm(cs,mean=starpars$h$par[3],sd=starpars$h$par[4]))
    lines(cs,funval,col="red")
    rm(cs); rm(funval)



###########################
## estimate galaxy locus ##
###########################


idx.gal.y<-(1:nrow(lasdatrv))[lasdatrv$miss.y.bin==0 & lasdatrv$prob_star<0.1] 
idx.gal.j<-(1:nrow(lasdatrv))[lasdatrv$miss.j.bin==0 & lasdatrv$prob_star<0.1] 
idx.gal.h<-(1:nrow(lasdatrv))[lasdatrv$miss.h.bin==0 & lasdatrv$prob_star<0.1] 
idx.gal.k<-(1:nrow(lasdatrv))[lasdatrv$miss.k.bin==0 & lasdatrv$prob_star<0.1] 

idx.gal.rv.y<-idx.gal.y
idx.gal.rv.j<-idx.gal.j
idx.gal.rv.h<-idx.gal.h
idx.gal.rv.k<-idx.gal.k


SmpSize<-50
Nbins_y<-trunc(length(idx.gal.y)/SmpSize)
Nbins_j<-trunc(length(idx.gal.j)/SmpSize)
Nbins_h<-trunc(length(idx.gal.h)/SmpSize)
Nbins_k<-trunc(length(idx.gal.k)/SmpSize)
rnk_y<-rank(lasdatrv$mapp_y[idx.gal.y],ties.method="random")
rnk_j<-rank(lasdatrv$mapp_j[idx.gal.j],ties.method="random")
rnk_h<-rank(lasdatrv$mapp_h[idx.gal.h],ties.method="random")
rnk_k<-rank(lasdatrv$mapp_k[idx.gal.k],ties.method="random")

meanmag<-vector("list",4); names(meanmag)<-c("y","j","h","k")
var.cs<-vector("list",4); names(var.cs)<-c("y","j","h","k")
var.mag<-vector("list",4); names(var.mag)<-c("y","j","h","k")
med<-vector("list",4); names(med)<-c("y","j","h","k")
meanmag$y<-rep(NA,Nbins_y+1); meanmag$j<-rep(NA,Nbins_j+1); meanmag$h<-rep(NA,Nbins_h+1); meanmag$k<-rep(NA,Nbins_k+1)
var.cs$y<-rep(NA,Nbins_y+1); var.cs$j<-rep(NA,Nbins_j+1); var.cs$h<-rep(NA,Nbins_h+1); var.cs$k<-rep(NA,Nbins_k+1)
var.mag$y<-rep(NA,Nbins_y+1); var.mag$j<-rep(NA,Nbins_j+1); var.mag$h<-rep(NA,Nbins_h+1); var.mag$k<-rep(NA,Nbins_k+1)
med$y<-rep(NA,Nbins_y+1); med$j<-rep(NA,Nbins_j+1); med$h<-rep(NA,Nbins_h+1); med$k<-rep(NA,Nbins_k+1)

for(band in 1:4){
  if(band==1){
    dtemp<-lasdatrv[idx.gal.y,]
    for(i in 1:Nbins_y){
      if(i==1){magdata.idx<-(1:length(idx.gal.y))[rnk_y<=(SmpSize*i)]}
      if(i>1 & i<Nbins_y){magdata.idx<-(1:length(idx.gal.y))[rnk_y>(SmpSize*(i-1)) & rnk_y<=(SmpSize*i)]}
      if(i==Nbins_y){magdata.idx<-(1:length(idx.gal.y))[rnk_y>(SmpSize*(i-1))]}
      meanmag$y[i]<-mean(dtemp$mapp_y[magdata.idx])
      var.cs$y[i]<-var(dtemp$cs_y[magdata.idx])
      var.mag$y[i]<-var(dtemp$mapp_y[magdata.idx])
      med$y[i]<-mean(dtemp$cs_y[magdata.idx])
    }
    idx.max<-(1:nrow(lasdat))[lasdat$miss.y.bin==0 & lasdat$mapp_y==max(lasdat$mapp_y)]
    idx.sim<-(1:nrow(lasdat))[lasdat$mapp_y==max(lasdat$mapp_y)]
    idx.max<-c(idx.max,(1:nrow(lasdat))[lasdat$miss.y.bin==0 & lasdat$mapp_y==max(lasdat$mapp_y[-idx.sim])])
    idx.sim<-unique(c(idx.sim,idx.max))
    idx.max<-c(idx.max,(1:nrow(lasdat))[lasdat$miss.y.bin==0 & lasdat$mapp_y==max(lasdat$mapp_y[-idx.sim])])
    dtemp<-lasdat[idx.max,]
    meanmag$y[Nbins_y+1]<-mean(dtemp$mapp_y)
    var.cs$y[Nbins_y+1]<-var(dtemp$cs_y)
    var.mag$y[Nbins_y+1]<-var(dtemp$mapp_y)
    med$y[Nbins_y+1]<-mean(dtemp$cs_y)
  }
  if(band==2){
    dtemp<-lasdatrv[idx.gal.j,]
    for(i in 1:Nbins_j){
      if(i==1){magdata.idx<-(1:length(idx.gal.j))[rnk_j<=(SmpSize*i)]}
      if(i>1 & i<Nbins_j){magdata.idx<-(1:length(idx.gal.j))[rnk_j>(SmpSize*(i-1)) & rnk_j<=(SmpSize*i)]}
      if(i==Nbins_j){magdata.idx<-(1:length(idx.gal.j))[rnk_j>(SmpSize*(i-1))]}
      meanmag$j[i]<-mean(dtemp$mapp_j[magdata.idx])
      var.cs$j[i]<-var(dtemp$cs_j[magdata.idx])
      var.mag$j[i]<-var(dtemp$mapp_j[magdata.idx])
      med$j[i]<-mean(dtemp$cs_j[magdata.idx])
    }
    idx.max<-(1:nrow(lasdat))[lasdat$miss.j.bin==0 & lasdat$mapp_j==max(lasdat$mapp_j)]
    idx.sim<-(1:nrow(lasdat))[lasdat$mapp_j==max(lasdat$mapp_j)]
    idx.max<-c(idx.max,(1:nrow(lasdat))[lasdat$miss.j.bin==0 & lasdat$mapp_j==max(lasdat$mapp_j[-idx.sim])])
    idx.sim<-unique(c(idx.sim,idx.max))
    idx.max<-c(idx.max,(1:nrow(lasdat))[lasdat$miss.j.bin==0 & lasdat$mapp_j==max(lasdat$mapp_j[-idx.sim])])
    dtemp<-lasdat[idx.max,]
    meanmag$j[Nbins_j+1]<-mean(dtemp$mapp_j)
    var.cs$j[Nbins_j+1]<-var(dtemp$cs_j)
    var.mag$j[Nbins_j+1]<-var(dtemp$mapp_j)
    med$j[Nbins_j+1]<-mean(dtemp$cs_j)
  }
  if(band==3){
    dtemp<-lasdatrv[idx.gal.h,]
    for(i in 1:Nbins_h){
      if(i==1){magdata.idx<-(1:length(idx.gal.h))[rnk_h<=(SmpSize*i)]}
      if(i>1 & i<Nbins_h){magdata.idx<-(1:length(idx.gal.h))[rnk_h>(SmpSize*(i-1)) & rnk_h<=(SmpSize*i)]}
      if(i==Nbins_h){magdata.idx<-(1:length(idx.gal.h))[rnk_h>(SmpSize*(i-1))]}
      meanmag$h[i]<-mean(dtemp$mapp_h[magdata.idx])
      var.cs$h[i]<-var(dtemp$cs_h[magdata.idx])
      var.mag$h[i]<-var(dtemp$mapp_h[magdata.idx])
      med$h[i]<-mean(dtemp$cs_h[magdata.idx])
    }
    idx.max<-(1:nrow(lasdat))[lasdat$miss.h.bin==0 & lasdat$mapp_h==max(lasdat$mapp_h)]
    idx.sim<-(1:nrow(lasdat))[lasdat$mapp_h==max(lasdat$mapp_h)]
    idx.max<-c(idx.max,(1:nrow(lasdat))[lasdat$miss.h.bin==0 & lasdat$mapp_h==max(lasdat$mapp_h[-idx.sim])])
    idx.sim<-unique(c(idx.sim,idx.max))
    idx.max<-c(idx.max,(1:nrow(lasdat))[lasdat$miss.h.bin==0 & lasdat$mapp_h==max(lasdat$mapp_h[-idx.sim])])
    dtemp<-lasdat[idx.max,]
    meanmag$h[Nbins_h+1]<-mean(dtemp$mapp_h)
    var.cs$h[Nbins_h+1]<-var(dtemp$cs_h)
    var.mag$h[Nbins_h+1]<-var(dtemp$mapp_h)
    med$h[Nbins_h+1]<-mean(dtemp$cs_h)
  }
  if(band==4){
    for(i in 1:Nbins_k){
    dtemp<-lasdatrv[idx.gal.k,]
      if(i==1){magdata.idx<-(1:length(idx.gal.k))[rnk_k<=(SmpSize*i)]}
      if(i>1 & i<Nbins_k){magdata.idx<-(1:length(idx.gal.k))[rnk_k>(SmpSize*(i-1)) & rnk_k<=(SmpSize*i)]}
      if(i==Nbins_k){magdata.idx<-(1:length(idx.gal.k))[rnk_k>(SmpSize*(i-1))]}
      meanmag$k[i]<-mean(dtemp$mapp_k[magdata.idx])
      var.cs$k[i]<-var(dtemp$cs_k[magdata.idx])
      var.mag$k[i]<-var(dtemp$mapp_k[magdata.idx])
      med$k[i]<-mean(dtemp$cs_k[magdata.idx])
    }
    idx.max<-(1:nrow(lasdat))[lasdat$miss.k.bin==0 & lasdat$mapp_k==max(lasdat$mapp_k)]
    idx.sim<-(1:nrow(lasdat))[lasdat$mapp_k==max(lasdat$mapp_k)]
    idx.max<-c(idx.max,(1:nrow(lasdat))[lasdat$miss.k.bin==0 & lasdat$mapp_k==max(lasdat$mapp_k[-idx.sim])])
    idx.sim<-unique(c(idx.sim,idx.max))
    idx.max<-c(idx.max,(1:nrow(lasdat))[lasdat$miss.k.bin==0 & lasdat$mapp_k==max(lasdat$mapp_k[-idx.sim])])
    dtemp<-lasdat[idx.max,]
    meanmag$k[Nbins_k+1]<-mean(dtemp$mapp_k)
    var.cs$k[Nbins_k+1]<-var(dtemp$cs_k)
    var.mag$k[Nbins_k+1]<-var(dtemp$mapp_k)
    med$k[Nbins_k+1]<-mean(dtemp$cs_k)
  }
}


## Estimating means or medians & variances for true galactic locus

  SSmedCSHyp<-function(para){
    f.idx<-round(Nbins_y/10)
    check<-0
    para[3]<-142
    if(check==0){
      mmm<-max( max(lasdat$mapp_y[lasdat$miss.y.bin==0]), max(lasdat$mapp_j[lasdat$miss.j.bin==0]) + para[6]-para[7], max(lasdat$mapp_h[lasdat$miss.h.bin==0]) + para[6]-para[8], max(lasdat$mapp_k[lasdat$miss.k.bin==0]) + para[6]-para[9])
      res1<-(med$y[f.idx:(Nbins_y+1)] - ( (1-(meanmag$y[f.idx:(Nbins_y+1)]+(para[6]-para[6]))/mmm)*( (para[1]*(meanmag$y[f.idx:(Nbins_y+1)]-para[6])^2 + para[2]*(meanmag$y[f.idx:(Nbins_y+1)]-para[6]) + para[3])^(para[5]) + para[4]) ) )^2
      res2<-(med$j[f.idx:(Nbins_j+1)] - ( (1-(meanmag$j[f.idx:(Nbins_j+1)]+(para[6]-para[7]))/mmm)*( (para[1]*(meanmag$j[f.idx:(Nbins_j+1)]-para[7])^2 + para[2]*(meanmag$j[f.idx:(Nbins_j+1)]-para[7]) + para[3])^(para[5]) + para[4]) ) )^2
      res3<-(med$h[f.idx:(Nbins_h+1)] - ( (1-(meanmag$h[f.idx:(Nbins_h+1)]+(para[6]-para[8]))/mmm)*( (para[1]*(meanmag$h[f.idx:(Nbins_h+1)]-para[8])^2 + para[2]*(meanmag$h[f.idx:(Nbins_h+1)]-para[8]) + para[3])^(para[5]) + para[4]) ) )^2
      res4<-(med$k[f.idx:(Nbins_k+1)] - ( (1-(meanmag$k[f.idx:(Nbins_k+1)]+(para[6]-para[9]))/mmm)*( (para[1]*(meanmag$k[f.idx:(Nbins_k+1)]-para[9])^2 + para[2]*(meanmag$k[f.idx:(Nbins_k+1)]-para[9]) + para[3])^(para[5]) + para[4]) ) )^2
      res<-sum(res1) + sum(res2) + sum(res3) + sum(res4)
    }
    res
  }

  gr.SSmedCSHyp<-function(para){
    res<-grad(SSmedCSHyp,para)
    res
  }


  ## NB para[3] has been fixed at 142. This parameter is in fact redundant in the model (making it technically non identifiable if non-fixed).
  ## anyway the above is an empirical function, simply chosen to give a good fit; feel free to use another empirical function.
  MedPars<-optim(c(10,8,142,-56,0.8,23,22,21,20),SSmedCSHyp,method="BFGS",gr=gr.SSmedCSHyp,control=list(maxit=250))
  MedPars$par; MedPars$val; MedPars$convergence #

  mmagshift<-vector("list",4)
  names(mmagshift)<-c("y","j","h","k")
  mmagshift$y<-MedPars$par[6]-MedPars$par[6]
  mmagshift$j<-MedPars$par[6]-MedPars$par[7]
  mmagshift$h<-MedPars$par[6]-MedPars$par[8]
  mmagshift$k<-MedPars$par[6]-MedPars$par[9]

  pa<-MedPars$par
  mmm<-max( max(lasdat$mapp_y[lasdat$miss.y.bin==0]), max(lasdat$mapp_j[lasdat$miss.j.bin==0]) + pa[6]-pa[7], max(lasdat$mapp_h[lasdat$miss.h.bin==0]) + pa[6]-pa[8], max(lasdat$mapp_k[lasdat$miss.k.bin==0]) + pa[6]-pa[9])
  m<-seq(11,23,length=250)
  plot(c(13,23),c(-2,25),type="n"); points(meanmag$y,med$y); points(meanmag$j,med$j,col="brown"); points(meanmag$h,med$h,col="orange"); points(meanmag$k,med$k,col="yellow")
  yyy<-(1-(m+pa[6]-pa[6])/mmm)*((pa[1]*(m-pa[6])^2 + pa[2]*(m-pa[6]) + pa[3])^pa[5] + pa[4]); lines(m,yyy,col="red")
  yyy<-(1-(m+pa[6]-pa[7])/mmm)*((pa[1]*(m-pa[7])^2 + pa[2]*(m-pa[7]) + pa[3])^pa[5] + pa[4]); lines(m,yyy,col="blue")
  yyy<-(1-(m+pa[6]-pa[8])/mmm)*((pa[1]*(m-pa[8])^2 + pa[2]*(m-pa[8]) + pa[3])^pa[5] + pa[4]); lines(m,yyy,col="green")
  yyy<-(1-(m+pa[6]-pa[9])/mmm)*((pa[1]*(m-pa[9])^2 + pa[2]*(m-pa[9]) + pa[3])^pa[5] + pa[4]); lines(m,yyy,col="violet")

  SSvariCS<-function(para){
    f.idx<-round(Nbins_y/70)
    check<-0
    if(para[2]<(-1)){res<-10^(12); check<-1}
    if(check==0){
      res1<-(var.cs$y[f.idx:(Nbins_y-1)] -1 - para[1]*(10^5)*((10)^(para[2]*(meanmag$y[f.idx:(Nbins_y-1)] + mmagshift$y - 11))) )^2
      res2<-(var.cs$j[f.idx:(Nbins_j-1)] -1 - para[1]*(10^5)*((10)^(para[2]*(meanmag$j[f.idx:(Nbins_j-1)] + mmagshift$j - 11))) )^2
      res3<-(var.cs$h[f.idx:(Nbins_h-1)] -1 - para[1]*(10^5)*((10)^(para[2]*(meanmag$h[f.idx:(Nbins_h-1)] + mmagshift$h - 11))) )^2
      res4<-(var.cs$k[f.idx:(Nbins_k-1)] -1 - para[1]*(10^5)*((10)^(para[2]*(meanmag$k[f.idx:(Nbins_k-1)] + mmagshift$k - 11))) )^2    
      res<-sum(res1) + sum(res2) + sum(res3) + sum(res4)
    }
    res
  }

  gr.SSvarCS<-function(para){
    res<-grad(SSvariCS,para)
    res
  }

  VarPars<-optim(c(1.5,-0.5),SSvariCS,method="BFGS",gr=gr.SSvarCS)

  ## NB the model pa[1]*(10^5)*10(pa[2]*(m-11)) is equivalent to the model 10^(pa[2]*(m-pa'[1]))
  ## where pa'[1]=11-(5+log10(pa[1]))/pa[2]
  ## we have given the simpler model (with pa[1]) in the thesis...
  ## NBB whether we specify that sigma=10^pa[2]*(m-pa[1]) or sigma^2=10^pa'[2]*(m-pa[1]) is unimportant: pa'[2]=2*pa[2]


  par(mfrow=c(1,1))
  plot(meanmag$y+mmagshift$y,var.cs$y-1)
  points(meanmag$j+mmagshift$j,var.cs$j-1,col="brown")
  points(meanmag$h+mmagshift$h,var.cs$h-1,col="orange")
  points(meanmag$k+mmagshift$k,var.cs$k-1,col="yellow")
  m<-seq(14,23,length=250)
  vcs<-VarPars$par[1]*(10^5)*((10)^(VarPars$par[2]*(m-11)))
  lines(m,vcs,col="red")


################################
## fit magnitude distribution ##
################################

  dtemp<-c(lasdatrv[lasdatrv$miss.y.bin==0,]$mapp_y+mmagshift$y,lasdatrv[lasdatrv$miss.j.bin==0,]$mapp_j+mmagshift$j,lasdatrv[lasdatrv$miss.h.bin==0,]$mapp_h+mmagshift$h,lasdatrv[lasdatrv$miss.k.bin==0,]$mapp_k+mmagshift$k)
  MagParsGl<-FitMagInit(dtemp,c(20,0.45,0.25,12))



###################
## number counts ##
###################

## NB as we want to get parameters that hold true simultaneously for all 4 bands, we use stars and galaxies for all bands
## this means that any object can be counted up to 4 times during the estimation procedure.
## however, this is more a weighting problem than anything else: as we normalise to get the final a(m), the number of objects (or their sum of weights) is not of relevance

idx.gal.y<-(1:nrow(lasdatrv))[lasdatrv$miss.y.bin==0 & lasdatrv$prob_star<0.5] 
idx.gal.j<-(1:nrow(lasdatrv))[lasdatrv$miss.j.bin==0 & lasdatrv$prob_star<0.5] 
idx.gal.h<-(1:nrow(lasdatrv))[lasdatrv$miss.h.bin==0 & lasdatrv$prob_star<0.5] 
idx.gal.k<-(1:nrow(lasdatrv))[lasdatrv$miss.k.bin==0 & lasdatrv$prob_star<0.5] 

m<-seq(12,23,length=2000)

  dtemp<-vector("list",2)
  names(dtemp)<-c("st","ga")
  sty<-lasdatrv$mapp_y[-idx.gal.y]; sty<-sty[lasdatrv$miss.y.bin[-idx.gal.y]==0 & (lasdatrv$mapp_y[-idx.gal.y]+mmagshift$y)>15 & (lasdatrv$mapp_y[-idx.gal.y]+mmagshift$y)<21]
  stj<-lasdatrv$mapp_j[-idx.gal.j]; stj<-stj[lasdatrv$miss.j.bin[-idx.gal.j]==0 & (lasdatrv$mapp_j[-idx.gal.j]+mmagshift$j)>15 & (lasdatrv$mapp_j[-idx.gal.j]+mmagshift$j)<21]
  sth<-lasdatrv$mapp_h[-idx.gal.h]; sth<-sth[lasdatrv$miss.h.bin[-idx.gal.h]==0 & (lasdatrv$mapp_h[-idx.gal.h]+mmagshift$h)>15 & (lasdatrv$mapp_h[-idx.gal.h]+mmagshift$h)<21]
  stk<-lasdatrv$mapp_k[-idx.gal.k]; stk<-stk[lasdatrv$miss.k.bin[-idx.gal.k]==0 & (lasdatrv$mapp_k[-idx.gal.k]+mmagshift$k)>15 & (lasdatrv$mapp_k[-idx.gal.k]+mmagshift$k)<21]
  dtemp$st<-c(sty+mmagshift$y,stj+mmagshift$j,sth+mmagshift$h,stk+mmagshift$k)
  gay<-lasdatrv$mapp_y[idx.gal.y]; gay<-gay[(lasdatrv$mapp_y[idx.gal.y]+mmagshift$y)>15 & (lasdatrv$mapp_y[idx.gal.y]+mmagshift$y)<21]
  gaj<-lasdatrv$mapp_j[idx.gal.j]; gaj<-gaj[(lasdatrv$mapp_j[idx.gal.j]+mmagshift$j)>15 & (lasdatrv$mapp_j[idx.gal.j]+mmagshift$j)<21]
  gah<-lasdatrv$mapp_h[idx.gal.h]; gah<-gah[(lasdatrv$mapp_h[idx.gal.h]+mmagshift$h)>15 & (lasdatrv$mapp_h[idx.gal.h]+mmagshift$h)<21]
  gak<-lasdatrv$mapp_k[idx.gal.k]; gak<-gak[(lasdatrv$mapp_k[idx.gal.k]+mmagshift$k)>15 & (lasdatrv$mapp_k[idx.gal.k]+mmagshift$k)<21]
  dtemp$ga<-c(gay+mmagshift$y,gaj+mmagshift$j,gah+mmagshift$h,gak+mmagshift$k)
  length(sty)+length(gay); sum(lasdatrv$miss.y.bin==0 & (lasdatrv$mapp_y+mmagshift$y)>15 & (lasdatrv$mapp_y+mmagshift$y)<21)
  length(stj)+length(gaj); sum(lasdatrv$miss.j.bin==0 & (lasdatrv$mapp_j+mmagshift$j)>15 & (lasdatrv$mapp_j+mmagshift$j)<21)
  length(sth)+length(gah); sum(lasdatrv$miss.h.bin==0 & (lasdatrv$mapp_h+mmagshift$h)>15 & (lasdatrv$mapp_h+mmagshift$h)<21)
  length(stk)+length(gak); sum(lasdatrv$miss.k.bin==0 & (lasdatrv$mapp_k+mmagshift$k)>15 & (lasdatrv$mapp_k+mmagshift$k)<21)
  rm(sty,stj,sth,stk,gay,gaj,gah,gak)

  hall<-hist(c(dtemp$st,dtemp$ga),breaks=250,plot=F)
  hst<-hist(dtemp$st,breaks=250,plot=F) ## stars
  hga<-hist(dtemp$ga,breaks=250,plot=F) ## galaxies

  plot(hall$mids,log(hall$counts),type="l",main="number counts (log scale)")
  lines(hst$mids,log(hst$counts),type="l",col="blue")
  lines(hga$mids,log(hga$counts),type="l",col="red")

  SSst<-function(para){
    check<-0
    if(para[1]<0.0001){check<-1; res<-10^(20)}
    if(check==0){
      m<-hst$mids[hst$counts>0]
      y<-log(hst$counts[hst$counts>0])
      res<-sum( ( y - log( 10^(para[1]*(m-para[2])) * 0.5*erfc((m-para[3])/para[4]) ) )^2 )
    }
      res 
  }

  gr.st<-function(para){
    res<-grad(SSst,para)
    res
  }

  SSga<-function(para){
    check<-0
    if(para[1]<0.0001){check<-1; res<-10^(20)}
    if(check==0){
      m<-hga$mids[hga$counts>0]
      y<-log(hga$counts[hga$counts>0])
      res<-sum( ( y - log( 10^(para[1]*(m-para[2])) * 0.5*erfc((m-para[3])/para[4]) ) )^2 )
    }
    res 
  }

  gr.ga<-function(para){
    res<-grad(SSga,para)
    res
  }

  stpars<-optim(c(0.25,5,20,0.5),SSst,method="BFGS")
  gapars<-optim(c(0.45,15,20,0.5),SSga,method="BFGS")
  
  lst<-10^(stpars$par[1]*(m-stpars$par[2])) * 0.5*erfc((m-stpars$par[3])/stpars$par[4])
  lga<-10^(gapars$par[1]*(m-gapars$par[2])) * 0.5*erfc((m-gapars$par[3])/gapars$par[4])
  lines(m,log(lst),col="turquoise")														
  lines(m,log(lga),col="orange")

  MxCfGl<-function(m){
    res <- ( 10^(stpars$par[1]*(m-stpars$par[2]))*0.5*erfc((m-stpars$par[3])/stpars$par[4]) ) / ( (10^(stpars$par[1]*(m-stpars$par[2])))*0.5*erfc((m-stpars$par[3])/stpars$par[4]) + 10^(gapars$par[1]*(m-gapars$par[2]))*0.5*erfc((m-gapars$par[3])/gapars$par[4]) )
    #res <- ( 10^(stpars$par[1]*(m-stpars$par[2])) ) / ( (10^(stpars$par[1]*(m-stpars$par[2]))) + 10^(gapars$par[1]*(m-gapars$par[2])) )
    res
  }

  plot(m,MxCfGl(m),type="l")
  lines(hall$mids,(hst$counts/(hst$counts+hga$counts)),col="red")



###################################
## number counts plot for Thesis ##
###################################

## Fig 1 with numbers counts as used to fit number counts model
  #postscript("Thesis/graphs/D_Other_NumberCounts_SDSScountsUKIDSSmodelparameters_amended.eps", fonts=c("serif"),horizontal=F,paper="special",width=10,height=10)
  par(mfrow=c(1,1),pty="s",mai=c(1,1.25,0.5,0))

  allnorm<-nrow(datcl)/14.38 # normalisation constant; calculated from data sample properties
    #(NB area of histogram gets normalised to an area of 1 by using $density rather than $counts below)
    #14.38 = deg^2 area of sample
  stagalnorm<-(sum(hst$counts)/sum(hst$density))/(sum(hall$counts)/sum(hall$density))
  galstanorm<-(sum(hga$counts)/sum(hga$density))/(sum(hall$counts)/sum(hall$density))
  plot(hall$mids,log10(hall$density*allnorm),type="l",main="",xlab=expression(italic(Y)),ylab=expression(paste(log[10]," ",bgroup("(",paste(frac(paste(d,italic(N)),paste(d,italic(Y))),~~deg^{-2}),")"))),cex.axis=1.75,cex.main=1.75,cex.lab=1.75,xlim=c(15.5,21),lwd=3,family=c("serif"))
  lines(hst$mids,log10(hst$density*allnorm*stagalnorm),type="l",col="blue",lwd=3)
  lines(hga$mids,log10(hga$density*allnorm*galstanorm),type="l",col="red",lwd=3)
  
  mm<-seq(15,21,length=250)
  lst<-10^(stpars$par[1]*(mm-stpars$par[2])) * 0.5*erfc((mm-stpars$par[3])/stpars$par[4])
  area.st<-sum(lst)*((21-15)/249)
  lst<-lst/area.st
  lga<-10^(gapars$par[1]*(mm-gapars$par[2])) * 0.5*erfc((mm-gapars$par[3])/gapars$par[4])
  area.ga<-sum(lga)*((21-15)/249)
  lga<-lga/area.ga
  lines(mm,log10(lst*allnorm*stagalnorm),col="turquoise",lty=2,lwd=4)														
  lines(mm,log10(lga*allnorm*galstanorm),col="orange",lty=2,lwd=4)
  lst<-10^(stpars$par[1]*(mm-stpars$par[2]))/area.st 
  lga<-10^(gapars$par[1]*(mm-gapars$par[2]))/area.ga
  lines(mm,log10(lst*allnorm*stagalnorm),col="turquoise",lty=2,lwd=4)														
  lines(mm,log10(lga*allnorm*galstanorm),col="orange",lty=2,lwd=4)

  #dev.off()























#############################
## ASSESSING THE MODEL FIT ##
#############################


## contours

  par(mfrow=c(2,2),pty="s",mai=c(0.5,0.3,0,0),mgp=c(1.5,0.6,0))
  px<-vector("list",8)
  names(px)<-c("y","j","jst","jga","h","hst","hga","k")
  pa<-MedPars$par; mmm<-max( max(lasdat$mapp_y[lasdat$miss.y.bin==0]), max(lasdat$mapp_j[lasdat$miss.j.bin==0]) + pa[6]-pa[7], max(lasdat$mapp_h[lasdat$miss.h.bin==0]) + pa[6]-pa[8], max(lasdat$mapp_k[lasdat$miss.k.bin==0]) + pa[6]-pa[9])
  for(band in 1:4){
    if(band==1){
      x<-lasdatrv[lasdatrv$miss.y.bin==0,c(18,6)]
      MagPars<-MagParsGl
      Mshift<-mmagshift$y
      starparas<-starpars$y$par
    }
    if(band==2){
      x<-lasdatrv[lasdatrv$miss.j.bin==0,c(19,8)]
      MagPars<-MagParsGl
      Mshift<-mmagshift$j
      starparas<-starpars$j$par
    }
    if(band==3){
      x<-lasdatrv[lasdatrv$miss.h.bin==0,c(20,10)]
      MagPars<-MagParsGl
      Mshift<-mmagshift$h
      starparas<-starpars$h$par
    }
    if(band==4){
      x<-lasdatrv[lasdatrv$miss.k.bin==0,c(21,12)]
      MagPars<-MagParsGl
      Mshift<-mmagshift$k
      starparas<-starpars$k$par
    }

    l<-250
    if(band==1){mm<-max(lasdat$mapp_y[lasdat$miss.y.bin==0])}
    if(band==2){mm<-max(lasdat$mapp_j[lasdat$miss.j.bin==0])}
    if(band==3){mm<-max(lasdat$mapp_h[lasdat$miss.h.bin==0])}
    if(band==4){mm<-max(lasdat$mapp_k[lasdat$miss.k.bin==0])}

    x1<-seq(-10,60,length=l)
    x2<-seq(10,mm,length=l)
    grid<-expand.grid(x1,x2)

    mu<-(1-(grid[,2]+Mshift)/mmm)*((pa[1]*(grid[,2]-(pa[6]-Mshift))^2 + pa[2]*(grid[,2]-(pa[6]-Mshift)) + pa[3])^pa[5] + pa[4])
    mu[mu<0.01]<-0.01
    rho<-VarPars$par[1]*(10^5)*((10)^(VarPars$par[2]*(grid[,2]+Mshift-11))) 

    galoc.fun<-function(x,cs,mu0,rho0){
      dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(cs-x,starpars=starparas)
    }
    galoc.intfun<-function(para){
      rel.tol<-.Machine$double.eps^0.25
      integrate(f=galoc.fun,lower=-10,upper=70,cs=para[1],mu0=para[3],rho0=para[4],stop.on.error=F,subdivisions=1000,rel.tol=rel.tol,abs.tol=rel.tol/1000)$value
    }
    galoc.val<-apply(X=cbind(grid[,1],grid[,2],mu,rho),MARGIN=1,FUN=galoc.intfun)

    if(band==1){px$y<-dmag_erf(grid[,2],mlim=(MagPars$par[1]-Mshift),sigma=MagPars$par[2],alpha=MagPars$par[3],minit=MagPars$par[4])*(MxCfGl(grid[,2]+mmagshift$y)*dstar(grid[,1],starpars=starparas) + (1-MxCfGl(grid[,2]+mmagshift$y))*galoc.val )}
    if(band==2){
      px$j<-dmag_erf(grid[,2],mlim=(MagPars$par[1]-Mshift),sigma=MagPars$par[2],alpha=MagPars$par[3],minit=MagPars$par[4])*(MxCfGl(grid[,2]+mmagshift$j)*dstar(grid[,1],starpars=starparas) + (1-MxCfGl(grid[,2]+mmagshift$j))*galoc.val )
      px$jst<-dmag_erf(grid[,2],mlim=(MagPars$par[1]-Mshift),sigma=MagPars$par[2],alpha=MagPars$par[3],minit=MagPars$par[4])*(MxCfGl(grid[,2]+mmagshift$j)*dstar(grid[,1],starpars=starparas) )
      px$jga<-dmag_erf(grid[,2],mlim=(MagPars$par[1]-Mshift),sigma=MagPars$par[2],alpha=MagPars$par[3],minit=MagPars$par[4])*( (1-MxCfGl(grid[,2]+mmagshift$j))*galoc.val )
    }
    if(band==3){
      px$h<-dmag_erf(grid[,2],mlim=(MagPars$par[1]-Mshift),sigma=MagPars$par[2],alpha=MagPars$par[3],minit=MagPars$par[4])*(MxCfGl(grid[,2]+mmagshift$h)*dstar(grid[,1],starpars=starparas) + (1-MxCfGl(grid[,2]+mmagshift$h))*galoc.val )}
      px$hst<-dmag_erf(grid[,2],mlim=(MagPars$par[1]-Mshift),sigma=MagPars$par[2],alpha=MagPars$par[3],minit=MagPars$par[4])*(MxCfGl(grid[,2]+mmagshift$h)*dstar(grid[,1],starpars=starparas) )
      px$hga<-dmag_erf(grid[,2],mlim=(MagPars$par[1]-Mshift),sigma=MagPars$par[2],alpha=MagPars$par[3],minit=MagPars$par[4])*( (1-MxCfGl(grid[,2]+mmagshift$h))*galoc.val )
    if(band==4){px$k<-dmag_erf(grid[,2],mlim=(MagPars$par[1]-Mshift),sigma=MagPars$par[2],alpha=MagPars$par[3],minit=MagPars$par[4])*(MxCfGl(grid[,2]+mmagshift$k)*dstar(grid[,1],starpars=starparas) + (1-MxCfGl(grid[,2]+mmagshift$k))*galoc.val )}

    if(band==1){
      plot(cbind(lasdat[lasdat$miss.y.bin==0,18],lasdat[lasdat$miss.y.bin==0,6]),col="yellow",xlim=c(-10,45),ylim=c(22,10),pch=".",xlab="z_Y",ylab="m_Y")
      contour(x1,x2,matrix(px$y,ncol=l),levels=c(0.0005,0.001,0.005,0.0075,0.01,0.02,0.04,0.08,0.5),add=T)
    }
    if(band==2){
      plot(cbind(lasdat[lasdat$miss.j.bin==0,19],lasdat[lasdat$miss.j.bin==0,8]),col="yellow",xlim=c(-10,45),ylim=c(22,10),pch=".",xlab="z_J",ylab="m_J")
      contour(x1,x2,matrix(px$j,ncol=l),levels=c(0.0005,0.001,0.005,0.0075,0.01,0.02,0.04,0.08,0.5),add=T)
      contour(x1,x2,matrix(px$jst,ncol=l),levels=c(0.0005,0.001,0.005,0.0075,0.01,0.02,0.04,0.08,0.5),add=T,col="blue")
      contour(x1,x2,matrix(px$jga,ncol=l),levels=c(0.0005,0.001,0.005,0.0075,0.01,0.02,0.04,0.08,0.5),add=T,col="red")
    }
    if(band==3){
      plot(cbind(lasdat[lasdat$miss.h.bin==0,20],lasdat[lasdat$miss.h.bin==0,10]),col="yellow",xlim=c(-10,45),ylim=c(22,10),pch=".",xlab="z_H",ylab="m_H")
      contour(x1,x2,matrix(px$h,ncol=l),levels=c(0.0005,0.001,0.005,0.0075,0.01,0.02,0.04,0.08,0.5),add=T)
      contour(x1,x2,matrix(px$hst,ncol=l),levels=c(0.0005,0.001,0.005,0.0075,0.01,0.02,0.04,0.08,0.5),add=T,col="blue")
      contour(x1,x2,matrix(px$hga,ncol=l),levels=c(0.0005,0.001,0.005,0.0075,0.01,0.02,0.04,0.08,0.5),add=T,col="red")
    }
    if(band==4){
      plot(cbind(lasdat[lasdat$miss.k.bin==0,21],lasdat[lasdat$miss.k.bin==0,12]),col="yellow",xlim=c(-10,45),ylim=c(22,10),pch=".",xlab="z_K",ylab="m_K") 
      contour(x1,x2,matrix(px$k,ncol=l),levels=c(0.0005,0.001,0.005,0.0075,0.01,0.02,0.04,0.08,0.5),add=T)
    }

  }



## one-dim distributions / histograms



  ## check
    mcheck<-18.5; band<-1
    if(band==1){
      dtemp<-lasdat[lasdat$miss.y.bin==0 & lasdat$mapp_y<mcheck+0.1 & lasdat$mapp_y>mcheck-0.1,]
      mcheck<-mean(dtemp$mapp_y)
      MagPars<-MagParsGl
      Mshift<-mmagshift$y
      mm<-max(lasdat$mapp_y[lasdat$miss.y.bin==0])
      starparas<-starpars$y$par
    }
    if(band==2){
      dtemp<-lasdat[lasdat$miss.j.bin==0 & lasdat$mapp_j<mcheck+0.1 & lasdat$mapp_j>mcheck-0.1,]
      mcheck<-mean(dtemp$mapp_j)
      MagPars<-MagParsGl
      Mshift<-mmagshift$j
      mm<-max(lasdat$mapp_j[lasdat$miss.j.bin==0])
      starparas<-starpars$j$par
    }
    if(band==3){
      dtemp<-lasdat[lasdat$miss.h.bin==0 & lasdat$mapp_h<mcheck+0.1 & lasdat$mapp_h>mcheck-0.1,]
      mcheck<-mean(dtemp$mapp_h)
      MagPars<-MagParsGl
      Mshift<-mmagshift$h
      mm<-max(lasdat$mapp_h[lasdat$miss.h.bin==0])
      starparas<-starpars$h$par
    }
    if(band==4){
      dtemp<-lasdat[lasdat$miss.k.bin==0 & lasdat$mapp_k<mcheck+0.1 & lasdat$mapp_k>mcheck-0.1,]
      mcheck<-mean(dtemp$mapp_k)
      MagPars<-MagParsGl
      Mshift<-mmagshift$k
      mm<-max(lasdat$mapp_k[lasdat$miss.k.bin==0])
      starparas<-starpars$k$par
    }
    cs<-seq(-20,50,length=250)
    pa<-MedPars$par; mmm<-max( max(lasdat$mapp_y[lasdat$miss.y.bin==0]), max(lasdat$mapp_j[lasdat$miss.j.bin==0]) + pa[6]-pa[7], max(lasdat$mapp_h[lasdat$miss.h.bin==0]) + pa[6]-pa[8], max(lasdat$mapp_k[lasdat$miss.k.bin==0]) + pa[6]-pa[9])
    mu<-(1-(mcheck+Mshift)/mmm)*((pa[1]*(mcheck-(pa[6]-Mshift))^2 + pa[2]*(mcheck-(pa[6]-Mshift)) + pa[3])^pa[5] + pa[4])
    mu[mu<0.01]<-0.01
    rho<-VarPars$par[1]*(10^5)*((10)^(VarPars$par[2]*(mcheck+Mshift-11))) 

    galoc.fun<-function(x,cs,mu0,rho0){
      dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(cs-x,starpars=starparas)
    }
    galoc.intfun<-function(para){
      rel.tol<-.Machine$double.eps^0.25
      integrate(f=galoc.fun,lower=-5,upper=60,cs=para[1],mu0=para[3],rho0=para[4],stop.on.error=F,subdivisions=1000,rel.tol=rel.tol,abs.tol=rel.tol/1000)$value
    }
    galoc.val<-apply(X=cbind(cs,-mcheck,mu,rho),MARGIN=1,FUN=galoc.intfun)

    pdfcs<-(MxCfGl(mcheck+Mshift))*dstar(cs,starpars=starparas) + (1-MxCfGl(mcheck+Mshift))*galoc.val
    pdfcsN<-(MxCfGl(mcheck+Mshift))*dstar(cs,starpars=starparas) 
    pdfcsLN<-(1-MxCfGl(mcheck+Mshift))*galoc.val
    if(band==1){hist(dtemp$cs_y,freq=F,breaks=100,xlab="z_Y",main=c("Y band; mean magnitude of sample:",round(mcheck*10)/10),xlim=c(-10,40),ylim=c(0,0.45))}
    if(band==2){hist(dtemp$cs_j,freq=F,breaks=100,xlab="cs_J",main=c("J band; mean magnitude of sample:",round(mcheck*10)/10))}
    if(band==3){hist(dtemp$cs_h,freq=F,breaks=100,xlab="cs_H",main=c("H band; mean magnitude of sample:",round(mcheck*10)/10))}
    if(band==4){hist(dtemp$cs_k,freq=F,breaks=100,xlab="cs_K",main=c("K band;mean magnitude of sample:",round(mcheck*10)/10))}
    lines(cs,pdfcs,col="darkgreen")
    lines(cs,pdfcsN,col="blue")
    lines(cs,pdfcsLN,col="red")



##########################################################
## SINGLE BAND CLASSIFICATIONS WITH TRUE LOCATION MODEL ##
##########################################################

pStar<-vector("list",4)
names(pStar)<-c("y","j","h,","k")
Ncols<-250
pr<-seq(0,1,length=(Ncols+1))
rnk<-vector("list",4)
names(rnk)<-c("y","j","h,","k")

par(mfrow=c(2,2),pty="s",mai=c(0.5,0.3,0,0),mgp=c(1.5,0.6,0))
for(band in 1:4){
  if(band==1){
    starparas<-starpars$y$par
    pa<-MedPars$par; mmm<-max( max(lasdat$mapp_y[lasdat$miss.y.bin==0]), max(lasdat$mapp_j[lasdat$miss.j.bin==0]) + pa[6]-pa[7], max(lasdat$mapp_h[lasdat$miss.h.bin==0]) + pa[6]-pa[8], max(lasdat$mapp_k[lasdat$miss.k.bin==0]) + pa[6]-pa[9])
    MagPars<-MagParsGl
    Mshift<-mmagshift$y
    mm<-max(lasdat$mapp_y[lasdat$miss.y.bin==0])
    mu<-(1-(lasdat$mapp_y+Mshift)/mmm)*((pa[1]*(lasdat$mapp_y-(pa[6]-Mshift))^2 + pa[2]*(lasdat$mapp_y-(pa[6]-Mshift)) + pa[3])^pa[5] + pa[4])
    mu[mu<0.01]<-0.01
    rho<-VarPars$par[1]*(10^5)*((10)^(VarPars$par[2]*(lasdat$mapp_y+Mshift-11)))

    galoc.fun<-function(x,cs,mu0,rho0){
      dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(cs-x,starpars=starparas)
    }
    galoc.intfun<-function(para){
      rel.tol<-.Machine$double.eps^0.25
      integrate(f=galoc.fun,lower=-10,upper=70,cs=para[1],mu0=para[3],rho0=para[4],stop.on.error=F,subdivisions=1000,rel.tol=rel.tol,abs.tol=rel.tol/1000)$value
    }
    galoc.val<-apply(X=cbind(lasdat$cs_y,lasdat$mapp_y,mu,rho),MARGIN=1,FUN=galoc.intfun)

    pS<-(MxCfGl(lasdat$mapp_y+mmagshift$y))*dstar(lasdat$cs_y,starpars=starparas)
    pG<-(1-MxCfGl(lasdat$mapp_y+mmagshift$y))*galoc.val
    pStar$y<-pS/(pS+pG)
    pStar$y[lasdat$miss.y.bin==1]<-rep(NA,sum(lasdat$miss.y.bin==1))

    ## plot

      rnk$y<-rep(NA,nrow(lasdat))
      for(i in 1:Ncols){
        rnk$y[lasdat$miss.y.bin==0 & pStar$y>=pr[i] & pStar$y<pr[i+1]]<-i
      }
      rnk$y[lasdat$miss.y.bin==0 & pStar$y==1]<-Ncols

      dtemp<-lasdat[lasdat$miss.y.bin==0,]
      t<-rainbow(Ncols,start=1,end=0.65)
      plot(dtemp$cs_y,dtemp$mapp_y,col=t[rnk$y[lasdat$miss.y.bin==0]],pch=".",xlim=c(-25,65),ylim=c(22.1,11.25),xlab="cs_Y",ylab="m_Y")
  }
  if(band==2){
    starparas<-starpars$j$par
    pa<-MedPars$par; mmm<-max( max(lasdat$mapp_y[lasdat$miss.y.bin==0]), max(lasdat$mapp_j[lasdat$miss.j.bin==0]) + pa[6]-pa[7], max(lasdat$mapp_h[lasdat$miss.h.bin==0]) + pa[6]-pa[8], max(lasdat$mapp_k[lasdat$miss.k.bin==0]) + pa[6]-pa[9])
    MagPars<-MagParsGl
    Mshift<-mmagshift$j
    mm<-max(lasdat$mapp_j[lasdat$miss.j.bin==0])
    mu<-(1-(lasdat$mapp_j+Mshift)/mmm)*((pa[1]*(lasdat$mapp_j-(pa[6]-Mshift))^2 + pa[2]*(lasdat$mapp_j-(pa[6]-Mshift)) + pa[3])^pa[5] + pa[4])
    mu[mu<0.01]<-0.01
    rho<-VarPars$par[1]*(10^5)*((10)^(VarPars$par[2]*(lasdat$mapp_j+Mshift-11))) 

    galoc.fun<-function(x,cs,mu0,rho0){
      dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(cs-x,starpars=starparas)
    }
    galoc.intfun<-function(para){
      rel.tol<-.Machine$double.eps^0.25
      integrate(f=galoc.fun,lower=-10,upper=70,cs=para[1],mu0=para[3],rho0=para[4],stop.on.error=F,subdivisions=1000,rel.tol=rel.tol,abs.tol=rel.tol/1000)$value
    }
    galoc.val<-apply(X=cbind(lasdat$cs_j,lasdat$mapp_j,mu,rho),MARGIN=1,FUN=galoc.intfun)

    pS<-(MxCfGl(lasdat$mapp_j+mmagshift$j))*dstar(lasdat$cs_j,starpars=starparas)
    pG<-(1-MxCfGl(lasdat$mapp_j+mmagshift$j))*galoc.val
    pStar$j<-pS/(pS+pG)
    pStar$j[lasdat$miss.j.bin==1]<-rep(NA,sum(lasdat$miss.j.bin==1))

    ## plot

      rnk$j<-rep(NA,nrow(lasdat))
      for(i in 1:Ncols){
        rnk$j[lasdat$miss.j.bin==0 & pStar$j>=pr[i] & pStar$j<pr[i+1]]<-i
      }
      rnk$j[lasdat$miss.j.bin==0 & pStar$j==1]<-Ncols

      dtemp<-lasdat[lasdat$miss.j.bin==0,]
      t<-rainbow(Ncols,start=1,end=0.65)
      plot(dtemp$cs_j,dtemp$mapp_j,col=t[rnk$j[lasdat$miss.j.bin==0]],pch=".",xlim=c(-25,65),ylim=c(22.1,11.25),xlab="cs_J",ylab="m_J")
  }
  if(band==3){
    starparas<-starpars$h$par
    pa<-MedPars$par; mmm<-max( max(lasdat$mapp_y[lasdat$miss.y.bin==0]), max(lasdat$mapp_j[lasdat$miss.j.bin==0]) + pa[6]-pa[7], max(lasdat$mapp_h[lasdat$miss.h.bin==0]) + pa[6]-pa[8], max(lasdat$mapp_k[lasdat$miss.k.bin==0]) + pa[6]-pa[9])
    MagPars<-MagParsGl
    Mshift<-mmagshift$h
    mm<-max(lasdat$mapp_h[lasdat$miss.h.bin==0])
    mu<-(1-(lasdat$mapp_h+Mshift)/mmm)*((pa[1]*(lasdat$mapp_h-(pa[6]-Mshift))^2 + pa[2]*(lasdat$mapp_h-(pa[6]-Mshift)) + pa[3])^pa[5] + pa[4])
    mu[mu<0.01]<-0.01
    rho<-VarPars$par[1]*(10^5)*((10)^(VarPars$par[2]*(lasdat$mapp_h+Mshift-11))) 

    galoc.fun<-function(x,cs,mu0,rho0){
      dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(cs-x,starpars=starparas)
    }
    galoc.intfun<-function(para){
      rel.tol<-.Machine$double.eps^0.25
      integrate(f=galoc.fun,lower=-10,upper=70,cs=para[1],mu0=para[3],rho0=para[4],stop.on.error=F,subdivisions=1000,rel.tol=rel.tol,abs.tol=rel.tol/1000)$value
    }
    galoc.val<-apply(X=cbind(lasdat$cs_h,lasdat$mapp_h,mu,rho),MARGIN=1,FUN=galoc.intfun)

    pS<-(MxCfGl(lasdat$mapp_h+mmagshift$h))*dstar(lasdat$cs_h,starpars=starparas)
    pG<-(1-MxCfGl(lasdat$mapp_h+mmagshift$h))*galoc.val
    pStar$h<-pS/(pS+pG)
    pStar$h[lasdat$miss.h.bin==1]<-rep(NA,sum(lasdat$miss.h.bin==1))

    ## plot

      rnk$h<-rep(NA,nrow(lasdat))
      for(i in 1:Ncols){
        rnk$h[lasdat$miss.h.bin==0 & pStar$h>=pr[i] & pStar$h<pr[i+1]]<-i
      }
      rnk$h[lasdat$miss.h.bin==0 & pStar$h==1]<-Ncols

      dtemp<-lasdat[lasdat$miss.h.bin==0,]
      t<-rainbow(Ncols,start=1,end=0.65)
      plot(dtemp$cs_h,dtemp$mapp_h,col=t[rnk$h[lasdat$miss.h.bin==0]],pch=".",xlim=c(-25,65),ylim=c(20.1,11.1),xlab="cs_H",ylab="m_H")
  }
  if(band==4){
    starparas<-starpars$k$par
    pa<-MedPars$par; mmm<-max( max(lasdat$mapp_y[lasdat$miss.y.bin==0]), max(lasdat$mapp_j[lasdat$miss.j.bin==0]) + pa[6]-pa[7], max(lasdat$mapp_h[lasdat$miss.h.bin==0]) + pa[6]-pa[8], max(lasdat$mapp_k[lasdat$miss.k.bin==0]) + pa[6]-pa[9])
    MagPars<-MagParsGl
    Mshift<-mmagshift$k
    mm<-max(lasdat$mapp_k[lasdat$miss.k.bin==0])
    mu<-(1-(lasdat$mapp_k+Mshift)/mmm)*((pa[1]*(lasdat$mapp_k-(pa[6]-Mshift))^2 + pa[2]*(lasdat$mapp_k-(pa[6]-Mshift)) + pa[3])^pa[5] + pa[4])
    mu[mu<0.01]<-0.01
    rho<-VarPars$par[1]*(10^5)*((10)^(VarPars$par[2]*(lasdat$mapp_k+Mshift-11))) 

    galoc.fun<-function(x,cs,mu0,rho0){
      dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(cs-x,starpars=starparas)
    }
    galoc.intfun<-function(para){
      rel.tol<-.Machine$double.eps^0.25
      integrate(f=galoc.fun,lower=-10,upper=70,cs=para[1],mu0=para[3],rho0=para[4],stop.on.error=F,subdivisions=1000,rel.tol=rel.tol,abs.tol=rel.tol/1000)$value
    }
    galoc.val<-apply(X=cbind(lasdat$cs_k,lasdat$mapp_k,mu,rho),MARGIN=1,FUN=galoc.intfun)

    pS<-(MxCfGl(lasdat$mapp_k+mmagshift$k))*dstar(lasdat$cs_k,starpars=starparas)
    pG<-(1-MxCfGl(lasdat$mapp_k+mmagshift$k))*galoc.val
    pStar$k<-pS/(pS+pG)
    pStar$k[lasdat$miss.k.bin==1]<-rep(NA,sum(lasdat$miss.k.bin==1))

    ## plot

      rnk$k<-rep(NA,nrow(lasdat))
      for(i in 1:Ncols){
        rnk$k[lasdat$miss.k.bin==0 & pStar$k>=pr[i] & pStar$k<pr[i+1]]<-i
      }                                                                                                                                                                                                                                                                                                                   
      rnk$k[lasdat$miss.k.bin==0 & pStar$k==1]<-Ncols

      dtemp<-lasdat[lasdat$miss.k.bin==0,]
      t<-rainbow(Ncols,start=1,end=0.65)
      plot(dtemp$cs_k,dtemp$mapp_k,col=t[rnk$k[lasdat$miss.k.bin==0]],pch=".",xlim=c(-25,65),ylim=c(19.8,10.75),xlab="cs_K",ylab="m_K")
  }
}

pStarSing<-pStar

sum(pStarSing$y[!is.na(pStarSing$y)]>0.4 & pStarSing$y[!is.na(pStarSing$y)]<0.6)/sum(!is.na(pStarSing$y))
sum(pStarSing$j[!is.na(pStarSing$j)]>0.4 & pStarSing$j[!is.na(pStarSing$j)]<0.6)/sum(!is.na(pStarSing$j))
sum(pStarSing$h[!is.na(pStarSing$h)]>0.4 & pStarSing$h[!is.na(pStarSing$h)]<0.6)/sum(!is.na(pStarSing$h))
sum(pStarSing$k[!is.na(pStarSing$k)]>0.4 & pStarSing$k[!is.na(pStarSing$k)]<0.6)/sum(!is.na(pStarSing$k))



###########################
## GLOBAL CLASSIFICATION ##
###########################

## computing mMean
  mY<-lasdat$mapp_y+mmagshift$y
  mJ<-lasdat$mapp_j+mmagshift$j
  mH<-lasdat$mapp_h+mmagshift$h
  mK<-lasdat$mapp_k+mmagshift$k
  mDat<-cbind(mY,mJ,mH,mK)
  midxDat<-cbind(lasdat$miss.y.bin,lasdat$miss.j.bin,lasdat$miss.h.bin,lasdat$miss.k.bin)
  mMean<-apply(X=cbind(mDat,midxDat),MARGIN=1,FUN=SelectiveMean)

## probability densities

  pa<-MedPars$par
  mmm<-max( max(lasdat$mapp_y[lasdat$miss.y.bin==0]), max(lasdat$mapp_j[lasdat$miss.j.bin==0]) + pa[6]-pa[7], max(lasdat$mapp_h[lasdat$miss.h.bin==0]) + pa[6]-pa[8], max(lasdat$mapp_k[lasdat$miss.k.bin==0]) + pa[6]-pa[9])
  mu<-(1-(mMean)/mmm)*((pa[1]*(mMean-(pa[6]))^2 + pa[2]*(mMean-(pa[6])) + pa[3])^pa[5] + pa[4])
  mu[mu<0.01]<-0.01
  rho<-VarPars$par[1]*(10^5)*((10)^(VarPars$par[2]*(mMean-11))) 

  DensComp<-function(dat){
    ## dat[1:n] = data vector (classstat)
    ## dat[n+1:2*n] = binary vector of same length than x(=dat[1:n]), indicating whether a coordinate is (1) or is not (0) missing, i.e. not taken (1) or taken (0) into account 
    ## dat[(2*n+1)] = value of mu
    ## dat[(2*n+2)] = value of rho

    galoc.fun<-function(x,csy,csj,csh,csk,idx.y,idx.j,idx.h,idx.k,mu0,rho0){
      if(idx.y==0 & idx.j==0 & idx.h==0 & idx.k==0){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csy,starpars=starpars$y$par)*dstar(x-csj,starpars=starpars$j$par)*dstar(x-csh,starpars=starpars$h$par)*dstar(x-csk,starpars=starpars$k$par)}
      if(idx.y==1 & idx.j==0 & idx.h==0 & idx.k==0){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csj,starpars=starpars$j$par)*dstar(x-csh,starpars=starpars$h$par)*dstar(x-csk,starpars=starpars$k$par)}
      if(idx.y==0 & idx.j==1 & idx.h==0 & idx.k==0){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csy,starpars=starpars$y$par)*dstar(x-csh,starpars=starpars$h$par)*dstar(x-csk,starpars=starpars$k$par)}
      if(idx.y==0 & idx.j==0 & idx.h==1 & idx.k==0){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csy,starpars=starpars$y$par)*dstar(x-csj,starpars=starpars$j$par)*dstar(x-csk,starpars=starpars$k$par)}
      if(idx.y==0 & idx.j==0 & idx.h==0 & idx.k==1){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csy,starpars=starpars$y$par)*dstar(x-csj,starpars=starpars$j$par)*dstar(x-csh,starpars=starpars$h$par)}
      if(idx.y==1 & idx.j==1 & idx.h==0 & idx.k==0){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csh,starpars=starpars$h$par)*dstar(x-csk,starpars=starpars$k$par)}
      if(idx.y==1 & idx.j==0 & idx.h==1 & idx.k==0){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csj,starpars=starpars$j$par)*dstar(x-csk,starpars=starpars$k$par)}
      if(idx.y==1 & idx.j==0 & idx.h==0 & idx.k==1){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csj,starpars=starpars$j$par)*dstar(x-csh,starpars=starpars$h$par)}
      if(idx.y==0 & idx.j==1 & idx.h==1 & idx.k==0){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csy,starpars=starpars$y$par)*dstar(x-csk,starpars=starpars$k$par)}
      if(idx.y==0 & idx.j==1 & idx.h==0 & idx.k==1){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csy,starpars=starpars$y$par)*dstar(x-csh,starpars=starpars$h$par)}
      if(idx.y==0 & idx.j==0 & idx.h==1 & idx.k==1){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csy,starpars=starpars$y$par)*dstar(x-csj,starpars=starpars$j$par)}
      if(idx.y==1 & idx.j==1 & idx.h==1 & idx.k==0){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csk,starpars=starpars$k$par)}
      if(idx.y==1 & idx.j==1 & idx.h==0 & idx.k==1){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csh,starpars=starpars$h$par)}
      if(idx.y==1 & idx.j==0 & idx.h==1 & idx.k==1){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csj,starpars=starpars$j$par)}
      if(idx.y==0 & idx.j==1 & idx.h==1 & idx.k==1){res<-dlnorm(x,meanlog=log(mu0)-( (1/2)*log( 1 + rho0/(mu0^2) ) ),sdlog=sqrt( log(1 + rho0/(mu0^2)) ) )*dstar(x-csy,starpars=starpars$y$par)}
      res
    }

    n<-(length(dat)-2)/2
    x<-dat[1:n]
    miss.idx<-dat[(n+1):(2*n)]
    mu0<-dat[(2*n+1)]
    rho0<-dat[(2*n+2)]
    if(n==0 | (n-round(n))!=0 | sum(miss.idx[miss.idx!=0]!=1)>0){stop("Input data in bad format!")}


    staloc.fun<-function(x,idx.y,idx.j,idx.h,idx.k){
      if(idx.y==0 & idx.j==0 & idx.h==0 & idx.k==0){res<-dstar(x[1],starpars=starpars$y$par)*dstar(x[2],starpars=starpars$j$par)*dstar(x[3],starpars=starpars$h$par)*dstar(x[4],starpars=starpars$k$par)}
      if(idx.y==1 & idx.j==0 & idx.h==0 & idx.k==0){res<-dstar(x[2],starpars=starpars$j$par)*dstar(x[3],starpars=starpars$h$par)*dstar(x[4],starpars=starpars$k$par)}
      if(idx.y==0 & idx.j==1 & idx.h==0 & idx.k==0){res<-dstar(x[1],starpars=starpars$y$par)*dstar(x[3],starpars=starpars$h$par)*dstar(x[4],starpars=starpars$k$par)}
      if(idx.y==0 & idx.j==0 & idx.h==1 & idx.k==0){res<-dstar(x[1],starpars=starpars$y$par)*dstar(x[2],starpars=starpars$j$par)*dstar(x[4],starpars=starpars$k$par)}
      if(idx.y==0 & idx.j==0 & idx.h==0 & idx.k==1){res<-dstar(x[1],starpars=starpars$y$par)*dstar(x[2],starpars=starpars$j$par)*dstar(x[3],starpars=starpars$h$par)}
      if(idx.y==1 & idx.j==1 & idx.h==0 & idx.k==0){res<-dstar(x[3],starpars=starpars$h$par)*dstar(x[4],starpars=starpars$k$par)}
      if(idx.y==1 & idx.j==0 & idx.h==1 & idx.k==0){res<-dstar(x[2],starpars=starpars$j$par)*dstar(x[4],starpars=starpars$k$par)}
      if(idx.y==1 & idx.j==0 & idx.h==0 & idx.k==1){res<-dstar(x[2],starpars=starpars$j$par)*dstar(x[3],starpars=starpars$h$par)}
      if(idx.y==0 & idx.j==1 & idx.h==1 & idx.k==0){res<-dstar(x[1],starpars=starpars$y$par)*dstar(x[4],starpars=starpars$k$par)}
      if(idx.y==0 & idx.j==1 & idx.h==0 & idx.k==1){res<-dstar(x[1],starpars=starpars$y$par)*dstar(x[3],starpars=starpars$h$par)}
      if(idx.y==0 & idx.j==0 & idx.h==1 & idx.k==1){res<-dstar(x[1],starpars=starpars$y$par)*dstar(x[2],starpars=starpars$j$par)}
      if(idx.y==1 & idx.j==1 & idx.h==1 & idx.k==0){res<-dstar(x[4],starpars=starpars$k$par)}
      if(idx.y==1 & idx.j==1 & idx.h==0 & idx.k==1){res<-dstar(x[3],starpars=starpars$h$par)}
      if(idx.y==1 & idx.j==0 & idx.h==1 & idx.k==1){res<-dstar(x[2],starpars=starpars$j$par)}
      if(idx.y==0 & idx.j==1 & idx.h==1 & idx.k==1){res<-dstar(x[1],starpars=starpars$y$par)}
      res
    }
    pStar<-staloc.fun(x,idx.y=miss.idx[1],idx.j=miss.idx[2],idx.h=miss.idx[3],idx.k=miss.idx[4])
    #pStar<-prod(dnorm_nozero(x[miss.idx==0]))

    rel.tol<-.Machine$double.eps^0.25
    pGal<-integrate(f=galoc.fun,lower=-10,upper=70,mu0=mu0,rho0=rho0,csy=x[1],csj=x[2],csh=x[3],csk=x[4],idx.y=miss.idx[1],idx.j=miss.idx[2],idx.h=miss.idx[3],idx.k=miss.idx[4],stop.on.error=F,subdivisions=1000,rel.tol=rel.tol,abs.tol=rel.tol/1000)$value

    res<-c(pStar,pGal)
    res    
  }

  X<-cbind(lasdat$cs_y,lasdat$cs_j,lasdat$cs_h,lasdat$cs_k,lasdat$miss.y.bin,lasdat$miss.j.bin,lasdat$miss.h.bin,lasdat$miss.k.bin,mu,rho)
  dens<-apply(X=X,MARGIN=1,FUN=DensComp)
  dens<-t(dens)
  colnames(dens)<-c("pStar","pGal")
  idx.miss<-(1:nrow(lasdat))[lasdat$miss.y.bin==1 & lasdat$miss.j.bin==1 &lasdat$miss.h.bin==1 &lasdat$miss.k.bin==1]
  idx.miss ## to check that all objects have been classified
  dens[,1]<-MxCfGl(mMean)*dens[,1]
  dens[,2]<-(1-MxCfGl(mMean))*dens[,2]
  pStar<-dens[,1]/(dens[,1]+dens[,2])
  if(sum(dens[,1]==0 & dens[,2]==0)>0){pStar[dens[,1]==0 & dens[,2]==0]<-rep(0.5,length=sum(dens[,1]==0 & dens[,2]==0))}
  par(mfrow=c(1,1),pty="s")
  hist(pStar,breaks=100,freq=F)

  Ncols<-250
  pr<-seq(0,1,length=(Ncols+1))
  rnkgl<-rep(NA,nrow(lasdat))
  for(i in 1:Ncols){
    rnkgl[pStar>=pr[i] & pStar<pr[i+1]]<-i
  }                                                                                                                                                                                                                                                                                                                   
  rnkgl[pStar==1]<-Ncols

  par(mfrow=c(2,2),mai=c(0.5,0.3,0,0),mgp=c(1.5,0.6,0),pty="s")
  ## plotting the global classification in the Y band
    dtemp<-lasdat[lasdat$miss.y.bin==0,]
    t<-rainbow(Ncols,start=1,end=0.65)
    plot(dtemp$cs_y,dtemp$mapp_y,col=t[rnkgl[lasdat$miss.y.bin==0]],pch=".",xlim=c(-25,65),ylim=c(22.1,11.25),xlab="z_Y",ylab="m_Y")

  ## plotting the global classification in the J band
    dtemp<-lasdat[lasdat$miss.j.bin==0,]
    t<-rainbow(Ncols,start=1,end=0.65)
    plot(dtemp$cs_j,dtemp$mapp_j,col=t[rnkgl[lasdat$miss.j.bin==0]],pch=".",xlim=c(-25,65),ylim=c(22.1,11.25),xlab="z_J",ylab="m_J")

  ## plotting the global classification in the H band
    dtemp<-lasdat[lasdat$miss.h.bin==0,]
    t<-rainbow(Ncols,start=1,end=0.65)
    plot(dtemp$cs_h,dtemp$mapp_h,col=t[rnkgl[lasdat$miss.h.bin==0]],pch=".",xlim=c(-25,65),ylim=c(20.1,11.1),xlab="z_H",ylab="m_H")

  ## plotting the global classification in the K band
    dtemp<-lasdat[lasdat$miss.k.bin==0,]
    t<-rainbow(Ncols,start=1,end=0.65)
    plot(dtemp$cs_k,dtemp$mapp_k,col=t[rnkgl[lasdat$miss.k.bin==0]],pch=".",xlim=c(-25,65),ylim=c(19.8,10.75),xlab="z_K",ylab="m_K")

  library(fields)
  par(mfrow=c(2,2))
    quilt.plot.mod(lasdat$cs_y[lasdat$miss.y.bin==0],lasdat$mapp_y[lasdat$miss.y.bin==0],pStar[lasdat$miss.y.bin==0],col=rainbow(64,start=1,end=0.65),nrow=1024,ncol=1024,xlim=c(-12,50),ylim=c(21.3,11),main="",xlab=expression(cs[Y]),ylab="Y",family="serif",cex.axis=2,cex.lab=2.25,axis.args=list(cex.axis=2,family="serif"))
    quilt.plot.mod(lasdat$cs_j[lasdat$miss.j.bin==0],lasdat$mapp_j[lasdat$miss.j.bin==0],pStar[lasdat$miss.j.bin==0],col=rainbow(64,start=1,end=0.65),nrow=1024,ncol=1024,xlim=c(-12,50),ylim=c(21.3,11),main="",xlab=expression(cs[J]),ylab="J",family="serif",cex.axis=2,cex.lab=2.25,axis.args=list(cex.axis=2,family="serif"))
    quilt.plot.mod(lasdat$cs_h[lasdat$miss.h.bin==0],lasdat$mapp_h[lasdat$miss.h.bin==0],pStar[lasdat$miss.h.bin==0],col=rainbow(64,start=1,end=0.65),nrow=1024,ncol=1024,xlim=c(-12,50),ylim=c(21.3,11),main="",xlab=expression(cs[H]),ylab="H",family="serif",cex.axis=2,cex.lab=2.25,axis.args=list(cex.axis=2,family="serif"))
    quilt.plot.mod(lasdat$cs_k[lasdat$miss.k.bin==0],lasdat$mapp_k[lasdat$miss.k.bin==0],pStar[lasdat$miss.k.bin==0],col=rainbow(64,start=1,end=0.65),nrow=1024,ncol=1024,xlim=c(-12,50),ylim=c(21.3,11),main="",xlab=expression(cs[K]),ylab="K",family="serif",cex.axis=2,cex.lab=2.25,axis.args=list(cex.axis=2,family="serif"))
  
   # sums
    sum(pStar[lasdat$miss.y.bin==0]>0.4 & pStar[lasdat$miss.y.bin==0]<0.6); sum(pStarSing$y[lasdat$miss.y.bin==0]>0.4 & pStarSing$y[lasdat$miss.y.bin==0]<0.6)
    sum(pStar[lasdat$miss.j.bin==0]>0.4 & pStar[lasdat$miss.j.bin==0]<0.6); sum(pStarSing$j[lasdat$miss.j.bin==0]>0.4 & pStarSing$j[lasdat$miss.j.bin==0]<0.6)
    sum(pStar[lasdat$miss.h.bin==0]>0.4 & pStar[lasdat$miss.h.bin==0]<0.6); sum(pStarSing$h[lasdat$miss.h.bin==0]>0.4 & pStarSing$h[lasdat$miss.h.bin==0]<0.6)
    sum(pStar[lasdat$miss.k.bin==0]>0.4 & pStar[lasdat$miss.k.bin==0]<0.6); sum(pStarSing$k[lasdat$miss.k.bin==0]>0.4 & pStarSing$k[lasdat$miss.k.bin==0]<0.6)
   # fractions
    sum(pStar[lasdat$miss.y.bin==0]>0.4 & pStar[lasdat$miss.y.bin==0]<0.6)/sum(lasdat$miss.y.bin==0); sum(pStarSing$y[lasdat$miss.y.bin==0]>0.4 & pStarSing$y[lasdat$miss.y.bin==0]<0.6)/sum(lasdat$miss.y.bin==0)
    sum(pStar[lasdat$miss.j.bin==0]>0.4 & pStar[lasdat$miss.j.bin==0]<0.6)/sum(lasdat$miss.j.bin==0); sum(pStarSing$j[lasdat$miss.j.bin==0]>0.4 & pStarSing$j[lasdat$miss.j.bin==0]<0.6)/sum(lasdat$miss.j.bin==0)
    sum(pStar[lasdat$miss.h.bin==0]>0.4 & pStar[lasdat$miss.h.bin==0]<0.6)/sum(lasdat$miss.h.bin==0); sum(pStarSing$h[lasdat$miss.h.bin==0]>0.4 & pStarSing$h[lasdat$miss.h.bin==0]<0.6)/sum(lasdat$miss.h.bin==0)
    sum(pStar[lasdat$miss.k.bin==0]>0.4 & pStar[lasdat$miss.k.bin==0]<0.6)/sum(lasdat$miss.k.bin==0); sum(pStarSing$k[lasdat$miss.k.bin==0]>0.4 & pStarSing$k[lasdat$miss.k.bin==0]<0.6)/sum(lasdat$miss.k.bin==0)
   














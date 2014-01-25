erfc<-function(x){
  ## this is the complementary error function, modelling observation decline near observation limit
  res<-2*pnorm(-sqrt(2)*x)
  res
}
dN<-function(x,alpha){
  ## this is modelling the number count of objects N(x)=10^(x-m0); ignores the existence of an observation limit
  ## N0 and m0 redundant as adds only multiplicative term, neutralised by normalisation in dmag
  res<-(10)^(alpha*x)
  res
}

dmag_erf<-function(x,minit,mlim,sigma,alpha){
  ## combining dN and erfc gives a model of the observed number counts
  temp<-function(x){
    n<-length(x)
    if(n==1){
      if(dN(x,alpha)!=Inf){res<-0.5*erfc((x-mlim)/sigma) * dN(x,alpha) * (1-0.5*erfc((x-minit)/sigma))}
      if(dN(x,alpha)==Inf){res<-0}
    }
    if(n>1){
      res<-rep(NA,length=n)
      res[dN(x,alpha)<Inf]<-0.5*erfc((x[dN(x,alpha)<Inf]-mlim)/sigma) * dN(x[dN(x,alpha)<Inf],alpha) * (1-0.5*erfc((x[dN(x,alpha)<Inf]-minit)/sigma))
      res[dN(x,alpha)==Inf]<-0 
    }
    res
  }
  rel.tol<-.Machine$double.eps^0.25
  norm<-integrate(temp,lower=-100,upper=100,subdivisions=1000,stop.on.error=F,rel.tol=rel.tol,abs.tol=rel.tol/100)$value
  res<-temp(x)/norm
  res
}

FitMagInit<-function(x,para){
  ## fits a model of observed number counts to magnitude data
  ## x=vector of magnitude
  ## para=c(mlim,sigma,alpha,minit) initial paramaters

  NegLog<-function(para){
    y<-dmag_erf(x,minit=para[4],mlim=para[1],sigma=para[2],alpha=para[3])
    if(length(y)==1){if(y<10^(-12)){y<-10^(-12)}}
    if(length(y)> 1){y[y<10^(-12)]<-10^(-12)}
    y<-(-sum(log( y )))
    y
  }
  pars<-optim(para,NegLog,method="BFGS")
  pars
}

SelectiveMean<-function(dat){
  ## dat[1:n] = data vector
  ## dat[n+1:2*n] = binary vector of same length than x, indicating whether a coordinate is (1) or is not (0) missing, i.e. not taken (1) or taken (0) into account 
  n<-length(dat)/2
  x<-dat[1:n]
  miss.idx<-dat[(n+1):(2*n)]
  if(n==0 | (n-round(n))!=0 | sum(miss.idx[miss.idx!=0]!=1)>0){stop("Input data in bad format!")}
  res<-mean(x[miss.idx==0])
  res
}

quilt.plot.mod<-function (x, y, z, nrow = 64, ncol = 64, grid = NULL, add.legend = TRUE, add = FALSE, col = tim.colors(64)[64:1], ...){
  ## identical to the quilt.plot function from fields package
  ## except for allowing to specify colour scheme
  ## default colour scheme is the inverted standard scheme used by quilt.plot

    x <- as.matrix(x)
    if (ncol(x) == 2) {
        z <- y
    }
    if (ncol(x) == 1) {
        x <- cbind(x, y)
    }
    if (ncol(x) == 3) {
        z <- x[, 3]
        x <- x[, 1:2]
    }
    out.p <- as.image(z, x = x, nrow = nrow, ncol = ncol, na.rm = TRUE)
    if (add.legend) {
        image.plot(out.p, col = col, add = add, ...)
    }
    else {
        image(out.p, col = col, add = add, ...)
    }
}

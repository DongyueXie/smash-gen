#'@title Smooth Poisson sequence, accounting for nugget effect
#'@param x observed Poisson sequence
#'@param nugget nugget effect
#'@param s Scale factor for Poisson observations: y~Pois(scale*lambda), can be a vector.
#'@param transformation transformation of Poisson data, either 'vst' or 'lik_expansion'; 'vst' for varaince stablizing transformation; 'lik_expansion' for likelihood expansion
#'@param robust whether perform robust wavelet regression
#'@param robust.q quantile to determine outliers
#'@param smooth_method smoothing method for Gaussian sequence, either 'smash' or 'ti.thresh'. When n is large, ti.thresh is much faster
#'@param nug.init init value of nugget effect, either a scalar or NULL
#'@param nug.est.method method for estimating nugget effect, either 'rmad' or 'diff'.
#'@param ash.pm If choose lik_expansion, whehter use ash posterior mean approxiamtion if x=0. If not x = x+eps.
#'@param eps If choose lik_expansion, if x=0, set x = x + eps. Either input a numerical value or 'estimate'. If estimate, eps = sum(x==1)/sum(x<=1)
#'@param filter.number,family wavelet basis, see wavethresh pakcage for more details
#'@param maxiter max iterations for estimating nugget effect
#'@param tol tolerance to stop iterations.
#'@param nug.est.limit when estimating nugget effect, the proportion of top s_t to use.
#'@return estimated smoothed lambda, estimated nugget effect.
#'@import smashr
#'@import ashr
#'@export

library(smashr)
library(wavethresh)
library(ashr)

smash.gen.poiss = function(x,nugget=NULL,s=1,
                           transformation = 'lik_expansion',
                           smooth_method='smash',
                           robust = T,
                           robust.q = 0.99,
                           nug.init = NULL,
                           nug.est.method = 'diff',
                           ash.pm=FALSE,
                           smash.pm=FALSE,
                           eps='estimate',
                           filter.number = 1,
                           family = "DaubExPhase",
                           maxiter=1,
                           tol=1e-2,
                           nug.est.limit=0.2){

  if(!ispowerof2(length(x))){
    reflect.x = reflect(x)
    x = reflect.x$x
    idx = reflect.x$idx
  }else{
    idx = 1:length(x)
  }

  n = length(x)

  if(transformation == 'vst'){
    y = sqrt(x+3/8)/sqrt(s)
    st = rep(sqrt(0.25/s),n)
  }

  if(transformation == 'lik_expansion'){

    if(smash.pm){
      lambda_tilde = smash.poiss(x)/s
    }
    lambda_tilde = x/s
    if(min(x)<1){
      if(ash.pm){
        x_pm = ash(rep(0,n),1,lik=lik_pois(x,scale=s,link='log'),
                   optmethod='mixSQP',pointmass=F)$result$PosteriorMean
        lambda_tilde[x<1] = x_pm[x<1]
      }else{
        if(eps == 'estimate'){
          eps = sum(round(x)==1)/sum(round(x)<=1)+0.1
        }
        lambda_tilde[x<1] = (x[x<1]+eps)/s
      }
    }
    # working data
    st=sqrt(1/(s*lambda_tilde))
    y=log(lambda_tilde)+(x-s*lambda_tilde)/(s*lambda_tilde)

  }


  # estimate nugget effect and estimate mu

  if(robust){

    #win.size = round(log(n,2)/2)*2+1
    #browser()
    #y.wd = wd(y,filter.number,family,'station')
    #y.wd.coefJ = accessD(y.wd,level = log(n,2)-1)
    y.wd.coefJ = (y-c(y[-1],y[1]))/sqrt(2)

    win.size = round(sqrt(n)/2)*2+1
    y.rmed = runmed(y,win.size)

    robust.z = qnorm(0.5+robust.q/2)

    if(is.null(nugget)){
      #nug.init = uniroot(normaleqn,c(-1e6,1e6),y=y,mu=y.rmed,st=st)$root
      #nug.init = max(c(0,nug.init))
      outlier.idx = which(abs(y-y.rmed)>=(robust.z*caTools::runmad(y.wd.coefJ,win.size,endrule = 'mad')))
    }else{
      outlier.idx = which(abs(y-y.rmed)>=robust.z*sqrt(st^2+nugget))
    }
    # ??
    if(length(outlier.idx)!=0){
      st[outlier.idx] = 1.4826*abs((y.wd.coefJ)[outlier.idx] - median(y.wd.coefJ))
    }

  }else{
    outlier.idx = NULL
  }


  if(is.null(nugget)){
    fit0 = NuggetEst(y,st,nug.init,nug.est.method,nug.est.limit,outlier.idx=outlier.idx,
                    smooth_method,filter.number = filter.number,family = family,maxiter,tol)
    nug.est = fit0$nug.est
  }else{
    nug.est = nugget
  }

  if(smooth_method=='smash'){
    fit = smash.gaus(y,sigma=sqrt(st^2+nug.est),
                     filter.number = filter.number,family = family,
                     post.var = TRUE)
  }
  if(smooth_method == 'ti.thresh'){
    fit = ti.thresh(y,sigma=sqrt(st^2+nug.est),filter.number = filter.number,family = family)
    fit = list(mu.est = fit, mu.est.var = rep(0,length(fit)))
  }
  mu.est = (fit$mu.est)[idx]
  mu.est.var = (fit$mu.est.var)[idx]

  if(transformation == 'vst'){
    lambda.est = mu.est^2-3/(8*s)
  }
  if(transformation == 'lik_expansion'){
    lambda.est = exp(mu.est+mu.est.var/2)
  }

  return(list(lambda.est=lambda.est,mu.est=mu.est,nugget.est=nug.est))
}



normaleqn=function(nug,y,mu,st){
  return(sum((y-mu)^2/(nug+st^2)^2)-sum(1/(nug+st^2)))
}

NuggetEst=function(y,st,nug.init=NULL,nug.est.method = 'diff',nug.est.limit,
                   outlier.idx=NULL,smooth_method,filter.number,family,maxiter,tol){
  #initialize nugget effect sigma^2
  #browser()
  n=length(y)

  st.temp = st
  if(!is.null(outlier.idx)){
    st.temp[outlier.idx] = Inf
  }
  top.idx = order(st.temp,decreasing = F)[1:round(n*nug.est.limit)]

  if(nug.est.method=='diff'){
    x.m=c(y[n],y,y[1])
    st.m=c(st[n],st,st[1])
    nug.init = ((x.m[2:(n+1)]-x.m[3:(n+2)])^2+(x.m[2:(n+1)]-x.m[1:(n)])^2-2*st.m[2:(n+1)]^2-st.m[1:(n)]^2-st.m[3:(n+2)]^2)/4
    # take out outliers
  }

  if(nug.est.method=='rmad'){
    win.size = round(n/10)
    odd.boo = (win.size%%2 == 1)
    win.size = win.size + (1 - odd.boo)
    #browser()
    #y.wd = wd(y,filter.number,family,'station')
    #y.wd.coefJ = accessD(y.wd,level = log(n,2)-1)
    y.wd.coefJ = (y-c(y[-1],y[1]))/sqrt(2)
    total_var = (caTools::runmad(y.wd.coefJ,win.size,endrule = 'mad'))^2
    nug.init = total_var - st^2
  }
  nug.est = mean(nug.init[top.idx])
  nug.est = max(0,nug.init)
#
#   if(is.null(nug.init)){
#     x.m=c(y[n],y,y[1])
#     st.m=c(st[n],st,st[1])
#     nug.init = ((x.m[2:(n+1)]-x.m[3:(n+2)])^2+(x.m[2:(n+1)]-x.m[1:(n)])^2-2*st.m[2:(n+1)]^2-st.m[1:(n)]^2-st.m[3:(n+2)]^2)/4
#
#     nug.init = nug.init[nug.init>0&nug.init<var(y)]
#     nug.init = median(nug.init)
#
#     #
#
#   }
  #given st and nug to estimate mean

  for(iter in 1:maxiter){

    #print(nug.est)

    # update mean
    if(smooth_method == 'smash'){
      est = smash.gaus(y,sigma=sqrt(st^2+nug.est),filter.number = filter.number,family = family,post.var = TRUE)
      mu.est = est$mu.est
      mu.est.var = est$mu.est.var
    }
    if(smooth_method == 'ti.thresh'){
      mu.est = ti.thresh(y,sigma=sqrt(st^2+nug.est),filter.number = filter.number,family = family)
      mu.est.var = rep(0,n)
    }

    # update nugget effect
    nug.est.new=uniroot(normaleqn,c(-1e6,1e6),y=y,mu=mu.est,st=st)$root
    nug.est.new = max(c(0,nug.est.new))


    if(abs(nug.est.new - nug.est)<=tol){
      break
    }else{
      nug.est = nug.est.new
    }
  }

  return(list(mu.est = mu.est,mu.est.var=mu.est.var,nug.est=nug.est))

}

#'@export
ispowerof2 <- function (x){
  x >= 1 & 2^ceiling(log2(x)) == x
}



smash.pois.gaus = function(x){
  l = smash.poiss(x)
  fit = smash.gaus(log(l),joint = TRUE,post.var = TRUE)
  return(list(mu.est = fit$mu.res$mu.est,
         lambda.est = exp(fit$mu.res$mu.est+fit$mu.res$mu.est.var/2),
         nug.est = fit$var.res$var.est))
}


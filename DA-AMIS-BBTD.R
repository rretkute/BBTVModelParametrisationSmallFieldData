# devtools::install_github("rretkute/AMISEpi")
library(AMISEpi)
library(matrixStats)
library(statip)
library(ggplot2)
library(gridExtra)
library(viridis)
library(zoo)
library(lubridate)


################################################################
##    Data for model parametrisation
################################################################

plants<-readRDS( "data_for_parametrisation_Benin.rds")
plants$planted<-as.Date(plants$planted)
plants$removed<-as.Date(plants$removed)

ggplot(plants[plants$planted==as.Date("2016-12-01"),]) +
  geom_point(aes(x=x, y=y, fill=Infected), pch=21, col="black", size=5) +
  scale_fill_manual(values=c("#80cdc1", "#a6611a"))+
  theme_bw()

# Length of observation window
Tobs<-max(plants$R[plants$R<Inf])+1

# Distances between plants
dist.plants<-as.matrix(dist(plants[, 2:3], method = "euclidean", 
                    diag = TRUE, upper = TRUE, p = 2))

################################################################
#  Functions
################################################################

# Leaf emergence rate
k0<-0.056; k1<-0.062; TT<-90
LER <- function(t, k0, k1, TT) {
  k0*sin(2*pi*(t-TT )/365) + k1
}

incubation.period<- function(d, param0, dT){
  ds<-d %% 365
  tmp<-cumsum(LER(seq(ds, min(ds+3*365, dT), by=1), param0[1], param0[2], param0[3]))
  wh<-which(tmp>=2)
  if(length(wh)>0){
    ans<-min(wh)
  }
  return(ans)
}

latent.period<- function(d, param0, dT){
  ds<-d %% 365
  tmp<-cumsum(LER(seq(ds, min(ds+3*365, dT), by=1), param0[1], param0[2], param0[3]))
  wh<-which(tmp>=3.7)
  if(length(wh)>0){
    ans<-min(wh)
  }
  return(ans)
}

#  Pre-calculate incubation & latent periods 
tt<-seq(-365, 365*10, 1)
p1<-sapply(1:length(tt), function(a) 
  incubation.period(tt[a], c(k0, k1, TT), 365*10))
p2<-sapply(1:length(tt), function(a) 
  latent.period(tt[a], c(k0, k1, TT), 365*10))
inc.lat.p<-data.frame(Time=tt, Incub.per=p1, Lat.per=p2)
inc.lat.p$Inf.day<-inc.lat.p$Time+inc.lat.p$Incub.per
inc.lat.p$Sympt.day<-inc.lat.p$Time+inc.lat.p$Lat.per
GetIncubPeriod<-approxfun(inc.lat.p$Time, inc.lat.p$Incub.per)
GetIncubPeriodRev<-approxfun(inc.lat.p$Inf.day[seq(1,nrow(inc.lat.p), 9)], 
                             inc.lat.p$Incub.per[seq(1,nrow(inc.lat.p), 9)])
GetLatPeriod<-approxfun(inc.lat.p$Time, inc.lat.p$Lat.per)
GetLatPeriodRev<-approxfun(inc.lat.p$Sympt.day[seq(1,nrow(inc.lat.p), 12)], 
                           inc.lat.p$Lat.per[seq(1,nrow(inc.lat.p), 12)])


# Transmission kernel
K<-function(d, alpha){
  exp(-d/alpha)
}

likelihood<- function(obs, dd, Tobs, parameters){
  ll<-0
  # Propose Infection times
  plants<-get.epi.times(obs, parameters)
  ## INFECTION & REMOVAL
  II<-which(plants$infected==1)
  for(ii in II){
    #  Infected suckers 
    if(plants$planted[ii]>as.Date("2016-12-01")){
      wh1<-which(plants$P<plants$E[ii] & plants$I<plants$E[ii])
      wh2<-which(plants$P<plants$E[ii] & plants$I<plants$E[ii] & plants$R>=plants$E[ii] )
      ll1<-parameters$epsilon*(plants$E[ii]-plants$P[ii])
      if(length(wh1)>0){
        for(jj in wh1){
          ll1<-ll1+parameters$beta*K(dd[ii,jj], parameters$alpha)*
            max(min(plants$E[ii], plants$R[jj])-max(plants$P[ii],plants$I[jj]),0)
        }
      }
      ll2<-parameters$p +(1-parameters$p)*(parameters$epsilon + parameters$beta* 
                                             sum(K(dd[ii,wh2], parameters$alpha)))*
        exp(-ll1)
      ll<- ll + log(ll2)
    } else {
      # Instantaneous FOI at time of infection 
      wh<-which(plants$P<plants$E[ii] & plants$I<plants$E[ii] & plants$R>=plants$E[ii] )
      ll<- ll+ log(parameters$epsilon + parameters$beta* 
                     sum(K(dd[ii,wh], parameters$alpha)))
      #FOI before infection
      ll<-ll-parameters$epsilon*(plants$E[ii]-plants$P[ii])
      wh<-which(plants$P<plants$E[ii] & plants$I<plants$E[ii])
      if(length(wh)>0){
        for(jj in wh){
          ll<-ll-parameters$beta*K(dd[ii,jj], parameters$alpha)*
            max(min(plants$E[ii], plants$R[jj])-max(plants$P[ii],plants$I[jj]),0)
        }
      }
    }
    ##  REMOVAL
    if(plants$R[ii]<Tobs){
      ll<-ll + log(parameters$gamma) 
    }
    ll<-ll -(parameters$gamma)*(min(plants$R[ii],Tobs) - plants$I[ii])
  }
  ##  NO INFECTION
  II<-which(plants$infected==0)
  for(ii in II){
    #  Healthy suckers 
    if(plants$planted[ii]>as.Date("2016-12-01")){
      ll<-ll+log(1-parameters$p)
    }
    ll<-ll-parameters$epsilon*(Tobs-plants$P[ii])
    wh<-which(plants$infected==1)
    if(length(wh)>0){
      for(jj in wh){
        ll<-ll-parameters$beta*K(dd[ii,jj], parameters$alpha)*
          max(min(plants$R[jj],Tobs)-max(plants$P[ii], plants$I[jj]), 0)
      }
    }
  }
  return(ll)
}

get.epi.times<-function(plants, parameters){
  obs<-plants
  obs$E<-Inf
  obs$I<-Inf
  for(i in 1:nrow(obs)){
    if(obs$infected[i]==1){
      obs$I[i]<-max(obs$P[i], obs$R[i]-1/parameters$gamma) # I->R
      pp<-GetLatPeriodRev(obs$I[i]) 
      obs$E[i]<-max(obs$P[i],obs$I[i]-pp) # E->I
    }
  }
  return(obs)
}

################################################################
#  DA-AMIS  
################################################################

start.time <- Sys.time()

# Set target effective sample size
ESS.R<-1000

# Minimal ESS for sampling from target distribution
ESS.min1<-5

# Max number of iterations
T.max<-1000

# Set number of particles for each iteration
NN<-rep(1000, T.max)

# Set ranges for parameters
n.param<-5
par.min<-c(0, 0, 5, 1/100, 0)
par.max<-c(1, 1, 30, 1/10, 1)

# Prior
dprop0<-function(pp){ 
  sum(sapply(1:length(pp), function(a) dunif(pp[a], min=par.min[a], max=par.max[a], log=TRUE)))
}

# Density value of prior 
prior.dns<-dprop0(c(0.5, 0.5, 10, 1/30, 0.5))

#  Set up for mixture function
proposal=mvtComp(df=3); mixture=mclustMix();
dprop <- proposal$d
rprop <- proposal$r

# Initialise parameter and proposalvariables
param<-matrix(NA, ncol=n.param+2, nrow=0)
Sigma <- list(NA)
Mean<-list(NA)
PP<-list(NA)
GG<-list(NA)

# Iteration 1
it<-1
cat(c("Started iteration ", it,"\n"))
ii<-nrow(param)
while(ii<NN[it]){
  new.par<-sapply(1:n.param, function(a) runif(1, min=par.min[a], max=par.max[a]))
  parameters<-data.frame(epsilon=new.par[1], beta=new.par[2], 
                         alpha=new.par[3], gamma=new.par[4], p=new.par[5])
  ll<- likelihood(plants, dist.plants, Tobs, parameters)
  if(!is.na(ll)){
    if(ll>-Inf & ll<Inf){
      param<-rbind(param, c(new.par, ll, it))
      ii<-ii+1
    }}
}

# Calculate ESS and weighst
q <- sapply(1:nrow(param), function(b)  dprop0(as.numeric(param[b,1:n.param])))
wl<- param[,n.param+1]-q
CC<-logSumExp(wl)
WW<-exp(wl-CC)
WW<-WW/sum(WW)
if(sum(WW)>0) {
  ess<-(sum((WW)^2))^(-1)
} else {
  ess<-0
}

print(ess)
ESS<-data.frame(iteration=1, ess=ess)

# Iterations 2+
stop<-0
while(stop==0){
  it<-it+1
  cat(c("\n Started iteration ", it,"\n"))
  if(ess<ESS.min1) {
    if(it<10){
      WW<-WW+10^(-3)
    } else {
      thr<-param[order(-param[,n.param+1]),n.param+1][100]
      WW<-0*WW
      WW[which(param[, n.param+1]>=thr)]<-1
    }
  }
  J<-sample(1:sum(NN[1:(it-1)]), NN[it], prob= WW, replace=TRUE)
  xx<-param[J,1:n.param]
  clustMix <- mixture(xx)
  G <- clustMix$G
  cluster <- clustMix$cluster
  ### Components of the mixture
  ppt <- clustMix$alpha
  muHatt <- clustMix$muHat
  varHatt <- clustMix$SigmaHat
  GG[[it-1]]<-G
  G1<-0; G2<-G
  if(it>2) {
    G1<-sum(sapply(1:(it-2), function(a) GG[[a]]))
    G2<-sum(sapply(1:(it-1), function(a) GG[[a]]))
  }
  for(i in 1:G){
    Sigma[[i+G1]] <- varHatt[,,i]
    Mean[[i+G1]] <- muHatt[i,]
    PP[[i+G1]]<-ppt[i]
  }
  
  # Draw new parameters and calculate log likelihood
  while(nrow(param)<sum(NN[1:it])){
    compo <- sample(1:G,1,prob=ppt)
    x1 <- t(rprop(1,muHatt[compo,], varHatt[,,compo]))
    new.param<-as.numeric(x1)
    if(dprop0(new.param)>-Inf){
      parameters<-data.frame(epsilon=new.param[1], beta=new.param[2], 
                             alpha=new.param[3], gamma=new.par[4], p=new.param[5])
      ll<- likelihood(plants, dist.plants, Tobs, parameters)
      if(!is.na(ll)){
        if(ll>-Inf & ll< Inf){
          param<-rbind(param, c(as.numeric(new.param), ll, it))
        }}}
  }
  
  q <- (NN[1]/sum(NN[1:it]))*exp(prior.dns) +
    (sum(NN[2:it])/sum(NN[1:it]))* rowSums(as.matrix(sapply(1:G2, 
      function(a) PP[[a]] * dprop(param[,1:n.param],mu= Mean[[a]], Sig=Sigma[[a]]))))
  # Calculate ESS and weights
  wl<- param[,n.param+1]-log(q)
  CC<-logSumExp(wl)
  WW<-exp(wl-CC)
  WW<-WW/sum(WW)
  if(sum(WW)>0) {
    ess<-(sum((WW)^2))^(-1)
  } else {
    ess<-0
  }
  ESS<-rbind(ESS, data.frame(iteration=it, ess=ess))
  if(min(ess)>=ESS.R) stop<-1
  if(it>= TT) stop<-1
  # Print effective sample size and maximum value of loglikelihood
  cat(c("\n", ess, "", max(param[, n.param+1]), "\n"))
}

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

nrow(param)

ggplot(ESS)+
  geom_path(aes(x=iteration, y=ess)) +
  geom_point(aes(x=iteration, y=ess)) +
  scale_y_continuous(trans='log10')+
  theme_bw()+xlab("Iteration") +ylab("Effective sample size")

yy<-data.frame(Iteration=param[,n.param+2], ll=param[, n.param+1])
ggplot(yy)+
  geom_boxplot(aes(x=as.factor(Iteration), y=ll)) +
  theme_bw() +xlab("Iteration")+
  ylab("Loglikelihood")

ggplot(yy[yy$Iteration>=10,])+
  geom_boxplot(aes(x=as.factor(Iteration), y=ll)) +
  theme_bw() +xlab("")+
  ylab("")

(n.max<-max(yy$Iteration))
range(yy$ll[yy$Iteration==1])
range(yy$ll[yy$Iteration==n.max])

whP<-sample(1:length(WW), 1000, prob=WW)
posterior_samples <-data.frame(epsilon=param[whP,1]*10^4,
                beta=param[whP,2]*10^4,
                alpha=param[whP,3],
                gamma=param[whP,4],
                p=param[whP,5])


labels <- c("epsilon*10^4", "beta*10^4","alpha", "gamma", "p")
# Pairwise plot
ggpairs(posterior_samples,
        columnLabels = labels,
        diag = list(continuous = wrap("densityDiag", alpha = 0.6, fill = "skyblue")),
        upper = list(continuous = wrap("cor", size = 4)),
        lower = list(continuous = wrap("points", alpha = 0.3, size=1))) +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


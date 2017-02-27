#' A function to simulate epidemics using the Reed-Frost model
#'
#' This function allows you to simulate epidemics.
#' @param Nsim.
#' @keywords Reed-Frost
#' @export
#' @examples
#' simulateReedFrost()

simulateReedFrost <- function(Nsim=100000, Npop=10000, p=0.00015, I0=1, Nit=50, fig=T)
{
  S0=Npop-I0
  I=t(matrix(rep(c(I0,rep(0,Nit-1)),Nsim),ncol = Nsim))
  S=t(matrix(rep(c(S0,rep(0,Nit-1)),Nsim),ncol = Nsim))
  for(t in 2:Nit)
  {
    I[,t]=rbinom(Nsim,size=S[,t-1],1-(1-p)^I[,t-1])
    S[,t]=S[,t-1]-I[,t]
  }

  if(fig)
    hist(apply(I,1,max),breaks = 100, freq = T)

  return(list(S=S,I=I))
}

RF_with_obs<-function(Npop=10000, p=0.00015)
{
  initial_pop=Npop-rpois(n = 1, lambda = 20)
  initial_I=rpois(n=1, lambda = 5)
  simu_epi=simulateReedFrost(Nsim=1,Npop=initial_pop, p=p, I0=initial_I, fig=F)
  X<-cbind(simu_epi$S[1,],simu_epi$I[1,]) #the "hidden" Markov chain
  Y<-rnbinom(n=length(X[,1]),mu=.05*X[,2], size=10) #the observered process

  return(list(X=X,Y=Y))
}


RF_SIR <- function(Y, Np=1000, Npop=10000)
{
  N=length(Y)

  #Sequential important sampling
  Xp<-array(rep(0,2*N*Np),dim=c(Np,N,2))      #particles
  gammap<-matrix(rep(0,N*Np),ncol=N)  #partial likelihood
  wp<-matrix(rep(0,N*Np),ncol=N)      #unormalised importance weights
  Wp<-matrix(rep(0,N*Np),ncol=N)      #normalised importance weights
  #Xrp<-matrix(rep(0,N*Np),ncol=N)     #resampled particles
  A<-matrix(rep(0,N*Np),ncol=N) #Ancestry of particle

  #Initialisation of the algorithm
  Xp[,1,1]<-Npop-rpois(n = Np, lambda = 20) #initialisation of the size of susceptible population for each particle
  Xp[,1,2]<-rpois(n=1, lambda = 5) #initialisation of the size of the infectious population for each particle
  gammap[,1]<-dnorm(Xp[,1],mean=0, sd=sqrt(sigma^2/(1-alpha^2)))*dnorm(Y[1],mean=0,sd=beta*exp(Xp[,1]/2))
  wp[,1]<-gammap[,1]/dnorm(Xp[,1],mean=0, sd=sqrt(sigma^2/(1-alpha^2)))
  Wp[,1]<-wp[,1]/sum(wp[,1])
}

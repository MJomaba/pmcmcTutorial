alpha=0.91
sigma=1
beta=0.5

N=500
Np=1000

X=rep(0,N)
Y=rep(0,N)

X[1]<-rnorm(1,mean=0, sd=sqrt(sigma^2/(1-alpha^2)))
Wn<-rnorm(1,mean=0,sd=1)
Y[1]<-beta*exp(X[1]/2)*Wn

for(i in 2:N)
{
  Vn<-rnorm(1,mean=0,sd=1)
  Wn<-rnorm(1,mean=0,sd=1)
  X[i]<-alpha*X[i-1]+sigma*Vn
  Y[i]<-beta*exp(X[i]/2)*Wn
}

plot(X,type="l",col="blue")
points(Y,pch="*",col="red")

#Sequential important sampling
Xp<-matrix(rep(0,N*Np),ncol=N)      #particles
gammap<-matrix(rep(0,N*Np),ncol=N)  #partial likelihood
wp<-matrix(rep(0,N*Np),ncol=N)      #unormalised importance weights
Wp<-matrix(rep(0,N*Np),ncol=N)      #normalised importance weights
Xrp<-matrix(rep(0,N*Np),ncol=N)     #resampled particles

#Initialisation of the algorithm
Xp[,1]<-rnorm(Np,mean=0, sd=sqrt(sigma^2/(1-alpha^2)))
gammap[,1]<-dnorm(Xp[,1],mean=0, sd=sqrt(sigma^2/(1-alpha^2)))*dnorm(Y[1],mean=0,sd=beta*exp(Xp[,1]/2))
wp[,1]<-gammap[,1]/dnorm(Xp[,1],mean=0, sd=sqrt(sigma^2/(1-alpha^2)))
Wp[,1]<-wp[,1]/sum(wp[,1])
index_T<-rep(0,N)
index_2<-rep(0,Np)

#Sequential calculation of particles and importance weights
for(i in 2:N)
{
  #IS step
  Xp[,i]<-rnorm(Np,mean=alpha*Xp[,i-1], sd=sigma)
  #gammap[,i]<-gammap[,i-1]*dnorm(Xp[,i],mean=alpha*Xp[,i-1], sd=1)*dnorm(Y[i],mean=0,sd=beta*exp(Xp[,i]/2))
  #wp[,i]<-gammap[,i]/dnorm(Xp[,i],mean=alpha*Xp[,i-1], sd=sigma)
  wp[,i]<-dnorm(Xp[,i],mean=alpha*Xp[,i-1], sd=1)*dnorm(Y[i],mean=0,sd=beta*exp(Xp[,i]/2))/dnorm(Xp[,i],mean=alpha*Xp[,i-1], sd=sigma)
  Wp[,i]<-wp[,i]/sum(wp[,i])
  
  #Resampling step -using systematic resampling
  U1<-runif(1,min = 0,max = 1/Np)
  if(i==2)
    U2<-U1
  cumWj<-0
  lU<-U1
  index<-1
  for(j in 1:Np)
  {
    cumWjminus<-cumWj
    cumWj<-cumWj+Wp[j,i]
    if(lU<cumWj) #test if at least one Uk is between two corresponding cumulative Wp
    {
      Nji<-1+floor((cumWj-lU)*Np) 
      Xrp[index:(index+Nji-1),i]<-Xp[j,i]
      lU<-lU+Nji/Np
      index<-index+Nji
    }
  }
  Xp[,i]<-Xrp[,i]
  Wp[,i]<-1/Np
  index_T[i]<-index-1
}
  
#Compute the mean particles path
mu_Xp<-rep(0,N)
ESS<-rep(0,N)
low_Xp<-rep(0,N)
up_Xp<-rep(0,N)
for(i in 1:N)
{
  ESS[i]<-1/sum(Wp[,i]^2)
  mu_Xp[i]<-sum(Xp[,i]*Wp[,i])
  points(i,mu_Xp[i],pch=19,col="#00aa0066")
}


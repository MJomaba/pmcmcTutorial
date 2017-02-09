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
#Use yule-walker equations to find estimates for phi for a AR(n) model, given vector of ACVF values up to lag-n.
#Specify confidence intervals with alpha

yw.p<-function(gamma, n, p, alpha){
  g = matrix(0, nrow=p, ncol=p)
  y = c(0, length=p)
  for(i in 1:p) {
    for(j in 1:p) {
      g[i, j] = gamma[abs(j-i)+1]
    }
  }
  for(i in 1:p) {
    y[i] = gamma[i+1]
  }
  
  phi = c(solve(g) %*% y)
  sigma = drop(gamma[1] - t(phi)%*%y)
  variance = solve(g)*(sigma[1]/n)
  
  a = (1-alpha)/2
  sd.phi = sqrt(variance[1,1])
  ci = list()
  for(k in 1:p) {
    ci[[paste("phi_", k)]] = c(qnorm(a, mean=phi[k], sd=sd.phi), 
              qnorm(1-a, mean=phi[k], sd=sd.phi))
  }
  
  final = list("phi"=phi, "sigma"=sigma, "variance"=variance[1,1], 
               "confidence interval" = ci)
  return(final)
}

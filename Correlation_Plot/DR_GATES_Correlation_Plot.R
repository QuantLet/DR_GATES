if(!require("clusterGeneration")) install.packages("clusterGeneration"); library("clusterGeneration")
if(!require("mvtnorm")) install.packages("mvtnorm"); library("mvtnorm")
if(!require("psych")) install.packages("psych"); library("psych")



correlation_plot<- function(N,p) {
  
  N = N # Number of Observations
  p = p # Number of Covariates 
  
  # = Generate covariance matrix of z = #
  sigma <- genPositiveDefMat(p, "unifcorrmat")$Sigma
  sigma <- cov2cor(sigma)
  
  
  z <- rmvnorm(N, sigma = sigma) # = Generate z = #
  z_plot <- as.data.frame(z)
  
  pairs.panels(z[,1:10], 
               method = "pearson", # correlation method
               hist.col = "blue",
               density = TRUE,  # show density plots
               ellipses = TRUE, # show correlation ellipses
               show.points = F)
}

correlation_plot(100,10)

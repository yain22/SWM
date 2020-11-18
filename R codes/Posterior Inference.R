# Load Spatial Weibull Model (SWM)
{
  setwd("C:/Users/seyoonlee/OneDrive/Documents/GitHub/SWM/R codes")
  source("SWM.R")
}

# Load Dataset
{
  setwd("C:/Users/seyoonlee/OneDrive/Documents/GitHub/SWM/Dataset/")
  load("Production_Results_Eagle_Ford.RData")
  
  # Set number of testing wells
  N.test = 36
  N = 360 - N.test # number of training wells   
  
  
  Y = log(Production_Results_Eagle_Ford$P[1:N,])
  X = scale(Production_Results_Eagle_Ford$X, center = TRUE, scale = TRUE)[1:N,]
  Loc = Production_Results_Eagle_Ford$Loc[1:N,]
}

# Implement Gibbs sampler
{
  Post.Inf.Res = SWM(Y = Y, X = X, Loc = Loc,
                     burn=100,nmc=100,thin=10)
}


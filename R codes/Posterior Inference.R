# Load Spatial Weibull Model (SWM)
{
  setwd("C:/Users/seyoonlee/OneDrive/Documents/GitHub/SWM/R codes")
  source("SWM.R")
}

# Load Dataset
{
  setwd("C:/Users/seyoonlee/OneDrive/Documents/GitHub/SWM/Dataset/")
  library("readr")
  time_series_F <- read_csv("time_series(F).csv", col_types = cols(API10 = col_number()))
  cov_data_F <- read_csv("cov_data(F).csv", col_types = cols(API10 = col_number()))
  time_series_F = as.matrix(time_series_F) 
  cov_data_F = as.matrix(cov_data_F)
}

# Divide the full dataset into test:training dataset
{
  #1. Set number of test set
  N.test = 36 
  #2. Set number of total number of set (maximum is 360)
  Total.no.well = 360
  #3. Then number of training set is automatically set
  N = Total.no.well - N.test  # number of training set
  
  # Response
  Y = log(time_series_F[1:N,-1])
  # Design matrix 
  X = scale(cov_data_F[,5:15], center = TRUE, scale = TRUE)[1:N,]
  
  # Spatial locations
  s1 = as.matrix(cov_data_F[,16]) # Latitude
  colnames(s1) = "Latitude"
  s2 = as.matrix(cov_data_F[,17]) # Longitude
  colnames(s2) = "Longitude"
  Loc = cbind(s1, s2)[1:N,] # (Latitude, Longitude)
}

# Implement Gibbs sampler
{
  Post.Inf.Res = SWM(Y = Y, X = X, Loc = Loc,
                     burn=1000,nmc=1000,thin=20)
}

Post.Inf.Res$thinned.beta.1[1,]


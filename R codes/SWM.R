# Spatial Weibull Model (SWM) 
# Update: 11/17/2020
# (c) 2020 Se Yoon Lee and Bani Mallick All Rights Reserved
# User should contact seyoonlee@stat.tamu.edu and bmallick@stat.tamu.edu 
# for any use and modification for the code for the publication purpose or industrial use

SWM = function(Y,X,Loc,
               seed.no=1,burn=100,nmc=100,thin=1,prop.var.theta.2=1/pi^5,prop.var.theta.3=1/pi^5,
               rho.1=-4,rho.2=-4,rho.3=-4){
  
  # Training dataset
    # Y: N-by-T matrix for log-scaled oil production rate time sereis data for N wells
    # X: N-by-p matrix from the N wells: X should be column-wised standardized "in advance"
    # Loc: N-bt-2 matrix for the location coordinates (Latitude, Longitude) from N wells
  
  # Simulation Settings
    # seed.no: number of seed for the Gibbs sampling algorithm 
    # burn: no of burn
    # nmc: no of iterations after burn
    # thin: no of thining for the nmc
    # pro.var.theta.2 and pro.var.theta.2: proposal variances for the MH algorithms to sample from theta.2.i and theta.3.i  
  
  # Hyper-parameters
    # rho.1, rho.2, and rho.3: range parameters for the Gaussian Kernel
  
  # Read packages
  {
    library("mvtnorm")  
    library("fields")
  }
  
  # Data work
  {
    # Define index
    N = nrow(Y) ; p = ncol(X)
    
    # Log-scaled oil production rate of the ith well
    y = function(i) {
      temp.y = as.numeric(Y[i,]) # LOG PRODUCTION
      start.i = sum(as.numeric(is.na(temp.y))) + 1
      y = c(temp.y[start.i:ncol(Y)])
      return(y)
    }
    
    # Distance matrix
    dist.mat <- rdist(Loc)
    
    # Time matrix
    t = c()
    for (i in 1:N){
      t[i] = length(y(i))
    }
    T.mat = diag(t) 
  }
  
  # MCMC setting
  {
    burn = burn
    nmc = nmc
    thin = thin
    S = burn + nmc
  }
  
  # Make a room to store MCMC samples and set initial values
    ## I. Nonspatial random variables
    {
    # Rooms
    alpha.1 = matrix(rep(0,S), nrow = 1) # 1 dim vector
    alpha.2 = matrix(rep(0,S), nrow = 1) # 1 dim vector 
    alpha.3 = matrix(rep(0,S), nrow = 1) # 1 dim vector 
    beta.1 = matrix(rep(0,p*S),nrow = p) # p dim vector
    beta.2 = matrix(rep(0,p*S),nrow = p) # p dim vector
    beta.3 = matrix(rep(0,p*S),nrow = p) # p dim vector
    lambda.sq.1 = matrix(rep(0,p*S),nrow = p) # p dim vector
    lambda.sq.2 = matrix(rep(0,p*S),nrow = p) # p dim vector
    lambda.sq.3 = matrix(rep(0,p*S),nrow = p) # p dim vector
    nu.1 = matrix(rep(0,p*S),nrow = p) # p dim vector
    nu.2 = matrix(rep(0,p*S),nrow = p) # p dim vector
    nu.3 = matrix(rep(0,p*S),nrow = p) # p dim vector
    tau.sq.1 = matrix(rep(0,S), nrow = 1) # 1 dim vector
    tau.sq.2 = matrix(rep(0,S), nrow = 1) # 1 dim vector
    tau.sq.3 = matrix(rep(0,S), nrow = 1) # 1 dim vector
    xi.1 = matrix(rep(0,S), nrow = 1) # 1 dim vector
    xi.2 = matrix(rep(0,S), nrow = 1) # 1 dim vector
    xi.3 = matrix(rep(0,S), nrow = 1) # 1 dim vector
    sigma.sq.1 = matrix(rep(0,S),nrow = 1) # 1 dim vector
    sigma.sq.2 = matrix(rep(0,S),nrow = 1) # 1 dim vector
    sigma.sq.3 = matrix(rep(0,S),nrow = 1) # 1 dim vector
    theta.1 = matrix(rep(0,N*S),nrow = N) # N dim vector
    theta.2 = matrix(rep(0,N*S),nrow = N) # N dim vector
    theta.3 = matrix(rep(0,N*S),nrow = N) # N dim vector
    sigma.sq = matrix(rep(0,S),nrow = 1) # 1 dim vector
    phi = matrix(rep(0,S), nrow = 1) # 1 dim vector
    
    # Initial values
    alpha.1[1] = 16.5 # 1 dim vector
    alpha.2[1] = -2.5 # 1 dim vector 
    alpha.3[1] = -0.8 # 1 dim vector 
    beta.1[,1] = rep(0,p) # p dim vector
    beta.2[,1] = rep(0,p) # p dim vector
    beta.3[,1] = rep(0,p) # p dim vector
    lambda.sq.1[,1] = rep(1,p) # p dim vector
    lambda.sq.2[,1] = rep(1,p) # p dim vector
    lambda.sq.3[,1] = rep(1,p) # p dim vector
    nu.1[,1] = rep(1,p) # p dim vector
    nu.2[,1] = rep(1,p) # p dim vector
    nu.3[,1] = rep(1,p) # p dim vector
    tau.sq.1[1] = 1 # 1 dim vector
    tau.sq.2[1] = 1 # 1 dim vector
    tau.sq.3[1] = 1 # 1 dim vector
    xi.1[1] = 1 # 1 dim vector
    xi.2[1] = 1 # 1 dim vector
    xi.3[1] = 1 # 1 dim vector
    sigma.sq.1[1] = 1 # 1 dim vector
    sigma.sq.2[1] = 1 # 1 dim vector
    sigma.sq.3[1] = 1 # 1 dim vector
    theta.1[,1] = rep(16,N) # N dim vector
    theta.2[,1] = rep(-2,N) # N dim vector
    theta.3[,1] = rep(0.5,N) # N dim vector
    sigma.sq[1] = 1 # 1 dim vector
    phi[1] = 1 # 1 dim vector
  }
  
    ## II. Spatial random variables
    {
    # Rooms
    gamma.sq.1 = matrix(rep(0,S), nrow = 1) # 1 dim vector
    gamma.sq.2 = matrix(rep(0,S), nrow = 1) # 1 dim vector
    gamma.sq.3 = matrix(rep(0,S), nrow = 1) # 1 dim vector
    omega.1 = matrix(rep(0,S), nrow = 1) # 1 dim vector
    omega.2 = matrix(rep(0,S), nrow = 1) # 1 dim vector
    omega.3 = matrix(rep(0,S), nrow = 1) # 1 dim vector
    
    # Initial values
    gamma.sq.1[1] = 1 # 1 dim vector
    gamma.sq.2[1] = 1 # 1 dim vector
    gamma.sq.3[1] = 1 # 1 dim vector
    omega.1[1] = 1 # 1 dim vector
    omega.2[1] = 1 # 1 dim vector
    omega.3[1] = 1 # 1 dim vector
  }
  
  # Gibbs sampler (Appendix)
  set.seed(seed.no)
  for (s in 1:(S-1)){
    
    if (s == 1){ # Define necessary functions
      
      # Log-scaled function 
      {
        h.it =  function(t, theta.2.i, theta.3.i) {
          y = theta.2.i + theta.3.i + (exp(theta.3.i) - 1)*log(t) - exp(theta.2.i)*t^exp(theta.3.i)
          return(y)
        } 
        mu.it = function(t, theta.1.i, theta.2.i, theta.3.i){
          y = theta.1.i + h.it(t = t, theta.2.i = theta.2.i, theta.3.i = theta.3.i)
          return(y)
        }
        h.i = function(i, theta.2.i, theta.3.i){
          
          h.i.vector = c(rep(0,as.numeric(T.mat[i,i]))) # Room for vector
          
          for (t in 1:as.numeric(T.mat[i,i])){
            h.i.vector[t] = h.it(t, theta.2.i = theta.2.i, theta.3.i = theta.3.i)
          }
          
          return(h.i.vector)  
        }
        mu.i = function(i, theta.1.i, theta.2.i, theta.3.i){
          
          mu.i.vector = c(rep(0,as.numeric(T.mat[i,i]))) # room for vector
          
          for (t in 1:as.numeric(T.mat[i,i])){
            mu.i.vector[t] = mu.it(t = t, theta.1.i = theta.1.i, theta.2.i = theta.2.i, theta.3.i = theta.3.i)
          }
          
          return(mu.i.vector)  
          
        }
      }
      
      # Norm function
      {
        norm_vec <- function(x,y) {
          norm_vector = sqrt(sum((x-y)^2))
          if (norm_vector == Inf){
            norm_vector = 1e+308
          } else {
            norm_vector = norm_vector
          } 
          return(norm_vector)
        }  
      }
      
      # Rue's algorithm to sample from p-dimensional coefficient vector
      {
        rmvt.Rue.sample = function(Q, b){
        # Goal:
        # Sample from beta = N[mu, Sigma]
        # such that
        # mu = solve(Q)%*%b # p dim vector
        # Sigma = solve(Q) # p by p matrix
        
        # Q : p by p matrix, Precision matrix
        # b : p dim vector
        
        # Useful 1. when n >>> p (i.e. Number of data is more than covariate)
        #        2. Possibly, We want to utilize precision structure
        #        3. Use cholesky and avoiding inverting precision matrix by solving three linear equations
        p = dim(Q)[1]
        
        # Step 1
        L = t(chol(Q))
        
        # Step 2
        z = rnorm(n = p ,mean = 0,sd = 1)
        
        # Step 3
        y = solve(t(L), z)
        
        # Step 4
        v = solve(L, b)
        
        # Step 5
        theta = solve(t(L), v)
        
        # Step 6
        beta = y + theta
        beta = round((beta),10)
        return(beta)
        
      }
      }
      
      # Define necessary functions for spatial modeling
      {
        # Gaussian correlation function with range parameter rho; rho is a real number
        b0 = function(rho, distance){
          g = exp(-exp(rho)*distance^2)
          return(g)
        }
        # density of inverse gamma distribution
        dIG = function(x, alpha, beta){
          res  = ((beta^alpha)/gamma(alpha))*(x^(-alpha - 1))*exp(-beta/x)
          return(res)
        }
        # log(dIG(nomi.x | alpha, nomi.beta)/dIG(deno.x | alpha, deno.beta) ) 
        log.ratio.dIG.common.shape = function(nomi.x, common.alpha, nomi.beta, deno.x, deno.beta){
          # log(dIG(nomi.x|alpha, nomi.beta)/dIG(deno.x|alpha, deno.beta) )
          # We don't need to use gamma function
          res = common.alpha*log(nomi.beta/deno.beta) - (common.alpha + 1)*log(nomi.x/deno.x) - nomi.beta/nomi.x + deno.beta/deno.x
          return(res)
        }
        # log(dIG(nomi.x | alpha, beta)/dIG(deno.x | alpha, beta) ) 
        log.ratio.dIG.common.shape.scale = function(nomi.x, common.alpha, common.beta, deno.x){
          # log(dIG(nomi.x | alpha, beta)/dIG(deno.x | alpha, beta) ) 
          # We don't need to use gamma function
          res = - (common.alpha + 1)*log(nomi.x/deno.x) - (1/nomi.x - 1/deno.x)*common.beta
          return(res)
        }
        # log(det(B.nomi)/det(B.deno))
        log_ratio_det = function(B.nomi, B.deno, method = c("chol", "eigen")){
          if (method == "chol"){
            L1 = t(chol(B.nomi))
            L2 = t(chol(B.deno))
            g = 2*sum(log(diag(L1)*diag(1/L2)))  
          }
          
          if (method == "eigen"){
            
            
            g = log(prod(abs(eigen(B.nomi)$values)/abs(eigen(B.deno)$values)))
            
            
          }
          
          return(g)
          
        }
        # solve(B) by eigen decomposition
        solve.eigen = function(B){
          # Suppose B is positive definte
          eigen.info = eigen(B)
          
          U = eigen.info$vectors
          inv.Lamda = diag(1/abs(eigen.info$values),nrow = dim(B), ncol = dim(B))
          G = U%*%inv.Lamda%*%t(U)
          return(G)
          
        }
        # sqrt(B) by eigen decomposition
        sqrt.eigen = function(B){
          # Suppose B is positive definte
          eigen.info = eigen(B)
          
          U = eigen.info$vectors
          sqrt.Lamda = diag(sqrt(abs(eigen.info$values)) ,nrow = dim(B), ncol = dim(B))
          G = U%*%sqrt.Lamda%*%t(U)
          return(G)
        }
        
        # Let B.tilda = a*diag(N) + b*B.mat 
        # t(c.vec)%*%solve(B.tilda)%*%d.vec
        quad.inverse.B.tilda.vec = function(a, b, B.mat, c.vec, d.vec){
          #B.tilda = a*diag(N) + b*B.mat
          #t(c.vec)%*%solve(B.tilda)%*%d.vec
          N = dim(B.mat)[1]
          
          eigen.result = eigen(B.mat) # eigen decomposition
          
          D = c(abs(eigen.result$values))
          U = eigen.result$vectors
          new.c.vec = t(U)%*%c.vec
          new.d.vec = t(U)%*%d.vec
          
          res = as.numeric(sum( (new.c.vec*new.d.vec)/(a + (b)*D) ))
          
          return(res)
          
        }
        #det(B.tilda)
        det.B.tilda = function(a, b, B.mat){
          
          #B.tilda = a*diag(N) + b*B.mat
          #det(B.tilda)
          
          eigen.result = eigen(B.mat) # eigen decomposition
          
          D = c(abs(eigen.result$values))
          res = as.numeric(prod((a + (b)*D))  )
          return(res)
          
        }
        #t(X)%*%solve(B.tilda)%*%Y
        quad.inverse.B.tilda.matrix = function(a, b, B.mat, X, Y){
          #B.tilda = a*diag(N) + b*B.mat
          #t(X)%*%solve(B.tilda)%*%Y
          X = as.matrix(X)
          Y = as.matrix(Y)
          
          M = dim(X)[2]
          L = dim(Y)[2]
          
          temp.mat = matrix(rep(0,M*L),nrow = M,ncol = L)
          
          N = dim(B.mat)[1]
          
          eigen.result = eigen(B.mat) # eigen decomposition
          
          D = c(abs(eigen.result$values))
          U = eigen.result$vectors
          
          # Possible Parrellel 
          for(m in 1:M){
            for(l in 1:L){
              
              new.c.vec = t(U)%*%X[,m]
              new.d.vec = t(U)%*%Y[,l]
              
              temp.mat[m,l] = as.numeric(sum( (new.c.vec*new.d.vec)/(a + (b)*D) ))
              
            }
          }
          
          return(temp.mat)
          
          
        }
        # solve(B.tilda)
        inverse.B.tilda = function(a, b, B.mat){
          #B.tilda = a*diag(N) + b*B.mat
          #solve(B.tilda)
          eigen.result = eigen(B.mat)
          D = c(abs(eigen.result$values))
          U = eigen.result$vectors
          
          res = U%*%diag(1/(a + (b)*D))%*%t(U)
          return(res)
          
        }
        # log(det(nomi.B.tilda)/det(deno.B.tilda))
        # such that nomi.B.tilda = nomi.a*diag(N) + nomi.b*nomi.B.mat
        #           deno.B.tilda = deno.a*diag(N) + deno.b*deno.B.mat
        log.det.ratio.B.tildas = function(nomi.a, deno.a , nomi.b, deno.b, nomi.B.mat, deno.B.mat){
          
          # nomi.B.tilda = nomi.a*diag(N) + nomi.b*nomi.B.mat
          # deno.B.tilda = deno.a*diag(N) + deno.b*deno.B.mat
          # log(det(nomi.B.tilda)/det(deno.B.tilda))
          
          
          # nomi.B.mat
          nomi.eigen.result = eigen(nomi.B.mat) # eigen decomposition
          nomi.D = c(abs(nomi.eigen.result$values))
          
          
          # deno.B.mat
          deno.eigen.result = eigen(deno.B.mat) # eigen decomposition
          deno.D = c(abs(deno.eigen.result$values))
          
          res = log(prod( c((nomi.a + (nomi.b)*nomi.D))/(deno.a + (deno.b)*deno.D) ))
          
          return(res)
          
        }
        # Correlation matrix for magnitude (1), scale (2), and shape (3) corresponding to traning wells
        B0.1 = b0(rho = rho.1, distance = dist.mat)
        B0.2 = b0(rho = rho.2, distance = dist.mat)
        B0.3 = b0(rho = rho.3, distance = dist.mat)
      }
    }
    
    # Step 1.
    {
      for (j in 1:p){
        nu.1[j,(s+1)] = 1/rgamma(n=1, shape = 1, rate = 1 + 1/lambda.sq.1[j,s])
        nu.2[j,(s+1)] = 1/rgamma(n=1, shape = 1, rate = 1 + 1/lambda.sq.2[j,s])
        nu.3[j,(s+1)] = 1/rgamma(n=1, shape = 1, rate = 1 + 1/lambda.sq.3[j,s])
      }
    }
    
    # Step 2.
    {
      xi.1[s+1] = 1/rgamma(n=1, shape = 1, rate = 1 + 1/tau.sq.1[s])  
      xi.2[s+1] = 1/rgamma(n=1, shape = 1, rate = 1 + 1/tau.sq.2[s])  
      xi.3[s+1] = 1/rgamma(n=1, shape = 1, rate = 1 + 1/tau.sq.3[s])  
    }
    
    # Step 3.
    {
      for (j in 1:p){
        lambda.sq.1[j,(s+1)] = 1/rgamma(n = 1, shape = 1, rate = 1/nu.1[j,(s+1)] + (((beta.1[j,s])^2)/(2*sigma.sq.1[s]*tau.sq.1[s])))
        lambda.sq.2[j,(s+1)] = 1/rgamma(n = 1, shape = 1, rate = 1/nu.2[j,(s+1)] + (((beta.2[j,s])^2)/(2*sigma.sq.2[s]*tau.sq.2[s])))
        lambda.sq.3[j,(s+1)] = 1/rgamma(n = 1, shape = 1, rate = 1/nu.3[j,(s+1)] + (((beta.3[j,s])^2)/(2*sigma.sq.3[s]*tau.sq.3[s])))
      }
      Lambda.1 = diag(c(lambda.sq.1[,(s+1)]))
      Lambda.2 = diag(c(lambda.sq.2[,(s+1)]))
      Lambda.3 = diag(c(lambda.sq.3[,(s+1)]))
    }
    
    # Step 4.
    {
      inv.Lambda.1 = diag(1/c(lambda.sq.1[,(s+1)]))
      tau.sq.1[s+1] = 1/rgamma(n = 1, shape = (p+1)/2, rate = 1/xi.1[s+1] + (1/(2*sigma.sq.1[s]))*(t(beta.1[,s])%*%inv.Lambda.1%*%beta.1[,s]))
      inv.Lambda.2 = diag(1/c(lambda.sq.2[,(s+1)]))
      tau.sq.2[s+1] = 1/rgamma(n = 1, shape = (p+1)/2, rate = 1/xi.2[s+1] + (1/(2*sigma.sq.2[s]))*(t(beta.2[,s])%*%inv.Lambda.2%*%beta.2[,s]))
      inv.Lambda.3 = diag(1/c(lambda.sq.3[,(s+1)]))
      tau.sq.3[s+1] = 1/rgamma(n = 1, shape = (p+1)/2, rate = 1/xi.3[s+1] + (1/(2*sigma.sq.3[s]))*(t(beta.3[,s])%*%inv.Lambda.3%*%beta.3[,s]))
    }
    
    # Step 5.
    {
    # magnitude (1)
    inv.Lamda.star.1 = (1/tau.sq.1[s+1])*inv.Lambda.1
    theta.1.tilda = theta.1[,s] - alpha.1[s]*c(rep(1,N)) - X%*%beta.1[,s]
    c.1 = gamma.sq.1[s]/sigma.sq.1[s]
    
    # proposal
    sigma.sq.1.new = 1/rgamma(n = 1, shape = (N + p)/2 + 2 -2, rate = (1/2)*(quad.inverse.B.tilda.vec(a = 1, b = c.1, B.mat = B0.1, c.vec = theta.1.tilda, d.vec = theta.1.tilda)  + t(beta.1[,s])%*%inv.Lamda.star.1%*%beta.1[,s]) )
    
    # MH ratio
    log.A = dmvnorm(x = theta.1[,s], mean = alpha.1[s]*c(rep(1,N)) + X%*%beta.1[,s], sigma = sigma.sq.1.new*diag(N) + gamma.sq.1[s]*B0.1, log = TRUE)
    log.B = dmvnorm(x = beta.1[,s], mean = c(rep(0,p)), sigma = sigma.sq.1.new*Lambda.1, log = TRUE)
    log.E = dmvnorm(x = theta.1[,s], mean = alpha.1[s]*c(rep(1,N)) + X%*%beta.1[,s], sigma = sigma.sq.1[s]*diag(N) + gamma.sq.1[s]*B0.1, log = TRUE)
    log.F = dmvnorm(x = beta.1[,s], mean = c(rep(0,p)), sigma = sigma.sq.1[s]*Lambda.1, log = TRUE)
    
    log.C.over.G = log(sigma.sq.1[s]/sigma.sq.1.new)
    log.D.over.H = log.ratio.dIG.common.shape(nomi.x = sigma.sq.1[s], deno.x = sigma.sq.1.new, common.alpha = (N+p)/2 +2 -2
                                              , nomi.beta = (1/2)*(quad.inverse.B.tilda.vec(a = 1, b = gamma.sq.1[s]/sigma.sq.1.new, B.mat = B0.1, c.vec = theta.1.tilda, d.vec = theta.1.tilda) + t(beta.1[,s])%*%inv.Lamda.star.1%*%beta.1[,s]) 
                                              , deno.beta = (1/2)*(quad.inverse.B.tilda.vec(a = 1, b = gamma.sq.1[s]/sigma.sq.1[s], B.mat = B0.1, c.vec = theta.1.tilda, d.vec = theta.1.tilda) + t(beta.1[,s])%*%inv.Lamda.star.1%*%beta.1[,s]) )
    
    R = log.A + log.B - (log.E + log.F) + log.C.over.G + log.D.over.H
    
    r = exp(R)
    
    # criterion
    rat = min(r,1)
    U = runif(n=1, min = 0, max = 1)
    
    if (U < rat){
      sigma.sq.1[(s+1)] = sigma.sq.1.new
      
    } else {
      sigma.sq.1[(s+1)] = sigma.sq.1[s]
    }
    
    # scale (2)
    inv.Lamda.star.2 = (1/tau.sq.2[s+1])*inv.Lambda.2
    theta.2.tilda = theta.2[,s] - alpha.2[s]*c(rep(1,N)) - X%*%beta.2[,s]
    c.2 = gamma.sq.2[s]/sigma.sq.2[s]
    
    # proposal
    sigma.sq.2.new = 1/rgamma(n = 1, shape = (N + p)/2 + 2 -2, 
                              rate = (1/2)*(quad.inverse.B.tilda.vec(a = 1, b = c.2, B.mat = B0.2, c.vec = theta.2.tilda, d.vec = theta.2.tilda)  + t(beta.2[,s])%*%inv.Lamda.star.2%*%beta.2[,s]) )
    
    # MH ratio
    log.A = dmvnorm(x = theta.2[,s], mean = alpha.2[s]*c(rep(1,N)) + X%*%beta.2[,s], sigma = sigma.sq.2.new*diag(N) + gamma.sq.2[s]*B0.2, log = TRUE)
    log.B = dmvnorm(x = beta.2[,s], mean = c(rep(0,p)), sigma = sigma.sq.2.new*Lambda.2, log = TRUE)
    log.E = dmvnorm(x = theta.2[,s], mean = alpha.2[s]*c(rep(1,N)) + X%*%beta.2[,s], sigma = sigma.sq.2[s]*diag(N) + gamma.sq.2[s]*B0.2, log = TRUE)
    log.F = dmvnorm(x = beta.2[,s], mean = c(rep(0,p)), sigma = sigma.sq.2[s]*Lambda.2, log = TRUE)
    
    log.C.over.G = log.C.over.G = log(sigma.sq.2[s]/sigma.sq.2.new)
    log.D.over.H = log.ratio.dIG.common.shape(nomi.x = sigma.sq.2[s], deno.x = sigma.sq.2.new, common.alpha = (N+p)/2 +2 -2
                                              , nomi.beta = (1/2)*(quad.inverse.B.tilda.vec(a = 1, b = gamma.sq.2[s]/sigma.sq.2.new, B.mat = B0.2, c.vec = theta.2.tilda, d.vec = theta.2.tilda) + t(beta.2[,s])%*%inv.Lamda.star.2%*%beta.2[,s]) 
                                              , deno.beta = (1/2)*(quad.inverse.B.tilda.vec(a = 1, b = gamma.sq.2[s]/sigma.sq.2[s], B.mat = B0.2, c.vec = theta.2.tilda, d.vec = theta.2.tilda) + t(beta.2[,s])%*%inv.Lamda.star.2%*%beta.2[,s]) )
    
    R = log.A + log.B - (log.E + log.F) + log.C.over.G + log.D.over.H
    
    r = exp(R)
    
    # criterion
    rat = min(r,1)
    U = runif(n=1, min = 0, max = 1)
    
    if (U < rat){
      sigma.sq.2[(s+1)] = sigma.sq.2.new
      
    } else {
      sigma.sq.2[(s+1)] = sigma.sq.2[s]
      
    }
    
    
    # shape (3)
    inv.Lamda.star.3 = (1/tau.sq.3[s+1])*inv.Lambda.3
    theta.3.tilda = theta.3[,s] - alpha.3[s]*c(rep(1,N)) - X%*%beta.3[,s]
    c.3 = gamma.sq.3[s]/sigma.sq.3[s]
    
    # proposal (p.166 - 4 HBM I)
    sigma.sq.3.new = 1/rgamma(n = 1, shape = (N + p)/2 + 2 -2, 
                              rate = (1/2)*(quad.inverse.B.tilda.vec(a = 1, b = c.3, B.mat = B0.3, c.vec = theta.3.tilda, d.vec = theta.3.tilda)  + t(beta.3[,s])%*%inv.Lamda.star.3%*%beta.3[,s]) )
    
    
    # MH ratio
    log.A = dmvnorm(x = theta.3[,s], mean = alpha.3[s]*c(rep(1,N)) + X%*%beta.3[,s], sigma = sigma.sq.3.new*diag(N) + gamma.sq.3[s]*B0.3, log = TRUE)
    log.B = dmvnorm(x = beta.3[,s], mean = c(rep(0,p)), sigma = sigma.sq.3.new*Lambda.3, log = TRUE)
    log.E = dmvnorm(x = theta.3[,s], mean = alpha.3[s]*c(rep(1,N)) + X%*%beta.3[,s], sigma = sigma.sq.3[s]*diag(N) + gamma.sq.3[s]*B0.3, log = TRUE)
    log.F = dmvnorm(x = beta.3[,s], mean = c(rep(0,p)), sigma = sigma.sq.3[s]*Lambda.3, log = TRUE)
    
    log.C.over.G = log.C.over.G = log(sigma.sq.3[s]/sigma.sq.3.new)
    log.D.over.H = log.ratio.dIG.common.shape(nomi.x = sigma.sq.3[s], deno.x = sigma.sq.3.new, common.alpha = (N+p)/2 +2 -2
                                              , nomi.beta = (1/2)*(quad.inverse.B.tilda.vec(a = 1, b = gamma.sq.3[s]/sigma.sq.3.new, B.mat = B0.3, c.vec = theta.3.tilda, d.vec = theta.3.tilda) + t(beta.3[,s])%*%inv.Lamda.star.3%*%beta.3[,s]) 
                                              , deno.beta = (1/2)*(quad.inverse.B.tilda.vec(a = 1, b = gamma.sq.3[s]/sigma.sq.3[s], B.mat = B0.3, c.vec = theta.3.tilda, d.vec = theta.3.tilda) + t(beta.3[,s])%*%inv.Lamda.star.3%*%beta.3[,s]))
    
    R = log.A + log.B - (log.E + log.F) + log.C.over.G + log.D.over.H
    
    r = exp(R)
    
    # criterion
    rat = min(r,1)
    U = runif(n=1, min = 0, max = 1)
    
    if (U < rat){
      sigma.sq.3[(s+1)] = sigma.sq.3.new
      
    } else {
      sigma.sq.3[(s+1)] = sigma.sq.3[s]
    }
    
    }
    
    # Step 6.
    {
      Q.beta.1 = quad.inverse.B.tilda.matrix(a = sigma.sq.1[s+1], b = gamma.sq.1[s], B.mat = B0.1, X = X, Y = X)  + (1/sigma.sq.1[s+1])*inv.Lamda.star.1
      b.beta.1 = quad.inverse.B.tilda.matrix(a = sigma.sq.1[s+1], b = gamma.sq.1[s], B.mat = B0.1, X = X, Y = c(theta.1[,s] - alpha.1[s]*c(rep(1,N))) )
      beta.1[,(s+1)] = rmvt.Rue.sample(Q = Q.beta.1, b = b.beta.1)
      
      Q.beta.2 = quad.inverse.B.tilda.matrix(a = sigma.sq.2[s+1], b = gamma.sq.2[s], B.mat = B0.2, X = X, Y = X)  + (1/sigma.sq.2[s+1])*inv.Lamda.star.2
      b.beta.2 = quad.inverse.B.tilda.matrix(a = sigma.sq.2[s+1], b = gamma.sq.2[s], B.mat = B0.2, X = X, Y = c(theta.2[,s] - alpha.2[s]*c(rep(1,N))) )
      beta.2[,(s+1)] = rmvt.Rue.sample(Q = Q.beta.2, b = b.beta.2)
      
      Q.beta.3 = quad.inverse.B.tilda.matrix(a = sigma.sq.3[s+1], b = gamma.sq.3[s], B.mat = B0.3, X = X, Y = X)  + (1/sigma.sq.3[s+1])*inv.Lamda.star.3
      b.beta.3 = quad.inverse.B.tilda.matrix(a = sigma.sq.3[s+1], b = gamma.sq.3[s], B.mat = B0.3, X = X, Y = c(theta.3[,s] - alpha.3[s]*c(rep(1,N))) )
      beta.3[,(s+1)] = rmvt.Rue.sample(Q = Q.beta.3, b = b.beta.3)
    }
    
    # Step 7.
    {
      deno.1 = quad.inverse.B.tilda.vec(a = sigma.sq.1[s+1], b = gamma.sq.1[s], B.mat = B0.1, c.vec = rep(1,N), d.vec = rep(1,N))
      nomi.1 = quad.inverse.B.tilda.vec(a = sigma.sq.1[s+1], b = gamma.sq.1[s], B.mat = B0.1, c.vec = rep(1,N), d.vec =  theta.1[,s] - X%*%beta.1[,(s+1)])
      alpha.1[s+1] = rnorm(n = 1, mean = nomi.1/deno.1, sd = sqrt(1/deno.1))
    
      deno.2 = quad.inverse.B.tilda.vec(a = sigma.sq.2[s+1], b = gamma.sq.2[s], B.mat = B0.2, c.vec = rep(1,N), d.vec = rep(1,N))
      nomi.2 = quad.inverse.B.tilda.vec(a = sigma.sq.2[s+1], b = gamma.sq.2[s], B.mat = B0.2, c.vec = rep(1,N), d.vec =  theta.2[,s] - X%*%beta.2[,(s+1)])
      alpha.2[s+1] = rnorm(n = 1, mean = nomi.2/deno.2, sd = sqrt(1/deno.2))
      
      deno.3 = quad.inverse.B.tilda.vec(a = sigma.sq.3[s+1], b = gamma.sq.3[s], B.mat = B0.3, c.vec = rep(1,N), d.vec = rep(1,N))
      nomi.3 = quad.inverse.B.tilda.vec(a = sigma.sq.3[s+1], b = gamma.sq.3[s], B.mat = B0.3, c.vec = rep(1,N), d.vec =  theta.3[,s] - X%*%beta.3[,(s+1)])
      alpha.3[s+1] = rnorm(n = 1, mean = nomi.3/deno.3, sd = sqrt(1/deno.3))
    }
    
    # Step 8.
    {
    r.elements = c()
    for (i in 1:N){
      r.elements[i] = sum(y(i) - h.i(i = i, theta.2.i = theta.2[i,s], theta.3.i = theta.3[i,s]))
    }
    r.vec = c(r.elements)
    
    inverse.matrix.B0.1.tilda = inverse.B.tilda(a = sigma.sq.1[s+1], b = gamma.sq.1[s], B.mat = B0.1)
    
    Q.theta.1 = (1/sigma.sq[s])*T.mat + inverse.matrix.B0.1.tilda
    b.theta.1 = (1/sigma.sq[s])*r.vec + inverse.matrix.B0.1.tilda%*%( X%*%beta.1[,(s+1)] + rep(1,N)*alpha.1[s+1] )
    
    theta.1[,(s+1)] = rmvt.Rue.sample(Q = Q.theta.1, b = b.theta.1)
    }
    
    # Step 9.
    {
    for (i in 1:N){
      
      # Goal: Sampling theta.2[i,(s+1)]
      
      # proposing a new value
      theta.2.i.new = rnorm(n=1, mean = theta.2[i,s], sd = sqrt(prop.var.theta.2))
      
      # MH ratio
      r = exp(-(1/(2*sigma.sq[s]))*(norm_vec(y(i), mu.i(i = i, theta.1.i = theta.1[i,(s+1)],theta.2.i = theta.2.i.new,theta.3.i = theta.3[i,s])))^2
              -(1/(2*( sigma.sq.2[s+1] + gamma.sq.2[s] ) ))*(theta.2.i.new - alpha.2[s+1] - t(X[i,])%*%beta.2[,(s+1)])^2
              
              +(1/(2*sigma.sq[s]))*(norm_vec(y(i), mu.i(i = i, theta.1.i = theta.1[i,(s+1)],theta.2.i = theta.2[i,s],theta.3.i = theta.3[i,s])))^2
              +(1/(2*( sigma.sq.2[s+1] + gamma.sq.2[s] ) ))*(theta.2[i,s] - alpha.2[s+1] - t(X[i,])%*%beta.2[,(s+1)])^2
      )
      
      rat = min(r,1)
      
      # Criterion
      U = runif(n=1, min = 0, max = 1)
      
      if (U < rat){
        theta.2[i,(s+1)] = theta.2.i.new
      } else {
        theta.2[i,(s+1)] = theta.2[i,s]
      }
      
      
      # Goal: Sampling theta.3[i,(s+1)]
      
      # Propose a new value
      
      theta.3.i.new = rnorm(n=1, mean = theta.3[i,s], sd = sqrt(prop.var.theta.3))
      
      # MH ratio
      r = exp(-(1/(2*sigma.sq[s]))*(norm_vec(y(i), mu.i(i = i, theta.1.i = theta.1[i,(s+1)],theta.2.i = theta.2[i,(s+1)],theta.3.i = theta.3.i.new)))^2
              -(1/(2*(sigma.sq.3[s+1] + gamma.sq.3[s]) ))*(theta.3.i.new - alpha.3[s+1] - t(X[i,])%*%beta.3[,(s+1)])^2
              
              +(1/(2*sigma.sq[s]))*(norm_vec(y(i), mu.i(i = i, theta.1.i = theta.1[i,(s+1)],theta.2.i = theta.2[i,(s+1)],theta.3.i = theta.3[i,s])))^2
              +(1/(2*(sigma.sq.3[s+1] + gamma.sq.3[s]) ))*(theta.3[i,s] - alpha.3[s+1] - t(X[i,])%*%beta.3[,(s+1)])^2
      )
      
      rat = min(r,1)
      
      # Criterion
      U = runif(n=1, min = 0, max = 1)
      
      if (U < rat){
        theta.3[i,(s+1)] = theta.3.i.new
        
      } else {
        theta.3[i,(s+1)] = theta.3[i,s]
      }
    }
    }
    
    # Step 10.
    {
      phi[s+1] = 1/rgamma(n = 1, shape = 1, rate = 1 + 1/sigma.sq[s])
    }
    
    # Step 11.
    {
      tempo.vec = c()
      for (i in 1:N){
        tempo.vec[i] = norm_vec(y(i), mu.i(i=i, theta.1.i = theta.1[i,(s+1)],theta.2.i = theta.2[i,(s+1)], theta.3.i = theta.3[i,(s+1)]))^2
      }
      sigma.sq[s+1] = 1/rgamma(n = 1, shape = (sum(T.mat)+1)/2, rate = 1/phi[s+1] + (1/2)*sum(tempo.vec))
    }
    
    # Step 12.
    {
      # magnitude (1)
      theta.1.tilda = theta.1[,(s+1)] - alpha.1[s+1]*c(rep(1,N)) - X%*%beta.1[,(s+1)]
      c.1 = gamma.sq.1[s]/sigma.sq.1[s+1]
      
      # proposal 
      gamma.sq.1.new = 1/rgamma(n = 1, shape = (N+1)/2 , rate = (1/2)*quad.inverse.B.tilda.vec(a = 1/c.1, b = 1, B.mat = B0.1, c.vec = theta.1.tilda, d.vec = theta.1.tilda) + 1/omega.1)
      
      # MH ratio
      log.A = dmvnorm(x = theta.1[,(s+1)], mean = alpha.1[s+1]*c(rep(1,N)) + X%*%beta.1[,(s+1)], sigma = sigma.sq.1[s+1]*diag(N) + gamma.sq.1.new*B0.1, log = TRUE)
      log.D = dmvnorm(x = theta.1[,(s+1)], mean = alpha.1[s+1]*c(rep(1,N)) + X%*%beta.1[,(s+1)], sigma = sigma.sq.1[s+1]*diag(N) + gamma.sq.1[s]*B0.1, log = TRUE)
      
      
      log.B.over.E = log.ratio.dIG.common.shape.scale(nomi.x = gamma.sq.1.new, deno.x = gamma.sq.1[s],common.alpha = 1/2, common.beta = 1/omega.1[s]   )   
      log.C.over.F = log.ratio.dIG.common.shape(nomi.x = gamma.sq.1[s], deno.x = gamma.sq.1.new, common.alpha = (N+1)/2, nomi.beta = (1/2)*quad.inverse.B.tilda.vec(a = sigma.sq.1[s+1]/gamma.sq.1.new, b = 1, B.mat = B0.1, c.vec = theta.1.tilda, d.vec = theta.1.tilda) + 1/omega.1[s], deno.beta = (1/2)*quad.inverse.B.tilda.vec(a = sigma.sq.1[s+1]/gamma.sq.1[s], b = 1, B.mat = B0.1, c.vec = theta.1.tilda, d.vec = theta.1.tilda)  + 1/omega.1[s])
      
      R = log.A  - log.D + log.B.over.E + log.C.over.F
      
      r = exp(R)
      
      # criterion
      rat = min(r,1)
      U = runif(n=1, min = 0, max = 1)
      
      if (U < rat){
        gamma.sq.1[(s+1)] = gamma.sq.1.new
        
      } else {
        gamma.sq.1[(s+1)] = gamma.sq.1[s]
        
      }
      
      
      # scale (2)
      theta.2.tilda = theta.2[,(s+1)] - alpha.2[s+1]*c(rep(1,N)) - X%*%beta.2[,(s+1)]
      c.2 = gamma.sq.2[s]/sigma.sq.2[s+1]
      
      # proposal 
      gamma.sq.2.new = 1/rgamma(n = 1, shape = (N+1)/2 , rate = (1/2)*quad.inverse.B.tilda.vec(a = 1/c.2, b = 1, B.mat = B0.2, c.vec = theta.2.tilda, d.vec = theta.2.tilda) + 1/omega.2)
      
      # MH ratio
      log.A = dmvnorm(x = theta.2[,(s+1)], mean = alpha.2[s+1]*c(rep(1,N)) + X%*%beta.2[,(s+1)], sigma = sigma.sq.2[s+1]*diag(N) + gamma.sq.2.new*B0.2, log = TRUE)
      log.D = dmvnorm(x = theta.2[,(s+1)], mean = alpha.2[s+1]*c(rep(1,N)) + X%*%beta.2[,(s+1)], sigma = sigma.sq.2[s+1]*diag(N) + gamma.sq.2[s]*B0.2, log = TRUE)
      
      
      log.B.over.E = log.ratio.dIG.common.shape.scale(nomi.x = gamma.sq.2.new, deno.x = gamma.sq.2[s],common.alpha = 1/2, common.beta = 1/omega.2[s]   )   
      log.C.over.F = log.ratio.dIG.common.shape(nomi.x = gamma.sq.2[s], deno.x = gamma.sq.2.new, common.alpha = (N+1)/2, nomi.beta = (1/2)*quad.inverse.B.tilda.vec(a = sigma.sq.2[s+1]/gamma.sq.2.new, b = 1, B.mat = B0.2, c.vec = theta.2.tilda, d.vec = theta.2.tilda) + 1/omega.2[s], deno.beta = (1/2)*quad.inverse.B.tilda.vec(a = sigma.sq.2[s+1]/gamma.sq.2[s], b = 1, B.mat = B0.2, c.vec = theta.2.tilda, d.vec = theta.2.tilda)  + 1/omega.2[s])
      
      R = log.A  - log.D + log.B.over.E + log.C.over.F
      
      r = exp(R)
      
      # criterion
      rat = min(r,1)
      U = runif(n=1, min = 0, max = 1)
      
      if (U < rat){
        gamma.sq.2[(s+1)] = gamma.sq.2.new
        
      } else {
        gamma.sq.2[(s+1)] = gamma.sq.2[s]
        
      }
      
      
      # shape (3)
      theta.3.tilda = theta.3[,(s+1)] - alpha.3[s+1]*c(rep(1,N)) - X%*%beta.3[,(s+1)]
      c.3 = gamma.sq.3[s]/sigma.sq.3[s+1]
      
      # proposal 
      gamma.sq.3.new = 1/rgamma(n = 1, shape = (N+1)/2 , rate = (1/2)*quad.inverse.B.tilda.vec(a = 1/c.3, b = 1, B.mat = B0.3, c.vec = theta.3.tilda, d.vec = theta.3.tilda) + 1/omega.3)
      
      # MH ratio
      log.A = dmvnorm(x = theta.3[,(s+1)], mean = alpha.3[s+1]*c(rep(1,N)) + X%*%beta.3[,(s+1)], sigma = sigma.sq.3[s+1]*diag(N) + gamma.sq.3.new*B0.3, log = TRUE)
      log.D = dmvnorm(x = theta.3[,(s+1)], mean = alpha.3[s+1]*c(rep(1,N)) + X%*%beta.3[,(s+1)], sigma = sigma.sq.3[s+1]*diag(N) + gamma.sq.3[s]*B0.3, log = TRUE)
      
      
      log.B.over.E = log.ratio.dIG.common.shape.scale(nomi.x = gamma.sq.3.new, deno.x = gamma.sq.3[s],common.alpha = 1/2, common.beta = 1/omega.3[s]   )   
      log.C.over.F = log.ratio.dIG.common.shape(nomi.x = gamma.sq.3[s], deno.x = gamma.sq.3.new, common.alpha = (N+1)/2, nomi.beta = (1/2)*quad.inverse.B.tilda.vec(a = sigma.sq.3[s+1]/gamma.sq.3.new, b = 1, B.mat = B0.3, c.vec = theta.3.tilda, d.vec = theta.3.tilda) + 1/omega.3[s], deno.beta = (1/2)*quad.inverse.B.tilda.vec(a = sigma.sq.3[s+1]/gamma.sq.3[s], b = 1, B.mat = B0.3, c.vec = theta.3.tilda, d.vec = theta.3.tilda)  + 1/omega.3[s])
      
      R = log.A  - log.D + log.B.over.E + log.C.over.F
      
      r = exp(R)
      
      # criterion
      rat = min(r,1)
      U = runif(n=1, min = 0, max = 1)
      
      if (U < rat){
        gamma.sq.3[(s+1)] = gamma.sq.3.new
        
      } else {
        gamma.sq.3[(s+1)] = gamma.sq.3[s]
        
      }
    }
    
    # Step 13.
    {
      omega.1[s+1] = 1/rgamma(n = 1, shape = 1, rate = 1 + 1/gamma.sq.1[s+1] )
      omega.2[s+1] = 1/rgamma(n = 1, shape = 1, rate = 1 + 1/gamma.sq.2[s+1] )
      omega.3[s+1] = 1/rgamma(n = 1, shape = 1, rate = 1 + 1/gamma.sq.3[s+1] )
    }
    
    print(s)
    
  }
  
  # Print-out Results
  {
  mc.index = seq(from = burn + 1, to = burn + nmc, by = thin)
  res = list(thinned.theta.1 = theta.1[,mc.index],
             thinned.theta.2 = theta.2[,mc.index],
             thinned.theta.3 = theta.3[,mc.index],
             thinned.sigma.sq = sigma.sq[mc.index],
             thinned.alpha.1 = alpha.1[mc.index],
             thinned.alpha.2 = alpha.2[mc.index],
             thinned.alpha.3 = alpha.3[mc.index],
             thinned.beta.1 = beta.1[,mc.index],
             thinned.beta.2 = beta.2[,mc.index],
             thinned.beta.3 = beta.3[,mc.index],
             thinned.sigma.sq.1 = sigma.sq.1[,mc.index],
             thinned.sigma.sq.2 = sigma.sq.2[,mc.index],
             thinned.sigma.sq.3 = sigma.sq.3[,mc.index],
             thinned.gamma.sq.1 = gamma.sq.1[,mc.index],
             thinned.gamma.sq.2 = gamma.sq.2[,mc.index],
             thinned.gamma.sq.3 = gamma.sq.3[,mc.index])
  }
  return(res)
}

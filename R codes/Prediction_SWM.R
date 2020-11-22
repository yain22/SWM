# Spatial Weibull Model (SWM) 
# Update: 11/17/2020
# @ 2020 Se Yoon Lee and Bani Mallick All Rights Reserved
# User should contact seyoonlee@stat.tamu.edu and bmallick@stat.tamu.edu 
# for any use and modification for the code for the publication purpose or industrial use

Prediction_SWM = function(index.of.test.well, log_scale = TRUE){
  
  #1. Select index for a single testing well (i.t)
  i.t = index.of.test.well # New well index for Prediction test; the index should not chosen from training well index
  #2. Formulate full location set corresponding to the (1) individual testing well and (2) training wells
  full.Loc = Production_Results_Eagle_Ford$Loc[c(i.t, 1:N),]
  full.dist.mat <- rdist(full.Loc) 
  
  #3. Full correlation matrix
  # Gaussian Correlation Function
  b0 = function(rho, distance){
    g = exp(-exp(rho)*distance^2)
    return(g)
  }
  dist.mat <- rdist(Production_Results_Eagle_Ford$Loc[1:N,]) 
  B0.1 = b0(rho = rho.1, distance = dist.mat)
  B0.2 = b0(rho = rho.2, distance = dist.mat)
  B0.3 = b0(rho = rho.3, distance = dist.mat)
  full.B0.1 = b0(rho = rho.1, distance = full.dist.mat) 
  full.B0.2 = b0(rho = rho.2, distance = full.dist.mat) 
  full.B0.3 = b0(rho = rho.3, distance = full.dist.mat) 
  
  #4. N-dimenstional correlation vector corresponding to (one test well, N training wells)
  b.1 = full.B0.1[c(2:(N+1) ),1]
  b.2 = full.B0.2[c(2:(N+1) ),1]
  b.3 = full.B0.3[c(2:(N+1) ),1]
  
  #5. Length of production period for test well 
  T.p = length(y(i.t))
  #6. Covariate data for testweill i.t
  x.p.trans = matrix(
    scale(Production_Results_Eagle_Ford$X, center = TRUE, scale = TRUE)[i.t,],
    nrow =1, ncol = 11) 
  
  # Prediction Sampling Setup
  # Number of realized posterior samples
  S = ncol(theta.1) 
  
  # Make room for Y_p Tp by S , T.p = length(y(i.t))
  y.p = matrix(rep(0,T.p*S), nrow = T.p, ncol = S) # room for prediction y.p 
  theta.1.p = rep(0,S)
  theta.2.p = rep(0,S)
  theta.3.p = rep(0,S)
  
  set.seed(2) # For reproducibility
  for (s in 1:(S-1)){
    
    if (s == 1){
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
      
      q.i = function(i, M.i, b.i, k.i){
        q.i.vector = c(rep(0,as.numeric(T.mat[i,i]))) # Room for vector
        
        for (t in 1:as.numeric(T.mat[i,i])){
          q.i.vector[t] = Weibull.model.i.fuction(t = t, M.i = M.i, b.i = b.i, k.i = k.i)
        }
        
        return(q.i.vector)  
      }
      # Log-scaled function 
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
      
      # For prediction prediction 
      mu.p = function(prod.period , theta.1.p, theta.2.p, theta.3.p){
        
        mu.i.vector = c(rep(0,as.numeric(prod.period))) # room for vector
        
        for (t in 1:as.numeric(prod.period)){
          
          mu.i.vector[t] = mu.it(t = t, theta.1.i = theta.1.p, theta.2.i = theta.2.p, theta.3.i = theta.3.p)
          
        }
        
        return(mu.i.vector)
        
      }
      
    }
    
    # Note: 
    # 1. Each prediction equation is same with the paper
    # 2. We tried to match notations in the paper to avoid confusion
    # 3. x.p.trans & b.1 involve information from the test well
    
    # Step 1. Obtain posterior sample from spatial Weibull model for the training wells.
    
    # Step 2. For each s = 1, ..., S, sample theta.1.p, theta.2.p, theta.3.p independently
    # Maginitude (1)
    theta.1.tilda = c(theta.1[,s] - alpha.1[s]*c(rep(1,N)) - X%*%beta.1[,s])
    mu.1. =  (alpha.1[s] + x.p.trans%*%beta.1[,s]) + (gamma.sq.1[s])*quad.inverse.B.tilda.vec(a = sigma.sq.1[s], b = gamma.sq.1[s], B.mat = B0.1, c.vec = b.1, d.vec = theta.1.tilda)
    sigma.sq.1. =  (sigma.sq.1[s] + gamma.sq.1[s]) - ((gamma.sq.1[s])^2)*quad.inverse.B.tilda.vec(a = sigma.sq.1[s], b = gamma.sq.1[s], B.mat = B0.1, c.vec = b.1, d.vec = b.1)
    theta.1.p[s] = rnorm(n = 1,mean = mu.1., sd = sqrt(sigma.sq.1.))
    
    # Scale (2)
    theta.2.tilda = c(theta.2[,s] - alpha.2[s]*c(rep(1,N)) - X%*%beta.2[,s])
    mu.2. =  (alpha.2[s] + x.p.trans%*%beta.2[,s]) + (gamma.sq.2[s])*quad.inverse.B.tilda.vec(a = sigma.sq.2[s], b = gamma.sq.2[s], B.mat = B0.2, c.vec = b.2, d.vec = theta.2.tilda)
    sigma.sq.2. =  (sigma.sq.2[s] + gamma.sq.2[s]) - ((gamma.sq.2[s])^2)*quad.inverse.B.tilda.vec(a = sigma.sq.2[s], b = gamma.sq.2[s], B.mat = B0.2, c.vec = b.2, d.vec = b.2)
    theta.2.p[s] = rnorm(n = 1,mean = mu.2., sd = sqrt(sigma.sq.2.))
    
    # Shape (3)
    theta.3.tilda = c(theta.3[,s] - alpha.3[s]*c(rep(1,N)) - X%*%beta.3[,s])
    mu.3. =  (alpha.3[s] + x.p.trans%*%beta.3[,s]) + (gamma.sq.3[s])*quad.inverse.B.tilda.vec(a = sigma.sq.3[s], b = gamma.sq.3[s], B.mat = B0.3, c.vec = b.3, d.vec = theta.3.tilda)
    sigma.sq.3. =  (sigma.sq.3[s] + gamma.sq.3[s]) - ((gamma.sq.3[s])^2)*quad.inverse.B.tilda.vec(a = sigma.sq.3[s], b = gamma.sq.3[s], B.mat = B0.3, c.vec = b.3, d.vec = b.3)
    theta.3.p[s] = rnorm(n = 1,mean = mu.3., sd = sqrt(sigma.sq.3.))
    
    
    # Step 3. For each s = 1, ..., S, sample y.p 
    
    y.p[ ,s] = rmvnorm(  n = 1, mean = mu.p(prod.period = T.p, theta.1.p = theta.1.p[s], theta.2.p = theta.2.p[s], theta.3.p = theta.3.p[s]), sigma = sigma.sq[s]*diag(T.p)) 
     
    if ( s%%500 == 0 ){
      print(s)
    }
    
  }
  
  # Plotting
  
  #1. Plotting in original-scale
  if (log_scale == FALSE){
    
    # Turning back to original scale
    q.p = exp(y.p) # T.p by S dimensional vector
    
    post.mean.q.p.given.y_1.N = rowMeans(q.p) 
    post.median.q.p.given.y_1.N = apply(X = q.p,MARGIN = 1,FUN = median) 
    cred.interval.q.p.given.y_1.N = matrix(rep(0,T.p*2),nrow = T.p, ncol = 2)
    lower.prop = 0.05
    upper.prop = 0.95
    for (t in 1:T.p){
      cred.interval.q.p.given.y_1.N[t , ] = quantile(q.p[t,], probs = c(lower.prop, upper.prop))  
    }
    
    plot(post.median.q.p.given.y_1.N, ylab = "Production (BBL/Month)", xlab = "Month",type = "l", ylim = c(min(cred.interval.q.p.given.y_1.N), max(cred.interval.q.p.given.y_1.N)), lwd =2, col = "red")
    polygon(x = c(c(1:T.p),rev(c(1:T.p))), y = c(cred.interval.q.p.given.y_1.N[,1],rev(cred.interval.q.p.given.y_1.N[,2])),border = FALSE, col = "grey")
    points(post.median.q.p.given.y_1.N, ylab = "Production (BBL/Month)", xlab = "Month",type = "l", ylim = c(min(cred.interval.q.p.given.y_1.N), max(cred.interval.q.p.given.y_1.N)), lwd =2, col = "red")
    lines(x = c(1:T.p), y = exp(y(i.t)), type ="p", col ="black",pch = 20)
    legend("topright", legend = c( "Observations","Posterior predictive median ", "95% credible interval")
           ,lty = c(0,1, 1), pch = c(20,-1,-1), col = c("black","red","grey"), lwd = c(1,2,7))
    title(paste("API 10 =", Production_Results_Eagle_Ford$API10[i.t,]))
    
    #res = list(post.median.q.p.given.y_1.N)
    #names(res) = c("Post.median.q.p (original)")
  }
  
  
  #2. Plotting in log-scale
  if (log_scale == TRUE){
    post.mean.y.p.given.y_1.N = rowMeans(y.p) 
    post.median.y.p.given.y_1.N = apply(X = y.p,MARGIN = 1,FUN = median) 
    cred.interval.y.p.given.y_1.N = matrix(rep(0,T.p*2),nrow = T.p, ncol = 2)
    lower.prop = 0.05
    upper.prop = 0.95
    for (t in 1:T.p){
      cred.interval.y.p.given.y_1.N[t , ] = quantile(y.p[t,], probs = c(lower.prop, upper.prop))  
    }
    
    y.margin.log.scaled = 1
    plot(post.median.y.p.given.y_1.N, ylab = "Log-production (BBL/Month)", xlab = "Month",type = "l", ylim = c(min(cred.interval.y.p.given.y_1.N), max(cred.interval.y.p.given.y_1.N) + y.margin.log.scaled ), lwd =2, col = "red")
    polygon(x = c(c(1:T.p),rev(c(1:T.p))), y = c(cred.interval.y.p.given.y_1.N[,1],rev(cred.interval.y.p.given.y_1.N[,2])),border = FALSE, col = "grey")
    points(post.median.y.p.given.y_1.N, ylab = "Log-production (BBL/Month)", xlab = "Month",type = "l", ylim = c(min(cred.interval.y.p.given.y_1.N), max(cred.interval.y.p.given.y_1.N) + y.margin.log.scaled ), lwd =2, col = "red")
    lines(x = c(1:T.p), y = y(i.t), type ="p", col ="black", pch = 20)
    legend("topright", legend = c( "Observations","Posterior predictive median ", "95% credible interval")
           ,lty = c(0,1, 1), pch = c(20,-1,-1), col = c("black","red","grey"), lwd = c(1,2,7))
    title(paste("API 10 =", Production_Results_Eagle_Ford$API10[i.t,]))
    
    #res = list(post.median.y.p.given.y_1.N)
    #names(res) = c("Post.median.y.p (log-scaled)")
  }
  
  #return(res)
  
}


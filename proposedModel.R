library(nimble)
  library(coda)
  library(MASS)
  library(sp)
  library(fossil)
  library(purrr)
  library(parallel)
  library(foreach)
  library(doParallel)
  library(DescTools)
  set.seed(1000)
  load("DATA.Rdata")
  
  load("china_adj.RData")
  colnames <- china_adj$chinapoly[[6]]
  
  adj_matrix <- china_adj$adj[-c(9,28,29),-c(9,28,29)]
  str(adj_matrix)
  rowSums(adj_matrix)
  as.vector(colSums(adj_matrix))
  
  library(Matrix)
  row_adj <- Matrix::Diagonal(x = 1 / sqrt(Matrix::rowSums(adj_matrix^2))) %*% adj_matrix
  
  unitind <- function(x) {
    #return(which(x == "1"))
    return(which(x != 0))
  }
  
  adj1 <- as.vector(apply(row_adj, 1, unitind))
  adjvector <- c()
  adjnum <- c()
  for (i in 1:28) {
    adjvector <- append(adjvector, as.vector(adj1[[i]]))
    adjnum[i] <- length(as.vector(adj1[[i]]))
  }
  
  eigenadj <- eigen(adj_matrix)$values
  maxeigenadj <- max(eigenadj); maxeigenadj
  mineigenadj <- min(eigenadj); mineigenadj

##### proposed model ######
NigCode <- nimbleCode({
  for (i in 1:N) {
    y[i] ~ dcat(py[i,1:5])
    
    probit(Q[i,1]) <- eta[i,1] - (beta[1] + beta[2]*z1[i] + beta[3]*z2[i] + beta[4]*z3[i] + beta[5]*z4[i] +
                                    beta[6]*z5[i] + betax[prov[i]]*x[i] + Wy[prov[i]])
    py[i,1] <- Q[i,1]
    for (j in 2:4) {
      probit(Q[i,j]) <- eta[i,j] - (beta[1] + beta[2]*z1[i] + beta[3]*z2[i] + beta[4]*z3[i] + beta[5]*z4[i] +
                                      beta[6]*z5[i] + betax[prov[i]]*x[i]+ Wy[prov[i]])
      py[i,j] <- Q[i,j] - Q[i, j-1]
    }
    py[i,5] <- 1 - Q[i,4]
    
    eta[i,1] <- 0
    for (j in 2:4) {  
      eta[i,j] <- eta[i,j-1] + exp(phi[j,1] + phi[j,2] * z1[i] + phi[j,3] * z2[i] + phi[j,4] * z3[i] + 
                                     phi[j,5] * z4[i] + phi[j,6] * z5[i])
    }
    
    
    
    x[i] ~ dnorm(mu_x[i], tau_x)
    mu_x[i] <- alpha[1] + alpha[2] * z1[i] + alpha[3] * z2[i] + alpha[4] * z3[i] + 
      alpha[5] * z4[i] + alpha[6] * z5[i] 
    
    
    R[i] ~ dbern(mu_R[i])
    logit(mu_R[i]) <- mu_R0[i]
    xm[i] <- (x[i] > medcri[prov[i]]) + 0
    mu_R0[i] <- mu_R1[i] + gammax0[prov[i]] * xm[i] 
    mu_R1[i] <- gamma[1] + gamma[2] * z1[i] + gamma[3] * z2[i]  + gamma[4] * z3[i] +
      gamma[5] * z4[i] + gamma[6] * z5[i] + gamma[7] * ((y[i]==1)+0) + 
      gamma[8] * ((y[i]==2)+0) + gamma[9] * ((y[i]==4)+0) + gamma[10] * ((y[i]==5)+0) + Wr[prov[i]]
    uDevR[i] <- -2 * (R[i] * log(mu_R[i]) + (1 - R[i]) * log(1 - mu_R[i]))
  }
  
  for (j in 1:4) {
    for (k in 1:6) {
      phi[j,k] ~ dnorm(0, 0.1)
    }
  }
  
  Wy[1:28] ~ dmnorm(mean = mu_wy[1:28], prec = L_wy[1:28, 1:28])
  L_wy[1:28, 1:28] <- (I[1:28, 1:28] - gamma_wy * A[1:28, 1:28]) * 1/sigma_wy^2
  gamma_wy ~ dunif(mine, maxe)
  
  Wr[1:28] ~ dmnorm(mean = mu_wr[1:28], prec = L_wr[1:28, 1:28])
  L_wr[1:28, 1:28] <- (I[1:28, 1:28] - gamma_wr * A[1:28, 1:28]) * 1/sigma_wr^2
  gamma_wr ~ dunif(mine, maxe)
  
  for (j in 1:6) {
    alpha[j] ~ dnorm(0, 0.01)
  }
  for (j in 1:6) {
    beta[j] ~ dnorm(0, 0.01)
  }
  for (j in 1:10) {
    gamma[j] ~ dnorm(0, 0.01)
  }
  tau_x  ~ dgamma(0.1, 0.1)
  for (j in 1:28) {
    medcri[j] <- medx[j] * xm_s
    betax[j] ~ dnorm(0, 0.1)
    gammax0[j] <- delta_g1[j] * (clu[j] == 1) + delta_g2[j] * (clu[j] == 2) + delta_g3[j] * (clu[j] == 3)
    delta_g1[j] ~ T(dnorm(mu_g1, tau_g1),,0)
    delta_g3[j] ~ T(dnorm(mu_g3, tau_g3),0,)
    clu[j] ~ dcat(pclu[1:3])
  }
  
  pclu[1:3] ~ ddirch(alpha_pclu[1:3])
  
  mu_g1 ~ T(dnorm(0, tau_mu1),,0)
  tau_mu1 <- tau_g1 * 0.1  
  mu_g3 ~ T(dnorm(0, tau_mu3),0,)
  tau_mu3 <- tau_g3 * 0.1  
  
  tau_g1 ~ dgamma(0.1,0.1)
  tau_g3 ~ dgamma(0.1,0.1)
  
  sigma_wy ~ dgamma(1, 1)
  sigma_wr ~ dgamma(1, 1)
  
  xm_s ~ T(dnorm(1, 0.1),0,)
})


R <- (is.na(DATA[[k]]$x)+0)   # kth replicate
    NigData <- list(y = DATA[[k]]$y, z1 = DATA[[k]]$z1, z2 = DATA[[k]]$z2, 
                    z3 = DATA[[k]]$z3, z4 = DATA[[k]]$z4, z5 = DATA[[k]]$z5, 
                    x = DATA[[k]]$x, medx = DATA[[k]]$xm,
                    alpha_pclu = c(1/3,1/3,1/3), 
                    R = R, delta_g2 = rep(0,28),
                    mine = 1/mineigenadj,
                    maxe = 1/maxeigenadj, I = diag(28), A = adj_matrix,
                    mu_wy = rep(0,28), mu_wr=rep(0,28))
    NigConsts <- list(N = length(NigData$y),  
                      prov = DATA[[k]]$prov, L = 124)
    x_init <- ifelse(is.na(DATA[[k]]$x),1,NA)
    eta_init <- matrix(rep(c(0,1,2,3), NigConsts$N), nrow = NigConsts$N, ncol = 4, byrow = T)
    NigInits <- list(beta = rnorm(6), alpha = rnorm(6),
                     gamma = rnorm(10), betax = rnorm(28),
                     sigma_wy = 1, sigma_wr = 1, 
                     Wy = rep(1, 28), Wr = rep(1, 28),
                     x = x_init, tau_x = 1, xm_s = 1,
                     mu_R = runif(NigConsts$N), eta = eta_init, 
                     gammax0 = rnorm(28), mu_x = runif(NigConsts$N), 
                     mu_g1 = -1, mu_g3 = 1, tau_g1 = 1, tau_g3 = 1, 
                     delta_g1 = runif(28, -1,0), delta_g3 = runif(28, 0,1),
                     pclu = c(0.1,0.3,0.6), phi = matrix(0, nrow = 4, ncol = 6),
                     clu = c(rep(1,10), rep(2,10), rep(3, 8)),
                     gamma_wy = 0, gamma_wr = 0,
                     L_wy = diag(28), L_wr = diag(28),
                     py = matrix(rep(c(0,0.1,0.2,0.3,0.4),NigConsts$N), nrow=NigConsts$N, ncol = 5, byrow = T)
    )
    Nig <- nimbleModel(code = NigCode, name = "Nig", constants = NigConsts,
                       data = NigData, inits = NigInits)
    cNig <- compileNimble(Nig)
    Nigconf <- configureMCMC(Nig, print = FALSE)
    Nigconf$addMonitors(c("beta","alpha","gamma", "gammax0",
                          "tau_x","phi",
                          "sigma_wy","sigma_wr",
                          "gamma_wy", "gamma_wr","clu", 
                          "xm_s"
    ))
    Nigmcmc <- buildMCMC(Nigconf)
    cNigmcmc <- compileNimble(Nigmcmc, project = Nig)
    cNig$setInits(NigInits)
    mcmc.out <- runMCMC(cNigmcmc, niter = 150000, setSeed = 123, thin = 15)
    pos_mcmc <- as.mcmc(mcmc.out[-c(1:4000),])







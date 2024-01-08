XMRi <- function(y1,X1,y2,X2,yt,Xt, beta=0, sig_alpha, a_Sig, B_Sig, a_alpha, b_alpha, a_tau, b_tau,
                      gibbsIter = 1000, burnin_prop = 0.2, null_prop = 0.1, fix.beta = FALSE){
  ### Get the number of samples (n) and variables (p).
  n1 <- nrow(X1)

  n2 <- nrow(X2)

  nt <- nrow(Xt)

  p <- ncol(X1)

  burninIter <- ceiling(burnin_prop * gibbsIter)

  ### Get the number of samples (n) and variables (p).


  s1 <- 1/sqrt(n1)
  s2 <- 1/sqrt(n2)
  st <- 1/sqrt(nt)

  bh1 <- t(X1) %*% y1/n1
  bh2 <- t(X2) %*% y2/n2
  bht <- t(Xt) %*% yt/nt


  # LD mat

  R1 <- t(X1) %*%X1/n1
  R2 <- t(X2) %*%X2/n2
  Rt <- t(Xt) %*%Xt/nt

  R1.tilde <- R1-diag(p)
  R2.tilde <- R2-diag(p)
  Rt.tilde <- Rt-diag(p)


  # initial values of hyper parameters
  tau <- 1/p
  if(missing(sig_alpha)) sig_alpha <- 0.01/p
  Sigma <- diag(rep(0.01,2))/(p*tau)
  # beta <- 0


  # beta prior of tau: with mean 1/p
  if(missing(a_tau))
    a_tau <- 0.5

  if(missing(b_tau))
    b_tau <- a_tau * (1-tau) / tau


  # inverse gamma prior of sigma_alpha: with mean sig_alpha
  # if(missing(b_alpha))
  #   b_alpha <- 0.001/p
  #
  # if(missing(a_alpha))
  #   a_alpha <- b_alpha/(sig_alpha) + 1

  if(missing(a_alpha))
    a_alpha <- 1

  if(missing(b_alpha))
    b_alpha <- sig_alpha * (a_alpha-1)

  # inverse Wishart prior of gamma: with mean diag(rep(0.01,2))/(p*tau)
  if(missing(B_Sig))
    B_Sig <- diag(rep(1,2))

  if(missing(a_Sig))
    a_Sig <- mean(diag(B_Sig) / diag(Sigma)) + 3

  # Gibbs samples
  sample_beta <- rep(0,gibbsIter)
  sample_tau <- rep(0,gibbsIter)
  sample_sig_alpha <- rep(0,gibbsIter)
  sample_sig_beta <- rep(0,gibbsIter)
  sample_mu_beta <- rep(0,gibbsIter)
  sample_Sig <- array(0,dim = c(2,2,gibbsIter))
  sample_w <- matrix(0,gibbsIter,p)


  # initial values of latent variables
  Z <- rep(0,p)
  gamma <- matrix(0,p,2)
  alpha <- rep(0,p)
  ### Gibbs sampling
  for (iter in 1 : gibbsIter){

    # sample latent variables
    inv_Sigma_j <- diag(c(1/s1^2+beta^2/st^2,1/s2^2)) + solve(Sigma)
    Sigma_j <- solve(inv_Sigma_j)
    logDS_j <- determinant(Sigma_j)$modulus[1]
    logDSigma <- determinant(Sigma)$modulus[1]
    for (j in 1 : p ){

      mu_j <- Sigma_j %*% c(bh1[j]/s1^2 - sum(R1.tilde[,j]*gamma[,1]*Z)/s1^2 + beta*(bht[j]/st^2 - beta*sum(Rt.tilde[,j]*gamma[,1]*Z)/st^2 - sum(alpha*Rt[,j])/st^2),
                            bh2[j]/s2^2 - sum(R2.tilde[,j]*gamma[,2]*Z)/s2^2)


      u <- log(tau / (1 - tau)) + logDS_j / 2 - logDSigma / 2 + t(mu_j) %*% inv_Sigma_j %*% mu_j / 2
      sample_w[iter,j] <- 1 / (1 + exp(-u))
      # part1 <- exp(log(tau) + logDS_j / 2 - determinant(Sigma)$modulus[1] / 2 + t(mu_j) %*% inv_Sigma_j %*% mu_j / 2)
      # part2 <- 1-tau
      # sample_w[iter,j] <- part1/(part1+part2)

      randu_num <- runif(1)
      if(rbinom(1,1,sample_w[iter,j])==1){
        Z[j] <- 1
        gamma[j,] <- mvtnorm::rmvnorm(1,mu_j,Sigma_j)
      } else {
        Z[j] <- 0
        gamma[j,] <- mvtnorm::rmvnorm(1,rep(0,2),Sigma)
      }

      sig_alpha_j <- 1 / (1/st^2 + 1/sig_alpha)
      mu_alpha_j <- sig_alpha_j * (bht[j]/st^2 - beta*sum(gamma[,1]*Z*Rt[,j])/st^2 - sum(Rt.tilde[,j]*alpha)/st^2)
      alpha[j] <- rnorm(1,mu_alpha_j,sqrt(sig_alpha_j))
    }

    # sample hyper-parameters
    shape1_tau <- sum(Z) + a_tau
    shape2_tau <- sum(1-Z) + b_tau
    tau <- rbeta(1,shape1_tau,shape2_tau)

    shape_sig_alpha <- p/2 + a_alpha
    scale_sig_alpha <- sum(alpha^2)/2 + b_alpha
    sig_alpha <- MCMCpack::rinvgamma(1,shape_sig_alpha,scale_sig_alpha)

    df_Sigma <- p + a_Sig
    scale_Sigma <- t(gamma) %*% gamma + B_Sig
    Sigma <- MCMCpack::riwish(df_Sigma,scale_Sigma)

    # df_Sigma <- sum(Z) + a_Sig
    # scale_Sigma <- t(gamma*Z) %*% (gamma*Z) + B_Sig
    # Sigma <- MCMCpack::riwish(df_Sigma,scale_Sigma)


    # sig_beta <- 1 / (sum(gamma[,1]*(Rt%*%gamma[,1]))/st^2)
    # mu_beta <- sig_beta * (sum(bht*gamma[,1])/st^2 - sum(alpha*(Rt%*%gamma[,1]))/st^2)
    # beta <- ifelse(iter>0,rnorm(1,mu_beta,sqrt(sig_beta)),0)


    sig_beta <- 1 / (sum((gamma[,1]*Z)*(Rt%*%(gamma[,1]*Z)))/st^2 + 1/1000000)
    mu_beta <- sig_beta * (sum(bht*gamma[,1]*Z)/st^2 - sum(alpha*(Rt%*%(gamma[,1]*Z)))/st^2)
    beta <- ifelse(!fix.beta&iter>ceiling(null_prop*gibbsIter),rnorm(1,mu_beta,sqrt(sig_beta)),beta)

    # record MCMC samples
    sample_beta[iter] <- beta
    sample_sig_alpha[iter] <- sig_alpha
    sample_sig_beta[iter] <- sig_beta
    sample_mu_beta[iter] <- mu_beta
    sample_tau[iter] <- tau
    sample_Sig[,,iter] <- Sigma

    # if(iter>burninIter){
    #   cat(burninIter," burn-in iterations finished.\n")
    #   if(iter%%500==0){
    #     cat(iter,"-th iteration: mean(beta)=",mean(sample_beta[(burninIter+1):iter]),"; mean(sig_alpha)=",mean(sample_sig_alpha[(burninIter+1):iter]),"; mean(tau)=",mean(sample_tau[(burninIter+1):iter]),".\n")
    #   }
    # }
    cat(iter,"-th iteration: mean(beta)=",mean(sample_beta[1:iter]),"; mean(sig_alpha)=",mean(sample_sig_alpha[1:iter]),"; mean(tau)=",mean(sample_tau[1:iter]),".\n")
  }

  mean_beta <- mean(sample_beta[(burninIter+1):gibbsIter])
  sd_beta <- sd(sample_beta[(burninIter+1):gibbsIter])

  return(  list(sample_beta=sample_beta,sample_sig_alpha=sample_sig_alpha,sample_Sig=sample_Sig,
                sample_tau=sample_tau,sample_w=sample_w,sample_sig_beta=sample_sig_beta,sample_mu_beta=sample_mu_beta,
                mean_beta=mean_beta,sd_beta=sd_beta,Z=Z,gamma=gamma,alpha=alpha)  )

}

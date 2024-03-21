### Codes in this script are downloaded from https://github.com/Minzhe/geneck/tree/master/program/RemoteR/geneck
### The original article is: Zhang, M. et al. GeNeCK: a web server for gene network construction and visualization. BMC Bioinformatics 20, 12 (2019).

###                         BayesianGLasso.R                          ###
### ================================================================= ###
# This R script is function modified from blockGLasso function in BayesianGLasso package.

#' Title
#'
#' @param X data.frame
#' @param iterations interger
#' @param burnIn interger
#' @param lambdaPriora numeric
#' @param lambdaPriorb numeric
#' @param verbose logic
#'
#' @return double
#' @export
#'
#' @examples
#' print('ok')
blockGLasso_s <- function(X, iterations = 2000, burnIn = 1000, lambdaPriora = 1, lambdaPriorb = 1 / 10, verbose = TRUE) {
  
  totIter <- iterations + burnIn
  S <- t(X) %*% X
  Sigma <- stats::cov(X)
  n <- nrow(X)
  Omega <- MASS::ginv(Sigma)
  p <- dim(Omega)[1]
  Rho <- matrix(0, nrow = p, ncol = p)
  indMat <- matrix(1:p ^ 2, ncol = p, nrow = p)
  perms <- matrix(NA, nrow = p - 1, ncol = p)
  permInt <- 1:p
  for (i in 1:ncol(perms)) {
    perms[, i] <- permInt[-i]
  }
  # SigmaMatList <- OmegaMatList <- list()
  # lambdas <- rep(NA, totIter)
  tau <- matrix(NA, nrow = p, ncol = p)
  lambdaPosta <- (lambdaPriora + (p * (p + 1) / 2))
  count <- 0
  
  for (iter in 1:totIter) {
    lambdaPostb <- (lambdaPriorb + sum(abs(c(Omega))) / 2)
    lambda <-
      stats::rgamma(1, shape = lambdaPosta, scale = 1 / lambdaPostb)
    OmegaTemp <- Omega[lower.tri(Omega)]
    rinvgaussFun <- function(x) {
      x <- ifelse(x < 1e-12, 1e-12, x)
      return(statmod::rinvgauss(n = 1, mean = x, shape = lambda ^ 2))
    }
    tau[lower.tri(tau)] <- 1 / sapply(sqrt(lambda ^ 2 / (OmegaTemp ^ 2)), rinvgaussFun)
    tau[upper.tri(tau)] <- t(tau)[upper.tri(t(tau))]
    for (i in 1:p) {
      tauI <- tau[perms[, i], i]
      Sigma11 <- Sigma[perms[, i], perms[, i]]
      Sigma12 <- Sigma[perms[, i], i]
      S21 <- S[i, perms[, i]]
      Omega11inv <- Sigma11 - Sigma12 %*% t(Sigma12) / Sigma[i,i]
      Ci <- (S[i, i] + lambda) * Omega11inv + diag(1 / tauI)
      CiChol <- chol(Ci)
      mui <- solve(-Ci, S[perms[, i], i])
      beta <- mui + solve(CiChol, stats::rnorm(p - 1))
      Omega[perms[, i], i] <- beta
      Omega[i, perms[, i]] <- beta
      gamm <- stats::rgamma(n = 1, shape = n / 2 + 1, rate = (S[1,1] + lambda) / 2)
      Omega[i, i] <- gamm + t(beta) %*% Omega11inv %*% beta
      OmegaInvTemp <- Omega11inv %*% beta
      Sigma[perms[, i], perms[, i]] <- Omega11inv + (OmegaInvTemp %*% t(OmegaInvTemp)) / gamm
      Sigma[perms[, i], i] <- Sigma[i, perms[, i]] <- (-OmegaInvTemp / gamm)
      Sigma[i, i] <- 1 / gamm
    }
    # if (iter%%100 == 0) {
    #   cat("Total iterations= ", iter, "Iterations since burn in= ",
    #       ifelse(iter - burnIn > 0, iter - burnIn, 0),
    #       "\n")
    # }
    # lambdas[iter] <- lambda
    # SigmaMatList[[iter]] <- Sigma
    # OmegaMatList[[iter]] <- Omega
    if (verbose) {
      if (floor(iter / (iterations + burnIn) * 100) == count) {
        print(paste0(count, "% has been done"))
        count <- count + 10
      }
    }
    
    if (iter > burnIn) {
      for (j in 1:p) {
        for (jj in 1:p) {
          Rho[j, jj] <- Rho[j, jj] + (-Omega[j, jj] / sqrt(Omega[j, j] * Omega[jj, jj]))
          
        }
      }
    }
  }
  Rho <- Rho / iterations
  
  # list(Sigmas = SigmaMatList, Omegas = OmegaMatList, lambdas = lambdas)
  return (Rho)
  
}

#' Title
#'
#' @param G 1
#' @param Gval 1
#' @param order 1
#' @param dat 1
#' @param t 1
#' @param lambda 1
#' @param nargin 1
#'
#' @return 1
#' @export
#'
#' @examples
#' print('')
pcacmi_edgereduce <- function(G, Gval, order, dat, t, lambda, nargin = 2) {
  if (order == 0) {
    for (i in 1:(nrow(G) - 1)) {
      for (j in (i + 1):nrow(G)) {
        if (G[i, j] != 0) {
          cmiv = cmi(dat[i,], dat[j,])
          Gval[i, j] = cmiv
          Gval[j, i] = cmiv
          
          if (cmiv < lambda) {
            G[i, j] = 0
            G[j, i] = 0
          }
        }
      }
    }
    t = t + 1
    
  } else {
    for (i in 1:(nrow(G) - 1)) {
      for (j in (i + 1):nrow(G)) {
        if (G[i, j] != 0) {
          adj = numeric()
          
          for (k in 1:nrow(G)) {
            if (G[i, k] != 0 && G[j, k] != 0) {
              adj = c(adj, k)
            }
          }
          if (length(adj) >= order) {
            if (length(adj) == 1) {
              combntnslist = adj
              combntnsrow = 1
            } else {
              combntnslist = t(utils::combn(adj, order))
              combntnsrow = nrow(combntnslist)
            }
            cmiv = 0
            v1 = dat[i,]
            v2 = dat[j,]
            
            if (combntnsrow == 1) {
              vcs = dat[combntnslist,]
              cmiv = cmi(v1, v2, 3, vcs)
            } else {
              for (k in 1:combntnsrow) {
                vcs = dat[combntnslist[k,],]
                
                a = cmi(v1,
                        v2,
                        3,
                        vcs)
                
                cmiv = max(cmiv, a)
                
              }
            }
            Gval[i, j] = cmiv
            Gval[j, i] = cmiv
            
            if (cmiv < lambda) {
              G[i, j] = 0
              G[j, i] = 0
            }
            t = t + 1
          }
        }
      }
    }
  }
  return(list(G = G, Gval = Gval, t = t))
}

#' Title
#'
#' @param dat data.frame
#' @param lambda double
#' @param order0 interger
#' @param nargin interger
#'
#' @return list
#' @export
#'
#' @examples
#' print('')
pca_cmi <- function(dat, lambda, order0 = 0, nargin = 2) {
  # function [G,Gval,order]=pca_cmi(data,lamda,order0)
  
  n_gene = nrow (dat)
  G = matrix(1, n_gene, n_gene)
  G = t(as.matrix(Matrix::tril(G, -1)))
  
  G = G + t(G)
  Gval = G
  order = -1
  t = 0
  
  while (t == 0) {
    order = order + 1
    
    if (nargin == 3) {
      if (order > order0) {
        G = t(as.matrix(Matrix::tril(G, -1)))
        
        Gval = t(as.matrix(Matrix::tril(G, -1)))
        
        order = order - 1
        # The value of order is the last order of pc algorith
        return(list(
          G = G,
          Gval = Gval,
          order = order
        ))
      }
    }
    
    res = pcacmi_edgereduce(G, Gval, order, dat, t, lambda)
    
    G = res$G
    Gval = res$Gval
    t = res$t
    
    if (t == 0) {
      # cat('\n No edge is reduce! Algorithm  finished!')
      break
    } else {
      t = 0
    }
  }
  
  G = t(as.matrix((Matrix::tril(G, -1))))
  
  Gval = t(as.matrix((Matrix::tril(Gval, -1))))
  
  order = order - 1
  # The value of order is the last order of pc algorith
  return(list(G = G,
              Gval = Gval,
              order = order))
}

###                    GeneNet.R                     ###
### ================================================ ###
# This R script is function to use GeneNet to constrcut gene network.

# setwd(paste(Sys.getenv("remoter_path"), "geneck/", sep = ""))

#' Title
#'
#' @param expr.data data.frame
#' @param fdr double
#'
#' @return data.frame
#' @export
#'
#' @examples
#' print('')
network.GeneNet <- function(expr.data, fdr) {
  
  if (fdr <= 0 | fdr >= 1) {
    stop('Input error: false discovery rate for GeneNet should bewteen 0 and 1.')
  }
  p <- ncol(expr.data)
  n <- nrow(expr.data)
  gene.index <- colnames(expr.data)
  
  expr.mat <- scale(as.matrix(expr.data), center = TRUE, scale = FALSE)
  
  pcor_est <- corpcor::pcor.shrink(expr.mat)
  test_pcor <- GeneNet::network.test.edges(pcor_est)
  
  est_res <- test_pcor[test_pcor[, 6] >= (1-fdr), ]
  est_edge <- est_res[, c(2, 3, 1)]
  est_edge[,3] <- signif(est_edge[,3], 4)
  names(est_edge)[3] <- "p.corr"
  est_edge[,1] <- gene.index[est_edge[,1]]
  est_edge[,2] <- gene.index[est_edge[,2]]
  
  file.remove("Rplots.pdf")
  
  return(est_edge)
}

# glasso-SF by Liu and Ihler

#' Title
#'
#' @param X data.frame
#' @param alpha double
#' @param niter interger
#'
#' @return data.frame
#' @export
#'
#' @examples
#' print('')
glasso_sf <- function(X, alpha, niter = 5) {
  
  n = nrow(X)
  p = ncol(X)
  
  S = t(X) %*% X / n
  
  eps = 1
  bt = 2 * alpha / eps
  
  lam_mat = matrix(2 * alpha, p, p)
  for (iter in 1:niter) {
    out = glasso::glasso(S, rho = lam_mat)
    # diag(out$wi) = 0
    th_m = apply(abs(out$wi), 2, sum)
    th_m = th_m - abs(diag(out$wi))
    # Update weights of tuning parameter
    for (i in 1:(p - 1))
      for (j in (i + 1):p)
        lam_mat[i, j] <- lam_mat[j, i] <- alpha * (1 / (th_m[i] + eps) + 1 / (th_m[j] + eps))
    # diag(lam_mat) <- 2*alpha
  }
  
  return(out)
}

#' Title
#'
#' @param Y.m matrix
#' @param lam1 double
#' @param lam2 double
#' @param sig NULL
#' @param weight NULL
#' @param iter integer
#'
#' @return list
#' @export
#'
#' @examples
#' print('')
space.joint <- function (Y.m, lam1, lam2 = 0, sig = NULL, weight = NULL, iter = 2)
{
  n = nrow(Y.m)
  p = ncol(Y.m)
  ITER = iter
  if (!is.null(sig)) {
    SIG.update = F
    SIG = sig
  } else {
    SIG.update = T
    SIG = rep(1, p)
  }
  if (length(weight) == 0 | (length(weight) > 1 & length(weight) <
                             p)) {
    WEIGHT.tag = 0
    WEIGHT = rep(1, p)
    WEIGHT.update = F
  }
  if (length(weight) == 1) {
    if (weight == 1) {
      WEIGHT.tag = 1
      WEIGHT = SIG
      WEIGHT.update = T
      ITER = max(2, iter)
    } else {
      WEIGHT.tag = 2
      WEIGHT = rep(1, p)
      WEIGHT.update = T
      ITER = max(2, iter)
    }
  }
  if (length(weight) == p) {
    WEIGHT.tag = 3
    WEIGHT = weight
    WEIGHT.update = F
  }
  for (i in 1:iter) {
    print(paste("iter=", i, sep = ""))
    Y.u <- Y.m * matrix(sqrt(WEIGHT), n, p, byrow = TRUE)
    sig.u <- SIG/WEIGHT
    ParCor.fit <- jsrm(Y.u, sig.u, n, p, lam1, lam2)
    diag(ParCor.fit) <- 1
    coef <- ParCor.fit[upper.tri(ParCor.fit)]
    beta.cur <- Beta.coef(coef, SIG)
    if (!WEIGHT.update & !SIG.update) {
      break
    } else {
      if (SIG.update) {
        SIG <- InvSig.diag.new(Y.m, beta.cur)
      }
      if (WEIGHT.update) {
        if (WEIGHT.tag == 1) {
          WEIGHT = SIG
        } else {
          temp.w <- apply(abs(ParCor.fit) > 1e-06, 1,
                          sum)
          temp.w <- temp.w + max(temp.w)
          WEIGHT <- temp.w/sum(temp.w) * p
        }
      }
    }
  }
  result <- list(ParCor = ParCor.fit, sig.fit = SIG)
  return(result)
}

#' Title
#'
#' @param Y matrix
#' @param sig.use double
#' @param n integer
#' @param p double
#' @param lam1 double
#' @param lam2 double
#' @param n_iter integer
#'
#' @return matrix
#' @export
#'
#' @examples
#' print('')
jsrm<-function(Y,sig.use,n,p,lam1,lam2, n_iter=1000)
{
  lambda1=lam1
  lambda2=lam2
  sigma_sr=sig.use^0.5
  
  Beta_output<-rep(0, p*p)
  Y_data<-as.vector(t(Y))
  n_iter=500
  
  #dyn.load("JSRM.so") ### compiled from "JSRM.c"
  
  iter_count=0
  
  junk<-.C("JSRM",
           as.integer(n),
           as.integer(p),
           as.single(lambda1),
           as.single(lambda2),
           as.single(Y_data),
           as.single(sigma_sr),
           as.integer(n_iter),
           iter.count=as.integer(iter_count),
           beta.estimate=as.single(Beta_output)
  )
  
  #print(junk$iter.count)
  beta.v<-junk$beta.estimate
  beta.m<-matrix(beta.v, p,p, byrow=T)
  return(beta.m)
}

#' Title
#'
#' @param Y matrix
#' @param Beta matrix
#'
#' @return matrix
#' @export
#'
#' @examples
#' print('')
InvSig.diag.new<-function(Y, Beta){
  ################### parameters
  ### Y:    n by p data matrix; should be standardized to sd 1;
  ### Beta: beta matrix for the regression model
  p=ncol(Y)
  Beta.temp=Beta
  diag(Beta.temp)=0
  esti.Y=Y%*%Beta.temp
  residue=Y-esti.Y
  result=apply(residue^2, 2, mean)
  return(1/result)
}

#' Title
#'
#' @param coef double
#' @param sig.fit numeric
#'
#' @return matrix
#' @export
#'
#' @examples
#' print('')
Beta.coef<-function(coef, sig.fit){
  ############## parameter
  ### coef: rho^{ij};
  ### sig.fit: sig^{ii}
  p<-length(sig.fit)
  result<-matrix(0,p,p)
  result[upper.tri(result)]<-coef
  result<-result+t(result)
  result<-diag(1/sqrt(sig.fit))%*%result%*%diag(sqrt(sig.fit))
  result<-t(result)
  return(result)
}

#' Title
#'
#' @param v1 character
#' @param v2 character
#' @param nargin numeric
#' @param vcs numeric
#'
#' @return numeric
#' @export
#'
#' @examples
#' print('')
cmi <- function(v1, v2, nargin = 2, vcs = NULL) {
  if (nargin == 2) {
    c1 = stats::var(v1)
    
    c2 = stats::var(v2)
    
    c3 = det(stats::cov(t(rbind(v1, v2))))
    
    cmiv = 0.5 * log(c1 * c2 / c3)
    
  } else if (nargin == 3) {
    c1 = det(stats::cov(t(rbind(v1, vcs))))
    
    c2 = det(stats::cov(t(rbind(v2, vcs))))
    
    if (is.null(dim(vcs))) {
      c3 = stats::var(vcs)
    } else {
      c3 = det(stats::cov(t(vcs)))
      
    }
    c4 = det(stats::cov(t(rbind(v1, v2, vcs))))
    
    cmiv = 0.5 * log((c1 * c2) / (c3 * c4))
    
  }
  
  #if(nargin==2)
  #cat('\n values: ',c1,c2,c3,cmiv,sep=' ')
  #if(nargin==3)
  #cat('\n values: ',c1,c2,c3,c4,cmiv,sep=' ')
  
  if (cmiv == Inf) {
    cmiv = 1e10
    
  }
  
  return(cmiv)
}

#' Title
#'
#' @param dat data.frame
#' @param lambda double
#' @param order0 interger
#' @param nargin interger
#'
#' @return list
#' @export
#'
#' @examples
#' print('')
cmi2ni <- function(dat, lambda, order0 = 0, nargin = 2) {
  #function [G,Gval,order]=pca_cmi(data,lamda,order0)
  
  n_gene = nrow (dat)
  G = matrix(1, n_gene, n_gene)
  G = t(as.matrix(Matrix::tril(G, -1)))
  
  G = G + t(G)
  
  Gval = G
  
  order = -1
  
  t = 0
  
  while (t == 0) {
    order = order + 1
    
    if (nargin == 3) {
      if (order > order0) {
        G = t(as.matrix(Matrix::tril(G, -1)))
        
        Gval = t(as.matrix(Matrix::tril(G, -1)))
        
        order = order - 1
        # The value of order is the last order of pc algorith
        return(list(
          G = G,
          Gval = Gval,
          order = order
        ))
      }
    }
    
    res = edgereduce(G, Gval, order, dat, t, lambda)
    
    G = res$G
    Gval = res$Gval
    t = res$t
    
    if (t == 0) {
      # cat('\n No edge is reduce! Algorithm  finished!')
      
      break
      
    } else {
      t = 0
      
    }
  }
  
  G = t(as.matrix(Matrix::tril(G, -1)))
  
  Gval = t(as.matrix(Matrix::tril(Gval, -1)))
  
  order = order - 1
  # The value of order is the last order of pc algorith
  return(list(
    G = G,
    Gval = Gval,
    order = order
  ))
}

#' Title
#'
#' @param x character
#' @param y character
#' @param z character
#'
#' @return numeric
#' @export
#'
#' @examples
#' print('')
cas <- function(x, y, z) {
  # x=rand(10,1)';y=rand(10,1)';z=rand(10,2)';
  if (is.null(dim(z))) {
    n1 = 1
  } else {
    n1 = nrow(z)
  }
  n = n1 + 2
  
  
  Cov = as.matrix(stats::var(x))
  
  Covm = stats::cov(t(rbind(x, y, z)))
  
  Covm1 = stats::cov(t(rbind(x, z)))
  
  
  InvCov = solve(Cov)
  
  InvCovm = solve(Covm)
  
  InvCovm1 = solve(Covm1)
  
  
  C11 = InvCovm1[1, 1]
  
  C12 = 0
  
  C13 = InvCovm1[1, 2:(1 + n1)]
  
  C23 = InvCovm[2, 3:(2 + n1)] - InvCovm[1, 2] * (1 / (InvCovm[1, 1] - InvCovm1[1, 1] +
                                                         InvCov[1, 1])) *
    (InvCovm[1, 3:(2 + n1)] - InvCovm1[1, 2:(1 + n1)])
  
  C22 = InvCovm[2, 2] - InvCovm[1, 2] ^ 2 * (1 / (InvCovm[1, 1] - InvCovm1[1, 1] +
                                                    InvCov[1, 1]))
  
  C33 = InvCovm[3:(2 + n1), 3:(2 + n1)] - (1 / (InvCovm[1, 1] - InvCovm1[1, 1] +
                                                  InvCov[1, 1])) *
    ((InvCovm[1, 3:(2 + n1)] - InvCovm1[1, 2:(1 + n1)]) %o% (InvCovm[1, 3:(2 +
                                                                             n1)] - InvCovm1[1, 2:(1 + n1)]))
  
  InvC = rbind(c(C11, C12, C13), c(C12, C22, C23), cbind(C13, C23, C33))
  
  #%C = inv(InvC);
  
  C0 = Cov[1, 1] * (InvCovm[1, 1] - InvCovm1[1, 1] + InvCov[1, 1])
  
  CS = 0.5 * (sum(diag(InvC %*% Covm)) + log(C0) - n)
  
  
  return(CS)
}

#' Title
#'
#' @param G matrix
#' @param Gval matrix
#' @param order numeric
#' @param dat matrix
#' @param t numeric
#' @param lambda numeric
#' @param nargin numeric
#'
#' @return list
#' @export
#'
#' @examples
#' print('')
edgereduce <- function(G, Gval, order, dat, t, lambda, nargin = 2) {
  if (order == 0) {
    for (i in 1:(nrow(G) - 1)) {
      for (j in (i + 1):nrow(G)) {
        if (G[i, j] != 0) {
          cmiv = cmi(dat[i, ], dat[j, ])
          
          Gval[i, j] = cmiv
          
          Gval[j, i] = cmiv
          
          if (cmiv < lambda) {
            G[i, j] = 0
            
            G[j, i] = 0
            
          }
        }
      }
    }
    t = t + 1
    
  } else {
    for (i in 1:(nrow(G) - 1)) {
      for (j in (i + 1):nrow(G)) {
        if (G[i, j] != 0) {
          adj = numeric()
          
          for (k in 1:nrow(G)) {
            if (G[i, k] != 0 && G[j, k] != 0) {
              adj = c(adj, k)
              
            }
          }
          if (length(adj) >= order) {
            if (length(adj) == 1) {
              combntnslist = adj
              combntnsrow = 1
            } else {
              combntnslist = t(utils::combn(adj, order))
              
              combntnsrow = nrow(combntnslist)
              
            }
            cmiv = 0
            
            v1 = dat[i, ]
            
            v2 = dat[j, ]
            
            if (combntnsrow == 1) {
              vcs = dat[combntnslist, ]
              #cat('\n ',adj,'adj',combntnsrow,k,i,j,' list ', combntnslist,' dim ', nrow(vcs), ncol(vcs),sep=' ')
              cmiv = MI2(v1, v2, vcs)
            } else {
              for (k in 1:combntnsrow) {
                vcs = dat[combntnslist[k, ], ]
                
                #cat('\n ',adj,'adj',combntnsrow,k,i,j,' list ', combntnslist[k,],' dim ', nrow(vcs), ncol(vcs),sep=' ')
                a = MI2(v1, v2, vcs)
                
                cmiv = max(cmiv, a)
                
              }
            }
            Gval[i, j] = cmiv
            
            Gval[j, i] = cmiv
            
            if (cmiv < lambda) {
              G[i, j] = 0
              
              G[j, i] = 0
              
            }
            t = t + 1
            
          }
        }
        
      }
    }
  }
  return(list(G = G, Gval = Gval, t = t))
}

#' Title
#'
#' @param x character
#' @param y character
#' @param z character
#'
#' @return numeric
#' @export
#'
#' @examples
#' print('')
MI2 <- function(x, y, z) {
  r_dmi = (cas(x, y, z) + cas(y, x, z)) / 2
  
  return(r_dmi)
}



###                       bayesianglasso.R                       ###
### ============================================================ ###
# This R script is function to use bayesianglasso to constrcut gene network.

# setwd(paste(Sys.getenv("remoter_path"), "geneck/", sep = ""))
# source("lib.BayesianGLasso.R")

#' Title
#'
#' @param expr.data matrix
#' @param prob numeric
#'
#' @return data.frame
#' @export
#'
#' @examples
#' print('ok')
network.bayesianglasso <- function(expr.data, prob) {
  if (prob <= 0 | prob >= 1) {
    stop('Input error: probability for BayesianGLasso should bewteen 0 and 1.')
  }
  p <- ncol(expr.data)
  n <- nrow(expr.data)
  gene.index <- colnames(expr.data)
  
  expr.mat <- scale(as.matrix(expr.data), center = TRUE, scale = FALSE)
  
  a <- 10^(-2); b <- 10^(-6); iter <- 2000; burn <- 1000
  out <- blockGLasso_s(expr.mat, iterations = iter, burnIn = burn, lambdaPriora = a, lambdaPriorb = b, verbose = FALSE)
  est_edge <- which(abs(out) > prob, T)
  est_edge <- est_edge[est_edge[, 1] < est_edge[, 2], ]
  if (length(est_edge) == 2) est_edge <- matrix(est_edge, 1, 2)
  
  est_edge <- as.data.frame(est_edge)
  colnames(est_edge) <- c("node1", "node2")
  est_edge$p.corr <- 0
  for (i in 1:nrow(est_edge)) {
    est_edge$p.corr[i] <- out[est_edge$node1[i], est_edge$node2[i]]
  }
  
  est_edge$p.corr <- signif(est_edge$p.corr, 4)
  est_edge <- est_edge[with(est_edge, order(-abs(p.corr))),]
  est_edge[,1] <- gene.index[est_edge[,1]]
  est_edge[,2] <- gene.index[est_edge[,2]]
  
  return(est_edge)
}


###                       cmi2ni.R                         ###
### ====================================================== ###
# This R script is function to use cmi2ni to constrcut gene network.

# setwd(paste(Sys.getenv("remoter_path"), "geneck/", sep = ""))
# source("lib.CMI2NI.R")

#' Title
#'
#' @param expr.data data.frame
#' @param lambda double
#'
#' @return data.frame
#' @export
#'
#' @examples
#' print('ok')
network.cmi2ni <- function(expr.data, lambda) {
  if (lambda <= 0) {
    stop('Input error: parameter lambda for cmi2ni should be larger than 0.')
  }
  p <- ncol(expr.data)
  n <- nrow(expr.data)
  gene.index <- colnames(expr.data)
  
  expr.mat <- scale(as.matrix(expr.data), center = TRUE, scale = FALSE)
  
  out <- cmi2ni(t(expr.mat), lambda)
  est_edge <- which(out$G == 1, T)
  est_edge <- est_edge[est_edge[, 1] < est_edge[, 2], ]
  if (length(est_edge) == 2) est_edge <- matrix(est_edge, 1, 2)
  
  est_edge <- as.data.frame(est_edge)
  colnames(est_edge) <- c("node1", "node2")
  est_edge$strength <- 0
  for (i in 1:nrow(est_edge)) {
    est_edge$strength[i] <- out$Gval[est_edge$node1[i], est_edge$node2[i]]
  }
  
  est_edge$strength <- signif(est_edge$strength, 4)
  est_edge <- est_edge[with(est_edge, order(-abs(strength))),]
  est_edge[,1] <- gene.index[est_edge[,1]]
  est_edge[,2] <- gene.index[est_edge[,2]]
  
  return(est_edge)
}

###                       glasso.R                       ###
### ==================================================== ###
# This R script is function to use glasso to constrcut gene network.

# setwd(paste(Sys.getenv("remoter_path"), "geneck/", sep = ""))

#' Title
#'
#' @param expr.data data.frame
#' @param lambda double
#'
#' @return data.frame
#' @export
#'
#' @examples
#' print('')
network.glasso <- function(expr.data, lambda) {
  if (lambda <= 0) {
    stop('Input error: parameter lambda for ns should be larger than 1.')
  }
  p <- ncol(expr.data)
  n <- nrow(expr.data)
  gene.index <- colnames(expr.data)
  
  expr.mat <- scale(as.matrix(expr.data), center = TRUE, scale = FALSE)
  
  S <- t(expr.mat) %*% expr.mat / n
  out <- glasso::glasso(S, rho = lambda)
  est_edge <- which(abs(out$wi) > 0, T)
  est_edge <- est_edge[est_edge[, 1] < est_edge[, 2], ]
  if (length(est_edge) == 2) est_edge <- matrix(est_edge, 1, 2)
  
  est_edge <- as.data.frame(est_edge)
  colnames(est_edge) <- c("node1", "node2")
  est_edge$p.corr <- 0
  for (i in 1:nrow(est_edge)) {
    est_edge$p.corr[i] <- out$wi[est_edge$node1[i], est_edge$node2[i]]
  }
  
  est_edge$p.corr <- signif(est_edge$p.corr, 4)
  est_edge <- est_edge[with(est_edge, order(-abs(p.corr))),]
  est_edge[,1] <- gene.index[est_edge[,1]]
  est_edge[,2] <- gene.index[est_edge[,2]]
  
  return(est_edge)
}

###                       glassosf.R                       ###
### ====================================================== ###
# This R script is function to use glassosf to constrcut gene network.

# setwd(paste(Sys.getenv("remoter_path"), "geneck/", sep = ""))

#' Title
#'
#' @param expr.data data.frame
#' @param alpha double
#'
#' @return data.frame
#' @export
#'
#' @examples
#' print('')
network.glassosf <- function(expr.data, alpha) {
  if (alpha <= 0) {
    stop('Input error: parameter alpha for glassosf should be larger than 0.')
  }
  p <- ncol(expr.data)
  n <- nrow(expr.data)
  gene.index <- colnames(expr.data)
  
  expr.mat <- scale(as.matrix(expr.data), center = TRUE, scale = FALSE)
  
  out <- glasso_sf(expr.mat, alpha)
  est_edge <- which(abs(out$wi) > 0, T)
  est_edge <- est_edge[est_edge[, 1] < est_edge[, 2], ]
  if (length(est_edge) == 2) est_edge <- matrix(est_edge, 1, 2)
  
  est_edge <- as.data.frame(est_edge)
  colnames(est_edge) <- c("node1", "node2")
  est_edge$p.corr <- 0
  for (i in 1:nrow(est_edge)) {
    est_edge$p.corr[i] <- out$wi[est_edge$node1[i], est_edge$node2[i]]
  }
  
  est_edge$p.corr <- signif(est_edge$p.corr, 4)
  est_edge <- est_edge[with(est_edge, order(-abs(p.corr))),]
  est_edge[,1] <- gene.index[est_edge[,1]]
  est_edge[,2] <- gene.index[est_edge[,2]]
  
  return(est_edge)
}


###                       ns.R                       ###
### ================================================ ###
# This R script is function to use neighborhood selection to constrcut gene network.

# setwd(paste(Sys.getenv("remoter_path"), "geneck/", sep = ""))

#' Title
#'
#' @param expr.data data.frame
#' @param alpha double
#'
#' @return data.frame
#' @export
#'
#' @examples
#' print('')
network.ns <- function(expr.data, alpha) {
  if (alpha <= 0) {
    stop('Input error: parameter alpha for ns should be larger than 1.')
  }
  p <- ncol(expr.data)
  n <- nrow(expr.data)
  gene.index <- colnames(expr.data)
  
  expr.mat <- scale(as.matrix(expr.data), center = TRUE, scale = FALSE)
  
  est_res <- matrix(0, p, p)
  for (k in 1:p) {
    rsp <- expr.mat[, k]
    prd <- t(expr.mat[, -k])
    lam <- sqrt(sum(rsp ^ 2)) * stats::qnorm(alpha / (2 * p ^ 2), lower.tail = F)
    out <- CDLasso::l2.reg(prd, rsp, lambda = lam)
    est_res[k, -k] <- out$estimate
  }
  
  est_edge <- which(abs(est_res) > 0 & abs(t(est_res)) > 0, arr.ind = TRUE)
  est_edge <- est_edge[est_edge[, 1] < est_edge[, 2], ]
  if (length(est_edge) == 2) est_edge <- matrix(est_edge, 1, 2)
  
  est_edge <- as.data.frame(est_edge)
  colnames(est_edge) <- c("node1", "node2")
  est_edge$p.corr <- 0
  for (i in 1:nrow(est_edge)) {
    est_edge$p.corr[i] <- mean(est_res[est_edge$node1[i], est_edge$node2[i]], est_res[est_edge$node2[i], est_edge$node1[i]])
  }
  
  est_edge$p.corr <- signif(est_edge$p.corr, 4)
  est_edge <- est_edge[with(est_edge, order(-abs(p.corr))),]
  est_edge[,1] <- gene.index[est_edge[,1]]
  est_edge[,2] <- gene.index[est_edge[,2]]
  
  return(est_edge)
}

###                       space.R                          ###
### ====================================================== ###
# This R script is function to use space to constrcut gene network.

# setwd(paste(Sys.getenv("remoter_path"), "geneck/", sep = ""))

#' Title
#'
#' @param expr.data data.frame
#' @param alpha double
#'
#' @return data.frame
#' @export
#'
#' @examples
#' print('')
network.space <- function(expr.data, alpha) {
  if (alpha <= 0) {
    stop('Input error: parameter alpha for space should be larger than 0.')
  }
  p <- ncol(expr.data)
  n <- nrow(expr.data)
  gene.index <- colnames(expr.data)
  
  expr.mat <- scale(as.matrix(expr.data), center = TRUE, scale = FALSE)
  
  out <- space::space.joint(expr.mat, lam1 = alpha * n, iter = 5)
  est_edge <- which(abs(out$ParCor) > 0, T)
  est_edge <- est_edge[est_edge[, 1] < est_edge[, 2], ]
  if (length(est_edge) == 2) est_edge <- matrix(est_edge, 1, 2)
  
  est_edge <- as.data.frame(est_edge)
  colnames(est_edge) <- c("node1", "node2")
  est_edge$p.corr <- 0
  for (i in 1:nrow(est_edge)) {
    est_edge$p.corr[i] <- out$ParCor[est_edge$node1[i], est_edge$node2[i]]
  }
  
  est_edge$p.corr <- signif(est_edge$p.corr, 4)
  est_edge <- est_edge[with(est_edge, order(-abs(p.corr))),]
  est_edge[,1] <- gene.index[est_edge[,1]]
  est_edge[,2] <- gene.index[est_edge[,2]]
  
  return(est_edge)
  
  bic_val = sum(log(n / out$sig.fit)) + log(n) / n * nrow(est_edge) * 2
}


###                       pcacmi.R                         ###
### ====================================================== ###
# This R script is function to use pcacmi to constrcut gene network.

# setwd(paste(Sys.getenv("remoter_path"), "geneck/", sep = ""))
# source("lib.PCA_CMI.R")

#' Title
#'
#' @param expr.data data.frame
#' @param lambda double
#'
#' @return data.frame
#' @export
#'
#' @examples
#' print('')
network.pcacmi <- function(expr.data, lambda) {
  if (lambda <= 0) {
    stop('Input error: parameter alpha for pcacmi should be larger than 0.')
  }
  p <- ncol(expr.data)
  n <- nrow(expr.data)
  gene.index <- colnames(expr.data)
  
  expr.mat <- scale(as.matrix(expr.data), center = TRUE, scale = FALSE)
  
  out <- pca_cmi(t(expr.mat), lambda)
  est_edge <- which(out$G == 1, T)
  est_edge <- est_edge[est_edge[, 1] < est_edge[, 2], ]
  if (length(est_edge) == 2) est_edge <- matrix(est_edge, 1, 2)
  
  est_edge <- as.data.frame(est_edge)
  colnames(est_edge) <- c("node1", "node2")
  est_edge$strength <- 0
  for (i in 1:nrow(est_edge)) {
    est_edge$strength[i] <- out$Gval[est_edge$node1[i], est_edge$node2[i]]
  }
  
  est_edge$strength <- signif(est_edge$strength, 4)
  est_edge <- est_edge[with(est_edge, order(-abs(strength))),]
  est_edge[,1] <- gene.index[est_edge[,1]]
  est_edge[,2] <- gene.index[est_edge[,2]]
  
  return(est_edge)
}



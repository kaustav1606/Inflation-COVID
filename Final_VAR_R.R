
library(svars)
## import data

library(readxl)
dataind_seasadj_endog <- read_excel("C:/Users/user/Downloads/haroon code/dataind_seasadj_endog.xlsx",
                                    range = "B1:I111")
View(dataind_seasadj_endog)
# plot data


plot(dataind_seasadj_endog, nc = 2)
summary(dataind_seasadj_endog)
adf1 <- summary(ur.df(dataind_seasadj_endog["wpi"], type = "trend", lags = 2))



dataind_seasadj_endog_last = subset(dataind_seasadj_endog_last, select = -c(asset) )
dataind_seasadj_endog_last <- dataind_seasadj_endog_last[-c(97,98), ]

#var
VARselect(dataind_seasadj_endog_last, lag.max = 20, type = "both")

###############################################################
drop <- c("date")
Data = Data[,!(names(Data) %in% drop)]

options(max.print=999999)


vareg<-VAR(New_Microsoft_Excel_Worksheet, p = 1, type = "const")
summary(vareg)


#vars
var_x<-VAR(dataind_seasadj_endog_last, p = 6, type = "both")
summary(var_x)
Eigenvalues <- roots(var_x, modulus = TRUE)
normality.test(var_x)
serial.test(var_x, lags.pt = 16, type = "PT.asymptotic")
arch.test(var_x, lags.multi = 5,multivariate.only = TRUE)

# structural matrix
#amat -  contemporaneous
amat <- diag(9)
amat[1, 3] <- NA
amat[1, 4] <- NA
amat[1, 5] <- NA
amat[2, 1] <- NA
amat[2, 5] <- NA
amat[2, 6] <- NA
amat[3, 2] <- NA
amat[4, 3] <- NA
amat[4, 5] <- NA
amat[5, 1] <- NA
amat[5, 2] <- NA
amat[5, 4] <- NA
amat[1,7] <- NA #1.e  2.pi always have effect 
amat[1,8] <- NA
amat[1,9] <- NA
amat[2,7] <- NA
amat[2,8] <- NA
amat[2,9] <- NA
amat[3,7] <- NA
amat[3,8] <- NA
amat[3,9] <- NA
amat[4,7] <- NA
amat[4,8] <- NA
amat[4,9] <- NA
amat[5,7] <- NA
amat[5,8] <- NA
amat[5,9] <- NA
amat[6,7] <- NA
amat[6,8] <- NA
amat[6,9] <- NA

#restrict
restrict <- matrix(c(1, 1, 1,  1, 1, 1, 1, 1 ,1,1,1,1,1,1,1,1,1,1,1,
                     1, 1, 1,  1, 1, 1, 1, 1 ,1,1,1,1,1,1,1,1,1,1,1,
                     1, 1, 1,  1, 1, 1, 1, 1 ,1,1,1,1,1,1,1,1,1,1,1,
                     1, 1, 1,  1, 1, 1, 1, 1 ,1,1,1,1,1,1,1,1,1,1,1,
                     1, 1, 1,  1, 1, 1, 1, 1 ,1,1,1,1,1,1,1,1,1,1,1,
                     1, 1, 1,  1, 1, 1, 1, 1 ,1,1,1,1,1,1,1,1,1,1,1,
                     0, 0, 0,  0, 0, 1, 0, 0 ,0,0,0,0,0,0,1,0,0,1,1,
                     0, 0, 0,  0, 0, 0, 1, 0 ,0,0,0,0,0,0,0,1,0,1,1,
                     0, 0, 0,  0, 0, 0, 0, 1 ,0,0,0,0,0,0,0,0,1,1,1
),
nrow=8, ncol=20, byrow=TRUE)



#svars
svar.a <-SVAR(res, estmethod =  "direct", Amat = amat , Bmat =NULL ,max.iter = 100,
              lrtest = TRUE)
summary(svar.a)
res<-restrict(var_x, method = "man", resmat = restrict)
summary(res)

par(mar = c(1, 1, 1, 1))

#predict
varprd<-predict(res, n.ahead = 14, ci = 0.95)
plot(varprd)

#irf
irffood<-irf(svar.a, impulse = "gscpi", response = "cpi", boot =
               TRUE)
plot(irffood)
irf(svar.a, impulse =c("reer","wpi","wacr","asset","cf_cycle","wage", "oil","gscpi","food"), response ="wpi" , boot = FALSE)

irf(svar.a, impulse ="oil", response ="wpi" , boot = FALSE)
irffood<-irf(svar.a, impulse ="oil", response ="wpi" , boot = FALSE)
plot(irffood)
irffood<-irf(svar.a, impulse ="food", response ="wpi" , boot = FALSE)
plot(irffood)
irfgscpi<-irf(svar.a, impulse ="gscpi", response ="wpi" , boot = FALSE)
plot(irfgscpi)
irfgscpi<-irf(svar.a, impulse ="gscpi", response ="wpi" , boot = TRUE)
plot(irfgscpi)
irfgscpi<-irf(svar.a, impulse ="wacr", response ="wpi" , boot = TRUE)
plot(irfgscpi)

### VAR with t distribution


library(gmvarkit)
vartgm <- fitGSMVAR(dataind_seasadj_endog_last, p=4, M=2, model="StMVAR", ncalls=16,parametrization ="mean",maxit = 1000, ncores=8,seeds=1:16)
plot(fit12)
summary(fit12)
print_std_errors(fit12)
profile_logliks(fit12)


irf2 <- linear_IRF(vartgm, regime=1, N=20, ci=0.90, bootstrap_reps=10,
                   ncalls=1, seeds=1:10, ncores=1)
plot(irf2)

### historical decomposition




#VARhd Ambrosio
VARhd <- function(Estimation){
  
  ## make X and Y
  nlag    <- Estimation$p   # number of lags
  DATA    <- Estimation$y   # data
  QQ      <- VARmakexy(DATA,nlag,1)
  
  
  ## Retrieve and initialize variables 
  invA    <- t(chol(as.matrix(summary(Estimation)$covres)))   # inverse of the A matrix
  Fcomp   <- companionmatrix(Estimation)                      # Companion matrix
  
  #det     <- c_case                                           # constant and/or trends
  F1      <- t(QQ$Ft)                                         # make comparable to notes
  eps     <- ginv(invA) %*% t(residuals(Estimation))          # structural errors 
  nvar    <- Estimation$K                                     # number of endogenous variables
  nvarXeq <- nvar * nlag                                      # number of lagged endogenous per equation
  nvar_ex <- 0                                                # number of exogenous (excluding constant and trend)
  Y       <- QQ$Y                                             # left-hand side
  #X       <- QQ$X[,(1+det):(nvarXeq+det)]                    # right-hand side (no exogenous)
  nobs    <- nrow(Y)                                          # number of observations
  
  
  ## Compute historical decompositions
  
  # Contribution of each shock
  invA_big <- matrix(0,nvarXeq,nvar)
  invA_big[1:nvar,] <- invA
  Icomp <- cbind(diag(nvar), matrix(0,nvar,(nlag-1)*nvar))
  HDshock_big <- array(0, dim=c(nlag*nvar,nobs+1,nvar))
  HDshock <- array(0, dim=c(nvar,(nobs+1),nvar))
  
  for (j in 1:nvar){  # for each variable
    eps_big <- matrix(0,nvar,(nobs+1)) # matrix of shocks conformable with companion
    eps_big[j,2:ncol(eps_big)] <- eps[j,]
    for (i in 2:(nobs+1)){
      HDshock_big[,i,j] <- invA_big %*% eps_big[,i] + Fcomp %*% HDshock_big[,(i-1),j]
      HDshock[,i,j] <-  Icomp %*% HDshock_big[,i,j]
    } 
    
  } 
  
  HD.shock <- array(0, dim=c((nobs+nlag),nvar,nvar))   # [nobs x shock x var]
  
  for (i in 1:nvar){
    
    for (j in 1:nvar){
      HD.shock[,j,i] <- c(rep(NA,nlag), HDshock[i,(2:dim(HDshock)[2]),j])
    }
  }
  
  return(HD.shock)
  
  
  
  
}


#VARmakexy Biacnhi toolbox


VARmakexy <- function(DATA,lags,c_case){
  
  nobs <- nrow(DATA)
  
  #Y matrix 
  Y <- DATA[(lags+1):nrow(DATA),]
  Y <- DATA[-c(1:lags),]
  
  #X-matrix 
  if (c_case==0){
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    } 
  } else if(c_case==1){ #constant
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    }
    X <- cbind(matrix(1,(nobs-lags),1), X) 
  } else if(c_case==2){ # time trend and constant
    X <- NA
    for (jj in 0:(lags-1)){
      X <- rbind(DATA[(jj+1):(nobs-lags+jj),])
    }
    trend <- c(1:nrow(X))
    X <-cbind(matrix(1,(nobs-lags),1), t(trend))
  }
  A <- (t(X) %*% as.matrix(X)) 
  B <- (as.matrix(t(X)) %*% as.matrix(Y))
  
  Ft <- ginv(A) %*% B
  
  retu <- list(X=X,Y=Y, Ft=Ft)
  return(retu)
}

companionmatrix <- function (x) 
{
  if (!(class(x) == "varest")) {
    stop("\nPlease provide an object of class 'varest', generated by 'VAR()'.\n")
  }
  K <- x$K
  p <- x$p
  A <- unlist(Acoef(x))
  companion <- matrix(0, nrow = K * p, ncol = K * p)
  companion[1:K, 1:(K * p)] <- A
  if (p > 1) {
    j <- 0
    for (i in (K + 1):(K * p)) {
      j <- j + 1
      companion[i, j] <- 1
    }
  }
  return(companion)
}



HD <- VARhd(Estimation=var_x)
HD[,,2]


#alt models - vecm

summary(ca.jo(dataind_seasadj_endog_last, type = "eigen", ecdet = "trend", K = 3,spec = "transitory"))
vecm <- ca.jo(dataind_seasadj_endog[, c("reer", "cpi", "wacr", "asset","cf_cycle","wage","gscpi", "food growth", "oil growth")], type = "trace",ecdet = "trend", K = 3, spec = "transitory")
vecm.r1 <- cajorls(vecm, r = 3)
vecvar_x<-vec2var(vecm, r = 3)
predict(vecvar_x, n.ahead = 14, ci = 0.95)

LR <- matrix(NA, nrow = 4, ncol = 4)
LR[1, 2:4] <- 0
LR[2:4, 4] <- 0
LR <- matrix(NA, nrow = 9, ncol = 9)
LR[1, 2:4] <- 0
LR[2:4, 4] <- 0
View(LR)
View(LR)
svec <- SVEC(vecm, LR = LR, SR = amat, r = 1, lrtest = FALSE, boot = TRUE,
             + runs = 100)
svec <- SVEC(vecm, LR = LR, SR = amat, r = 3, lrtest = FALSE, boot = TRUE, runs = 100)
LR <- matrix(NA, nrow = 8, ncol = 8)
LR[1, 2:4] <- 0
LR[2:4, 4] <- 0
svec <- SVEC(vecm, LR = LR, SR = amat, r = 3, lrtest = FALSE, boot = TRUE, runs = 100)
vecm <- ca.jo(dataind_seasadj_endog[, c("reer", "cpi", "wacr", "asset","cf_cycle","wage","gscpi", "food growth", "oil growth")], type = "trace",ecdet = "trend", K = 3, spec = "transitory")
vecm.r1 <- cajorls(vecm, r = 3)
LR <- matrix(NA, nrow = 9, ncol = 9)
LR[1, 2:4] <- 0
LR[2:4, 4] <- 0
svec <- SVEC(vecm, LR = LR, SR = amat, r = 3, lrtest = FALSE, boot = TRUE, runs = 100)
summary(svec)
vec2var(vecm.r1, r = 3)
vec2var(vecm, r = 3)
vecvar_x<-vec2var(vecm, r = 3)
predict(vecvar_x, n.ahead = 8, ci = 0.95)
predict(res, n.ahead = 8, ci = 0.95)

x1 <- id.ngml(var_x)
x2 <- irf(x1, n.ahead = 10)
plot(x2)
x2 <- hd(x1, series = 2)
plot(x2)
x1 <- id.garch(v1)
summary(x1)
x1 <- id.garch(var_x)
summary(x1)
i1 <- irf(x1, n.ahead = 30)
plot(i1, scales = 'free_y')

###


### other models - phillips curve



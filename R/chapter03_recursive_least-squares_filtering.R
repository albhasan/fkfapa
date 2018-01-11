# chapter03

#---- MAKING ZEROTH-ORDER LEAST SQUARES FILTER RECURSIVE ----




#---- Properties of zeroth-order or one-state filter ----




# table 3.1
#---- Table 3.1 Sample measurement data ----
smd <- data.frame(k = 1:4, T_s = 0:3, xstar_k = c(1.2, 0.2, 2.9, 2.1))
K <- 1/(smd$k)
Res <- rep(NA, length(K))
xhat.vec <- rep(NA, length(K))
xhat <- 0
for(k in smd$k){
  Res[k] <- smd$xstar_k[k] - xhat
  xhat.vec[k] <- xhat + K[k] * Res[k]
  xhat <- xhat.vec[k]
}

xhat <- 100
k <- 1
(Res1 <- smd$xstar_k[k] - xhat)
(xhat1 <- xhat + K[k] * Res1)


#---- Listing 3.1 Simulation for testing zeroth-order recursive least-squares filter ----
# ACT     - Constants signal
# XS      - Noise
# XHERR   - Actual error in the estimate
# SP11    - P_k


# Zeroth-order recursive least-squares filter 
#
# @param TEE    Numeric. Vector of time of observations
# @param XS     Numeric. Vector of noisy observations
# @return       A data frame
zeroth_order_filter <- function(TEE, XS){
  XH <- 0
  XN <- 1:length(TEE)
  XK <- 1 / XN
  for(i in 1:length(XS)){
    RES <- XS[i] - XH[length(XH)]
    XH[length(XH) + 1] <- XH[length(XH)] + XK[i] * RES
  }
  df.res <- data.frame(TEE, XS, XH = XH[-length(XH)])
  df.res[nrow(df.res) + 1, ] <- c(NA, NA, XH[length(XH)])
  return(df.res)
}
  



TS <- 0.1 # sampling time in seconds
SIGNOISE <- 1
A0 <- 1
A1 <- 0
TEE <- seq(0, 10, TS)
XNOISE <- rnorm(length(TEE), 0, SIGNOISE)
ACT <- A0 + A1 * TEE
XS <- ACT + XNOISE
XN <- 1:length(TEE)
XK <- 1 / XN
SP11 <- SIGNOISE / sqrt(XN)
EPS <- 0.5 * A1 * TS * (XN - 1)
#XH <- 0
#for(i in 1:length(XS)){
#  RES <- XS[i] - XH[length(XH)]
#  XH[length(XH) + 1] <- XH[length(XH)] + XK[i] * RES
#}
res <- zeroth_order_filter(TEE, XS)
res$ACT <- c(ACT, NA)
res$XHERR <- res$ACT - res$XH
res$EPS <- c(EPS, NA)
res$SP11 <- c(SP11, NA)
res$SP11_neg <- -1 * res$SP11


#---- Fig. 3.2. Zeroth-order recursive least-squares filter is able to track zeroth-order polynomial plus noise ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$ACT, col = "True signal")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XS, col = "Measurement")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XH, col = "Estimate")) + 
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 3.2. Zeroth-order recursive least-squares filter is able to track zeroth-order polynomial plus noise")

#---- Fig. 3.3. Single-run simulation results agree with theoretical formula ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XHERR, col = "Simulation")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$SP11, col = "Theory")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$SP11_neg, col = "Theory")) + 
  ggplot2::ylab("Error in estimate") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 3.3. Single-run simulation results agree with theoretical formula")







SIGNOISE <- 0
A0 <- 1
A1 <- 2
TS <- 0.1
TEE <- seq(0, 10, TS)
XNOISE <- rnorm(length(TEE), 0, SIGNOISE)
ACT <- A0 + A1 * TEE
XS <- ACT + XNOISE
XN <- 1:length(TEE)
XK <- 1 / XN
SP11 <- SIGNOISE / sqrt(XN)
EPS <- 0.5 * A1 * TS * (XN - 1)
res <- zeroth_order_filter(TEE, XS)
res$ACT <- c(ACT, NA)
res$XHERR <- res$ACT - res$XH
res$EPS <- c(EPS, NA)
res$SP11 <- c(SP11, NA)
res$SP11_neg <- -1 * res$SP11


#---- Fig. 3.4. Zeroth-order recursive least-squares filter is unable to track first-order polynomial ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$ACT, col = "True signal")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XH, col = "Estimate")) + 
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 3.4. Zeroth-order recursive least-squares filter is unable to track first-order polynomial")



#---- Fig. 3.5. Simulation results and truncation error formula are in excellent agreement ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XHERR, col = "Simulation")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$SP11, col = "Theory")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$SP11_neg, col = "Theory")) + 
  ggplot2::ylab("Error in estimate") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 3.5. Simulation results and truncation error formula are in excellent agreement")



#---- Listing 3.2 Montecarlo simulation for testing zeroth-order recursive least-squares filter ----
TS <- 0.1
SIGNOISE <- 1
A0 <- 1
A1 <- 0
TEE <- seq(0, 10, TS)
ACT <- A0 + A1 * TEE
XN <- 1:length(TEE)
XK <- 1 / XN
SP11 <- SIGNOISE / sqrt(XN)
SP11_neg <- -1 * SP11
EPS <- 0.5 * A1 * TS * (XN - 1)
res.lt <- list()
for(k in 1:5){
  XNOISE <- rnorm(length(TEE), 0, SIGNOISE)
  XS <- ACT + XNOISE
  res <- zeroth_order_filter(TEE, XS)
  res$XHERR <- c(ACT, NA) - res$XH
  res.lt[[k]] <- res
}


#---- Fig. 3.6. Monte Carlo results lie within the theoretical bounds approximately 68% of the time ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = TEE, y = SP11, col = "Theory")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = TEE, y = SP11_neg, col = "Theory")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = TEE, y = res.lt[[1]]$XHERR[-length(TEE)], col = "Simulation")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = TEE, y = res.lt[[2]]$XHERR[-length(TEE)], col = "Simulation")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = TEE, y = res.lt[[3]]$XHERR[-length(TEE)], col = "Simulation")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = TEE, y = res.lt[[4]]$XHERR[-length(TEE)], col = "Simulation")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = TEE, y = res.lt[[5]]$XHERR[-length(TEE)], col = "Simulation")) + 
  ggplot2::ylab("Error in estimate") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 3.6. Monte Carlo results lie within the theoretical bounds approximately 68% of the time")




#---- Properties of first-order or one-state filter ----






#---- Table 3.2 Sample measurement data ----
smd <- data.frame(k = 1:4, T_s = 0:3, xstar_k = c(1.2, 0.2, 2.9, 2.1))
T_s <- 1
K1.vec <- rep(NA, nrow(smd))
K2.vec <- rep(NA, nrow(smd))
Res.vec <- rep(NA, nrow(smd))
xhat.vec <- rep(NA, nrow(smd))
xdhat.vec <- rep(NA, nrow(smd))
xhat <- 0
xdhat <- 0
for(k in smd$k){
  K1.vec[k] <- 2 * (2 * k - 1) / (k * (k + 1))
  K2.vec[k] <- 6  / (k * (k + 1) * T_s)
  Res.vec[k] <- smd$xstar[k] - xhat - xdhat * T_s
  xhat.vec[k] <- xhat + xdhat * T_s + K1.vec[k] * Res.vec[k]
  xdhat.vec[k] <- xdhat + K2.vec[k] * Res.vec[k]
  xhat <- xhat.vec[k]
  xdhat <- xdhat.vec[k]
}
res.df <- data.frame(smd$T_s)

# Re-computing Listing 2.1
AINV <- solve(matrix(c(nrow(smd), sum(smd$T_s), sum(smd$T_s), sum(smd$T_s^2)), ncol = 2))
B <- c(sum(smd$xstar_k), sum(smd$T_s * smd$xstar_k))
ANS <- AINV %*% B
bproc <- ANS[1] + ANS[2] * smd$T_s

#---- Fig. 3.7. First-order least-squares recursive and batch-processing least-squares filter yeilds the same answers after all measurements are taken ----
ggplot2::ggplot() +
  ggplot2::geom_point(mapping = ggplot2::aes(x = smd$T_s, y = smd$xstar_k, col = "Measurement")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = smd$T_s, y = xhat.vec, col = "Recursive"))  +
  ggplot2::geom_line(mapping = ggplot2::aes(x = smd$T_s, y = bproc, col = "Batch processing"))  +
  ggplot2::ylab("xhat") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 3.7. First-order least-squares recursive and batch-processing least-squares filter yeilds the same answers after all measurements are taken")



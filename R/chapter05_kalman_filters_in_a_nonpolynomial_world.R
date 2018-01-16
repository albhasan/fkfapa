#---- Listing 5.1 First-order polynomial kalman filter and sinusoidal measurement ----

ORDER <- 2
PHIS <- 0
TS <- 0.1
XH <- 0
XDH <- 0
SIGNOISE <- 1
PHI <- matrix(c(1, TS, 0, 1), nrow = 2, byrow = TRUE)
P <- matrix(c(99999999999, 0, 0, 99999999999), nrow = 2, byrow = TRUE)
IDNP <- diag(ORDER)
Q <- matrix(0, nrow = 2, ncol = 2)
Q[1,1] <- PHIS * TS^3 / 3
Q[1,2] <- PHIS * TS^2 / 2
Q[2,1] <- Q[1,2]
Q[2,2] <- PHIS * TS
RMAT <- SIGNOISE^2
HMAT <- t(c(1, 0))
#
TEE <- seq(from = 0, to = 20, by = TS)
XNOISE <- SIGNOISE * rnorm(length(TEE), 0, 1)
X <- sin(TEE)
XD <- cos(TEE)
XS <- X + XNOISE

XH.vec <- rep(NA, length(TEE))
XDH.vec <- rep(NA, length(TEE))

for(k in 1:length(TEE)){
  M <- PHI %*% P %*% t(PHI) + Q
  K <- M %*% t(HMAT) %*% solve((HMAT %*% M %*% t(HMAT) + RMAT))
  P <- (IDNP - K %*% HMAT) %*% M
  #XNOISE <- SIGNOISE * rnorm(1, mean = 0, sd = SIGNOISE)
  RES <- XS[k] - XH - TS * XDH
  XH <- XH + XDH %*% TS + K[1,1] %*% RES
  XDH <- XDH + K[2,1] %*% RES
  # TEE
  # X
  # XS
  XH.vec[k] <- XH
  # XD
  XDH.vec[k] <- XDH
}
res <- data.frame(TEE, X, XS, XH.vec, XD, XDH.vec)


#---- Fig. 5.1. Sinusoidal measurement is very noisy ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$X, col = "True signal")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XS, col = "Measurement")) + 
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 5.1. Sinusoidal measurement is very noisy")



#---- Fig. 5.2. First-order polynomial Kalman-filter has difficulty in tracking the sinusoidal signal ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$X, col = "True signal")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XH, col = "Estimate")) + 
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 5.2. First-order polynomial Kalman-filter has difficulty in tracking the sinusoidal signal")



#---- Fig. 5.3. First-order polynomial Kalman-filter does poorly in estimating the derivative of the true signal ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XD, col = "True x dot")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XDH, col = "Estimate")) + 
  ggplot2::ylab("xdot") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 5.3. First-order polynomial Kalman-filter does poorly in estimating the derivative of the true signal")





PHIS <- 10
Q <- matrix(0, nrow = 2, ncol = 2)
Q[1,1] <- PHIS * TS^3 / 3
Q[1,2] <- PHIS * TS^2 / 2
Q[2,1] <- Q[1,2]
Q[2,2] <- PHIS * TS
for(k in 1:length(TEE)){
  M <- PHI %*% P %*% t(PHI) + Q
  K <- M %*% t(HMAT) %*% solve((HMAT %*% M %*% t(HMAT) + RMAT))
  P <- (IDNP - K %*% HMAT) %*% M
  RES <- XS[k] - XH - TS * XDH
  XH <- XH + XDH %*% TS + K[1,1] %*% RES
  XDH <- XDH + K[2,1] %*% RES
  XH.vec[k] <- XH
  XDH.vec[k] <- XDH
}
res <- data.frame(TEE, X, XS, XH.vec, XD, XDH.vec)




#---- Fig. 5.4. Adding process noise yields better tracking at expense of noiser estimate ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$X, col = "True signal")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XH, col = "Estimate")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 5.4. Adding process noise yields better tracking at expense of noiser estimate")


#---- Fig. 5.5. First-order filter with process noise is now able to provide noisy estimate of derivative of true signal ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XD, col = "True x dot")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XDH, col = "Estimate")) + 
  ggplot2::ylab("xdot") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 5.5. First-order filter with process noise is now able to provide noisy estimate of derivative of true signal")

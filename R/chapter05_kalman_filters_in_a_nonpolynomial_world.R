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





#---- Listing 5.2 Second-order polynomial Kalman filter and sinusoidal measurement ----
ORDER <- 3
PHIS <- 0
TS <- 0.1
XH <- 0
XDH <- 0
XDDH <- 0
SIGNOISE <- 1
PHI <- matrix(c(1, TS, 0.5 * TS^2, 0, 1, TS, 0, 0, 1), nrow = 3, byrow = TRUE)
P <- matrix(c(99999999999, 0, 0, 0, 99999999999, 0, 0, 0, 99999999999), nrow = 3, byrow = TRUE)
IDNP <- diag(ORDER)
Q <- matrix(0, nrow = 3, ncol = 3)
Q[1,1] <- PHIS * TS^5 / 20
Q[1,2] <- PHIS * TS^4 / 8
Q[1,3] <- PHIS * TS^3 / 6
Q[2,1] <- Q[1,2]
Q[2,2] <- PHIS * TS^3 / 3
Q[2,3] <- PHIS * TS^2 / 2
Q[3,1] <- Q[1,3]
Q[3,2] <- Q[2,3]
Q[3,3] <- PHIS * TS
RMAT <- SIGNOISE^2
HMAT <- t(c(1, 0, 0))
#
TEE <- seq(from = 0, to = 20, by = TS)
XNOISE <- SIGNOISE * rnorm(length(TEE), 0, 1)
X <- sin(TEE)
XD <- cos(TEE)
XS <- X + XNOISE
#
XH.vec <- rep(NA, length(TEE))
XDH.vec <- rep(NA, length(TEE))
XDDH.vec <- rep(NA, length(TEE))
for(k in 1:length(TEE)){
  M <- PHI %*% P %*% t(PHI) + Q
  K <- M %*% t(HMAT) %*% solve((HMAT %*% M %*% t(HMAT) + RMAT))
  P <- (IDNP - K %*% HMAT) %*% M
  RES <- XS[k] - XH - TS * XDH - 0.5 * TS^2 * XDDH
  XH <- XH + XDH %*% TS + 0.5 * TS^2 * XDDH + K[1,1] %*% RES
  XDH <- XDH + XDDH * TS + K[2,1] %*% RES
  XDDH <- XDDH + K[3,1] %*% RES
  XH.vec[k] <- XH
  XDH.vec[k] <- XDH
  XDDH.vec[k] <- XDDH
}
res <- data.frame(TEE, X, XS, XH.vec, XD, XDH.vec, XDDH.vec)

#---- Fig. 5.6. Higher-order polynomial Kalman filter with zero process noise yields better but noiser estimates ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$X, col = "True signal")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XH, col = "Estimate")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 5.6. Higher-order polynomial Kalman filter with zero process noise yields better but noiser estimates")


#---- Fig. 5.7. Higher-order polynomial Kalman filter does a better job of tracking derivative of true signal ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XD, col = "True x dot")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XDH, col = "Estimate")) + 
  ggplot2::ylab("xdot") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 5.7. Higher-order polynomial Kalman filter does a better job of tracking derivative of true signal")



PHIS <- 10
Q <- matrix(0, nrow = 3, ncol = 3)
Q[1,1] <- PHIS * TS^5 / 20
Q[1,2] <- PHIS * TS^4 / 8
Q[1,3] <- PHIS * TS^3 / 6
Q[2,1] <- Q[1,2]
Q[2,2] <- PHIS * TS^3 / 3
Q[2,3] <- PHIS * TS^2 / 2
Q[3,1] <- Q[1,3]
Q[3,2] <- Q[2,3]
Q[3,3] <- PHIS * TS
XH.vec <- rep(NA, length(TEE))
XDH.vec <- rep(NA, length(TEE))
XDDH.vec <- rep(NA, length(TEE))
for(k in 1:length(TEE)){
  M <- PHI %*% P %*% t(PHI) + Q
  K <- M %*% t(HMAT) %*% solve((HMAT %*% M %*% t(HMAT) + RMAT))
  P <- (IDNP - K %*% HMAT) %*% M
  RES <- XS[k] - XH - TS * XDH - 0.5 * TS^2 * XDDH
  XH <- XH + XDH %*% TS + 0.5 * TS^2 * XDDH + K[1,1] %*% RES
  XDH <- XDH + XDDH * TS + K[2,1] %*% RES
  XDDH <- XDDH + K[3,1] %*% RES
  XH.vec[k] <- XH
  XDH.vec[k] <- XDH
  XDDH.vec[k] <- XDDH
}
res <- data.frame(TEE, X, XS, XH.vec, XD, XDH.vec, XDDH.vec)



#---- Fig. 5.8. Filter lag has been removed by the addition of process noise ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$X, col = "True signal")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XH, col = "Estimate")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 5.8. Filter lag has been removed by the addition of process noise")



#---- Fig. 5.9. Estimate of derivative has been improved by the addition of process noise ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XD, col = "True x dot")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XDH, col = "Estimate")) + 
  ggplot2::ylab("xdot") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 5.9. Estimate of derivative has been improved by the addition of process noise")




#---- Listing 5.3 New first-orer Kalman filter and sinusoidal measurement ----
ORDER <- 2
PHIS <- 0
W <- 1
A <- 1
TS <- 0.1
XH <- 0
XDH <- 0
SIGNOISE <- 1
PHI <- matrix(c(cos(W * TS), sin(W * TS)/W, -W * sin(W * TS), cos(W* TS)), nrow = 2, byrow = TRUE)
P <- matrix(c(99999999999, 0, 0, 99999999999), nrow = 2, byrow = TRUE)
IDNP <- diag(ORDER)
Q <- matrix(0, nrow = 2, ncol = 2)
Q[1,1] <- PHIS * (0.5 * W * TS - 0.25 * sin(2 * W * TS))
Q[1,2] <- 0.5 * PHIS * sin(W * TS)^2/W^2
Q[2,1] <- Q[2,1]
Q[2,2] <- PHIS * (0.5 * W * TS + 0.25 * sin(2 * W * TS))/W
RMAT <- SIGNOISE^2
HMAT <- t(c(1, 0))
#
TEE <- seq(from = 0, to = 20, by = TS)
XNOISE <- SIGNOISE * rnorm(length(TEE), 0, 1)
X <- A * sin(W * TEE)
XD <- A * W * cos(W * TEE)
XS <- X + XNOISE
XH.vec <- rep(NA, length(TEE))
XDH.vec <- rep(NA, length(TEE))

for(k in 1:length(TEE)){
  M <- PHI %*% P %*% t(PHI) + Q
  K <- M %*% t(HMAT) %*% solve((HMAT %*% M %*% t(HMAT) + RMAT))
  P <- (IDNP - K %*% HMAT) %*% M
  XHOLD <- XH
  RES <- XS[k] - XH * cos(W * TS) - sin(W * TS) * XDH/W
  XH <- cos(W * TS) * XH + XDH * sin(W * TS) / W + K[1,1] * RES
  XDH <- -W * sin(W*TS) * XHOLD + XDH * cos(W * TS) + K[2,1] * RES
  XH.vec[k] <- XH
  XDH.vec[k] <- XDH
}
res <- data.frame(TEE, X, XS, XH.vec, XD, XDH.vec)


#---- Fig. 5.10. New filter dramatically improves estimate of signal ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$X, col = "True signal")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XH, col = "Estimate")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 5.10. New filter dramatically improves estimate of signal")


#---- Fig. 5.11. New filter dramatically improves estimate of derivative signal ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XD, col = "True x dot")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XDH, col = "Estimate")) + 
  ggplot2::ylab("xdot") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 5.11. New filter dramatically improves estimate of derivative signal")




TS <- 0.1
W <- 1
PHI <- matrix(c(1, TS, -W^2 * TS, 1), nrow = 2, byrow = TRUE)
ORDER <- 2
PHIS <- 0
A <- 1
XH <- 0
XDH <- 0
SIGNOISE <- 1
P <- matrix(c(99999999999, 0, 0, 99999999999), nrow = 2, byrow = TRUE)
IDNP <- diag(ORDER)
Q <- matrix(0, nrow = 2, ncol = 2)
Q[1,1] <- PHIS * (0.5 * W * TS - 0.25 * sin(2 * W * TS))
Q[1,2] <- 0.5 * PHIS * sin(W * TS)^2/W^2
Q[2,1] <- Q[2,1]
Q[2,2] <- PHIS * (0.5 * W * TS + 0.25 * sin(2 * W * TS))/W
RMAT <- SIGNOISE^2
HMAT <- t(c(1, 0))
#
TEE <- seq(from = 0, to = 20, by = TS)
#XNOISE <- SIGNOISE * rnorm(length(TEE), 0, 1)
X <- A * sin(W * TEE)
XD <- A * W * cos(W * TEE)
XS <- X + XNOISE
XH.vec <- rep(NA, length(TEE))
XDH.vec <- rep(NA, length(TEE))

for(k in 1:length(TEE)){
  M <- PHI %*% P %*% t(PHI) + Q
  K <- M %*% t(HMAT) %*% solve((HMAT %*% M %*% t(HMAT) + RMAT))
  P <- (IDNP - K %*% HMAT) %*% M
  XHOLD <- XH
  RES <- XS[k] - XH * cos(W * TS) - sin(W * TS) * XDH/W
  XH <- cos(W * TS) * XH + XDH * sin(W * TS) / W + K[1,1] * RES
  XDH <- -W * sin(W*TS) * XHOLD + XDH * cos(W * TS) + K[2,1] * RES
  XH.vec[k] <- XH
  XDH.vec[k] <- XDH
}
res <- data.frame(TEE, X, XS, XH.vec, XD, XDH.vec)



#---- Fig. 5.12. Estimater of signal is worse when fundamental matrix is approximate ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$X, col = "True signal")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XH, col = "Estimate")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 5.12. Estimater of signal is worse when fundamental matrix is approximate")



#---- Fig. 5.13. Estimate of signal derivative is also worse when fundamental matrix is approximate
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XD, col = "True x dot")) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = res$TEE, y = res$XDH, col = "Estimate")) + 
  ggplot2::ylab("xdot") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 5.13. Estimate of signal derivative is also worse when fundamental matrix is approximate")



#TODO: missing figures

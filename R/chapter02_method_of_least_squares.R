# chapter02


#---- ZEROTH-ORDER OR ONE-STATE FILTER ----
R <- vector()

x_star <- c(1.2, 0.2, 2.9, 2.1)
x_caret <- sum(x_star)/length(x_star)
k <- 1:length(x_star)
T_s <- 1 # sampling time in seconds
t_vec <- (k - 1) * T_s

#---- Table 2.1 Measurement data for example ----
data.frame(t_vec, x_star)

#---- Table 2.2 Measurement data from table 2.1 expressed more mathematically ----
data.frame(k, t_vec, x_star)

#---- Fig. 2.1. Constant 1.6 is an unreasonable fit to measurement data ----
ggplot2::ggplot() +
  ggplot2::geom_point(mapping = ggplot2::aes(x = x_star, y = t_vec, col = "Measurement")) +
  ggplot2::geom_hline(yintercept = x_caret, mapping = ggplot2::aes(col = "Estimate")) +
  ggplot2::ylab("xhat") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.1. Constant 1.6 is an unreasonable fit to measurement data")
R <- append(R, sum((x_caret - x_star)^2))
R[length(R)]



#---- FIRST_ORDER OR TWO-STATE FILTER ----



#---- Listing 2.1 Solving for least-squares coefficients in the first-order fiter ----
TEE <- 0:3
X <- c(1.2, 0.2, 2.9, 2.1)
N <- 4
SUM <- rep(0, length(X))
for(I in 1:4){
  SUM[1] <- SUM[1] + TEE[I]
  SUM[2] <- SUM[2] + TEE[I] * TEE[I]
  SUM[3] <- SUM[3] + X[I]
  SUM[4] <- SUM[4] + TEE[I] * X[I]
}
A <- matrix(c(N, SUM[1], SUM[1], SUM[2]), ncol = 2)
#DET <- det(A)
AINV <- solve(A)
B <- c(SUM[3], SUM[4])
ANS <- AINV %*% B


T_s <- 1 # sampling time in seconds
k <- 1:length(X)
t_vec <- (k - 1) * T_s
est <- ANS[1, 1] + (ANS[2, 1] * t_vec)
#---- Fig. 2.2. Straight-line fit to data is better than constant fit ----
ggplot2::ggplot() +
  ggplot2::geom_point(mapping = ggplot2::aes(x = X, y = t_vec, col = "Measurement")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t_vec, y = est, col = "Estimate")) +
  ggplot2::ylab("xhat") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.2. Straight-line fit to data is better than constant fit")
R <- append(R, sum((est - X)^2))
R[length(R)]


#---- SECOND-ORDER OR THREE-STATE LEAST-SQUARES FILTER ----


#---- Listing 2.2 Solving for least-squared coeficients with three-state least-squares filter ----
T_s <- 1 # sampling time in seconds
TEE <- 0:3
X <- c(1.2, 0.2, 2.9, 2.1)
N <- 4
k <- 1:length(X)
t_vec <- (k - 1) * T_s
A <- matrix(rep(0, 9), nrow = 3)
A[1, 1] <- N
A[1, 2] <- sum((k - 1) * T_s)
A[1, 3] <- sum(((k - 1) * T_s)^2)
A[2, 2] <- A[1, 3]
A[2, 3] <- sum(((k - 1) * T_s)^3)
A[3, 3] <- sum(((k - 1) * T_s)^4)
A[2, 1] <- A[1, 2]
A[3, 1] <- A[1, 3]
A[3, 2] <- A[2, 3]
#DET <- det(A)
AINV <- solve(A)
B <- c(sum(X), sum((t_vec * X)), sum((t_vec^2 * X)))
ANS <- AINV %*% B

x_caret <- ANS[1] + (ANS[2] * (k - 1) * T_s) + (ANS[3] * ((k - 1) * T_s)^2)
#---- Fig. 2.3. Parabolic fit to data is pretty good, too ----
ggplot2::ggplot() +
  ggplot2::geom_point(mapping = ggplot2::aes(x = X, y = t_vec, col = "Measurement")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t_vec, y = x_caret, col = "Estimate")) +
  ggplot2::ylab("xhat") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.3. Parabolic fit to data is pretty good, too")
R <- append(R, sum((x_caret - X)^2))
R[length(R)]


#---- THIRD-ORDER SYSTEM ----


T_s <- 1 # sampling time in seconds
TEE <- 0:3
X <- c(1.2, 0.2, 2.9, 2.1)
N <- 4
k <- 1:length(X)
t_vec <- (k - 1) * T_s
A <- matrix(rep(0, 16), nrow = 4)
A[1, 1] <- N
A[1, 2] <- sum((k - 1) * T_s)
A[1, 3] <- sum(((k - 1) * T_s)^2)
A[1, 4] <- sum(((k - 1) * T_s)^3)
A[2, 2] <- A[1, 3]
A[2, 3] <- A[1, 4]
A[2, 4] <- sum(((k - 1) * T_s)^4)
A[3, 3] <- A[2, 4]
A[3, 4] <- sum(((k - 1) * T_s)^5)
A[4, 4] <- sum(((k - 1) * T_s)^6)
A[2, 1] <- A[1, 2]
A[3, 1] <- A[1, 3]
A[4, 1] <- A[1, 4]
A[3, 2] <- A[2, 3]
A[4, 2] <- A[2, 4]
A[4, 3] <- A[3, 4]
B <- c(sum(X), sum((t_vec * X)), sum((t_vec^2 * X)), sum((t_vec^3 * X)))
AINV <- solve(A)
ANS <- AINV %*% B


x_caret <- ANS[1] + (ANS[2] * (k - 1) * T_s) + (ANS[3] * ((k - 1) * T_s)^2) + (ANS[4] * ((k - 1) * T_s)^3)
t_dense = seq(from = t_vec[1], to = t_vec[length(t_vec)], by = 0.1)
x_dense = ANS[1] + (ANS[2] * t_dense) + (ANS[3] * (t_dense)^2) + (ANS[4] * (t_dense^3))

#---- Fig. 2.3. Parabolic fit to data is pretty good, too ----
ggplot2::ggplot() +
  ggplot2::geom_point(mapping = ggplot2::aes(x = t_vec, y = X, col = "Measurement")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t_vec, y = x_caret, col = "Estimate")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t_dense, y = x_dense, col = "Dense Estimate")) +
  ggplot2::ylab("xhat") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.3. Parabolic fit to data is pretty good, too")
R <- append(R, sum((x_caret - X)^2))
R[length(R)]

#---- Table 2.3 Residual decreases as order of least-squaredpolynomial increses ----
data.frame(System_order = 0:(length(R)-1), R = round(R, 2))


#---- Experiments with Zeroth-order or One-state filter ----



#---- Listing 2.3 One-state filter for extracting signal from measurement ----

# One-state filter for extracting signal from measurement
#
# @param t_vec  Numeric. Vector of time of observations
# @param X      Numeric. Vector of noisy observations
# @param X1     Numeric. Vector of true observations
# @return       A data frame
one_state_filter <- function(t_vec, X, X1){
  A <- matrix(length(t_vec) -1)
  B <- matrix(sum(X))
  AINV <- 1/A
  ANS <- AINV %*% B
  XHAT <- rep(ANS[1, 1], length(t_vec))
  ERRX2 <- sum((X1 - XHAT)^2)
  ERRXP2 <- sum((X - XHAT)^2)
  return(data.frame(
    TEE = t_vec, 
    X1 = X1, 
    X = X, 
    XHAT = XHAT, 
    ERRX = X1 - XHAT, 
    ERRXP = X - XHAT,
    ERRX2 <- (X1 - XHAT)^2,
    ERRXP2 <- (X - XHAT)^2,
    SUMPZ1 <- sum(ERRX2),
    SUMPZ2 <- sum(ERRXP2)
  ))
}

# N <- 0
TS <- 0.1 # sampling time in seconds
t_vec <- seq(0, 10, TS)
SIGNOISE <- 1
XNOISE <- rnorm(length(t_vec), 0, SIGNOISE)
X1 <- rep(1, length(t_vec))
X <- X1 + XNOISE
# SUM3 <- sum(X)
# NMAX <- length(t_vec) -1
# N <- NMAX
# 
# A <- matrix(N)
# B <- matrix(SUM3)
# AINV <- 1/A
# ANS <- AINV %*% B
# 
# i <- 1
# SUMPZ1 <- 0
# SUMPZ2 <- 0
# res <- list()
# for(i in 1:NMAX){
#   TEE <- TS * (i-1)
#   XHAT <- ANS[1,1]
#   ERRX <- X1[i] - XHAT
#   ERRXP <- X[i] - XHAT
#   ERRX2 <- (X1[i] - XHAT)^2
#   ERRXP2 <- (X[i] - XHAT)^2
#   SUMPZ1 <- ERRX2 + SUMPZ1
#   SUMPZ2 <- ERRXP2 + SUMPZ2
#   res[[i]] <- c(TEE, X1[i], X[i], XHAT, ERRX, ERRXP, SUMPZ1, SUMPZ2)
# }
# res.df <- as.data.frame(do.call(rbind, res))
# colnames(res.df) <- c("TEE", "X1", "X", "XHAT", "ERRX", "ERRXP", "SUMPZ1", "SUMPZ2")
res.df <- one_state_filter(t_vec, X, X1)

compare_res <- data.frame(
  time = res.df$TEE, 
  trueSignal = res.df$X1, 
  obsSignal = res.df$X, 
  zerothOrder = res.df$XHAT
)

#---- Fig. 2.5. One-state, least-squares filter smoothesnoise measurements ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$X, col = "Measurements")) +
  ggplot2::geom_hline(yintercept = res.df$XHAT[1]) + 
  ggplot2::annotate("text", 1, res.df$XHAT[1], vjust = -1, label = "Estimate") + 
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.5. One-state, least-squares filter smoothesnoise measurements") + 
  ggplot2::coord_cartesian(ylim = c(-2, 4))

#---- Fig. 2.6. One-state filteryields near perfect estimate of constant signal ----
ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = res.df$XHAT[1]) +
  ggplot2::geom_hline(yintercept = res.df$X1[1]) + 
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.6. One-state filter yields near perfect estimate of constant signal") + 
  ggplot2::coord_cartesian(ylim = c(0.0, 1.4))

res.df["s_e"] <- res.df$X1 - res.df$XHAT
res.df["m_e"] <- res.df$X - res.df$XHAT

#---- Fig. 2.7 Estimation errors are nearly zero for one-state least-squared filter ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$s_e, col = "Signal & estiamte")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$m_e, col = "Measurement & estiamte")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.7 Estimation errors are nearly zero for one-state least-squared filter") + 
  ggplot2::coord_cartesian(ylim = c(-4, 4))



#---- Listing 2.3 (modified) One-state filter for extracting signal from measurement ----
# N <- 0
TS <- 0.1 # sampling time in seconds
t_vec <- seq(0, 10, TS)
SIGNOISE <- 5                     # modified
XNOISE <- rnorm(length(t_vec), 0, SIGNOISE)
X1 <- rep(t_vec + 3)              # modified
X <- X1 + XNOISE
# SUM3 <- sum(X)
# NMAX <- length(t_vec) -1
# N <- NMAX
# 
# A <- matrix(N)
# B <- matrix(SUM3)
# AINV <- 1/A
# ANS <- AINV %*% B
# 
# i <- 1
# SUMPZ1 <- 0
# SUMPZ2 <- 0
# res <- list()
# for(i in 1:NMAX){
#   TEE <- TS * (i-1)
#   XHAT <- ANS[1,1]
#   ERRX <- X1[i] - XHAT
#   ERRXP <- X[i] - XHAT
#   ERRX2 <- (X1[i] - XHAT)^2
#   ERRXP2 <- (X[i] - XHAT)^2
#   SUMPZ1 <- ERRX2 + SUMPZ1
#   SUMPZ2 <- ERRXP2 + SUMPZ2
#   res[[i]] <- c(TEE, X1[i], X[i], XHAT, ERRX, ERRXP, SUMPZ1, SUMPZ2)
# }
# res.df <- as.data.frame(do.call(rbind, res))
# colnames(res.df) <- c("TEE", "X1", "X", "XHAT", "ERRX", "ERRXP", "SUMPZ1", "SUMPZ2")
res.df <- one_state_filter(t_vec, X, X1)


#---- Fig. 2.8. Zeroth-order least-squares filter does not capture upward trend ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$XHAT, col = "Estiamte")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$X, col = "Measurement")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.8. Zeroth-order least-squares filter does not capture upward trend") + 
  ggplot2::coord_cartesian(ylim = c(-10, 30))

#---- Fig. 2.9. Zeroth-order least-squares filter cannot estimate slope of actual sign ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$XHAT, col = "Estiamte")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$X1, col = "Actual")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.9. Zeroth-order least-squares filter cannot estimate slope of actual sign") + 
  ggplot2::coord_cartesian(ylim = c(-10, 30))

res.df["s_e"] <- res.df$X1 - res.df$XHAT
res.df["m_e"] <- res.df$X - res.df$XHAT

#---- Fig. 2.10 Errors in the estimate of the signal appear to grow with time ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$s_e, col = "Signal & estiamte")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$m_e, col = "Measurement & estiamte")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.10 Errors in the estimate of the signal appear to grow with time") + 
  ggplot2::coord_cartesian(ylim = c(-20, 20))

sum(res.df["s_e"]^2)
sum(res.df["m_e"]^2)



#---- Experiments with First-order or Two-state filter ----


#---- Listing 2.4 Two-state least-squares statements invoke the absoft random number generator on the macintosch ----
# NOTE: Run several times to see the book's results


# Two-state filter for extracting signal from measurement
#
# @param t_vec  Numeric. Vector of time of observations
# @param X      Numeric. Vector of noisy observations
# @param X1     Numeric. Vector of true observations
# @param XD     Numeric. First derivative of the model
# @return       A data frame
two_state_filter <- function(t_vec, X, X1, XD){
  A <- matrix(c(length(t_vec), sum(t_vec), sum(t_vec), sum(t_vec^2)), nrow = 2, byrow = TRUE)
  B <- matrix(c(sum(X), sum(t_vec * X)), ncol = 1)  
  AINV <- solve(A)
  ANS <- AINV %*% B
  XHAT <- ANS[1, 1] + ANS[2, 1] * t_vec
  XDHAT <- ANS[2, 1]
  ERRX <- X1 - XHAT
  ERRXD <- XD - XDHAT
  ERRXP <- X - XHAT
  ERRX2 <- (X1 - XHAT)^2
  ERRXP2 <-(X -XHAT)^2
  SUMPZ1 <- sum(ERRX2)
  SUMPZ2 <- sum(ERRXP2)
  return(data.frame(TEE = t_vec, X1, X, XHAT, ERRX, ERRXD, SUMPZ1, SUMPZ2))
}



#TS <- 0.1
#t_vec <- seq(0, 10, TS)
#SIGNOISE <- 1
#TEE <- t_vec # 0:(length(t_vec) - 1)
#N <- 1:length(t_vec)
#XNOISE <- rnorm(length(N), 0, SIGNOISE)
#X1 <- rep(1, length(N))
XD <- rep(0, length(compare_res$obsSignal)) # first derivative???
#X <- compare_res$obsSignal # X1 + XNOISE
# SUM1 <- sum(TEE)
# SUM2 <- sum(TEE^2)
# SUM3 <- sum(X)
# SUM4 <- sum(TEE * X)
# SUMPZ1 <- 0
# SUMPZ2 <- 0
# NMAX <- max(N)
# A <- matrix(c(N[length(N)], SUM1, SUM1, SUM2), nrow = 2, byrow = TRUE)
# B <- matrix(c(SUM3, SUM4), ncol = 1)
# #DET <- det(A)
# AINV <- solve(A)
# ANS <- AINV %*% B
# #TEE <- 0.1 * (1:NMAX)
# XHAT <- ANS[1, 1] + ANS[2, 1] * TEE
# XDHAT <- ANS[2, 1]
# ERRX <- X1 - XHAT
# ERRXD <- XD - XDHAT
# ERRXP <- X - XHAT
# ERRX2 <- (X1 - XHAT)^2
# ERRXP2 <-(X -XHAT)^2
# SUMPZ1 <- cumsum(ERRX2)
# SUMPZ2 <- cumsum(ERRXP2)
# res.df <- data.frame(TEE, X1, X, XHAT, ERRX, ERRXD, SUMPZ1, SUMPZ2)
#C <- A %*% B

res.df <- two_state_filter(compare_res$time, compare_res$obsSignal, compare_res$obsSignal, XD)
compare_res$firstOrder <- res.df$XHAT

#----- Fig. 2.11 First-order filter has trouble in estimating zeroth-order signal ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$X1, col = "Actual")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$XHAT, col = "Estimate")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.11 First-order filter has trouble in estimating zeroth-order signal") #+ 
  #ggplot2::coord_cartesian(ylim = c(0, 1.4))



#---- Fig. 2.12 Errors in the estimate of the signal and its derivative are not too large -----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRX, col = "Actual x\nminus estimate")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRXD, col = "Actual x Dot\nminus estimate ")) +
  ggplot2::ylab("Differences") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.12 Errors in the estimate of the signal and its derivative are not too large") #+ 
  #ggplot2::coord_cartesian(ylim = c(-0.2, 0.2))
# sum of the squares of the difference between the measurement and estiamte
sum((res.df$X - res.df$XHAT)^2)



#---- Listing 2.4 (modified) Two-state least-squares statements invoke the absoft random number generator on the macintosch ----
# NOTE: Run several times to see the book's results
# TS <- 0.1
# t_vec <- seq(0, 10, TS)
SIGNOISE <- 5                                     # modified
# TEE <- t_vec # 0:(length(t_vec) - 1)
# 
# N <- 1:length(t_vec)
XNOISE <- rnorm(length(compare_res$time), 0, SIGNOISE)
X1 <- compare_res$time + 3 # TEE + 3             # modified
XD <- rep(1, length(compare_res$time))           # modified
X <- X1 + XNOISE
# SUM1 <- sum(TEE)
# SUM2 <- sum(TEE^2)
# SUM3 <- sum(X)
# SUM4 <- sum(TEE * X)
# SUMPZ1 <- 0
# SUMPZ2 <- 0
# NMAX <- max(N)
# A <- matrix(c(N[length(N)], SUM1, SUM1, SUM2), nrow = 2, byrow = TRUE)
# B <- matrix(c(SUM3, SUM4), ncol = 1)
# #DET <- det(A)
# AINV <- solve(A)
# ANS <- AINV %*% B
# 
# #TEE <- 0.1 * (1:NMAX)
# XHAT <- ANS[1, 1] + ANS[2, 1] * TEE
# XDHAT <- ANS[2, 1]
# ERRX <- X1 - XHAT
# ERRXD <- XD - XDHAT
# ERRXP <- X - XHAT
# ERRX2 <- (X1 - XHAT)^2
# ERRXP2 <-(X -XHAT)^2
# SUMPZ1 <- cumsum(ERRX2)
# SUMPZ2 <- cumsum(ERRXP2)
# res.df <- data.frame(TEE, X1, X, XHAT, ERRX, ERRXD, SUMPZ1, SUMPZ2)
# C <- A %*% B
res.df <- two_state_filter(compare_res$time, X, X1, XD)


#----- Fig. 2.13 First-order filter does a much better job of estimating first-order signal than does a zeroth-order filter ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$X1, col = "Actual")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$XHAT, col = "Estimate")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.13 First-order filter does a much better job of estimating first-order signal than does a zeroth-order filter") + 
  ggplot2::coord_cartesian(ylim = c(0, 14))

#---- Fig. 2.14 First-order filter is able to estimate derivative and first-order signal accurately ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRX, col = "Actual x\nminus estimate")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRXD, col = "Actual x Dot\nminus estimate ")) +
  ggplot2::ylab("Differences") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.14 First-order filter is able to estimate derivative and first-order signal accurately") + 
  ggplot2::coord_cartesian(ylim = c(-1, 1))
sum((res.df$X - res.df$XHAT)^2)



#---- Listing 2.4 (modified) Two-state least-squares statements invoke the absoft random number generator on the macintosch ----
# NOTE: Run several times to see the book's results
# TS <- 0.1
# t_vec <- seq(0, 10, TS)
SIGNOISE <- 50                                     # modified
TEE <- compare_res$time # t_vec
# 
N <- 1:length(TEE)
XNOISE <- rnorm(length(N), 0, SIGNOISE)
X1 <- 5 * TEE^2 - 2 * TEE + 2                    # modified
XD <- 10 * TEE - 2                               # modified
X <- X1 + XNOISE
# SUM1 <- sum(TEE)
# SUM2 <- sum(TEE^2)
# SUM3 <- sum(X)
# SUM4 <- sum(TEE * X)
# SUMPZ1 <- 0
# SUMPZ2 <- 0
# NMAX <- max(N)
# A <- matrix(c(N[length(N)], SUM1, SUM1, SUM2), nrow = 2, byrow = TRUE)
# B <- matrix(c(SUM3, SUM4), ncol = 1)
# #DET <- det(A)
# AINV <- solve(A)
# ANS <- AINV %*% B
# 
# #TEE <- 0.1 * (1:NMAX)
# XHAT <- ANS[1, 1] + ANS[2, 1] * TEE
# XDHAT <- ANS[2, 1]
# ERRX <- X1 - XHAT
# ERRXD <- XD - XDHAT
# ERRXP <- X - XHAT
# ERRX2 <- (X1 - XHAT)^2
# ERRXP2 <-(X -XHAT)^2
# SUMPZ1 <- cumsum(ERRX2)
# SUMPZ2 <- cumsum(ERRXP2)
# res.df <- data.frame(TEE, X1, X, XHAT, ERRX, ERRXD, SUMPZ1, SUMPZ2)
# C <- A %*% B
res.df <- two_state_filter(compare_res$time, X, X1, XD)


#---- Fig. 2.15 First-order filter attempts to track second-order measurements ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$X, col = "Measurement")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$XHAT, col = "Estimate")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.15 First-order filter attempts to track second-order measurements") + 
  ggplot2::coord_cartesian(ylim = c(-100, 500))

#---- Fig. 2.16 On the average first-order filter estimates second-order signal ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$X1, col = "Actual")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$XHAT, col = "Estimate")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.16 On the average first-order filter estimates second-order signal") + 
  ggplot2::coord_cartesian(ylim = c(-100, 500))

#---- Fig. 2.17 Large estimation error result when first-order filter attempts to track second-order signal ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRX, col = "Actual x\nminus estimate")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRXD, col = "Actual x Dot\nminus estimate ")) +
  ggplot2::ylab("Differences") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.17 Large estimation error result when first-order filter attempts to track second-order signal") + 
  ggplot2::coord_cartesian(ylim = c(-40, 100))
sum((res.df$X1 - res.df$XHAT)^2)
sum((res.df$X - res.df$XHAT)^2)




#---- Experiments with Second-order or Three-state filter ----



#---- Listing 2.5 Three-state least-squares filter for extracting signal from measurement ----


# Three-state least-squares filter for extracting signal from measurement
#
# @param t_vec  Numeric. Vector of time of observations
# @param X      Numeric. Vector of noisy observations
# @param X1     Numeric. Vector of true observations
# @param XD     Numeric. First derivative of the model
# @param XD     Numeric. Second derivative of the model
# @return       A data frame
three_state_filter <- function(t_vec, X, X1, XD, XDD){
  A <- matrix(c(length(t_vec), sum(t_vec), sum(t_vec^2), sum(t_vec), 
                sum(t_vec^2), sum(t_vec^3), sum(t_vec^2), sum(t_vec^3), 
                sum(t_vec^4)),  ncol = 3, byrow = 3)
  B <- matrix(c(sum(X), sum(t_vec * X), sum(t_vec^2 * X)), ncol = 1, byrow = 3)
  AINV <- solve(A)
  ANS <- AINV %*% B
  XHAT <- ANS[1, 1] + ANS[2, 1] * t_vec + ANS[3, 1] * t_vec^2
  XDHAT <- ANS[2, 1] + 2 * ANS[3, 1] * t_vec
  XDDHAT <- 2 * ANS[3, 1]
  ERRX <- X1 - XHAT
  ERRXD <- XD - XDHAT
  ERRXDD <- XDD - XDDHAT
  ERRXP <- X - XHAT
  ERRX2 <- (X1 - XHAT)^2
  ERRXP2 <- (X - XHAT)^2
  SUMPZ1 <- sum(ERRX2)
  SUMPZ2 <- sum(ERRXP2)
  data.frame(TEE = t_vec, X, X1, XHAT, ERRX, ERRXD, ERRXDD, SUMPZ1, SUMPZ2)
}


# SIGNOISE <- 1
# TS <- 0.1
# N <- 0
# SUM1 <- SUM2 <- SUM3 <- SUM4 <- SUM5 <- SUM6 <- SUM7 <- SUMPZ1 <- SUMPZ2 <- 0
# TEE <- seq(0, 10, TS)
# N <- length(TEE)
# XNOISE <- rnorm(N, mean = 0, sd = SIGNOISE)
# X1 <- 1
XD <- rep(0, length(compare_res$time))
XDD <- rep(0, length(compare_res$time))
# X <- compare_res$obsSignal # X1 + XNOISE
# SUM1 <- sum(TEE)
# SUM2 <- sum(TEE^2)  
# SUM3 <- sum(X)
# SUM4 <- sum(TEE * X)
# SUM5 <- sum(TEE^3)
# SUM6 <- sum(TEE^4)
# SUM7 <- sum(TEE^2 * X)
# NMAX <- max(N)
# A <- matrix(c(N, SUM1, SUM2, SUM1, SUM2, SUM5, SUM2, SUM5, SUM6),  ncol = 3, byrow = 3)
# B <- matrix(c(SUM3, SUM4, SUM7), ncol = 1, byrow = 3)
# #DET <- det(A)
# AINV <- solve(A)
# ANS <- AINV %*% B
# #TEE <- 0.1 * (0:NMAX)
# XHAT <- ANS[1, 1] + ANS[2, 1] * TEE + ANS[3, 1] * TEE^2
# XDHAT <- ANS[2, 1] + 2 * ANS[3, 1] * TEE
# XDDHAT <- 2 * ANS[3, 1]
# ERRX <- X1 - XHAT
# ERRXD <- XD - XDHAT
# ERRXDD <- XDD - XDDHAT
# ERRXP <- X - XHAT
# ERRX2 <- (X1 - XHAT)^2
# ERRXP2 <- (X - XHAT)^2
# SUMPZ1 <- sum(ERRX2)
# SUMPZ2 <- sum(ERRXP2)
# res.df <- data.frame(TEE, X1, X, XHAT, ERRX, ERRXD, ERRXDD, SUMPZ1, SUMPZ2)
res.df <- three_state_filter(compare_res$time, compare_res$obsSignal, compare_res$trueSignal, XD, XDD)

compare_res$secondOrder <- res.df$XHAT

#---- Fig. 2.18 Second-order filter estimates signal is parabola even though it is a constant ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$X1, col = "Actual")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$XHAT, col = "Estimate")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.18 Second-order filter estimates signal is parabola even though it is a constant") + 
  ggplot2::coord_cartesian(ylim = c(0, 1.4))

#---- Fig. 2.19 Estimation errors between estimates and states of signal are not terrible when the order of filter is too high ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRX, col = "Actual x\nminus estimate")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRXD, col = "Actual x Dot\nminus estimate ")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRXDD, col = "Actual x Double Dot\nminus estimate ")) +
  ggplot2::ylab("Differences") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.19 Estimation errors between estimates and states of signal are not terrible when the order of filter is too high") + 
  ggplot2::coord_cartesian(ylim = c(-0.1, 0.4))
sum((res.df$X1 - res.df$XHAT)^2)
sum((res.df$X - res.df$XHAT)^2)


#---- Listing 2.5 (modified) Three-state least-squares filter for extracting signal from measurement ----
SIGNOISE <- 5                   # modified
# TS <- 0.1
# N <- 0
# SUM1 <- SUM2 <- SUM3 <- SUM4 <- SUM5 <- SUM6 <- SUM7 <- SUMPZ1 <- SUMPZ2 <- 0
TEE <- compare_res$time # seq(0, 10, TS)
N <- length(TEE)
XNOISE <- rnorm(N, mean = 0, sd = SIGNOISE)
X1 <- TEE + 3                   # modified
XD <- rep(1, length(compare_res$time))                         # modified
XDD <- rep(0, length(compare_res$time))
# X <- X1 + XNOISE
# SUM1 <- sum(TEE)
# SUM2 <- sum(TEE^2)  
# SUM3 <- sum(X)
# SUM4 <- sum(TEE * X)
# SUM5 <- sum(TEE^3)
# SUM6 <- sum(TEE^4)
# SUM7 <- sum(TEE^2 * X)
# NMAX <- max(N)
# A <- matrix(c(N, SUM1, SUM2, SUM1, SUM2, SUM5, SUM2, SUM5, SUM6),  ncol = 3, byrow = 3)
# B <- matrix(c(SUM3, SUM4, SUM7), ncol = 1, byrow = 3)
# #DET <- det(A)
# AINV <- solve(A)
# ANS <- AINV %*% B
# #TEE <- 0.1 * (0:NMAX)
# XHAT <- ANS[1, 1] + ANS[2, 1] * TEE + ANS[3, 1] * TEE^2
# XDHAT <- ANS[2, 1] + 2 * ANS[3, 1] * TEE
# XDDHAT <- 2 * ANS[3, 1]
# ERRX <- X1 - XHAT
# ERRXD <- XD - XDHAT
# ERRXDD <- XDD - XDDHAT
# ERRXP <- X - XHAT
# ERRX2 <- (X1 - XHAT)^2
# ERRXP2 <- (X - XHAT)^2
# SUMPZ1 <- sum(ERRX2)
# SUMPZ2 <- sum(ERRXP2)
# res.df <- data.frame(TEE, X1, X, XHAT, ERRX, ERRXD, ERRXDD, SUMPZ1, SUMPZ2)
res.df <- three_state_filter(compare_res$time, X, X1, XD, XDD)


#---- Fig. 2.20 Second-order filter attempts to fit first-order signal with a parabola ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$X1, col = "Actual")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$XHAT, col = "Estimate")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.20 Second-order filter attempts to fit first-order signal with a parabola") #+ 
  #ggplot2::coord_cartesian(ylim = c(0, 14))

# Fig. 2.21 Second fit to first-order signal yields larger errors than first-order fit
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRX, col = "Actual x\nminus estimate")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRXD, col = "Actual x Dot\nminus estimate ")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRXDD, col = "Actual x Double Dot\nminus estimate ")) +
  ggplot2::ylab("Differences") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.21 Second fit to first-order signal yields larger errors than first-order fit") #+ 
  #ggplot2::coord_cartesian(ylim = c(-0.5, 2))
sum((res.df$X1 - res.df$XHAT)^2)
sum((res.df$X - res.df$XHAT)^2)





#---- Listing 2.5 (modified) Three-state least-squares filter for extracting signal from measurement ----
SIGNOISE <- 50                              # modified
# TS <- 0.1
# N <- 0
# SUM1 <- SUM2 <- SUM3 <- SUM4 <- SUM5 <- SUM6 <- SUM7 <- SUMPZ1 <- SUMPZ2 <- 0
TEE <- compare_res$time # seq(0, 10, TS)
N <- length(TEE)
XNOISE <- rnorm(N, mean = 0, sd = SIGNOISE)
X1 <- 5 * TEE^2 - 2 * TEE + 2               # modified
XD <- 10 * TEE - 2                          # modified
XDD <- rep(10, length(compare_res$time))                                   # modified
# X <- X1 + XNOISE
# SUM1 <- sum(TEE)
# SUM2 <- sum(TEE^2)  
# SUM3 <- sum(X)
# SUM4 <- sum(TEE * X)
# SUM5 <- sum(TEE^3)
# SUM6 <- sum(TEE^4)
# SUM7 <- sum(TEE^2 * X)
# NMAX <- max(N)
# A <- matrix(c(N, SUM1, SUM2, SUM1, SUM2, SUM5, SUM2, SUM5, SUM6),  ncol = 3, byrow = 3)
# B <- matrix(c(SUM3, SUM4, SUM7), ncol = 1, byrow = 3)
# #DET <- det(A)
# AINV <- solve(A)
# ANS <- AINV %*% B
# #TEE <- 0.1 * (0:NMAX)
# XHAT <- ANS[1, 1] + ANS[2, 1] * TEE + ANS[3, 1] * TEE^2
# XDHAT <- ANS[2, 1] + 2 * ANS[3, 1] * TEE
# XDDHAT <- 2 * ANS[3, 1]
# ERRX <- X1 - XHAT
# ERRXD <- XD - XDHAT
# ERRXDD <- XDD - XDDHAT
# ERRXP <- X - XHAT
# ERRX2 <- (X1 - XHAT)^2
# ERRXP2 <- (X - XHAT)^2
# SUMPZ1 <- sum(ERRX2)
# SUMPZ2 <- sum(ERRXP2)
# res.df <- data.frame(TEE, X1, X, XHAT, ERRX, ERRXD, ERRXDD, SUMPZ1, SUMPZ2)
res.df <- three_state_filter(compare_res$time, X, X1, XD, XDD)


#---- Fig. 2.22 Second-order filter provides near perfect estimates of second-order signal ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$X1, col = "Actual")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$XHAT, col = "Estimate")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.22 Second-order filter provides near perfect estimates of second-order signal") #+ 
  #ggplot2::coord_cartesian(ylim = c(0, 500))

#---- Fig. 2.23 Error in the estimates of all states of second-order filter against second-order signal are better than all other filter fits ----
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRX, col = "Actual x\nminus estimate")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRXD, col = "Actual x Dot\nminus estimate ")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$ERRXDD, col = "Actual x Double Dot\nminus estimate ")) +
  ggplot2::ylab("Differences") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.23 Error in the estimates of all states of second-order filter against second-order signal are better than all other filter fits") #+ 
  #ggplot2::coord_cartesian(ylim = c(-0.5, 20))
sum((res.df$X1 - res.df$XHAT)^2)
sum((res.df$X - res.df$XHAT)^2)



#---- Comparison of filters



#---- Fig. 2.24 Zeroth-order least-squares filter best tracks zeroth-order measurement
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = compare_res$time, y = compare_res$trueSignal, col = "True signal")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = compare_res$time, y = compare_res$zerothOrder, col = "Zeroth-Order")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = compare_res$time, y = compare_res$firstOrder, col = "First-Order")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = compare_res$time, y = compare_res$secondOrder, col = "Second-Order")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.24 Zeroth-order least-squares filter best tracks zeroth-order measurement") + 
  ggplot2::coord_cartesian(ylim = c(0, 1.4))


#---- Prepare Fig. 2.25 
TS <- 0.1 # sampling time in seconds
t_vec <- seq(0, 10, TS)
SIGNOISE <- 1
XNOISE <- rnorm(length(t_vec), 0, SIGNOISE)
X1 <- t_vec + 3
X <- X1 + 5 * XNOISE
XD <- rep(1, length(t_vec))
XDD <- rep(0, length(t_vec))
res1s.df <- one_state_filter(t_vec, X, X1)
res2s.df <- two_state_filter(t_vec, X, X1, XD)
res3s.df <- three_state_filter(t_vec, X, X1, XD, XDD)

#---- Fig. 2.25 First-order least-squares fits best tracks first-order measurement
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res1s.df$TEE, y = res1s.df$X1, col = "True signal")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res1s.df$TEE, y = res1s.df$XHAT, col = "Zeroth-Order")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res2s.df$TEE, y = res2s.df$XHAT, col = "First-Order")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res3s.df$TEE, y = res3s.df$XHAT, col = "Second-Order")) +
  ggplot2::ylab("Estimate of x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.25 First-order least-squares fits best tracks first-order measurement") #+ 
#ggplot2::coord_cartesian(ylim = c(0, 500))



#---- Prepare Fig. 2.26
TS <- 0.1 # sampling time in seconds
t_vec <- seq(0, 10, TS)
SIGNOISE <- 1
XNOISE <- rnorm(length(t_vec), 0, SIGNOISE)
X1 <- 5 * t_vec^2 - 2 * t_vec + 2
X <- X1 + 50 * XNOISE
XD <- 10 * t_vec - 2
XDD <- rep(10, length(t_vec))
res1s.df <- one_state_filter(t_vec, X, X1)
res2s.df <- two_state_filter(t_vec, X, X1, XD)
res3s.df <- three_state_filter(t_vec, X, X1, XD, XDD)


#---- Fig. 2.26 Second-order filter tracks parabolic signal quite well
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res1s.df$TEE, y = res1s.df$X1, col = "True signal")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res1s.df$TEE, y = res1s.df$XHAT, col = "Zeroth-Order")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res2s.df$TEE, y = res2s.df$XHAT, col = "First-Order")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res3s.df$TEE, y = res3s.df$XHAT, col = "Second-Order")) +
  ggplot2::ylab("Estimate of x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.26 Second-order filter tracks parabolic signal quite well") #+ 
#ggplot2::coord_cartesian(ylim = c(0, 500))



#--- Table 2.4 TODO Where are the data comming from?
sum((res1s.df$X1 - res1s.df$XHAT)^2)
sum((res2s.df$X1 - res2s.df$XHAT)^2)
sum((res3s.df$X1 - res3s.df$XHAT)^2)
#--- Table 2.5 TODO Where are the data comming from?
sum((res1s.df$X - res1s.df$XHAT)^2)
sum((res2s.df$X - res2s.df$XHAT)^2)
sum((res3s.df$X - res3s.df$XHAT)^2)





#---- Accelerometer testing example ----
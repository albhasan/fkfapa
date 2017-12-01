# chapter02


#---- Zeroth-Order or One-State filter ----
R <- vector()

x_star <- c(1.2, 0.2, 2.9, 2.1)
x_caret <- sum(x_star)/length(x_star)
k <- 1:length(x_star)
T_s <- 1 # sampling time in seconds
t_vec <- (k - 1) * T_s

# Table 2.1 Measurement data for example
data.frame(t_vec, x_star)

# Table 2.2 Measurement data from table 2.1 expressed more mathematically
data.frame(k, t_vec, x_star)

ggplot2::ggplot() +
  ggplot2::geom_point(mapping = ggplot2::aes(x = x_star, y = t_vec, col = "Measurement")) +
  ggplot2::geom_hline(yintercept = x_caret, mapping = ggplot2::aes(col = "Estimate")) +
  ggplot2::ylab("xhat") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.1. Constant 1.6 is an unreasonable fit to measurement data")
R <- append(R, sum((x_caret - x_star)^2))
R[length(R)]



#---- First-Order or Two-State filter ----


# Listing 2.1 Solving for least-squares coefficients in the first-order fiter
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
ggplot2::ggplot() +
  ggplot2::geom_point(mapping = ggplot2::aes(x = X, y = t_vec, col = "Measurement")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t_vec, y = est, col = "Estimate")) +
  ggplot2::ylab("xhat") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.2. Straight-line fit to data is better than constant fit")
R <- append(R, sum((est - X)^2))
R[length(R)]

#---- Second-Order or Three-State Least-Squares Filter ----


# Listing 2.2 Solving for least-squared coeficients with three-state least-squares filter
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
ggplot2::ggplot() +
  ggplot2::geom_point(mapping = ggplot2::aes(x = X, y = t_vec, col = "Measurement")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t_vec, y = x_caret, col = "Estimate")) +
  ggplot2::ylab("xhat") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.3. Parabolic fit to data is pretty good, too")
R <- append(R, sum((x_caret - X)^2))
R[length(R)]

#---- Third-Order System ----

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


ggplot2::ggplot() +
  ggplot2::geom_point(mapping = ggplot2::aes(x = t_vec, y = X, col = "Measurement")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t_vec, y = x_caret, col = "Estimate")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t_dense, y = x_dense, col = "Dense Estimate")) +
  ggplot2::ylab("xhat") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.3. Parabolic fit to data is pretty good, too")
R <- append(R, sum((x_caret - X)^2))
R[length(R)]

# Table 2.3 Residual decreases as order of least-squaredpolynomial increses
data.frame(System_order = 0:(length(R)-1), R = round(R, 2))

# Listing 2.3 One-state filter for extracting signal from measurement
N <- 0
TS <- 0.1 # sampling time in seconds
t_vec <- seq(0, 10, TS)
SIGNOISE <- 1
XNOISE <- rnorm(length(t_vec), 0, SIGNOISE)
X1 <- rep(1, length(t_vec))
X <- X1 + XNOISE
SUM3 <- sum(X)
NMAX <- length(t_vec) -1
N <- NMAX

A <- matrix(N)
B <- matrix(SUM3)
AINV <- 1/A
ANS <- AINV %*% B

i <- 1
SUMPZ1 <- 0
SUMPZ2 <- 0
res <- list()
for(i in 1:NMAX){
  TEE <- TS * (i-1)
  XHAT <- ANS[1,1]
  ERRX <- X1[i] - XHAT
  ERRXP <- X[i] - XHAT
  ERRX2 <- (X1[i] - XHAT)^2
  ERRXP2 <- (X[i] - XHAT)^2
  SUMPZ1 <- ERRX2 + SUMPZ1
  SUMPZ2 <- ERRXP2 + SUMPZ2
  res[[i]] <- c(TEE, X1[i], X[i], XHAT, ERRX, ERRXP, SUMPZ1, SUMPZ2)
}
res.df <- as.data.frame(do.call(rbind, res))
colnames(res.df) <- c("TEE", "X1", "X", "XHAT", "ERRX", "ERRXP", "SUMPZ1", "SUMPZ2")

ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t_vec, y = X, col = "Measurements")) +
  ggplot2::geom_hline(yintercept = XHAT) + 
  ggplot2::annotate("text", 1, XHAT, vjust = -1, label = "Estimate") + 
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.5. One-state, least-squares filter smoothesnoise measurements") + 
  ggplot2::coord_cartesian(ylim = c(-2, 4))

ggplot2::ggplot() +
  ggplot2::geom_hline(yintercept = XHAT) +
  ggplot2::geom_hline(yintercept = X1[1]) + 
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.6. One-state filteryields near perfect estimate of constant signal") + 
  ggplot2::coord_cartesian(ylim = c(0.0, 1.4))

res.df["s_e"] <- res.df$X1 - res.df$XHAT
res.df["m_e"] <- res.df$X - res.df$XHAT

ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$s_e, col = "Signal & estiamte")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$m_e, col = "Measurement & estiamte")) +
  ggplot2::ylab("x") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.7 Estimation errors are nearly zero for one-state least-squared filter") + 
  ggplot2::coord_cartesian(ylim = c(-4, 4))

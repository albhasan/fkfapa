# chapter02


#---- Zeroth-Order or One-State filter ----


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
(R <- sum((x_caret - x_star)^2))



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
(R <- sum((est - X)^2))


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
(R <- sum((x_caret - X)^2))
  
  
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


# TODO: 
# - finish
# - check former formula, I think it is wrong
x_caret <- ANS[1] + (ANS[2] * (k - 1) * T_s) + (ANS[3] * ((k - 1) * T_s)^2)
#ANS * c()


ggplot2::ggplot() +
  ggplot2::geom_point(mapping = ggplot2::aes(x = X, y = t_vec, col = "Measurement")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = t_vec, y = x_caret, col = "Estimate")) +
  ggplot2::ylab("xhat") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 2.3. Parabolic fit to data is pretty good, too")
(R <- sum((x_caret - X)^2))



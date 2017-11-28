# chapter01

#---- Euler numerical integration ----
# listing 1.6
# page 16/26
W <- 2
T <- 0
S <- 0
X <- 0
XD <- W
H <- 0.01
tmp <- list()
i <- 1
while(T <= 10){
  S <- S + H
  XDD <- -W * W * X
  XD <- XD + H * XDD
  X <- X + H * XD
  T <- T + H
  if(S >= 0.09999){
    S <- 0
    XTHEORY <- sin(W*T)
    tmp[[i]] <- c(T, X, XTHEORY)
    i <- i + 1
  }
}
res <- as.data.frame(do.call("rbind", tmp), stringsAsFactors = FALSE)
colnames(res) <- c("T", "X", "XTHEORY")
# fig 1.2
ggplot2::ggplot(data = res, ggplot2::aes(x=group, y=weight)) + 
  ggplot2::geom_path(mapping = ggplot2::aes(x = T, y = XTHEORY, col = "Theoretical")) + 
  ggplot2::geom_point(mapping = ggplot2::aes(x = T, y = X, col = "Numerical")) + 
  ggplot2::xlab("Time (Sec)") + ggplot2::ylab("X") + 
  ggplot2::ggtitle("Fig. 1.2. Euler integration")

#---- Runge-Kutta (second order) numerical integration ----
# listing 1.7
# page 19/29
W <- 2
T <- 0
S <- 0
X <- 0
XD <- W
H <- 0.01
tmp <- list()
i <- 1
while(T <= 10){
  S <- S + H
  XOLD <- X
  XDOLD <- XD
  XDD <- -W * W * X # differential equation
  X <- X + H * XD
  XD <- XD + H * XDD
  T <- T + H
  XDD <- -W * W * X # differential equation
  X <- 0.5 * (XOLD + X + H * XD)
  XD <- 0.5 * (XDOLD + XD + H * XDD)
  if(S >= 0.09999){
    S <- 0
    XTHEORY <- sin(W*T)
    tmp[[i]] <- c(T, X, XTHEORY)
    i <- i + 1
  }
}
res <- as.data.frame(do.call("rbind", tmp), stringsAsFactors = FALSE)
colnames(res) <- c("T", "X", "XTHEORY")
# fig 1.2
ggplot2::ggplot(data = res, ggplot2::aes(x=group, y=weight)) + 
  ggplot2::geom_path(mapping = ggplot2::aes(x = T, y = XTHEORY, col = "Theoretical")) + 
  ggplot2::geom_point(mapping = ggplot2::aes(x = T, y = X, col = "Numerical")) + 
  ggplot2::xlab("Time (Sec)") + ggplot2::ylab("X") + 
  ggplot2::ggtitle("Fig. 1.4. Second order Runge-Kutta integration")


#---- Gaussian noise example ----
set.seed(666)
rn <- rnorm(1000)
hx <- hist(rn, breaks = 50, plot = FALSE)
xseq <- seq(range(hx$breaks[-1])[1], range(hx$breaks[-1])[2], length.out = length(hx$breaks))
densities <- dnorm(xseq, mean(rn), sd(rn))
#plot(x = hx$breaks[-1], y = hx$density, type = "l", xlim = c(-6, 6))
#lines(x = hx$breaks, y = densities)
res <- data.frame(x = hx$breaks[-1], y = hx$density, type = "Calculated")
res <- rbind(res, data.frame(x = hx$breaks, y = densities, type = "Theory"))

ggplot2::ggplot() + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = 1:length(rn), y = rn)) + 
  ggplot2::xlab("Number of random numbers") + ggplot2::ylab("x") +
  ggplot2::ggtitle("Fig. 1.7. One thousand random numbers with Gaussian distribution")

ggplot2::ggplot(data = res) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y, group = type, col = type)) + 
  ggplot2::xlab("x") + ggplot2::ylab("p(x)") + 
  ggplot2::ggtitle("Fig. 1.8. Sampled Gaussian distribution")

set.seed(666)
rn <- rnorm(5000)
hx <- hist(rn, breaks = 250, plot = FALSE)
xseq <- seq(range(hx$breaks[-1])[1], range(hx$breaks[-1])[2], length.out = length(hx$breaks))
densities <- dnorm(xseq, mean(rn), sd(rn))
res <- data.frame(x = hx$breaks[-1], y = hx$density, type = "Calculated")
res <- rbind(res, data.frame(x = hx$breaks, y = densities, type = "Theory"))
ggplot2::ggplot(data = res) + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = x, y = y, group = type, col = type)) + 
  ggplot2::xlab("x") + ggplot2::ylab("p(x)") + 
  ggplot2::ggtitle("Fig. 1.9. Sampled Gaussian distribution is even better for 5000 random numbers")



#--- Calculating standard deviation ----




n <- 1000
for(i in 1:n){
  
}



# listing 1.10
set.seed(666)
sd_actual = 1
mean_actual = 0
rms_actual = 1
number_of_samples <- 1000
samples <- rnorm(number_of_samples, mean_actual, sd_actual)
m <- rep(NA, number_of_samples)
std <- rep(NA, number_of_samples)
rms <- rep(NA, number_of_samples)
for(n in 1:number_of_samples){
  m[n] <- mean(samples[1:n])
  std[n] <- sd(samples[1:n])
  rms[n] <- sqrt(mean(samples[1:n]^2))
}
ggplot2::ggplot() + 
  ggplot2::geom_line(mapping = ggplot2::aes(x = 1:number_of_samples, y = m, col = "mean")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = 1:number_of_samples, y = std, col = "std")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = 1:number_of_samples, y = rms, col = "rms")) + 
  ggplot2::xlab("Number of samples") + ggplot2::ylab("Calculated statistics") + 
  ggplot2::ggtitle("Fig. 1.10. Many samples must be taken before computer-calculated statistics are accurate")


# listing 1.11 Simulation of low pass filter driven by white noise
TAU <- 0.2
PHI <- 1
TEE <- 0
H <- 0.01
SIG <- sqrt(PHI/H)
Y <- 0
res <- list()
n <- 1
while(TEE <= 4.999){
  X <- rnorm(1, sd = SIG)
  YOLD <- Y
  YD <- (X - Y)/TAU
  Y <- Y + H * YD
  TEE <- TEE + H
  YD <- (X - Y)/TAU
  Y <- (YOLD + Y)/2 + 0.5 * H * YD
  SIGPLUS <- sqrt(PHI * (1 - exp(-2 * TEE/TAU))/(2 * TAU))
  SIGMINUS = -SIGPLUS
  res[[n]] <- c(TEE, Y, SIGPLUS, SIGMINUS)
  n <- n + 1 
}
res.df <- as.data.frame(do.call(rbind, res), stringsAsFactors = FALSE)
colnames(res.df) <- c("TEE", "Y", "SIGPLUS", "SIGMINUS")
ggplot2::ggplot() +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$SIGPLUS, col = "Theory")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$SIGMINUS, col = "Theory")) +
  ggplot2::geom_line(mapping = ggplot2::aes(x = res.df$TEE, y = res.df$Y, col = "Simulation")) +
  ggplot2::ylab("y") + ggplot2::xlab("Time (Sec)") + 
  ggplot2::ggtitle("Fig. 1.12. Low-pass filter output agrees with theory")

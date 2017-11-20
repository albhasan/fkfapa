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
ggplot2::ggplot(data = res, aes(x=group, y=weight)) + 
  ggplot2::geom_path(mapping = aes(x = T, y = XTHEORY, col = "Theoretical")) + 
  ggplot2::geom_point(mapping = aes(x = T, y = X, col = "Numerical")) + 
  xlab("Time (Sec)") + ylab("X") + 
  ggtitle("Fig. 1.2. Euler integration")

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
ggplot2::ggplot(data = res, aes(x=group, y=weight)) + 
  ggplot2::geom_path(mapping = aes(x = T, y = XTHEORY, col = "Theoretical")) + 
  ggplot2::geom_point(mapping = aes(x = T, y = X, col = "Numerical")) + 
  xlab("Time (Sec)") + ylab("X") + 
  ggtitle("Fig. 1.4. Second order Runge-Kutta integration")


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



ggplot2::ggplot(data = res) + 
  ggplot2::geom_line(mapping = aes(x = x, y = y, group = type, col = type)) + 
  xlab("x") + ylab("p(x)") + 
  ggtitle("Fig. 1.8. Sampled Gaussian distribution")








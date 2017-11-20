# euler numerical integration
# listing 1.6
# page 16/26
W <- 2
T <- 0
S <- 0
X <- 0
XD <- W
H <- 0.01
res <- list()

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
    res[[i]] <- c(T, X, XTHEORY)
    i <- i + 1
  }
}
res <- do.call("rbind", res)
colnames(res) <- c("T", "X", "XTHEORY")
# fig 1.2
plot(res$XTHEORY)

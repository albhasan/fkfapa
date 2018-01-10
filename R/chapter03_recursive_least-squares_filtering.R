# chapter03

#---- MAKING ZEROTH-ORDER LEAST SQUARES FILTER RECURSIVE ----

#---- Properties of zeroth-order or one-state filter ----

# table 3.1
sample_measurement_data <- data.frame(k = 1:4, T_s = 0:3, xstar_k = c(1.2, 0.2, 2.9, 2.1))
smd <- sample_measurement_data

K <- 1/(smd$k)
Res <- rep(NA, length(K))
xhat.vec <- rep(NA, length(K))
xhat <- 0

for(i in 1:(nrow(smd))){
  Res <- smd[i, 'xstar_k'] - xhat[i]
  xhat.vec
}

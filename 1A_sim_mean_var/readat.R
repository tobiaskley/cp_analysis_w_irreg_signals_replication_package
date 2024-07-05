library(tidyverse)

# Loading data
Z_mean_var <- data.frame(theta    = numeric(),
                         mean     = numeric(),
                         mean_sd  = numeric(),
                         mean_CI  = numeric(),
                         var1     = numeric(),
                         var1_sd  = numeric(),
                         var1_CI  = numeric(),
                         var2     = numeric(),
                         var2_sd  = numeric(),
                         var2_CI  = numeric(),
                         var      = numeric(),
                         var_sd   = numeric(),
                         var_CI   = numeric())

THETA <- c(0.2, 0.3, 0.4) # length = 6

for (j in 1:length(THETA)) {
  data <- NULL
  dataSq <- NULL
  data2 <- c()
  
  for (i in 1:300) {
    path <- paste0("res_files/res_X_bar_", i, ".Rdata")
    load(path)
    if (theta == THETA[j]) {
      data <- rbind(data, res)
      dataSq <- rbind(dataSq, resSq)
    }
  }
  
  for (l in 1:2) {
    data2 <- c(data2, mean( 100 * (data[,l] - mean(data))^2 ))
  }
  
  #data <- data[1:5000,]
  #dataSq <- dataSq[1:5000,]
  
  x1 <- (sqrt(100000) * (data[,1] - mean(data)))^2
  x2 <- (sqrt(100000) * (data[,2] - mean(data)))^2
  x  <- (sqrt(200000) * (rowMeans(data) - mean(data)))^2
  sqrtR <- sqrt(nrow(data))
  Z_mean_var <- rbind(Z_mean_var,
                      data.frame(theta   = THETA[j],
                                 mean    = mean(data),
                                 mean_sd = sqrt( (mean(dataSq) - mean(data)^2)/ nrow(data) ),
                                 var1     = mean( x1 ),
                                 var1_sd  = sd( x1 ) / sqrtR,
                                 var2     = mean( x2 ),
                                 var2_sd  = sd( x2 ) / sqrtR,
                                 var    = mean( x ),
                                 var_sd = sd( x ) / sqrtR
                      ))
}

Z_mean_var$mean_CI <- paste("[", round(Z_mean_var$mean - qnorm(0.95) * Z_mean_var$mean_sd,3),
                            ",", round(Z_mean_var$mean + qnorm(0.95) * Z_mean_var$mean_sd,3), "]", sep = "")
Z_mean_var$var1_CI <- paste("[", round(Z_mean_var$var1 - qnorm(0.95) * Z_mean_var$var1_sd,3),
                            ",", round(Z_mean_var$var1 + qnorm(0.95) * Z_mean_var$var1_sd,3), "]", sep = "")
Z_mean_var$var2_CI <- paste("[", round(Z_mean_var$var2 - qnorm(0.95) * Z_mean_var$var2_sd,3),
                            ",", round(Z_mean_var$var2 + qnorm(0.95) * Z_mean_var$var2_sd,3), "]", sep = "")
Z_mean_var$var_CI <- paste("[", round(Z_mean_var$var - qnorm(0.95) * Z_mean_var$var_sd,3),
                           ",", round(Z_mean_var$var + qnorm(0.95) * Z_mean_var$var_sd,3), "]", sep = "")

# for special case of theta == 0, expectation is 0 and variance is 1
Z_mean_var[4, 1] <- 0 ## theta
Z_mean_var[4, 2] <- 0 ## mean
Z_mean_var[4, 8] <- 1 ## mean

save(Z_mean_var, file = "Z_mean_var.Rdata")

load("Z_mean_var.Rdata")
round(Z_mean_var %>%
  select(theta, mean, var) %>%
  filter(theta >= 0), 3)

require(changepoint)

# Classical CUSUM method
cusum_method <- function(data) {
  
  cumsum_tau <- which.min(cumsum(data - mean(data)))

  return(cumsum_tau + 1)
}

# AMOC method
amoc_method <- function(data) {
  
  amoc_tau <- changepoint::cpts(cpt.mean(data = data, method = "AMOC"))
  
  return(amoc_tau + 1)
}

# Binary Segmentation method
binseg_method <- function(data, first = 1) {
  
  out <- wbs::changepoints(wbs::sbs(data))
  return(sort(out$cpt.th[[1]])[1:first] + 1)

}

# Simple Change Point Regression on slope
slope_change_method <- function(data) {
  
  n <- length(data)
  
  aux_RSS <- function(k) {
    time <- c(rep(0, k), 1:(n-k)) # k-1 times 0, then 0, 1, ... n-k
    lm0 <- lm(data ~ 1 + time)
    return(sum(resid(lm0)^2))
  }
  aux_RSS <- Vectorize(aux_RSS)
  
  return(which.min(aux_RSS(1:(n-1))))
  
}



# Calculate the test statistic using precomputed cumulative sums
test_statistic_dp <- function(s, e, b, cum_sum) {
  n <- e - s + 1
  mean1 <- ifelse(s == 1, cum_sum[b] / b, (cum_sum[b] - cum_sum[s-1]) / (b - s + 1))
  mean2 <- (cum_sum[e] - cum_sum[b]) / (e - b)
  test_stat <- sqrt((e - b) * (b - s + 1) / n) * (mean1 - mean2)
  return(test_stat)
}

# The recursive Binary Segmentation function with dynamic programming
binseg_rec_dp <- function(s, e, X, zeta, cum_sum, memo) {
  if (e - s < 1) {
    return(NULL)
  } else if (!is.null(memo[[paste(s, e, sep="-")]])) {
    return(memo[[paste(s, e, sep="-")]])
  } else {
    # Find the b0 that maximizes the absolute test statistic using dp approach
    b0_values <- sapply(s:(e-1), function(b) abs(test_statistic_dp(s, e, b, cum_sum)))
    b0 <- which.max(b0_values) + s - 1
    test_stat <- test_statistic_dp(s, e, b0, cum_sum)
    if (abs(test_stat) > zeta) {
      change_points_left <- binseg_rec_dp(s, b0, X, zeta, cum_sum, memo)
      change_points_right <- binseg_rec_dp(b0 + 1, e, X, zeta, cum_sum, memo)
      change_points <- c(b0, change_points_left, change_points_right)
      memo[[paste(s, e, sep="-")]] <- change_points
      return(change_points)
    } else {
      memo[[paste(s, e, sep="-")]] <- NULL
      return(NULL)
    }
  }
}


binseg_zeta <- function(X, zeta) {
  s <- 1
  e <- length(X)
  cum_sum <- cumsum(X) 
  memo <- list() # Initialize memoization list
  result <- binseg_rec_dp(s, e, X, zeta, cum_sum, memo)
  if(is.null(result)){
    print(paste0("No Change Point from BS, Return NA, zeta = ", zeta))
    return(NA)
  }else{
    return(min(result) + 1)
  }
}

# # Set seed for reproducibility
# set.seed(42)
# 
# n = 100
# # Generate synthetic data with known change-points
# segment1 <- rnorm(n/4, mean = 0)   # 100 points with mean 0
# segment2 <- rnorm(n/2, mean = 100)   # 100 points with mean 5 (change-point at 100)
# segment3 <- rnorm(n/4, mean = 5)  
# X <- c(segment1, segment2, segment3) # Concatenate segments
# 
# # Visualize the data
# plot(X, type = 'l', main = 'Synthetic Data with Change-Points', xlab = 'Index', ylab = 'Value')
# 
# # Define a threshold zeta. This value may need to be adjusted based on the noise level of the data.
# zeta <- 1 * sqrt(2*log(n))
# 
# # Call the binseg function with the synthetic data
# change_points <- binseg_zeta(X, zeta)
# print(change_points)
# 
# # Adding vertical lines to the plot for the detected change-points
# if (!is.null(change_points)) {
#   abline(v = change_points, col = 'red', lwd = 2)
# }


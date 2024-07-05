## (1) to run on cluster:
task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

library(parallel)
library(changepoint)
library(zoo)

# Initialise L'Ecuyer's ?RngStreams?
RNGkind("L'Ecuyer-CMRG")
set.seed(12345)

# fetch seed for this task_id's substream
s <- .Random.seed
for (i in 1:task_id) {
  s <- nextRNGStream(s)  # send s to worker i as .Random.seed
}
.Random.seed <- s

source("../R/functions_for_CPA.R")
source("../R/functions_for_data_generation.R")
source("../R/functions_for_classical_locating.R")

# helper function for box plots
box_stats <- function(data) {
  if (is.null(data) || length(data) == 0) {
    stop("Data must not be NULL or empty.")
  }
  
  # Ensure data is a numeric vector
  if (!is.numeric(data)) {
    stop("Data must be numeric.")
  }
  
  # Calculate the quartiles and IQR
  Q1 <- as.numeric(quantile(data, 0.25, na.rm = TRUE))
  Q3 <- as.numeric(quantile(data, 0.75, na.rm = TRUE))
  IQR <- Q3 - Q1
  median_value <- median(data, na.rm = TRUE)
  
  # Calculate whiskers
  upper_whisker <- min(max(data[data <= Q3 + 1.5 * IQR], na.rm = TRUE), na.rm = TRUE)
  lower_whisker <- max(min(data[data >= Q1 - 1.5 * IQR], na.rm = TRUE), na.rm = TRUE)
  
  # Return the results as a list
  return(list(
    median = median_value,
    first_quantile = Q1,
    third_quantile = Q3,
    lower_whisker = lower_whisker,
    upper_whisker = upper_whisker
  ))
}

# Define a function to calculate box statistics and create a data frame
create_data_frame <- function(data, n, s, theta, method) {
  stats <- box_stats(data)
  data.frame(
    n = n, s = s, theta = theta, method = method,
    median = stats$median,
    first_quantile = stats$first_quantile,
    third_quantile = stats$third_quantile,
    lower_whisker = stats$lower_whisker,
    upper_whisker = stats$upper_whisker
  )
}

# Set simulation parameters
N <- c(50, 100, 300, 500, 2000) # length = 5
THETA <- c(0, 0.2, 0.3, 0.4) # length = 4
d <- 1
S <- c(0.4, 0.5, 0.6, 0.7, 0.8, 1.0) # length = 6
WS <- c(1/3, 1/2, 2/3) # length = 3
sd <- 0.5
alpha <- 0.05
nexp <- 100000

# Generate configurations for simulation
cfg <- list()
for(n in N){
  for(theta in THETA){
    for(s in S){
      cfg <- c(cfg, list(list(n  = n, theta = theta, s = s)))
    }
  }
}

# detach parameters from configuration
n  <- cfg[[task_id]]$n
s <- cfg[[task_id]]$s
theta <- cfg[[task_id]]$theta
tau <- floor(n * 0.4)

# mean of z_prime
load("../1A_sim_mean_var/Z_mean_var.Rdata")
z_prime_mean <- sign(theta) * sd * Z_mean_var[Z_mean_var$theta == abs(theta), ]$mean
z_prime_sigma <- sd * sqrt(Z_mean_var[Z_mean_var$theta == abs(theta), ]$var)

load("../1B_sim_Q_bridge/Q_Bridge.Rdata")
quantile_fin <- Q_Bridge[Q_Bridge$n == n & Q_Bridge$perc == alpha, ]$q_mean
c_1 <- 1

# Initialize variables for storing results
serror_tw <- array(dim = c(nexp, length(WS), 2)) # [,,1] true d, [,,2] \hat d (rest^(1/2) best one)
serror_tw_25 <- array(dim = c(nexp, length(WS), 2))
serror_tw_75 <- array(dim = c(nexp, length(WS), 2))
serror_cs <- numeric()
serror_ac <- numeric()
serror_bs <- numeric()
serror_bsz1 <- numeric()
num_reject <- 0

# Generate trend for simulation
trend <- trend_simulations(n = n, tau = tau, tau_prime = floor(0.6 * n),
                           tau_double_prime = floor(0.8 * n), mu_1 = 0, s = s)

# ... then do all the computations 
start_time <- Sys.time()
prc <- 0
for(i in 1:nexp){
  # Generate noise and data
  noise <- gen_z(n = n, theta = theta, estimated_mean = z_prime_mean, sd = sd)
  data <- trend + noise
  
  # Test for change point
  tw_test <- test_CP(data, asymptotic = TRUE,
                     quantile_val = quantile_fin, J = 3)$REJECT
  
  if(tw_test) {
    num_reject <- num_reject + 1
  } else {
    print("Unable to reject H0.")
  }
  
  # Estimate change point location using different methods
  for(ws in 1:length(WS)){
    k <- ceiling(n^WS[ws])
    
    # with true sigma and true d
    ts_tau <- ifelse(tw_test, IR_CPE(data = data, k = k,sigma = z_prime_sigma,
                                     hat_d = s * d)$hat_tau, NA)
    ts_tau_25 <- ifelse(tw_test, IR_CPE(data = data, k = k, sigma = z_prime_sigma,
                                        hat_d = s * d, rho = 0.25)$hat_tau, NA)
    ts_tau_75 <- ifelse(tw_test, IR_CPE(data = data, k = k, sigma = z_prime_sigma,
                                        hat_d = s * d, rho = 0.75)$hat_tau, NA)
    
    # with estimated sigma and estimated d
    ts_dhat_tau <- ifelse(tw_test, IR_CPE(data = data, k = k, J = 3)$hat_tau, NA)
    ts_dhat_tau_25 <- ifelse(tw_test, IR_CPE(data = data, k = k, J = 3,
                                             rho = 0.25)$hat_tau, NA)
    ts_dhat_tau_75 <- ifelse(tw_test, IR_CPE(data = data, k = k, J = 3,
                                             rho = 0.75)$hat_tau, NA)
    
    hat_sigma <- IR_CPE(data = data, k = k, J = 3)$hat_sigma
    
    # Store scaled errors for each method
    serror_tw[i,ws,1] <- abs(ts_tau - tau)
    serror_tw[i,ws,2] <- abs(ts_dhat_tau - tau)
    serror_tw_25[i,ws,1] <- abs(ts_tau_25 - tau)
    serror_tw_25[i,ws,2] <- abs(ts_dhat_tau_25 - tau)
    serror_tw_75[i,ws,1] <- abs(ts_tau_75 - tau)
    serror_tw_75[i,ws,2] <- abs(ts_dhat_tau_75 - tau)
  }
  
  # Estimate change point location using classical methods
  zeta_1 <- c_1 * hat_sigma * sqrt(2*log(n))
  
  cumsum_tau <- cusum_method(data)
  amoc_tau <- amoc_method(data = data)
  binseg_tau <- binseg_method(data = data)
  bsz1_tau <- binseg_zeta(data, zeta_1)
  
  # Store absolute errors for classical methods
  serror_cs <- c(serror_cs, abs(cumsum_tau - tau))
  serror_ac <- c(serror_ac, abs(amoc_tau - tau))
  serror_bs <- c(serror_bs, abs(binseg_tau - tau))
  serror_bsz1 <- c(serror_bsz1, abs(bsz1_tau - tau))
  
  end_time <- Sys.time()
  gc(reset = TRUE)
  
  # Print progress
  if (floor(100*i/nexp) > prc) {
    prc <- floor(100*i/nexp)
    cat("replication ", i, " of ", nexp, " completed. ",
        prc, "done.",
        difftime(end_time, start_time, units="mins"), "min elapsed.\n")
  }
}

# Store results
result_alt <- c(n = n,
                tau = tau,
                s  = s,
                theta = theta,
                sigma = z_prime_sigma,
                serror_tw_50_ = apply(serror_tw, c(2,3), function(x) mean(x, na.rm = TRUE)),
                serror_cs = mean(serror_cs, na.rm = TRUE),
                serror_ac = mean(serror_ac, na.rm = TRUE),
                serror_bs = mean(serror_bs, na.rm = TRUE),
                serror_bsz1 = mean(serror_bsz1, na.rm = TRUE),
                serror_tw_25_ = apply(serror_tw_25, c(2,3), function(x) mean(x, na.rm = TRUE)),
                serror_tw_75_ = apply(serror_tw_75, c(2,3), function(x) mean(x, na.rm = TRUE)),
                reject_ratio_tw = num_reject / nexp,
                na_ratio_tw = apply(serror_tw, c(2,3), function(x) mean(is.na(x))),
                na_ratio_tw_25_ = apply(serror_tw_25, c(2,3), function(x) mean(is.na(x))),
                na_ratio_tw_75_ = apply(serror_tw_75, c(2,3), function(x) mean(is.na(x))))

# Using the function to generate data frames for each case
data_tw_trued_50 <- create_data_frame(serror_tw[,1,1], n, s, theta,
                                      "Our Method w/ true d, rho = .50")
data_tw_estd_50 <- create_data_frame(serror_tw[,1,2], n, s, theta,
                                     "Our Method w/ est. d, rho = .50")
data_cs <- create_data_frame(serror_cs, n, s, theta, "CUSUM")
data_ac <- create_data_frame(serror_ac, n, s, theta, "LIKELIHOOD")
data_bs <- create_data_frame(serror_bs, n, s, theta, "Binary Segmentation")
data_bsz1 <- create_data_frame(serror_bsz1, n, s, theta, "Binary Segmentation, long_run")
data_tw_trued_25 <- create_data_frame(serror_tw_25[,1,1], n, s, theta,
                                      "Our Method w/ true d, rho = .25")
data_tw_trued_75 <- create_data_frame(serror_tw_75[,1,1], n, s, theta,
                                      "Our Method w/ true d, rho = .75")
data_tw_estd_25 <- create_data_frame(serror_tw_25[,1,2], n, s, theta,
                                     "Our Method w/ est. d, rho = .25")
data_tw_estd_75 <- create_data_frame(serror_tw_75[,1,2], n, s, theta,
                                     "Our Method w/ est. d, rho = .75")

# Combine all data frames into one
data_box <- rbind(
  data_tw_trued_50, data_tw_estd_50, data_cs, data_ac, data_bs, data_bsz1,
  data_tw_trued_25, data_tw_trued_75, data_tw_estd_25, data_tw_estd_75
)

# Save results
save(result_alt, file = paste("result_alt_", task_id, ".Rdata", sep = ""))
save(data_box, file = paste("result_box_", task_id, ".Rdata", sep = ""))

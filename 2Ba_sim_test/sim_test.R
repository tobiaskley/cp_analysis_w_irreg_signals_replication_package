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
  s <- nextRNGStream(s)
  # send s to worker i as .Random.seed
}
.Random.seed <- s

source("../R/functions_for_CPA.R")
source("../R/functions_for_data_generation.R")
source("../R/functions_for_classical_locating.R")

N <- c(50, 100, 300, 500, 2000) # length = 5
THETA <- c(-0.4, -0.2, 0, 0.2, 0.4) # length = 5
d <- 1
S <- seq(0, 0.045, length.out = 80) # length = 80

sd <- 0.5

nexp <- 100000
alpha = 0.05

cfg <- list()

for(n in N){
    for(theta in THETA){
      for(s in S){
        cfg <- c(cfg, list(list(n  = n, theta = theta, s = s
        )))
      
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

rej_alt      <- numeric()
rej_null_asy <- numeric()
rej_null_fin <- numeric()


trend <- trend_simulations(n = n, tau = tau, tau_prime = floor(0.6 * n),
                           tau_double_prime = floor(0.8 * n), mu_1 = 0, s = s)

load("../1B_sim_Q_bridge/Q_Bridge.Rdata")
quantile_fin <- Q_Bridge[Q_Bridge$n == n & Q_Bridge$perc == alpha, ]$q_mean
quantile_asy <- -1*sqrt(-0.5*log(alpha))

# ... then do all the computations 
start_time <- Sys.time()
prc <- 0

for(i in 1:nexp){
  
  # noise <- generate_z(n = n, theta = theta, estimated_mean = z_prime_mean, sd = sd)
  noise <- gen_z(n = n, theta = theta, estimated_mean = z_prime_mean, sd = sd)
  data <- trend + noise
  
  rej_alt <- c(rej_alt, test_CP(data, sigma = z_prime_sigma, asymptotic = TRUE,
                                quantile_val = quantile_asy)$REJECT)
  
  rej_null_asy <- c(rej_null_asy, test_CP(noise, sigma = z_prime_sigma,
                                          asymptotic = TRUE,
                                          quantile_val = quantile_asy)$REJECT)
  
  rej_null_fin <- c(rej_null_fin, test_CP(noise, sigma = z_prime_sigma,
                                          asymptotic = TRUE,
                                          quantile_val = quantile_fin)$REJECT)
  
  end_time <- Sys.time()
  gc(reset = TRUE)
  
  if (floor(100*i/nexp) > prc) {
    prc <- floor(100*i/nexp)
    cat("replication ", i, " of ", nexp, " completed. ",
        prc, "done.",
        difftime(end_time, start_time, units="mins"), "min elapsed.\n")  
  }
  
}

result_alt <- c(n = n,
                tau = tau,
                s  = s,
                theta = theta,
                sigma = z_prime_sigma,
                power = mean(rej_alt),
                t1e_asy = mean(rej_null_asy),
                t1e_fin = mean(rej_null_fin)
                
)


save(result_alt, file = paste("result_alt_", task_id, ".Rdata", sep = ""))

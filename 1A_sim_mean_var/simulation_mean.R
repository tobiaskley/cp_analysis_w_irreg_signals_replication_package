# This script simulates Z_1, ..., Z_{200000} and then
# computes ten block averages R_1 with k = 1000

## (1) to run on cluster:
task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
library(parallel)

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

source("../R/functions_for_data_generation.R")
source("../R/functions_for_CPA.R")

THETA <- c(0.2, 0.3, 0.4) # length = 3

sd <- 1

nexp <- 100000
ncopies <- 100

cfg <- list()

for(theta in THETA){
  for(cpy in 1:ncopies) {
    cfg <- c(cfg, list(list(theta = theta, cpy = cpy)))
  }
}

# detach parameters from configuration
theta  <- cfg[[task_id]]$theta

res   <- matrix(nrow = nexp, ncol = 2)
resSq <- matrix(nrow = nexp, ncol = 2)

# ... then do all the computations 
start_time <- Sys.time()
prc <- 0

for(r in 1:nexp){
  
  Z_prime <- gen_z(200000, theta, estimated_mean = 0, sd, burn_in = 10000)
  
  res[r,] <- get_data_blocks(Z_prime, k = 100000)
  resSq[r,] <- get_data_blocks(Z_prime^2, k = 100000)
  
  end_time <- Sys.time()
  gc(reset = TRUE)
  
  if (floor(100*r/nexp) > prc) {
    prc <- floor(100*r/nexp)
    cat("replication ", r, " of ", nexp, " completed. ",
        prc, "done.",
        difftime(end_time, start_time, units="mins"), "min elapsed.\n")  
  }
  
}

save(res, resSq, theta, file = paste("res_X_bar_", task_id, ".Rdata", sep = ""))

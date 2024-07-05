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

source("../R/functions_for_CPA.R")

N <- c(50, 100, 300, 500, 2000) # length = 5

nexp <- 100
ncopies <- 10
alpha = c(0.2, 0.1, 0.05, 0.01)

cfg <- list()

for(n in N){
  for(cpy in 1:ncopies) {
    cfg <- c(cfg, list(list(n = n, cpy = cpy)))
  }
}

# detach parameters from configuration
n  <- cfg[[task_id]]$n

res <- matrix(nrow = nexp, ncol = length(alpha))

# ... then do all the computations 
start_time <- Sys.time()

for(r in 1:nexp){
  
  res[r, ] <- as.numeric(Q_bridge(alpha, n = n, N = 1e6))
  
  end_time <- Sys.time()
  gc(reset = TRUE)
  cat("replication ", r, " of ", nexp, " completed. ", difftime(end_time, start_time, units="mins"), "min elapsed.\n")
  
}

save(res, n, file = paste("res_Q_bridge_", task_id, ".Rdata", sep = ""))

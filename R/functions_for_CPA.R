library(zoo)

# Utility function to compute data blocks
get_data_blocks <- function(data, k) {
  n <- length(data)
  indices <- seq(1, n, by = k)
  m <- floor(n / k)
  indices <- indices[1:m]
  data_block <- sapply(indices, function(idx) {
    mean(data[idx:min(idx + k - 1, n)])
  })
  return(data_block)
}

# LRV estimate from Section 2.3.2
compute_sigma <- function(data, ell, k = ceiling(length(data)^(1/3))) {
  data_block <- zoo::rollmean(data[1:ell], k) - mean(data[1:ell])
  return( sqrt( sum(data_block^2) * k / (ell - k + 1) ))
}


# plot for I_seq
plot_I <- function(I_seq, hat_eta, breaks = 1:length(I_seq)) {
  
  m <- length(I_seq)
  df <- rbind(data.frame(j = 1:m,
                   I = I_seq))
  df$I <- factor(df$I, levels = c(0,1), labels = c("0", "1")) 
  
  p <- ggplot(df) +
    geom_point(aes(x = j, y = as.numeric(I) - 1, color = I)) +
    geom_line(aes(x = j, y = 0), data = df[df$j <= hat_eta, ]) +
    geom_line(aes(x = j, y = 1), data = df[df$j > hat_eta, ]) +
    scale_color_manual(values = c("0" = "red", "1" = "blue")) +
    scale_y_continuous(breaks = c(0,1), limits = c(-0.1, 1.1)) +
    scale_x_continuous(breaks = breaks) +
    theme(legend.position="none") +
    xlab(expression(j)) + ylab(expression(hat(I)[j]))
 
  return(p)
}

# J specified the how many of the smallest blocks are to be taken into account
#   when obtaining \hat L

# Two-Step method for Change point estimation
IR_CPE <- function(data,
                   sigma = NULL,
                   k = NULL,
                   k2 = NULL,
                   J = 1,
                   hat_d = NULL,
                   hat_mu1 = NULL,
                   rho = 1/2,
                   eta = NULL) {
  n <- length(data)
  
  k <- ifelse(is.null(k), ceiling(n^(1/3)), k)
  m <- floor(n / k)

  data_block <- get_data_blocks(data, k)
  
  hat_L <- max(order(data_block)[1:J])
  l <- min(hat_L * k, n)
  hat_mu0 <- mean(data[1:l])
  
  hat_sigma <- ifelse(is.null(sigma), compute_sigma(data, l), sigma)
  I_seq <- NULL
  hat_eta <- NULL
  
  if (is.null(hat_mu1)) {
    D_seq <- sqrt(k) * (data_block - hat_mu0) / hat_sigma
    
    z <- qnorm(1 - 1 / m)
    I_seq <- 1 * (D_seq > z)
    cumsum_I <- cumsum(I_seq)
    t_values <- 1:(m - 1)
    sums <-
      cumsum_I[t_values] + (m - t_values) - (cumsum_I[m] - cumsum_I[t_values])
    
    
    hat_eta <- which(min(sums) == sums)

    hat_mu1 <- mean(data[1:(k * min(hat_eta))])
    
    if (is.null(hat_d)) {

      dat_d <- data[(k * max(hat_eta) + k + 1):n]
      
      rest_len <- n - (k * max(hat_eta) + k)
      
      k2 <- ifelse(is.null(k2), floor(rest_len^(1/2)), k2)
      
      averages <- zoo::rollmean(dat_d, k2) - hat_mu1
      
      if (length(averages)) {
        hat_d <- min(averages)  
      } else {
        hat_d <- 0
      }
      
    }
  }
  
  if(hat_d <= 0){
    print(paste0("Return NA because hat_d = ", hat_d))
    hat_tau <- NA
  }else{
    # Step 2: Compute hat_tau
    deviations <- cumsum(data - hat_mu1 - rho * hat_d)[-n]
    hat_tau <- which.min(deviations) + 1
  }

  return(
    list(
      n = n,
      k = k,
      hat_L = hat_L,
      hat_mu0 = hat_mu0,
      I_seq = I_seq,
      hat_sigma = hat_sigma,
      hat_eta = hat_eta,
      hat_mu1 = hat_mu1,
      hat_d = hat_d,
      hat_tau = hat_tau
    )
  )
}

# for simulation of Brownian Bridge

# Function to compute the cumulative sum statistic
compute_cusum_stat <- function(n) {
  x <- rnorm(n)
  min(cumsum(x - mean(x))) / sqrt(n)
}

cusum_simulation <- function(n, size = 300) {
  results <- replicate(size, compute_cusum_stat(n))
  return(list(sample = results, size = size, n = n))
}

Q_bridge <- function(alpha, n = 1e5, N = 1e6) {
  sample <- cusum_simulation(n, N)$sample
  return(quantile(sample, alpha))
}

# Testing for change point

test_CP <- function(data,
                    alpha = 0.05,
                    sigma = NULL,
                    k = NULL,
                    l_min = 1,
                    J = 1,
                    simulation_size = 500,
                    simulation_sample = NULL,
                    asymptotic = FALSE,
                    quantile_val = -1*sqrt(-0.5*log(alpha))) {

  n <- length(data)
  if (is.null(k)) {
    k <- ceiling(n^(1/3))
  }
  
  if (is.null(sigma)) {
    data_block <- get_data_blocks(data, k)
    hat_L <- max(order(data_block)[1:J])
    l <- max(hat_L * k, l_min)
    sigma <- compute_sigma(data, l)
  }
  
  # \hat T_n
  sta <- min(cumsum(data - mean(data))) / (sigma * sqrt(n))
  
  
  # Asymptotic test
  if (asymptotic) {
    Reject <- sta < quantile_val
    p_val <- exp(-2*sta^2)
  } else {
    # Compute p-value using either provided simulation sample or new simulation
    
    if (is.null(simulation_sample)) {
      simulated_values <-
        replicate(simulation_size, compute_cusum_stat(n))
      p_val <- mean(simulated_values <= sta)
    } else {
      p_val <- mean(simulation_sample <= sta)
    }
    
    Reject <- p_val < alpha
  }
  

  return(list(
    p = p_val,
    REJECT = Reject,
    T = sta,
    sigma = sigma
  ))
}

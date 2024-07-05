
# Generates the trend mu_i used in Section 4
trend_simulations <- function(n, tau, tau_prime, tau_double_prime, mu_1, s) {
  
  # Initialize trend vector
  trend <- numeric(n)
  
  # First segment of the trend (t = 1, ..., tau-1)
  trend[1:(tau-1)] <- mu_1
  
  # Second segment of the trend (t = tau, ..., tau_prime)
  trend[tau:tau_prime] <- mu_1 + s * ((2 * (tau:tau_prime) - 3 * tau + tau_prime) / (tau_prime - tau))
  
  # Third segment of the trend (t = tau_prime + 1, ..., tau_double_prime)
  trend[(tau_prime + 1):tau_double_prime] <- mu_1 + s * (2 + exp(2 * (((tau_prime + 1):tau_double_prime) - tau_prime) / (tau_double_prime - tau_prime)))
  
  # Fourth segment of the trend (t = tau_double_prime + 1, ..., n)
  trend[(tau_double_prime + 1):n] <- mu_1 + s * (2 + exp(2) * ((2 * n - tau_double_prime - ((tau_double_prime + 1):n)) / (2 * n - 2 * tau_double_prime)))
  
  return(trend)
}


## Helper Function
# Generate the causal non-linear process Z'

Rcpp::cppFunction("NumericVector gen_z_prime(double theta,
                                             NumericVector epsilon) {

    int n = epsilon.length();
    NumericVector z = NumericVector(n);
    
    z[0] = 0;
    z[1] = 0;
    for (int i = 2; i < n; ++i) {
      z[i] = theta * (fabs(z[i - 1]) + fabs(z[i - 2])) + epsilon[i];
    }
    return z;
    
  }")

gen_z <- function(n, theta, estimated_mean, sd, burn_in = 10000) {
  z <- gen_z_prime(theta, rnorm(n + burn_in, 0, sd))[-(1:burn_in)] - estimated_mean
  return(z)
}

library(tidyverse)
source("../R/functions_for_data_generation.R")

### (1) Plot trend for Figure 3

n <- 800     # Total number of time points
tau <- 320   # Time point where second segment starts
taup <- 500
taupp <- 640 # Time point where third segment starts
mu1 <- 0     # Initial mean value for the first segment
s <- 0.5     # Given constant d

# Generate the trend
df <- data.frame(i = 1:n,
                 mu = trend_simulations(n, tau, taup, taupp, mu1, s))
df <- rbind(df, c(tau - 0.5, NA))

# Create the plot
p <- df %>%
  ggplot(aes(x = i, y = mu)) +
  geom_line() +
  labs(title = "",
       x = expression(t),
       y = expression(mu[t]))

# Save the plot as a PDF file
pdf(file = "../outPlots/illustration_trend.pdf", width = 8, height = 4)
print(p)
dev.off()


### (1.5) Plot an example of data and trend
# Define parameters
n <- 800       # Total number of time points
tau <- 320     # Start of second segment
taup <- 500
taupp <- 640   # Start of third segment
theta <- 0.4
sd <- 0.5
mu1 <- 0       # Initial mean value for the first segment
s <- 0.5       # Constant value d

# Load data
load("../1A_sim_mean_var/Z_mean_var.Rdata")

# Calculate adjusted mean and variance
theta_abs <- abs(theta)
z_prime_mean <- sign(theta) * sd * Z_mean_var[Z_mean_var$theta == theta_abs, ]$mean
z_prime_sigma <- sd * sqrt(Z_mean_var[Z_mean_var$theta == theta_abs, ]$var)

# Generate the trend
mu <- trend_simulations(n, tau, taup, taupp, mu1, s)
X_t <- mu + gen_z(n = n, theta = theta, estimated_mean = z_prime_mean, sd = sd)

# Prepare data frame
df <- data.frame(i = 1:n, mu = mu, X_t = X_t)

# Optional: adding NA row at tau-0.5, might need further clarification on its purpose
df <- rbind(df, c(tau - 0.5, NA, NA))

# Create the plot
library(ggplot2)
p <- ggplot(df, aes(x = i)) +
  geom_line(aes(y = X_t, color = "Data (X_t)")) +
  geom_line(aes(y = mu, color = "Trend (μ_t)"), linewidth = 1.5) +
  labs(title = "Data and Trend Visualization",
       x = expression(t),
       y = "") +
  scale_color_manual(values = c("Trend (μ_t)" = "black", "Data (X_t)" = "red"),
                     name = "Legend",
                     labels = c("Trend (μ_t)" = expression(mu[t]), "Data (X_t)" = expression(X[t]))) +
  scale_y_continuous(name = "", sec.axis = sec_axis(~., name = expression(mu[t])))



# Save the plot as a PDF file
pdf(file = "../outPlots/illustration_data.pdf", width = 8, height = 4)
print(p)
dev.off()

### (2) Load Results from Simulation experiments and format it

library(tidyverse)
library(ggpubr)

# Initialize empty data frames
data <- data.frame()
data2 <- data.frame()

# Load and process data from result files
for (i in 1:120) {
  path_1 <- paste0("res_files/result_alt_", i, ".Rdata")
  load(path_1)
  path_2 <- paste0("res_files/result_box_", i, ".Rdata")
  load(path_2)
  
  data2 <- rbind(data2, data_box)
  
  for (k in 1:3) {
    data <- rbind(data, c(result_alt[1:5], k,
                          result_alt[c(k+5, k+8, 12:15, k+15, k+18, k+21, k+24,
                                       28, k+28, k+31, k+34, k+37, k+40, k+43)]))  
  }
  
}

row.names(data) <- 1:nrow(data)

# Assign column names
colnames(data) <- c(
  "n",
  "tau",
  "s",
  "theta",
  "sigma",
  "k",
  "Our Method/ rho = 0.5",
  "Our Method w/ est. d/ rho = 0.5",
  "CUSUM",
  "LIKELIHOOD",
  "Binary Segmentation",
  "Binary Segmentation, w/ long_run",
  "Our Method/ rho = 0.25",
  "Our Method w/ est. d/ rho = 0.25",
  "Our Method/ rho = 0.75",
  "Our Method w/ est. d/ rho = 0.75",
  "Rejection Ratio",
  "NA Ratio Our Method/ rho = 0.5",
  "NA Ratio Our Method w/ est. d/ rho = 0.5",
  "NA Ratio Our Method w/ est. d/ rho = 0.25",
  "NA Ratio Our Method w/ est. d/ rho = 0.75"
)

# Convert 's' to a factor variable
data$s <- factor(data$s)

# Normalize the error values by dividing by 'n'
data$`Our Method/ rho = 0.5` = data$`Our Method/ rho = 0.5`/ data$n
data$`Our Method/ rho = 0.25` = data$`Our Method/ rho = 0.25`/ data$n
data$`Our Method/ rho = 0.75` = data$`Our Method/ rho = 0.75`/ data$n
data$`Our Method w/ est. d/ rho = 0.5` = data$`Our Method w/ est. d/ rho = 0.5`/ data$n
data$`Our Method w/ est. d/ rho = 0.25` = data$`Our Method w/ est. d/ rho = 0.25`/ data$n
data$`Our Method w/ est. d/ rho = 0.75` = data$`Our Method w/ est. d/ rho = 0.75`/ data$n

data$`CUSUM` = data$`CUSUM`/ data$n
data$`LIKELIHOOD` = data$`LIKELIHOOD`/ data$n
data$`Binary Segmentation, w/ long_run` = data$`Binary Segmentation, w/ long_run`/ data$n
data$`Binary Segmentation` = data$`Binary Segmentation`/ data$n

# Convert 'theta' to a factor variable with expression labels
data$theta <- factor(data$theta, 
                     labels = c(expression(theta == 0),
                                expression(theta == 0.2),
                                expression(theta == 0.3),
                                expression(theta == 0.4)))

# Convert 'n' to a factor variable with expression labels
data$n <- factor(data$n, 
                 labels = c(expression(n == 50),
                            expression(n == 100),
                            expression(n == 300),
                            expression(n == 500),
                            expression(n == 2000)))

# Convert 'k' to a factor variable with expression labels
data$k <- factor(data$k, labels = c(expression(n^(1/3)),
                                    expression(n^(1/2)),
                                    expression(n^(2/3))))

# Convert 's' to a factor variable in data2
data2$s <- factor(data2$s)

# Normalize the error values by dividing by 'n' in data2
data2$`median` = data2$`median`/ data2$n
data2$`first_quantile` = data2$`first_quantile`/ data2$n
data2$`third_quantile` = data2$`third_quantile`/ data2$n
data2$`lower_whisker` = data2$`lower_whisker`/ data2$n
data2$`upper_whisker` = data2$`upper_whisker`/ data2$n

# Convert 'theta' to a factor variable with expression labels in data2
data2$theta <- factor(data2$theta, 
                      labels = c(expression(theta == 0),
                                 expression(theta == 0.2),
                                 expression(theta == 0.3),
                                 expression(theta == 0.4)))

# Convert 'n' to a factor variable with expression labels in data2
data2$n <- factor(data2$n, 
                  labels = c(expression(n == 50),
                             expression(n == 100),
                             expression(n == 300),
                             expression(n == 500),
                             expression(n == 2000)))

library(ggplot2)

### Bar plots across methods

# Select and filter the data for bar plots
data3 <- data %>%
  select("n", "tau", "s", "theta", "k",
         "Our Method w/ est. d/ rho = 0.5",
         "CUSUM",
         "LIKELIHOOD",
         "Binary Segmentation",
         "Binary Segmentation, w/ long_run"
  ) %>%
  filter(k == "n^(1/3)") %>%
  filter(theta %in% c(expression(theta == 0),
                      expression(theta == 0.4)),
         s %in% c(0.4, 0.8)
  ) %>%
  pivot_longer(
    cols = 6:10,
    names_to = "method",
    values_to = "error")

# Rename factor levels for 'method'
data3$method <- factor(data3$method, levels = c(
  "Our Method w/ est. d/ rho = 0.5",
  "CUSUM",
  "LIKELIHOOD",
  "Binary Segmentation",
  "Binary Segmentation, w/ long_run"
), c(
  "Our Method w/ est. d",
  "CUSUM",
  "LIKELIHOOD",
  "Binary Segmentation, w/ marginal_variance",
  "Binary Segmentation, w/ long_run_variance"
))

# Create the bar plot
p1 <- ggplot(data3, aes(x = as.factor(s), y = error, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(theta ~ n,
             scale = "free_x",
             labeller = label_parsed) +
  labs(y = "MAE / n",
       x = "s",
       title = "") +
  scale_fill_discrete(name = "method", labels = c(expression(Proposed ~ method ~ hat(tau)),
                                                  "CUSUM",
                                                  "LIKELIHOOD",
                                                  "1SBS, marginal variance",
                                                  "1SBS, long run variance"
  )) +
  theme(legend.position="bottom")

# Save the bar plot as a PDF file
pdf(file = "../outPlots/bar_plot_five_methods.pdf", width = 8, height = 3.5)
print(p1)
dev.off()

### Bar plots across rhos

# Select and filter the data for bar plots across rhos
data6 <- data %>%
  select("n", "tau", "s", "theta", "k",
         "Our Method w/ est. d/ rho = 0.25",
         "Our Method w/ est. d/ rho = 0.5",
         "Our Method w/ est. d/ rho = 0.75"
  ) %>%
  filter(k == "n^(1/3)") %>%
  filter(theta %in% c(expression(theta == 0),
                      expression(theta == 0.4)),
         s %in% c(0.4, 0.8)
  ) %>%
  pivot_longer(
    cols = 6:8,
    names_to = "method",
    values_to = "error")

# Rename factor levels for 'method'
data6$method <- factor(data6$method, levels = c(
  "Our Method w/ est. d/ rho = 0.25",
  "Our Method w/ est. d/ rho = 0.5",
  "Our Method w/ est. d/ rho = 0.75"
), c(
  "rho = .25",
  "rho = .50",
  "rho = .75"
))

# Create the bar plot across rhos
p2 <- ggplot(data6, aes(x = as.factor(s), y = error, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_grid(theta ~ n,
             scale = "free_x",
             labeller = label_parsed) +
  labs(y = "MAE / n",
       x = "s",
       title = "") +
  scale_fill_discrete(name = "method", labels = c(expression(rho == .25),
                                                  expression(rho == .50),
                                                  expression(rho == .75)
  )) +
  theme(legend.position="bottom")

# Save the bar plot across rhos as a PDF file
pdf(file = "../outPlots/bar_plot_rho.pdf", width = 8, height = 3.5)
print(p2)
dev.off()

### Box plot across methods

# Filter the data for box plots
data4 <- data2 %>%
  filter(method %in% c("Our Method w/ est. d, rho = .50", "Binary Segmentation", 
                       "Binary Segmentation, long_run", "LIKELIHOOD", "CUSUM"),
         s %in% c(0.4, 0.8),
         theta %in% c(expression(theta == 0),
                      expression(theta == 0.4)
         )
  )

# Rename factor levels for 'method'
data4$method <- factor(data4$method, levels = c(
  "Our Method w/ est. d, rho = .50",
  "CUSUM",
  "LIKELIHOOD",
  "Binary Segmentation",
  "Binary Segmentation, long_run"
), c(
  "Our Method w/ est. d, rho = .50",
  "CUSUM",
  "LIKELIHOOD",
  "Binary Segmentation",
  "Binary Segmentation, long_run"
))

# Create the box plot
p3 <- ggplot(data4, aes(x = as.factor(s), fill = method)) +
  geom_boxplot(aes(ymin=lower_whisker, lower=first_quantile, middle=median,
                   upper=third_quantile, ymax=upper_whisker),
               stat = "identity") +
  facet_grid(theta ~ n,
             scale = "free_x",
             labeller = label_parsed) +
  labs(y = "AE / n",
       x = "s",
       title = "") +
  scale_fill_discrete(name = "method", labels = c(expression(Proposed ~ method ~ hat(tau)),
                                                  "CUSUM",
                                                  "LIKELIHOOD",
                                                  "1SBS, marginal variance",
                                                  "1SBS, long run variance"
  )) +
  theme(legend.position="bottom")

# Save the box plot as a PDF file
pdf(file = "../outPlots/box_plot_five_methods.pdf", width = 8, height = 6)
print(p3)
dev.off()

### Box plot across rhos

# Filter the data for box plots across rhos
data5 <- data2 %>%
  filter(method %in% c("Our Method w/ est. d, rho = .25",
                       "Our Method w/ est. d, rho = .50",
                       "Our Method w/ est. d, rho = .75"), 
         s %in% c(0.4, 0.8), theta %in% c(expression(theta == 0),
                                          expression(theta == 0.4)))

# Rename factor levels for 'method'
data5$method <- factor(data5$method, levels = c(
  "Our Method w/ est. d, rho = .25",
  "Our Method w/ est. d, rho = .50",
  "Our Method w/ est. d, rho = .75"
), c(
  "rho = .25",
  "rho = .50",
  "rho = .75"
))

# Create the box plot across rhos
p4 <- ggplot(data5, aes(x = as.factor(s), fill = method)) +
  geom_boxplot(aes(ymin=lower_whisker, lower=first_quantile, middle=median,
                   upper=third_quantile, ymax=upper_whisker),
               stat = "identity") +
  facet_grid(theta ~ n,
             scale = "free_x",
             labeller = label_parsed) +
  labs(y = "AE / n",
       x = "s",
       title = "") +
  scale_fill_discrete(name = "method", labels = c(expression(rho == .25),
                                                  expression(rho == .50),
                                                  expression(rho == .75)
  )) +
  theme(legend.position="bottom")

# Save the box plot across rhos as a PDF file
pdf(file = "../outPlots/box_plot_rho.pdf", width = 8, height = 6)
print(p4)
dev.off()

# Box plot across rhos with true d

# Filter the data for box plots across rhos
data6 <- data2 %>%
  filter(method %in% c("Our Method w/ true d, rho = .25",
                       "Our Method w/ true d, rho = .50",
                       "Our Method w/ true d, rho = .75"), 
         s %in% c(0.4, 0.8), theta %in% c(expression(theta == 0),
                                          expression(theta == 0.4)))

# Rename factor levels for 'method'
data6$method <- factor(data6$method, levels = c(
  "Our Method w/ true d, rho = .25",
  "Our Method w/ true d, rho = .50",
  "Our Method w/ true d, rho = .75"
), c(
  "rho = .25",
  "rho = .50",
  "rho = .75"
))

# Create the box plot across rhos
p5 <- ggplot(data6, aes(x = as.factor(s), fill = method)) +
  geom_boxplot(aes(ymin=lower_whisker, lower=first_quantile, middle=median,
                   upper=third_quantile, ymax=upper_whisker),
               stat = "identity") +
  facet_grid(theta ~ n,
             scale = "free_x",
             labeller = label_parsed) +
  labs(y = "AE / n",
       x = "s",
       title = "") +
  scale_fill_discrete(name = "method", labels = c(expression(rho == .25),
                                                  expression(rho == .50),
                                                  expression(rho == .75)
  )) +
  theme(legend.position="bottom")
  
# Save the box plot across rhos as a PDF file
pdf(file = "../outPlots/box_plot_rho_true_d.pdf", width = 8, height = 6)
print(p5)
dev.off()

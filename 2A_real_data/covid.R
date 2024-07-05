library(changepoint)
library(tidyverse)
library(zoo)

source("../R/functions_for_CPA.R")
source("../R/functions_for_classical_locating.R")

# start and latest end date to be used throughout
start_date <- as.Date("2019-10-01")
end_date <- as.Date("2020-01-31")

### (0) Prepare data and plot Figure 1

# Read the CSV files
coug_data <- read.csv("data/coug_dat.csv", header = FALSE)[-1, c(1, 4)]
fev_data <- read.csv("data/fev_dat.csv", header = FALSE)[-1, c(1, 4)]
colnames(coug_data) <- c("time", "index")
colnames(fev_data) <- c("time", "index")

# Convert 'time' columns to Date format
coug_data$time <- as.Date(coug_data$time, format = "%Y-%m-%d")
fev_data$time <- as.Date(fev_data$time, format = "%Y-%m-%d")

Baidu_data <- full_join(coug_data, fev_data, by = "time")
colnames(Baidu_data) <- c("time", "cough", "fever")

Baidu_plot <-
  Baidu_data %>%
  filter(between(time, start_date, end_date)) %>%
  pivot_longer(cols = c("cough", "fever"),
               names_to = "word", values_to = "index") %>%
  group_by(word)

p <-
  ggplot(Baidu_plot, aes(x = time, y = index)) +
  geom_line() +
  facet_wrap("word") +
  labs(title = "",
       x = "Date",
       y = "Index value")

pdf(file = "../outPlots/Baidu_index_data.pdf", width = 8, height = 4)
  print(p)
dev.off()

### (1) Analysis and Figure 6 and Figure 


analyse_and_plot <- function(start_date = "2019-10-01",
                             end_date = "2020-01-31") {
  
  # Filter out the data based on the condition
  coug_data <- coug_data[(coug_data$time >= start_date)&(coug_data$time <= end_date),]
  fev_data <- fev_data[(fev_data$time >= start_date)&(fev_data$time <= end_date),]
  
  # Get start and end dates for each data set
  coug_start <- min(coug_data$time)
  coug_end <- max(coug_data$time)
  fev_start <- min(fev_data$time)
  fev_end <- max(fev_data$time)
  
  # Testing
  test_fev <- test_CP(fev_data$index, asymptotic = TRUE)
  test_coug <- test_CP(coug_data$index, asymptotic = TRUE)

  # Our Method - fever
  fev_twostep <- IR_CPE(fev_data$index, J=3)
  fev_twostep_date <- fev_data$time[fev_twostep$hat_tau]
  
  # Our Method - cough
  coug_twostep <- IR_CPE(coug_data$index, J=3)
  coug_twostep_date <- coug_data$time[coug_twostep$hat_tau]

  
  # CUSUM
  fev_cusum_date <- fev_data$time[cusum_method(fev_data$index)]
  coug_cusum_date <- coug_data$time[cusum_method(coug_data$index)]
  
  # Likelihood based method
  fev_amoc_date <- fev_data$time[amoc_method(fev_data$index)]
  coug_amoc_date <- coug_data$time[amoc_method(coug_data$index)]
  
  # Binary Segmentation based method
  fev_binseg_date <- fev_data$time[binseg_method(fev_data$index, first = 1)]
  coug_binseg_date <- coug_data$time[binseg_method(coug_data$index, first = 1)]
  
  # Binary Segmentation long-run variance
  fev_zeta <- IR_CPE(fev_data$index, J=3)$hat_sigma * sqrt(2*log(length(fev_data$index)))
  coug_zeta <- IR_CPE(coug_data$index, J=3)$hat_sigma * sqrt(2*log(length(coug_data$index)))
  fev_binseg_lrv_date <- fev_data$time[binseg_zeta(fev_data$index, fev_zeta)]
  coug_binseg_lrv_date <- coug_data$time[binseg_zeta(coug_data$index, coug_zeta)]
  
  # Simple regression change point estimate
  fev_amoc_slope_date <- fev_data$time[slope_change_method(fev_data$index)]
  coug_amoc_slope_date <- coug_data$time[slope_change_method(coug_data$index)]
  
  
  
  library(ggplot2)
  # Using ggplot for cough data
  coug_plot <- ggplot(coug_data, aes(x = time, y = index)) +
    geom_line() +
    labs(title = "",
         x = "Date",
         y = "Index for Cough") +
    geom_vline(xintercept = coug_cusum_date, color = "orange", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = coug_twostep_date, color = "red", linetype = "solid", linewidth = 1) +
    geom_vline(xintercept = coug_binseg_date, color = "purple", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = coug_binseg_lrv_date, color = "green", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = coug_amoc_date, color = "brown", linetype = "dashed", linewidth = 1) +
    geom_text(aes(x = coug_cusum_date, y = max(index) * 0.9, label = paste("CUSUM:", coug_cusum_date)), color = "orange", angle = 45, size = 4) +
    geom_text(aes(x = coug_twostep_date, y = max(index) * 0.7, label = paste("Proposed Method:", coug_twostep_date)), color = "red", angle = 45, size = 4) +
    geom_text(aes(x = coug_binseg_date, y = max(index) * 0.6, label = paste("1SBS:", coug_binseg_date)), color = "purple", angle = 45, size = 4) +
    geom_text(aes(x = coug_binseg_lrv_date, y = max(index) * 0.4, label = paste("1SBS with LRV:", coug_binseg_lrv_date)), color = "green", angle = 45, size = 4) +
    geom_text(aes(x = coug_amoc_date, y = max(index) * 0.5, label = paste("AMOC:", coug_amoc_date)), color = "brown", angle = 45, size = 4)


  
  # Using ggplot for fever data
  fev_plot <- ggplot(fev_data, aes(x = time, y = index)) +
    geom_line() +
    labs(title = "",
         x = "Date",
         y = "Index for Fever") +
    geom_vline(xintercept = fev_cusum_date, color = "orange", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = fev_twostep_date, color = "red", linetype = "solid", linewidth = 1) +
    geom_vline(xintercept = fev_binseg_date, color = "purple", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = fev_binseg_lrv_date, color = "green", linetype = "dashed", linewidth = 1) +
    geom_vline(xintercept = fev_amoc_date, color = "brown", linetype = "dashed", linewidth = 1) +
    geom_text(aes(x = fev_cusum_date, y = max(index) * 0.9, label = paste("CUSUM:", fev_cusum_date)), color = "orange", angle = 45, size = 4) +
    geom_text(aes(x = fev_twostep_date, y = max(index) * 0.8, label = paste("Proposed Method:", fev_twostep_date)), color = "red", angle = 45, size = 4) +
    geom_text(aes(x = fev_binseg_date, y = max(index) * 0.5, label = paste("1SBS:", fev_binseg_date)), color = "purple", angle = 45, size = 4) +
    geom_text(aes(x = fev_binseg_lrv_date, y = max(index) * 0.4, label = paste("1SBS with LRV:", fev_binseg_lrv_date)), color = "green", angle = 45, size = 4) +
    geom_text(aes(x = fev_amoc_date, y = max(index) * 0.5, label = paste("AMOC:", fev_amoc_date)), color = "brown", angle = 45, size = 4)
    
  return(list(T_fev = test_fev$T,
              T_coug = test_coug$T,
              p_asy_fev = test_fev$p,
              p_asy_coug = test_coug$p,
              fev_plot = fev_plot,
              coug_plot = coug_plot,
              fev_twostep = fev_twostep,
              coug_twostep = coug_twostep))
}
# run analysis with our choice of start and end date
res_paper <- analyse_and_plot(start_date = "2019-10-01",
                 end_date = "2020-01-31")

pdf(file = "../outPlots/coug_I_seq.pdf", width = 4, height = 1.5)
  plot_I(res_paper$coug_twostep$I_seq, res_paper$coug_twostep$hat_eta,
         breaks = c(1,5,10,15,20,25))
dev.off()

pdf(file = "../outPlots/fev_I_seq.pdf", width = 4, height = 1.5)
  plot_I(res_paper$fev_twostep$I_seq, res_paper$fev_twostep$hat_eta,
         breaks = c(1,5,10,15,20,25))
dev.off()

library(gridExtra)
# Print the combined plot
pdf(file = "../outPlots/covid_comp.pdf", width = 9, height = 9)
# Combine cough_plot and fever_plot
  grid.arrange(res_paper$coug_plot, res_paper$fev_plot, ncol = 1)
dev.off()

###########################################################

fev_data <- fev_data %>%
  filter(between(time, start_date, end_date))

coug_data <- coug_data %>%
  filter(between(time, start_date, end_date))

# Then, the numbers in the test are:

res_paper$T_coug
res_paper$T_fev
res_paper$p_asy_coug
res_paper$p_asy_fev

# Results for cough, as in paper
res_paper$coug_twostep$n
res_paper$coug_twostep$k
res_paper$coug_twostep$hat_L
res_paper$coug_twostep$hat_mu0
res_paper$coug_twostep$hat_sigma
res_paper$coug_twostep$hat_eta
res_paper$coug_twostep$hat_mu1
res_paper$coug_twostep$hat_d
res_paper$coug_twostep$hat_tau

res_paper$fev_twostep$n
res_paper$fev_twostep$k
res_paper$fev_twostep$hat_L
res_paper$fev_twostep$hat_mu0
res_paper$fev_twostep$hat_sigma
res_paper$fev_twostep$hat_eta
res_paper$fev_twostep$hat_mu1
res_paper$fev_twostep$hat_d
res_paper$fev_twostep$hat_tau

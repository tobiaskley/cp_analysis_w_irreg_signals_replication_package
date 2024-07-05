library(tidyverse)
library(ggpubr)

# Loading data
data <- data.frame()

for (i in 1:2000) {
  path <- paste0("res_files/result_alt_", i, ".Rdata")
  load(path)
  data <- rbind(data, result_alt)
}

row.names(data) <- 1:nrow(data)

# Assigning column names
colnames(data) <- c(
  "n",
  "tau",
  "s",
  "theta",
  "sigma",
  "power_asy",
  "power_fin"
)


# null 

null <- data %>%
  filter(data$s == 0) %>%
  select(n, theta, power_asy, power_fin)

null %>%
  select(n, theta, power_asy) %>%
  group_by(n, theta) %>%
  spread("theta", "power_asy")

null %>%
  select(n, theta, power_fin) %>%
  group_by(n, theta) %>%
  spread("theta", "power_fin")

# power

library(ggplot2)
library(dplyr)

data$theta <- factor(data$theta)
data <- data %>%
  select(n, s, theta, power_asy, power_fin) %>%
  filter(data$n %in% c(50, 300, 500, 2000))

# Check the unique levels of `n`
unique_n <- unique(data$n)
print(unique_n)  # This will show you the unique values of `n`

# Adjust the labels based on the unique levels of `n`
data <- data %>%
  mutate(n_label = factor(n, levels = unique_n, labels = c("n = 50", "n = 300", "n = 500", "n = 2000")))

# power plot using asym quantile
p <- ggplot(data, aes(x=s, y=power_asy, color=theta)) + 
  geom_line(aes(linetype=theta)) +
  facet_wrap(~n_label, scales = "free") +
  ylab("power") + xlab("s") + ylim(0,1) +
  labs(color = expression(theta), linetype = expression(theta)) +
  theme(legend.position="bottom")

# Print the plot
pdf(file = "../outPlots/pow_plot_hat_sigma.pdf", width = 8, height = 6)
  print(p)
dev.off()

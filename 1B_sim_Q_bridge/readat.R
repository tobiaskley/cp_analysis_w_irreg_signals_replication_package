library(tidyverse)

# Loading data
data <- data.frame(n = numeric(),
                   perc = numeric(),
                   val = numeric())

alpha = c(0.2, 0.1, 0.05, 0.01)

for (i in 1:50) {
  path <- paste0("res_files/res_Q_bridge_", i, ".Rdata")
  load(path)
  for (j in 1:nrow(res)) {
    data <- rbind(data,
                  data.frame(n = rep(n, length(alpha)),
                                     perc = alpha,
                                     val = res[j,]
                             )
                  )
  }
}

Q_Bridge <- data %>%
  group_by(n, perc) %>%
  summarize(q_mean = mean(val),
            q_sd   = sd(val))

Q_Bridge$q_sd <- Q_Bridge$q_sd/sqrt(1000)

save(Q_Bridge, file = "Q_Bridge.Rdata")

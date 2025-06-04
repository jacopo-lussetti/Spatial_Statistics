#install.packages("patchwork")
# Load required packages
library(MASS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# Number of nodes
n <- 100

# Function to create RW1 precision matrix
rw1_precision <- function(n) {
  Q <- matrix(0, n, n)
  diag(Q) <- 2
  Q[cbind(1:(n-1), 2:n)] <- -1
  Q[cbind(2:n, 1:(n-1))] <- -1
  Q[1, 1] <- Q[n, n] <- 1  # Fix edges
  return(Q)
}

# Function to create RW2 precision matrix (with second-order penalties)
rw2_precision <- function(n) {
  Q <- matrix(0, n, n)
  for (i in 3:n) {
    Q[i, i] <- Q[i, i] + 1
    Q[i, i-1] <- Q[i, i-1] - 2
    Q[i, i-2] <- Q[i, i-2] + 1
    
    Q[i-1, i] <- Q[i-1, i] - 2
    Q[i-1, i-1] <- Q[i-1, i-1] + 4
    Q[i-1, i-2] <- Q[i-1, i-2] - 2
    
    Q[i-2, i] <- Q[i-2, i] + 1
    Q[i-2, i-1] <- Q[i-2, i-1] - 2
    Q[i-2, i-2] <- Q[i-2, i-2] + 1
  }
  return(Q)
}

# Build precision matrices
Q_rw1 <- rw1_precision(n)
Q_rw2 <- rw2_precision(n)

# Compute generalized inverses
Sigma_rw1 <- ginv(Q_rw1)
Sigma_rw2 <- ginv(Q_rw2)

# Marginal standard deviations
std_rw1 <- sqrt(diag(Sigma_rw1))
std_rw2 <- sqrt(diag(Sigma_rw2))

# Create tidy data for ggplot2
df <- data.frame(
  index = 1:n,
  RW1 = std_rw1,
  RW2 = std_rw2
) %>%
  pivot_longer(cols = c("RW1", "RW2"), names_to = "Model", values_to = "StdDev")

#split the data frame for each model

df_rw1 <- df %>% filter(Model == "RW1")
df_rw2 <- df %>% filter(Model == "RW2")

# RW1 plot
p1 <- ggplot(df_rw1, aes(x = index, y = StdDev)) +
  geom_line(color = "blue", size = 1) +
  labs(title = "RW1: First-order IGMRF", x = "Node Index", y = "Marginal Std Dev") +
  theme_minimal()

# RW2 plot
p2 <- ggplot(df_rw2, aes(x = index, y = StdDev)) +
  geom_line(color = "darkorange", size = 1) +
  labs(title = "RW2: Second-order IGMRF", x = "Node Index", y = "Marginal Std Dev") +
  theme_minimal()

print(p1 + p2)  #required for patchwork to combine plots

# Compute geometric means (reference std devs)
sigma_ref_rw1 <- exp(mean(log(std_rw1)))
sigma_ref_rw2 <- exp(mean(log(std_rw2)))

cat(sprintf("Reference Std Dev RW1: %.2f\n", sigma_ref_rw1))
cat(sprintf("Reference Std Dev RW2: %.2f\n", sigma_ref_rw2))

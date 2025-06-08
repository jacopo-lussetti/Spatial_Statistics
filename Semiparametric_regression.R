library(INLA)
library(ggplot2)
library(dplyr)
library(Matrix)
library(MASS)

#--------------------------
# 1. Simulate Data
#--------------------------
set.seed(123)
n <- 100
x <- seq(0, 1, length.out = n)
f_true <- sin(2 * pi * x)
y <- f_true + rnorm(n, sd = 0.2)
data <- data.frame(y = y, x = x, id = 1:n)

#--------------------------
# 2. Compute Generalized Variance for Unscaled RW1
#--------------------------
# Build unscaled RW1 structure matrix R
R <- bandSparse(n, k = c(0, 1), diagonals = list(rep(1, n), rep(-1, n-1)))
R <- t(R) %*% R

# Generalized inverse to get marginal variances
R_inv <- ginv(as.matrix(R))
g_i <- diag(R_inv)
g_bar <- mean(g_i)  # Use average marginal variance

#--------------------------
# 3. Define Gamma Prior for Unscaled Model Based on U and alpha
#--------------------------
U <- 0.1    # Desired upper bound on marginal SD
alpha <- 0.001 # P(SD > U)
a <- 1        # Gamma shape

# Solve for rate b such that P(tau < g_bar / U^2) = alpha
tau_thresh <- g_bar / U^2
b <- qgamma(alpha, shape = a) / tau_thresh

# Define hyperprior
hyper_prec_unscaled <- list(prec = list(prior = "gamma", param = c(a, b)))

#--------------------------
# 4. Define PC Prior for Scaled Model (as before)
#--------------------------
hyper_prec_scaled <- list(prec = list(prior = "pc.prec", param = c(1, 0.001)))

#--------------------------
# 5. Fit RW1 Models
#--------------------------
result_unscaled <- inla(
  y ~ f(id, model = "rw1", hyper = hyper_prec_unscaled),
  data = data, family = "gaussian"
)

result_scaled <- inla(
  y ~ f(id, model = "rw1", scale.model = TRUE, hyper = hyper_prec_scaled),
  data = data, family = "gaussian"
)

#--------------------------
# 6. Combine Data for Plotting
#--------------------------
plot_data <- data %>%
  mutate(f_true = f_true,
         rw1_unscaled = result_unscaled$summary.random$id$mean,
         rw1_scaled   = result_scaled$summary.random$id$mean)

#--------------------------
# 7. Plot Fits Using ggplot2
#--------------------------
ggplot(plot_data, aes(x = x)) +
  geom_point(aes(y = y), color = "gray60", size = 1.5) +
  geom_line(aes(y = f_true), color = "blue", linetype = "dashed", size = 1.1) +
  geom_line(aes(y = rw1_unscaled), color = "red", size = 1) +
  geom_line(aes(y = rw1_scaled), color = "darkgreen", size = 1) +
  labs(title = "RW1 Fit: Custom Prior (Unscaled) vs PC Prior (Scaled)",
       y = "y", x = "x") +
  theme_minimal()

#--------------------------
# 8. Plot Posterior of Precision (kappa)
#--------------------------
df_prec <- bind_rows(
  as.data.frame(result_unscaled$marginals.hyperpar[["Precision for id"]]) %>%
    mutate(Model = "Unscaled"),
  as.data.frame(result_scaled$marginals.hyperpar[["Precision for id"]]) %>%
    mutate(Model = "Scaled")
)

colnames(df_prec) <- c("Precision", "Density", "Model")

ggplot(df_prec, aes(x = Precision, y = Density, color = Model)) +
  geom_line(size = 1.2) +
  labs(title = "Posterior of Precision (kappa)",
       x = "Precision", y = "Density") +
  theme_minimal() +
  scale_color_manual(values = c("red", "darkgreen"))

library(INLA)

# Simulate data
set.seed(123)
n <- 100
x <- seq(0, 1, length.out = n)
f_true <- sin(2 * pi * x)
y <- f_true + rnorm(n, sd = 0.2)

# Create data
data <- data.frame(y = y, x = x, id = 1:n)

# RW1 model
formula <- y ~ f(id, model = "rw1", scale.model = TRUE)

# Fit model
result <- inla(formula, data = data, family = "gaussian")

# Plot
plot(x, y, main = "RW1 Fit to Data", col = "grey")
lines(x, f_true, col = "blue", lwd = 2, lty = 2)
lines(x, result$summary.random$id$mean, col = "red", lwd = 2)
legend("topright", legend = c("True f(x)", "Estimated f(x)"), col = c("blue", "red"), lty = c(2, 1), bty = "n")

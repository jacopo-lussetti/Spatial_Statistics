library(INLA)
library(sf)
library(spData)
library(spdep)
library(Matrix)
library(MASS)   # For ginv
library(parallel)
library(ggplot2)

# Load Boston housing tracts
boston.tr <- st_read(system.file("shapes/boston_tracts.gpkg", package = "spData")[1], quiet = TRUE)

# Create adjacency matrix
boston.adj <- poly2nb(boston.tr)
W.boston <- nb2mat(boston.adj, style = "B")
D <- Diagonal(ncol(W.boston), x = rowSums(W.boston))
Q <- as(D - W.boston, "dgCMatrix")

# ---- Compute reference standard deviation σ_ref ----
Q_dense <- as.matrix(Q)
Q_star <- ginv(Q_dense)
diag_vars <- diag(Q_star)
sigma_ref <- exp(mean(log(sqrt(diag_vars))))
cat("σ_ref =", round(sigma_ref, 4), "\n")

# ---- Rescale hyperprior ----
a <- 1         # Gamma shape
alpha <- 0.001 # Tail probability
U <- 0.5       # Desired upper bound of marginal std dev

q_alpha <- qgamma(alpha, shape = a, rate = 1)
b_scaled <- (U^2 * q_alpha) / sigma_ref^2
cat("Scaled b =", round(b_scaled, 4), "\n")

# ---- Prepare data ----
boston.tr$CMEDV2 <- boston.tr$CMEDV
boston.tr$CMEDV2[boston.tr$CMEDV2 == 50.0] <- NA
boston.tr$ID <- 1:nrow(boston.tr)

boston.form <- log(CMEDV2) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) +
  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT)

# ---- Fit unscaled Besag model ----
model1 <- inla(update(boston.form, . ~ . + f(ID, model = "besag", graph = W.boston)),
               data = as.data.frame(boston.tr),
               family = "gaussian",
               control.predictor = list(compute = TRUE),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))

# ---- Fit scaled Besag model ----
model2 <- inla(update(boston.form, . ~ . +
                        f(ID, model = "besag", graph = W.boston,
                          scale.model = TRUE,
                          hyper = list(prec = list(prior = "gamma", param = c(a, b_scaled))))),
               data = as.data.frame(boston.tr),
               family = "gaussian",
               control.predictor = list(compute = TRUE),
               control.compute = list(dic = TRUE, waic = TRUE, cpo = TRUE))

# ---- Transform precision -> standard deviation for the spatial effect ----
sd1 <- inla.tmarginal(function(tau) tau^(-0.5), model1$marginals.hyperpar[[2]])  # Precision for ID
sd2 <- inla.tmarginal(function(tau) tau^(-0.5), model2$marginals.hyperpar[[2]])  # Precision for ID

# ---- Plot posterior marginals of spatial standard deviation ----
ggplot() +
  geom_line(data = as.data.frame(sd1), aes(x = x, y = y, color = "Unscaled"), size = 1) +
  geom_line(data = as.data.frame(sd2), aes(x = x, y = y, color = "Scaled"), linetype = "dashed", size = 1) +
  scale_color_manual(values = c("Unscaled" = "blue", "Scaled" = "red")) +
  labs(title = "Posterior of Spatial Standard Deviation",
       y = "Density", x = "SD (1/sqrt(precision))", color = "Model") +
  theme_minimal()

# ---- Plot posterior marginals of log(LSTAT) coefficient ----
df1 <- as.data.frame(model1$marginals.fixed$`log(LSTAT)`)
df2 <- as.data.frame(model2$marginals.fixed$`log(LSTAT)`)

ggplot() +
  geom_line(data = df1, aes(x = x, y = y, color = "Unscaled"), size = 1) +
  geom_line(data = df2, aes(x = x, y = y, color = "Scaled"), linetype = "dashed", size = 1) +
  scale_color_manual(values = c("Unscaled" = "blue", "Scaled" = "red")) +
  labs(title = "Posterior of log(LSTAT) Coefficient",
       y = "Density", x = "Coefficient", color = "Model") +
  theme_minimal()

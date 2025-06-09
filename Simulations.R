
#Replicate original dataset
library(INLA)
library(ggplot2)
library(MASS)
library(spam)
# Load the dataset from INLA
data(Leuk)
###################Discretising wbc and tpi####################################
###############################################################################

n<-50 #number of cuts
#wbc histogram
wbc<-ggplot(data=Leuk, aes(x=wbc))+
  geom_histogram(bins=n, fill="blue", color="black") +
  labs(title="Histogram of WBC", x="WBC", y="Frequency")
wbc

#k equally sized sub-interval for wbc
Leuk$wbc_dst<-cut(Leuk$wbc, breaks=n, labels = FALSE)
wbc_range<-max(Leuk$wbc) - min(Leuk$wbc)
k_wbc<-wbc_range / n #size of each sub-interval

#tpi histogram
tpi<-ggplot(data=Leuk, aes(x=tpi))+
  geom_histogram(bins=n, fill="orange", color="black") +
  labs(title="Histogram of WBC", x="WBC", y="Frequency")
tpi

Leuk$tpi_dst<-cut(Leuk$tpi, breaks=50, labels = FALSE)

tpi_range<-max(Leuk$tpi) - min(Leuk$tpi)
k_tpi<-tpi_range / n #size of each sub-interval

##################### compute the reference s.d.###############################
###############################################################################

# ref s.d. for wbc RW1

R_wbc_1<-as.matrix(spam::precmat.RW1(n)) # U_A, sparse matrix that does not contain 0

Q_wbc_1<-(k_wbc)^-1 *R_wbc_1 # Q<-k^-1*R 
#Spectral Decomposition
eig_wbc_1 <- eigen(Q_wbc_1)#last vector

# Remove the 0 eigenvalue (RW1 has 1 zero eigenvalue)
V_wbc_1 <- eig_wbc_1$vectors[, 1:(n - 1)]
D_inv_wbc_1 <- diag(1 / eig_wbc_1$values[1:(n - 1)])

# Reconstruct the generalized inverse from spectral decomposition
Sigma_spec_wbc_1 <- V_wbc_1 %*% D_inv_wbc_1 %*% t(V_wbc_1)
gen_inv<-ginv(Q_wbc_1) #generalized inverse of R
max(abs(Sigma_spec_wbc_1-gen_inv)) #check that the two are equal
sd_ref_wbc_1<-exp(mean(log(sqrt(diag(Sigma_spec_wbc_1))))) #reference s.d. for wbc
print(paste("Reference SD wbc(rw1):", round(sd_ref_wbc_1, 2)))

#------------------ RW2: Reference SD Computation ----------------------

R_wbc_2 <- as.matrix(spam::precmat.RW2(n)) # Structure matrix R for RW2

Q_wbc_2 <- (k_wbc)^-3 * R_wbc_2  #equation (3) page 42

# Spectral Decomposition of Q
eig_wbc_2 <- eigen(Q_wbc_2)

# RW2 has 2 zero eigenvalues → remove them
V_wbc_2 <- eig_wbc_2$vectors[, 1:(n - 2)]
D_inv_wbc_2 <- diag(1 / eig_wbc_2$values[1:(n - 2)])

# Reconstruct the generalized inverse from spectral decomposition
Sigma_spec_wbc_2 <- V_wbc_2 %*% D_inv_wbc_2 %*% t(V_wbc_2)

# Optional: Validate against numerical generalized inverse
gen_inv <- ginv(Q_wbc_2)
max(abs(Sigma_spec_wbc_2 - gen_inv))  # Should be close to zero

# Compute reference standard deviation (geometric mean of marginal sds)
sd_ref_wbc_2 <- exp(mean(log(sqrt(diag(Sigma_spec_wbc_2)))))
print(paste("Reference SD wbc(rw2):", round(sd_ref_wbc_2, 2)))


###################### RW1 for tpi (using ginv) #############################

R_tpi_1 <- as.matrix(spam::precmat.RW1(n))  # Structure matrix R for RW1
Q_tpi_1 <- (k_tpi)^-1 * R_tpi_1             # Q = k⁻¹ * R, RW1 scaling

# Generalized inverse via MASS
Sigma_tpi_1 <- ginv(Q_tpi_1)

# Reference standard deviation (geometric mean of marginal std devs)
sd_ref_tpi_1 <- exp(mean(log(sqrt(diag(Sigma_tpi_1)))))
print(paste("Reference SD tpi (RW1):", round(sd_ref_tpi_1, 2)))


###################### RW2 for tpi (using ginv) #############################

R_tpi_2 <- as.matrix(spam::precmat.RW2(n))  # Structure matrix R for RW2
Q_tpi_2 <- (k_tpi)^-3 * R_tpi_2             # Q = k⁻³ * R, RW2 scaling

# Generalized inverse via MASS
Sigma_tpi_2 <- ginv(Q_tpi_2)

# Reference standard deviation (geometric mean of marginal std devs)
sd_ref_tpi_2 <- exp(mean(log(sqrt(diag(Sigma_tpi_2)))))
print(paste("Reference SD tpi (RW2):", round(sd_ref_tpi_2, 2)))



# Hyperprior parameters
a <- 1           # shape
b <- 5e-5        # inverse-scale (rate)
alpha <- 0.001   # tail probability

# Function to compute U
compute_U <- function(sigma_ref, a, b, alpha) {
  q_gamma <- qgamma(alpha, shape = a, scale = 1)
  sqrt((b * sigma_ref^2) / q_gamma)
}

# Compute U values
U_wbc_RW1 <- compute_U(sd_ref_wbc_1, a, b, alpha)
U_wbc_RW2 <- compute_U(sd_ref_wbc_2, a, b, alpha)
U_tpi_RW1 <- compute_U(sd_ref_tpi_1, a, b, alpha)
U_tpi_RW2 <- compute_U(sd_ref_tpi_2, a, b, alpha)

# Print results
cat("\n--- Upper Limits U ---\n")
cat(sprintf("WBC RW1: %.4f\n", U_wbc_RW1))
cat(sprintf("WBC RW2: %.4f\n", U_wbc_RW2))
cat(sprintf("TPI RW1: %.4f\n", U_tpi_RW1))
cat(sprintf("TPI RW2: %.4f\n", U_tpi_RW2))




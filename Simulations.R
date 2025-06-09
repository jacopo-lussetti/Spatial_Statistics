
#Replicate original dataset
library(INLA)
library(ggplot2)
library(MASS)
library(spam)
# Load the dataset from INLA
data(Leuk)
load("LeukSurv.rda") # Load the Leuk dataset
Leuk<-LeukSurv
write.csv(Leuk, "LeukSurv.csv", row.names = FALSE) # Save the dataset as CSV
#--------------  first step: compute the reference s.d.-----------------------
##----- #discretised wbc and tpi

wbc<-ggplot(data=Leuk, aes(x=wbc))+
  geom_histogram(bins=50, fill="blue", color="black") +
  labs(title="Histogram of WBC", x="WBC", y="Frequency")
wbc

Leuk$wbc_dst<-cut(Leuk$wbc, breaks=50, labels = FALSE)

n<-50
R_RW1_wbc<-spam::precmat.RW1(n) # U_A, sparse matrix that does not contain 0
R_wbc <- as.matrix(R_RW1_wbc)
Q<-0.1*R_wbc #as we discretised the wbc, we will have new kth intervals, therefore


#--------------  first step: compute the reference s.d.-----------------------
##----- #discretised wbc and tpi

wbc<-ggplot(data=Leuk, aes(x=wbc))+
  geom_histogram(bins=50, fill="blue", color="black") +
  labs(title="Histogram of WBC", x="WBC", y="Frequency")
wbc

Leuk$wbc_dst<-cut(Leuk$wbc, breaks=50, labels = FALSE)

#-----tpi 

tpi<-ggplot(data=Leuk, aes(x=tpi))+
  geom_histogram(bins=50, fill="orange", color="black") +
  labs(title="Histogram of WBC", x="WBC", y="Frequency")
tpi
Leuk$tpi_dst<-cut(Leuk$tpi, breaks=50, labels = FALSE)

#precision matrix Q= tau (1) x R. so we need to compute then the structure matrix R
#with rank n-1, where n is the number of nodes (here 50)
n<-50
R <- matrix(0, n, n)
for (i in 1:(n - 1)) {
  R[i, i] <- R[i, i] + 1
  R[i, i + 1] <- R[i, i + 1] - 1
  R[i + 1, i] <- R[i + 1, i] - 1
  R[i + 1, i + 1] <- R[i + 1, i + 1] + 1
}
#compute Q
Q<- R 
#Replicate original dataset
library(INLA)
library(ggplot2)
library(MASS)
# Load the dataset from INLA
data(Leuk)

#--------------  first step: compute the reference s.d.-----------------------
##----- #discretised wbc and tpi

wbc<-ggplot(data=Leuk, aes(x=wbc))+
  geom_histogram(bins=50, fill="blue", color="black") +
  labs(title="Histogram of WBC", x="WBC", y="Frequency")
wbc

Leuk$wbc_dst<-cut(Leuk$wbc, breaks=50, labels = FALSE)

#-----tpi 

tpi<-ggplot(data=Leuk, aes(x=tpi))+
  geom_histogram(bins=50, fill="orange", color="black") +
  labs(title="Histogram of WBC", x="WBC", y="Frequency")
tpi
Leuk$tpi_dst<-cut(Leuk$tpi, breaks=50, labels = FALSE)
range(Leuk$wbc)
#precision matrix Q= tau (1) x R. so we need to compute then the structure matrix R
#with rank n-1, where n is the number of nodes (here 50)
n<-50
R <- matrix(0, n, n)
for (i in 1:(n - 1)) {
  R[i, i] <- R[i, i] + 1
  R[i, i + 1] <- R[i, i + 1] - 1
  R[i + 1, i] <- R[i + 1, i] - 1
  R[i + 1, i + 1] <- R[i + 1, i + 1] + 1
}
#compute Q
Q<- R * (1/500)

#now we are able to compute the reference s.d., thata are defined in formula (7) 
#from Sorbye, and Rue (2014)

#Spectral Decomposition

eig <- eigen(Q)#last vector

# Remove the 0 eigenvalue (RW1 has 1 zero eigenvalue)
V <- eig$vectors[, 1:(n - 1)]
D_inv <- diag(1 / eig$values[1:(n - 1)])

# Reconstruct the generalized inverse from spectral decomposition
Sigma_spec_wbc <- V %*% D_inv %*% t(V)
gen_inv<-ginv(Q) #generalized inverse of R
max(abs(Sigma_spec_wbc-gen_inv)) #check that the two are equal
sd_ref_wbc<-exp(mean(log(sqrt(diag(Sigma_spec_wbc))))) #reference s.d. for wbc
sd_ref_wbc

#RW2
R_RW2_wbc <- spam::precmat.RW2(n)
R_2wbc_RW2 <- as.matrix(R_RW2_wbc)

# Define Q = τ * R (assume τ = 0.1 here)
tau <- 0.001
Q_RW2 <- tau * R_2wbc_RW2

# Spectral decomposition
eig_RW2 <- eigen(Q_RW2)

# RW2 has 2 zero eigenvalues → keep 1:(n-2)
V_RW2 <- eig_RW2$vectors[, 1:(n - 2)]
D_inv_RW2 <- diag(1 / eig_RW2$values[1:(n - 2)])

# Generalized inverse from spectral decomposition
Sigma_spec_wbc_RW2 <- V_RW2 %*% D_inv_RW2 %*% t(V_RW2)

# Optional check against MASS::ginv (not sparse or efficient)
gen_inv <- ginv(Q_RW2)
max_diff <- max(abs(Sigma_spec_wbc_RW2 - gen_inv))  # should be close to 0
print(paste("Max difference:", round(max_diff, 10)))

# Compute reference standard deviation (geometric mean of marginal sds)
sd_ref_wbc_RW2 <- exp(mean(log(sqrt(diag(Sigma_spec_wbc_RW2)))))
print(paste("Reference SD (RW2):", round(sd_ref_wbc_RW2, 4)))



#-----tpi ----------------

tpi<-ggplot(data=Leuk, aes(x=tpi))+
  geom_histogram(bins=50, fill="orange", color="black") +
  labs(title="Histogram of WBC", x="WBC", y="Frequency")
tpi
Leuk$tpi_dst<-cut(Leuk$tpi, breaks=50, labels = FALSE)

#precision matrix Q= tau (1) x R. so we need to compute then the structure matrix R
#with rank n-1, where n is the number of nodes (here 50)
n<-50

R_tpi_1 <- spam::precmat.RW1(n = n)# U_A, sparse matrix that does not contain 0
R_tpi_1 <- as.matrix(R_tpi_1)
# as we discretised the wbc and tpi, we will have new kth intervals, therefore 
#we need to adjust the precision matrix Q accordingly accordigly to equation from
#equation (2)
rg <- range(Leuk$tpi)
diff <- rg[2]-rg[1]
k <- diff / 50 
Q_tpi_1<- (k*1)^-1 * R_tpi_1
## First compute the intervals between the cuts


#now we are able to compute the reference s.d., thata are defined in formula (7) 
#from Sorbye, and Rue (2014)
#We rist derive covariance matrix using Spectral decomposition of Q_RW1
eig_tpi_1 <- eigen(Q_tpi_1)#last vector

# Remove the 0 eigenvalue (RW1 has 1 zero eigenvalue)
V_tpi_1 <- eig_tpi_1$vectors[, 1:(n - 1)]
D_inv_tpi_1 <- diag(1 / eig_tpi_1$values[1:(n - 1)])

# Reconstruct the generalized inverse from spectral decomposition
Sigma_spec_tpi_1 <- V %*% D_inv_tpi_1 %*% t(V_tpi_1)
gen_inv_tpi_1<-ginv(Q_tpi_1) #generalized inverse of R
max(abs(Sigma_spec_wbc-gen_inv_tpi_1)) #check that the two are equal
sd_ref_tpi_1<-exp(mean(log(sqrt(diag(gen_inv_tpi_1))))) #reference s.d. for wbc
sd_ref_tpi_1

#tpi----RW2


#precision matrix Q= tau (1) x R. so we need to compute then the structure matrix R
#with rank n-1, where n is the number of nodes (here 50)
n<-50

R_tpi_2 <- spam::precmat.RW2(n = n)# U_A, sparse matrix that does not contain 0
R_tpi_2 <- as.matrix(R_tpi_2)
# as we discretised the wbc and tpi, we will have new kth intervals, therefore 
#we need to adjust the precision matrix Q accordingly accordigly to equation from
#equation (2)
rg <- range(Leuk$tpi)
diff <- rg[2]-rg[1]
k <- diff / 50 
Q_tpi_2<- (k*1)^-3 * R_tpi_2

gen_inv_tpi_2<-ginv(Q_tpi_2) #generalized inverse of R

sd_ref_tpi_2<-exp(mean(log(sqrt(diag(gen_inv_tpi_2))))) #reference s.d. for wbc
sd_ref_tpi_2

#define U for all the cases

#Chat _GPT----------------------------------------------------------
# Given reference SDs
'
# sd_ref_wbc_RW1  <- 8.68
# sd_ref_wbc_RW2  <- 458.08
# sd_ref_tpi_RW1  <- 1.55
# sd_ref_tpi_RW2  <- 2.68
'
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
U_wbc_RW1 <- compute_U(sd_ref_wbc, a, b, alpha)
U_wbc_RW2 <- compute_U(sd_ref_wbc_RW2, a, b, alpha)
U_tpi_RW1 <- compute_U(sd_ref_tpi_1, a, b, alpha)
U_tpi_RW2 <- compute_U(sd_ref_tpi_2, a, b, alpha)

# Print results
cat("\n--- Upper Limits U ---\n")
cat(sprintf("WBC RW1: %.4f\n", U_wbc_RW1))
cat(sprintf("WBC RW2: %.4f\n", U_wbc_RW2))
cat(sprintf("TPI RW1: %.4f\n", U_tpi_RW1))
cat(sprintf("TPI RW2: %.4f\n", U_tpi_RW2))



# Print the computed U values
library(INLA)
library(Matrix)
library(MASS)

# Load data
#Leuk <- read.csv("LeukSurv.csv")

# Step 1: Build spatial adjacency structure
# Extract district ID
district_ids <- sort(unique(Leuk$district))
n_districts <- length(district_ids)

# Create spatial neighborhood matrix: simple assumption of adjacent district IDs
# In real cases, use a shapefile or `spdep::poly2nb`
adjacency <- spam::specmat.RW2(n_districts)  # Using RW1 for simplicity
# Step 2: Structure matrix R for Besag model
R_besag <- -adjacency
diag(R_besag) <- rowSums(adjacency)

# Step 3: Q = τ * R with τ = 1
Q_besag <- R_besag

# Step 4: Compute generalized inverse using eigenvalue decomposition (remove 1 zero eig.)
eig_besag <- eigen(Q_besag)
V_besag <- eig_besag$vectors[, 1:(n_districts - 1)]
D_inv_besag <- diag(1 / eig_besag$values[1:(n_districts - 1)])
Sigma_besag <- V_besag %*% D_inv_besag %*% t(V_besag)

# Step 5: Compute reference sigma
sd_ref_besag <- exp(mean(log(sqrt(diag(Sigma_besag)))))
print(paste("Reference SD (Besag):", round(sd_ref_besag, 4)))

# Step 6: Compute U (using same a, b, alpha as before)
compute_U <- function(sigma_ref, a, b, alpha) {
  q_gamma <- qgamma(alpha, shape = a, scale = 1)
  sqrt((b * sigma_ref^2) / q_gamma)
}
a <- 1
b <- 5e-5
alpha <- 0.001

U_besag <- compute_U(sd_ref_besag, a, b, alpha)
print(paste("Upper limit U (Besag):", round(U_besag, 4)))





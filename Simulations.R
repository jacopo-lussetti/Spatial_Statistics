
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

<<<<<<< HEAD
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
#We rist derive covariance matrix using Spectral decomposition of Q_RW1
=======
#Spectral Decomposition
>>>>>>> 3a9753dab58562e55e74c3b13b5953eb851d243a
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

# Extract marginal standard deviations
marginal_sd_rw2 <- sqrt(diag(Sigma_spec_wbc_RW2))

# Plot: Figure 2 for RW2 applied to WBC
plot(Leuk$wbc, marginal_sd_rw2, type = "l", lwd = 2, col = "blue",
     xlab = expression(i), ylab = expression(sigma[i]),
     main = "Marginal SD for RW2 IGMRF (Discretized WBC)")




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
Q_tpi_1<-3.33*R_tpi_1
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
sd_ref_tpi_1<-exp(mean(log(sqrt(diag(Sigma_spec_tpi_1))))) #reference s.d. for wbc
sd_ref_tpi_1

#tpi----




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




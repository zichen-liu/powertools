## code to prepare `multisite.data` dataset goes here
m <- 20
J <- 10
N <- m*J

gamma00 <- 0
gamma10 <- 3
sigma_u0_sq <- 4
sigma_u1_sq <- 8
sigma_eps_sq <- 36

# generate REs for sites

set.seed(51124)

u0 <- rnorm(n=J,mean=0,sd=sqrt(sigma_u0_sq))
u1 <- rnorm(n=J,mean=0,sd=sqrt(sigma_u1_sq))

# make dataframe

id <- as.factor(seq(1:N))
i <- as.factor(rep(seq(1:m),J))
j <- as.factor(rep(seq(1:J),each=m))
xij <- rep(c(rep(-0.5,m/2),rep(0.5,m/2) ), J)
u0j <- rep(u0, each=m)
u1j <- rep(u1, each=m)
eij <- rnorm(n=N, mean=0, sd=sqrt(sigma_eps_sq))

# this Y includes site by group interaction:
Yij <- u0j + (gamma10+u1j)*xij + eij
# this Y has no interaction:
Wij <- u0j + gamma10*xij + eij

group <- xij+1.5

multisite.data <- data.frame(cbind(id,j,i,group,xij,Yij,Wij))

usethis::use_data(multisite.data, overwrite = TRUE)

library(CircStats)
library(rstan)

rstan_options(auto_write = TRUE)

graphics.off()

T = 300
theta.rho = 0.6
theta.mu = 0
step.shape  <- 5
step.scale <- 1
p.obs <- 1


step.length <- c()
theta <-c()
xy <- data.frame(x = 30, y = 20, obs = rbinom(1,1,p.obs))

dev.new()
xx <- seq(from=0,to=20,by=0.1)
yy <- dgamma(xx,shape = step.shape, scale=step.scale)
plot(xx,yy,type='l')

dev.new()
aa <- seq(from=-pi,to=pi,by=0.01)
yy <- dwrpcauchy(theta=aa,mu=theta.mu,rho=theta.rho)
plot(aa,yy,type='l')

for (n in c(1:T))
{
  theta[n+1] <- rwrpcauchy(1, location=theta.mu, rho=theta.rho)
  step.length[n+1] <- rgamma(1,shape = step.shape, scale = step.scale)
  xy[n+1,1:2] <- xy[n,1:2] + c(cos(theta[n+1])*step.length[n+1],sin(theta[n+1])*step.length[n+1])
  xy[n+1,3] <- rbinom(1,1,p.obs)
}


dev.new()
hist(step.length)

dev.new()
plot(xy[,1],xy[,2],type="o")
j <- (xy$obs == 1)
points(xy[j,1],xy[j,2], type="o", pch=22, col="red")


## Run STAN model

stan.code <- "

data {
  int T;
  real x[T];
  real y[T];
}

parameters {
  real <lower=-pi(), upper=pi()> theta_mu;
  real <lower=0> theta_kappa;
  real <lower=0> step_alpha;
  real <lower=0> step_beta;
  real <lower=0> stp[T];
  real <lower=0> theta[T];
  real <lower=0> sigma;
}


model {
  
  for (t in 2:T) {
    stp[t] ~ gamma(step_alpha,step_beta);
    theta[t] ~ von_mises(theta_mu, theta_kappa);
    x[t] ~ normal(x[t-1] + (stp[t] * cos(theta[t])),sigma) ;
    y[t] ~ normal(y[t-1] + (stp[t] * sin(theta[t])),sigma) ;
  }
}
"

stan.control <- list(adapt_delta = 0.9, max_treedepth = 20)

stan.data <- list(T=nrow(xy), x = xy[,1], y = xy[,2] )
fit <- stan(
  model_code = stan.code,
  data = stan.data,
  control = stan.control,
  chains = 1,
  warmup = 1000,
  thin = 3,
  iter = 2000,
  cores = 4,
  refresh=500
)

arr <- extract(fit, permuted = FALSE) 

dev.new(); plot(density(arr$theta_mu))
dev.new(); plot(density(arr$theta_kappa))
dev.new(); plot(density(arr$step_alpha))
dev.new(); plot(density(arr$step_beta))



print(fit,pars = c('theta_kappa','theta_mu','step_alpha','step_beta','sigma'))
plot(fit,pars = c('theta_kappa','theta_mu','step_alpha','step_beta','sigma'))


mcmc_acf(x, pars = c("alpha", "beta[1]"))


dev.new()
yy = dvm(aa,mu,kappa)
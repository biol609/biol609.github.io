#adapted from http://seananderson.ca/2014/04/08/gamma-glms.html and
# http://seananderson.ca/2014/05/18/gamma-hurdle.html
set.seed(1)
x <- 1:150
y <- rbinom(length(x), size = 1, prob = 0.7)
seas_m <- runif(length(x), 0.1, 5)
shape = 0.4
y_true = exp(2 - 0.5 * seas_m)
y <- y*rgamma(length(x), rate = shape / y_true, shape = shape)
non_zero <- ifelse(y > 0, 1, 0)
d <- data.frame(days_at_sea = y, seas_m = seas_m, non_zero = non_zero)

head(d)
plot(days_at_sea ~ seas_m, data=d)



#Fits ok....
library(rethinking)
mod_fish <- alist(
  #zero inflated gamma likelihood
  days_at_sea ~ dzagamma2(p, mu, scale),
  
  #days at sea
  log(mu) <- a + b*seas_m,
  
  #priors
  a ~ dnorm(0,100),
  b ~ dnorm(0,100),
  scale ~ dexp(2),
  p ~ dunif(0,1)
)

fit_fish <- map(mod_fish, data=d)



##OR....
dzagamma <- function (x, prob, shape, rate, log = FALSE) 
{
  K <- as.data.frame(cbind(x = x, prob = prob, shape = shape, rate = rate))
  llg <- dgamma(x, shape = shape, rate = rate, log = TRUE)
  ll <- ifelse(K$x == 0, log(K$prob), log(1 - K$prob) + llg)
  if (log == FALSE) 
    ll <- exp(ll)
  ll
}



mod_fish_1 <- alist(
  #zero inflated gamma likelihood
  days_at_sea ~ dzagamma(p, shape, rate),
  
  #days at sea
  rate <- shape/exp(log_mu),
  log_mu <- a + b*seas_m,
  
  #priors
  a ~ dnorm(0,100),
  b ~ dnorm(0,100),
  shape ~ dexp(2),
  p ~ dunif(0,1)
)

fit_fish_1 <- map(mod_fish_1, data=d)
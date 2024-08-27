library(brms)
library(rstan)
data <- readRDS("stan.dl_pred.RDS")
prior1 <- prior(normal(0,10), class = Intercept) +
  prior(cauchy(0,2), class = sd)
data$ind <- as.factor(data$ind)

ind <- data$ind_z
z <- data$z
data <- data.frame(ind,z)

data$ind <- as.factor(data$ind)
CompDHGLM <- brm(z ~ (1|ind), sigma ~ (1|ind),
                     data = data,
                     family=gaussian,
                     prior = prior1,
                     iter  = 1000, thin = 2,
                     chains = 2, cores = 4, seed = 12345)

get_stancode(CompDHGLM)

summary(CompDHGLM)
0.85^2


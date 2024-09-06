
library(tidyverse)
library(TMB)
library(rema)
options(scipen = 999)
theme_set(theme_minimal(base_size = 17))

genlogit <- function(x, low=0, upp=1){
  p = (x-low)/(upp-low)
  return(log(p/(1-p)))
}

geninvlogit <- function(x, low=0, upp=1){
  pprime = exp(x)/(1+exp(x))
  return(pprime*(upp-low)+low)
}

ci = function(par, se, p=0.975, lo = 0, hi = 1, type = "I", k = 1){
  ci = par + c(-1,1) * qnorm(0.975) * se
  if(type == "I") {
    return(c(par, se, ci))
  }
  if(type == "exp") {
    return(c(exp(par), exp(par)*se, exp(ci)))
  }
  if(type == "expit") { #Delta-method: V(lo + (hi-lo)/(1 + exp(-x))) ~ ((hi-lo) * p * (1-p))^2 * V(x)
    p = 1/(1 + exp(- k * par))
    dm.se = k * abs(hi-lo)*p*(1-p)*se
    return(c(geninvlogit(par, lo, hi), dm.se, lo + (hi-lo)/(1 + exp(-ci))))
  }
}

compile("src/ar1.cpp")
dyn.load(dynlib("src/ar1"))

rpn <- read_csv('data/bsai_turbot_rpns.csv') %>% 
  mutate(cv = sqrt(rpn_var) / rpn) %>% 
  select(-rpn_var) %>% 
  rename(strata = fmp, obs = rpn)

model_yrs = unique(rpn$year)

rpn <- rpn %>% 
  tidyr::expand(strata, year = model_yrs) %>% 
  left_join(rpn)

# data <- list(model_yrs = model_yrs,
#              obs = rpn %>% filter(strata == area) %>% pull(obs) %>% as.matrix(),
#              obs_cv = rpn %>% filter(strata == area) %>% pull(cv) %>% as.matrix())
data <- list(model_yrs = model_yrs,
             obs = rpn %>% 
               pivot_wider(id_cols = year, names_from = strata, values_from = obs) %>% 
               select(-year) %>% 
               as.matrix(),
             obs_cv = rpn %>% 
               pivot_wider(id_cols = year, names_from = strata, values_from = cv) %>% 
               select(-year) %>% 
               as.matrix())

par <- list(log_PE = rep(log(1), 2), 
            logit_rho = rep(genlogit(0.5, -1, 1), 2),
            log_pred = log(apply(X = data$obs, MARGIN = 2, FUN = zoo::na.approx, maxgap = 100, rule = 2)))

map <- par
map$log_PE <- as.factor(1:length(map$log_PE))
map$log_PE <- as.factor(rep(1, length(map$log_PE)))
map$logit_rho <- as.factor(1:length(map$logit_rho))
map$logit_rho <- as.factor(rep(1, length(map$logit_rho)))
map$log_pred <- as.factor(1:length(map$log_pred))

mod <- MakeADFun(data, par, map, random = "log_pred", DLL = "ar1")
# input$data = mod$simulate(complete=TRUE)
opt = nlminb(mod$par, mod$fn, mod$gr)

mod$rep = mod$report()
mod$rep$pred
mod$sdrep = sdreport(mod)
summary(mod$sdrep)
mod$env$parList()$log_pred

#check Laplace Approximation
# check <- checkConsistency(mod)
# summary(check)

pars = as.list(mod$sdrep, "Est")
sd = as.list(mod$sdrep, "Std")

ci(pars$log_PE[1], sd$log_PE, type = "exp") 
ci(pars$logit_rho, sd$logit_rho, lo = -1, hi = 1, type = "expit", k = 2) 
# rho is ~1, likely this model is overparameterized

# AI
pred <- matrix(data = NA, nrow = nrow(pars$log_pred), ncol = 4)
for(i in 1:nrow(pars$log_pred)) pred[i,] <- ci(pars$log_pred[i,1], sd$log_pred[i,1], type = "exp")
pred <- as.data.frame(pred); names(pred) <- c('pred', 'se', 'lci', 'uci')
preddf <- pred %>% mutate(year = model_yrs, strata = 'Aleutians')
# Bering
pred <- matrix(data = NA, nrow = nrow(pars$log_pred), ncol = 4)
for(i in 1:nrow(pars$log_pred)) pred[i,] <- ci(pars$log_pred[i,2], sd$log_pred[i,2], type = "exp")
pred <- as.data.frame(pred); names(pred) <- c('pred', 'se', 'lci', 'uci')
preddf <- preddf %>% bind_rows(pred %>% mutate(year = model_yrs, strata = 'Bering Sea'))

out <- rpn %>% 
  mutate(log_obs = ifelse(obs > 0, log(obs), NA),
         sd_log_obs = ifelse(obs > 0, sqrt(log(cv^2 + 1)), NA),
         obs_lci = exp(log_obs - qnorm(0.975) * sd_log_obs),
         obs_uci = exp(log_obs + qnorm(0.975) * sd_log_obs)) %>% 
  left_join(preddf) %>% 
  mutate(final_rpn_estimate = ifelse(is.na(obs), pred, obs))

out %>% 
  ggplot(aes(x = year, y = obs)) +
  geom_point() +
  geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci)) +
  geom_line(aes(y = pred), col = 'red') +
  geom_ribbon(aes(ymin = lci, ymax = uci), col = NA, fill = 'red', alpha = 0.2) +
  facet_wrap(~strata)


out %>% 
  mutate(next_obs = lead(obs, 1),
         pdiff = 100 * abs(next_obs - pred)/((next_obs + pred)/2)) %>% 
  group_by(strata) %>% 
  mutate(bsai_statistic = mean(pdiff, na.rm = TRUE)) %>% 
  print(n=Inf)

compare <- out %>% 
  select(strata, year, obs, obs_lci, obs_uci) %>% 
  left_join(
    read_csv('results/compare_methods_v1.csv') %>% 
      filter(method != 'statusquo_meanratio_olddata') %>%
      select(year, strata = fmp, final_rpn_estimate = rpn, method)
  )  %>% 
  mutate(method = ifelse(method == 'statusquo_meanratio_newdat', 'statquo_meanratio', 'linear_approx')) %>% 
  bind_rows(out %>% 
              mutate(method = 'ar1') %>% 
              select(strata, year, obs, obs_lci, obs_uci, final_rpn_estimate, method)) %>% 
  arrange(method, strata, year) %>% 
  group_by(method, strata) %>% 
  mutate(next_obs = lead(obs, 1),
         pdiff = 100 * abs(next_obs - final_rpn_estimate)/((next_obs + final_rpn_estimate)/2))

compare %>% summarize(pdiff_statistic = mean(pdiff, na.rm = TRUE)) %>% 
  select(strata, method, pdiff_statistic) %>% 
  arrange(strata, pdiff_statistic)

compare %>% 
  ggplot(aes(x = year, y = obs)) +
  geom_line(aes(y = final_rpn_estimate, col = method, lty = method), size = 1) +
  geom_point(aes(y = final_rpn_estimate, col = method, shape = method), size = 2) +
  geom_point() +
  geom_errorbar(aes(ymin = obs_lci, ymax = obs_uci)) +
  facet_wrap(~strata, ncol = 1) +
  ggthemes::scale_color_colorblind()


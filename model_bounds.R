
rm(list = ls())
gc()
.libPaths("~/Rlib")
library(cmdstanr)
library(jsonlite)
library(posterior)
library(bridgestan)

library(ggplot2)
library(scales)
library(latex2exp)

setwd("~/Code/vi_variance")
source("tools.r")
bridge_stan_path <- "../bridgestan/"

###############################################################################
## Functions to run experiment

bridgestan_model <- function(stan_file, data_file, bridge_stan_path) {
  current_wd <- getwd()
  model_so <- paste0(substr(stan_file, start = 0, 
                     stop = nchar(stan_file) - 5), "_model.so")
  
  setwd(bridge_stan_path)
  system(paste0("make ", current_wd, "/", model_so))
  
  setwd(current_wd)
  bridge_model <- StanModel$new(model_so, data_file, 1234, 0)
  
  return(bridge_model)
}

run_experiment <- function(stan_file, data_file, # n_dim, 
                           name, n_chains = 8,
                           n_sample = 2000, seed = 123) {
  stan_data <- read_json(data_file, simplifyVector = TRUE)
  mod <- cmdstan_model(stan_file)

  fit_mcmc <- mod$sample(data = stan_data, chains = n_chains, parallel_chains = n_chains,
                         iter_warmup = 1000, iter_sampling = n_sample, seed = seed,
                         adapt_delta = 0.95)
  
  fit_mcmc$save_object(file = file.path("deliv", paste0(name, ".mcmc.fit.RDS")))

  fit_advi <- mod$variational(data = stan_data, seed = seed, tol_rel_obj = 1e-3,
                              iter = 1e4, output_samples = 10000) # default max num of iter is 1e3

  fit_advi$save_object(file = file.path("deliv", paste0(name, ".advi.fit.RDS")))
  
  # Construct bridgestan model to get transformation to unconstrained scale
  bridge_model <- bridgestan_model(stan_file, data_file, bridge_stan_path)
  n_dim <- bridge_model$param_num()
  
  # Get MCMC draws on the unconstrained scale
  draws <- fit_mcmc$draws()[, , 2:(1 + n_dim)]
  draw_array_const <- array(draws, dim = c(n_sample * n_chains, n_dim))
  draw_array <- array(NA, dim = c(n_sample * n_chains, n_dim))
  for (i in 1:dim(draw_array)[1]) {
    draw_array[i, ] <- bridge_model$param_unconstrain(draw_array_const[i, ])
  }
  
  # Now put ADVI draws on the unconstrained scale
  draw_advi_const <- fit_advi$draws()[, 3:(2 + n_dim)]
  draw_advi <- array(NA, dim = dim(draw_advi_const))
  for (i in 1:dim(draw_advi)[1]) {
    draw_advi[i, ] <- bridge_model$param_unconstrain(draw_advi_const[i, ])
  }
  
  cov_mcmc = var(draw_array)
  cov_advi = var(draw_advi)
  var_mcmc = diag(cov_mcmc)
  var_advi = diag(cov_advi)

  S_advi = var_mcmc / var_advi   
  cov_mcmc_inv = solve(cov_mcmc)
  psi_gauss = 1 / diag(cov_mcmc_inv)
  S_gauss = diag(cov_mcmc) / psi_gauss
  log_det_sigma = log(det(cov_mcmc))
  
  return_list = list(var_mcmc = var_mcmc,
                     var_advi = var_advi,
                     psi_gauss = psi_gauss,
                     S_advi = S_advi,
                     S_gauss = S_gauss,
                     log_det_sigma = log_det_sigma)
  
  write(toJSON(return_list), file = file.path("deliv", paste0(name, ".json")))
  
  return(return_list)
}

if (FALSE) {

###############################################################################
## Eight Schools

data_file <- "data/eight_schools.json"

list_schools_nc = 
  run_experiment(stan_file = "stan_models/eight_schools_noncentered.stan",
                 data_file = data_file,
                 name = "8schools_nc")

list_schools_pool =
  run_experiment(stan_file = "stan_models/eight_schools_pool.stan",
                 data_file = data_file,
                 name = "8schools_pool")


###############################################################################
## Disease map of Finland

# disease_data <- read_json("data/disease_100.json", simplifyVector = TRUE)
data_file <- "data/disease_100.json"
list_disease_map =
  run_experiment(stan_file = "stan_models/disease_map.stan",
                 data_file = data_file, # n_dim = 102,
                 name = "disease_map")
 

###############################################################################
## Sparse kernel interaction model on prostate cancer model

data_file <- "data/prostate_200.json"
list_skim =
  run_experiment(stan_file = "stan_models/skim_logit.stan",
                 data_file = data_file, name = "SKIM")
                 # n_dim = 305)


###############################################################################
## GLM binomial
list_glm =
  run_experiment(stan_file = "stan_models/GLM_Binomial_model.stan",
                 data_file = "data/GLM_Binomial_data.json",
                 name = "glm_binomial")


###############################################################################
## Mixture of normals
list_mixture =
  run_experiment(stan_file = "stan_models/mixture.stan",
                 data_file = "data/low_dim_gauss_mix.json",
                 name = "Mixture")

# Correct variance will not be estimated by MCMC since target is bimodal.
# Compute true variance analytically...
mu = 1.5; sigma = 1;
list_mixture$var_mcmc <- rep(sigma + 0.25 * (2 * mu)^2, 2)
list_mixture$S_gauss <- c(1, 1)
list_mixture$S_advi <- list_mixture$var_mcmc / list_mixture$var_advi

write(toJSON(list_mixture),
      file = file.path("deliv", paste0("Mixture", ".json")))
}

###############################################################################
## Load saved lists

list_schools_nc <- read_json(file.path("deliv", "8schools_nc.json"), simplifyVector = TRUE)
list_schools_pool <- read_json(file.path("deliv", "8schools_pool.json"), simplifyVector = TRUE)
list_disease_map <- read_json(file.path("deliv", "disease_map.json"), simplifyVector = TRUE)
list_skim <- read_json(file.path("deliv", "SKIM.json"), simplifyVector = TRUE)
list_glm <- read_json(file.path("deliv", "glm_binomial.json"), simplifyVector = TRUE)
list_mixture <- read_json(file.path("deliv", "Mixture.json"), simplifyVector = TRUE)

###############################################################################
## Plot elements of shrinkage matrix for a single model

n_dim = 10
param <- factor(rep(paste0("z[", 1:n_dim, "]"), 2),
                levels = paste0("z[", 1:n_dim, "]"))

tex_levels <- rep(NA, n_dim)
for (i in 1:n_dim) {
  tex_levels[i] <- TeX(paste0("$z_{", i, "}$"))
}


S <- unname(c(S_advi = list_schools_nc$S_advi,
       S_gauss = list_schools_nc$S_gauss))
target <- rep(c("Posterior", "Gaussian"), each = n_dim)

plot.data <- data.frame(S = S, param = param, target = target)

p <- ggplot(data = plot.data, aes(x = param, y = S, fill = target)) +
  geom_bar(stat = "identity", color="black",
           position=position_dodge(), width = 0.75) +
  geom_hline(yintercept = 1.0, size = 0.5, linetype = "dashed") +
  # coord_flip() +
  ylab(TeX("$S_{ii} = \\Sigma_{ii} / \\Psi_{ii}$")) +
  xlab("Parameter") + theme_bw() +
  # theme(axis.text.x = element_text(angle = 15)) +
  # theme(legend.title=element_blank()) +
  theme(legend.position = c(0.15, 0.75)) +
  theme(text = element_text(size = 20)) +
  coord_cartesian(ylim = c(0.5, 1.8)) +
  scale_x_discrete(labels = tex_levels)
p
  


###############################################################################
## Plots comparing Shrinkage matrices

mean_trace = c(mean(list_glm$S_advi),
               mean(list_schools_pool$S_advi),
               mean(list_schools_nc$S_advi),
               mean(list_disease_map$S_advi),
               mean(list_skim$S_advi),
               mean(list_mixture$S_advi),
               mean(list_glm$S_gauss),
               mean(list_schools_pool$S_gauss),
               mean(list_schools_nc$S_gauss),
               mean(list_disease_map$S_gauss),
               mean(list_skim$S_gauss),
               mean(list_mixture$S_gauss))

plot_data <- data.frame(mean_trace = mean_trace,
                         model = factor(rep(c("glm_binomial",
                                                 "8schools_pool",
                                                 "8schools_nc",
                                                 "disease_map",
                                                 "SKIM",
                                                 "mixture"), 2),
                          level = c("glm_binomial",
                                     "8schools_pool", "8schools_nc",
                                     "disease_map", "SKIM", "mixture")),
                          target = rep(c("Posterior", "Gaussian"), 
                                       each = 6))


p <- ggplot(plot_data, aes(x = model, y = mean_trace, color = target,
                           fill = target)) + 
  geom_bar(stat="identity", color="black",
           position=position_dodge(), width = 0.5) + # coord_flip() +
  scale_y_sqrt() + 
  ylab("Trace(S) / n") + xlab(" ") +
  geom_hline(yintercept = 1.0, size = 0.5, linetype = "dashed") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 15)) +
  theme(legend.position = c(0.5, 0.7)) +
  theme(text = element_text(size = 20)) +
  coord_cartesian(ylim = c(0.5, 6.5))
p


##########################################################################
## Relative entropy between FG-VI targeting a Gaussian and a posterior

relative_difference <- function(list_model) {
  (sum(log(list_model$psi_gauss)) - sum(log(list_model$var_advi))) /
    length(list_model$var_mcmc)
}

entropy_empirical_bound <- function(list_model) {
  (list_model$log_det_sigma - (sum(log(list_model$psi)))) /
     length(list_model$var_mcmc)
}

relative_difference(list_glm)
relative_difference(list_schools_nc)
relative_difference(list_schools_pool)
relative_difference(list_disease_map)
relative_difference(list_skim)
relative_difference(list_mixture)

entropy_empirical_bound(list_glm)
entropy_empirical_bound(list_schools_nc)
entropy_empirical_bound(list_schools_pool)
entropy_empirical_bound(list_disease_map)
entropy_empirical_bound(list_skim)


# For mixture, don't trust MCMC (use var_mcmc computed analytically!)
sum(log(list_mixture$var_mcmc)) - sum(log(list_mixture$psi)) /
  length(list_mixture$var_mcmc)

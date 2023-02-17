
m(list = ls())
gc()
.libPaths("~/Rlib")
library(ggplot2)
library(scales)
library(nlme)
library(latex2exp)

# devtools::install_github("zeehio/facetscales")

setwd("~/Code/vi_variance")
source("tools.r")

custom_cols <- c("#888888", "#619CFF", "#AA3929")

###############################################################################
## Ellipse plots (Figure 1)

n <- c(2, 64, 64, 64)
eps <- c(0.75, 0.75, 0.9, 0.999)

n_sim = 1000
x <- array(NA, dim = c(length(eps) * n_sim, 2))
x_vi <- array(NA, dim = c(length(eps) * n_sim, 2))

for (j in 1:length(eps)) {
  # Construct target covariance
  C <- matrix(eps[j], nrow = n[j], ncol = n[j])
  diag(C) <- 1
  L <- chol(C)
  
  # VI approximation
  Psi <- 1 / diag(solve(C))
  
  # generate samples
  for (i in 1:n_sim) {
    z <- rnorm(n = n[j], mean = 0, sd = 1)
    w <- t(L) %*% z
    x[(j - 1) * n_sim + i,] <- w[1:2]
    x_vi[(j - 1) * n_sim + i, ] <- rnorm(n = 2, mean = 0, sd = sqrt(Psi[1:2]))
  }
}

plot.data <- data.frame(x = c(x[, 1], x_vi[, 1]),
                        y = c(x[, 2], x_vi[, 2]),
                        t = rep(c("target", "VI"), each = n_sim * length(eps)),
                        label = factor(rep(rep(paste0("eps = ", eps, ", n = ", n), 
                                               each = n_sim), 2),
                                       levels = c("eps = 0.75, n = 2",
                                                  "eps = 0.75, n = 64",
                                                  "eps = 0.9, n = 64",
                                                  "eps = 0.999, n = 64"))
)

levels(plot.data$label) <- c(TeX("$\\epsilon = 0.75$, $n = 2$"),
                             TeX("$\\epsilon = 0.75$, $n = 64$"),
                             TeX("$\\epsilon = 0.9$, $n = 64$"),
                             TeX("$\\epsilon = 0.999$, $n = 64$"))


p <- ggplot(data = plot.data, aes(x, y, fill = t)) +
  stat_ellipse(geom = "polygon", alpha = 0.5) + 
  theme_bw() + facet_wrap(~label, ncol = 4, labeller = label_parsed) +
  theme(text = element_text(size = 25)) +
  theme(legend.title=element_blank()) +
  theme(legend.position = c(0.95, 0.2)) +
  xlab(TeX("$z_1$")) + ylab(TeX("$z_2$"))
p

#####################################################################
## Utility functions for next figures

eps <- seq(from = 0, to = 0.95, by = 0.01)
n <- rep(c(10, 100, 200), each = length(eps))  # should be 10, 100, 200
eps <- rep(eps, 3)
set.seed(1954)
x = sort(runif(max(n), 0, 200))

bounds_comp <- function(eps, n, matrix_type, x) {
  # options for matrix type: Constant, Circulant, Kernel
  
  log_det_C <- rep(NA, length(eps))
  log_det_S <- rep(NA, length(eps))
  trace_S <- rep(NA, length(eps))
  
  joint_entropy_ub <- rep(NA, length(eps))
  entropy_ub <- rep(NA, length(eps))
  shrinkage_ub <- rep(NA, length(eps))
  delinkage_lb <- rep(NA, length(eps))
  condition_R <- rep(NA, length(eps))
  
  trace_lb <- rep(NA, length(eps))
  trace_ub <- rep(NA, length(eps))
  
  for (i in 1:length(eps)) {
    ## Generate covariance matrix (scaled, so acts as corr matrix)
    if (matrix_type == "Constant") {
      C <- matrix(eps[i], nrow = n[i], ncol = n[i])
    } else if (matrix_type == "Circulant") {
      C <- matrix(NA, nrow = n[i], ncol = n[i])
      
      for (j in 2:n[i]) {
        for (k in 1:(j - 1)) {
          C[k, j] = eps[i]^((j - k)^2)
          C[j, k] = C[k, j]
        }
      }
    } else if (matrix_type == "Kernel") {
      C <- matrix(NA, nrow = n[i], ncol = n[i])
      
      for (j in 2:n[i]) {
        for (k in 1:(j - 1)) {
          C[k, j] = exp(-(x[j] - x[k])^2 / eps[i]^2)
          C[j, k] = C[k, j]
        }
      }
    }
    diag(C) <- 1
    
    ## Obtain exact VI solution
    C_inv = solve(C)
    Psi = 1 / diag(C_inv)
    S = diag(C) / Psi
    
    log_det_C[i] <- determinant(C, logarithm = TRUE)$modulus[1]
    log_det_S[i] <- sum(log(S))
    trace_S[i] <- sum(S)
    
    ## Get the condition number
    lambda <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
    R <- max(lambda) / min(lambda)
    condition_R[i] <- R
    
    ub_list <- entropy_bound(R, n[i])
    joint_entropy_ub[i] <- ub_list$joint_entropy_ub
    entropy_ub[i] <- ub_list$entropy_ub
    shrinkage_ub[i] <- ub_list$shrinkage_ub
    delinkage_lb[i] <- ub_list$delinkage_lb
    
    trace_list <- trace_bound(R, n[i])
    trace_lb[i] <- trace_list$LbFun
    trace_ub[i] <- trace_list$UbFun
  }
  
  
  ## Save results in a data set for plots
  legend_levels = c("shrinkage", "delinkage", "entropy loss", "trace(S)",
                    "joint_entropy_ub", "entropy_ub",
                    "shrinkage_ub", "delinkage_lb",
                    "trace(S)_lb", "trace(S)_ub")
  n_lines = length(legend_levels)
  
  plot.data <- data.frame(eps = rep(eps, n_lines),
                          y = c(0.5 * log_det_S,                 # shrinkage
                                - 0.5 * log_det_C,               # delinkage
                                0.5 * (log_det_C + log_det_S),   # entropy loss
                                trace_S,                         # trace(S)
                                joint_entropy_ub,                # joint entropy bound
                                entropy_ub,                      # entropy upper bound
                                shrinkage_ub,                    # shrinkage upper bound
                                - delinkage_lb,                  # delinkage lower bound
                                trace_lb,                        # trace(S) lower bound
                                trace_ub),                       # trace(S) upper bound 
                          n = paste0("n = ", n),
                          R = rep(condition_R, n_lines),
                          legend = factor(rep(legend_levels, each = length(eps)),
                                          levels = legend_levels))

  return(plot.data)
}

plot_tradeoff <- function(plot.data, n_dim, matrix_type) {
  if (matrix_type == "Kernel") line_color = custom_cols[3]
  if (matrix_type == "Constant") line_color = custom_cols[2]
  
  n_dim_text <- paste0("n = ", n_dim)
  plot.data.subset <- plot.data[(plot.data$n == n_dim_text) & 
                                  (plot.data$legend == "shrinkage" | 
                                     plot.data$legend == "delinkage"), ]
  
  p <- ggplot(data = plot.data.subset, aes(x = eps, y = y, linetype = legend)) +
    geom_line(size = 1.5, color = line_color) + theme_bw() +
    ylab("") +
    theme(text = element_text(size = 20)) +
    theme(legend.title=element_blank()) +
    facet_wrap(~matrix_label) +
    theme(legend.text.align = 0)
  
  if (matrix_type == "Kernel") {
    p <- p + xlab(TeX("$\\rho$")) +
      theme(legend.position = c(0.275, 0.84)) +
      geom_segment(aes(x = 0.8, 
                       y = plot.data.subset[plot.data.subset$legend == "delinkage",]$y[80],
                       xend = 0.8,
                       yend = plot.data.subset[plot.data.subset$legend == "shrinkage", ]$y[80]
      ), color = "black", arrow = arrow(length = unit(0.03, "npc"), 
                                        ends = "both")) +
      annotate(geom = "text", x = 0.91, y = 7.5, label = "H(p) - H(q)",
               size = 4, parse = TRUE) +
      scale_linetype_discrete(breaks = c("shrinkage", "delinkage"),
                              labels = unname(TeX(c("Shrinkage: $log|S| / 2$",
                                                    "Delinkage: $log|C|^{-1} / 2"))))
    # scale_color_discrete(breaks = c("shrinkage", "delinkage"),
    #                      labels = unname(TeX(c("Shrinkage: $log|S| / 2$",
    #                                            "Delinkage: $log|C|^{-1} / 2"))))
  } else {
    p <- p + xlab(TeX("$\\epsilon")) + theme(legend.position="none")
  }
  
  p
}

###############################################################################
## Plot the shrinkage and delinkage terms (Figure 2)

plot.data.const <- bounds_comp(eps, n, matrix_type = "Constant", x)
plot.data.const$matrix_label <- rep("Constant off-diagonal", 
                                    length(plot.data.const$y))

p <- plot_tradeoff(plot.data.const, n_dim = 10, matrix_type = "Constant")

pdf(file = file.path(getwd(), "figures/entropy_loss_constant.pdf"), 
    width = 6, height = 4)
p
dev.off()

plot.data.kernel <- bounds_comp(eps, n, matrix_type = "Kernel", x)
plot.data.kernel$matrix_label <- rep("Squared exponential kernel")

p <- plot_tradeoff(plot.data.kernel, n_dim = 10, matrix_type = "Kernel")
pdf(file = file.path(getwd(), "figures/entropy_loss_kernel.pdf"), 
    width = 6, height = 4)
p
dev.off()



###############################################################################
## Plot bounds on shrinkage and delinkage (Figure 3)

n_sub = "n = 10"
if (n_sub == "n = 10") {
  xmax = 200
  ymax = 20
} else {
  xmax = 2000
  ymax = 300
}

plot.data.forces.10 <- plot.data.forces[plot.data.forces$n == n_sub, ]
plot.data.const.10 <- plot.data.forces.const[plot.data.forces.const$n == n_sub, ]
plot.data.bounds.10 <- plot.data.bounds[plot.data.bounds$n == n_sub, ]

p3 <- ggplot(
  data = plot.data.forces.10[plot.data.forces.10$legend == "shrinkage", ]) +
  geom_line(aes(x = R, y = y,
                color = "Squared exponential kernel"), 
            size = line_size, linetype = "dashed") +
  theme_bw() + ylab("") +
  theme(text = element_text(size = 20)) +
  theme(legend.title=element_blank()) +
  facet_wrap(~n, scales = "free") +
  scale_x_continuous(trans='log10', limits = c(1, xmax)) +
  ylim(0, ymax) +
  geom_line(aes(
    x = plot.data.forces.10[plot.data.forces.10$legend == "delinkage", ]$R,
    y = plot.data.forces.10[plot.data.forces.10$legend == "delinkage", ]$y,
    color = "Squared exponential kernel"), linetype = "solid", size = line_size) +
  geom_line(aes(
    x = plot.data.const.10[plot.data.forces.10$legend == "shrinkage", ]$R,
    y = plot.data.const.10[plot.data.forces.10$legend == "shrinkage", ]$y,
    color = "Constant off-diagonal"), linetype = "dashed", size = line_size) +
  geom_line(aes(
    x = plot.data.const.10[plot.data.forces.10$legend == "delinkage", ]$R,
    y = plot.data.const.10[plot.data.forces.10$legend == "delinkage", ]$y,
    color = "Constant off-diagonal"), linetype = "solid", size = line_size) +
  geom_line(aes(
    x = plot.data.bounds.10[plot.data.forces.10$legend == "shrinkage", ]$R,
    y = plot.data.bounds.10[plot.data.forces.10$legend == "shrinkage", ]$y,
    color = "Bounds"), linetype = "dashed", size = line_size) +
  geom_line(aes(
    x = plot.data.bounds.10[plot.data.forces.10$legend == "delinkage", ]$R,
    y = plot.data.bounds.10[plot.data.forces.10$legend == "delinkage", ]$y,
    color = "Bounds"), linetype = "solid", size = line_size) +
  scale_color_manual(values = custom_cols)
p3

# show difference in entropy
data.segm <- data.frame(x = 100,   # segment for entropy loss
                        xend = 100,
                        y = 2.5, 
                        yend = 5, 
                        legend = "shrinkage", n = "n = 10")

data.segm_ub <- data.frame(x = 150,   # segment on bounds of entropy loss
                           xend = 150,
                           y = 2, 
                           yend = 18, 
                           legend = "shrinkage", n = "n = 10")

data.anotate <- data.frame(x = 90, y = 18, label = "Bound on \n H(p) - H(q)",
                           legend = "shrinkage", n = "n = 10")

if (n_sub == "n = 10") {
  p3 <- p3 + geom_segment(data = data.segm, 
                          aes(x = x, y = y, xend = xend, yend = yend),
                          inherit.aes = FALSE,
                          arrow = arrow(length = unit(0.03, "npc"), 
                                        ends = "both"))  #+
  # geom_segment(data = data.segm_ub, 
  #              aes(x = x, y = y, xend = xend, yend = yend),
  #              inherit.aes = FALSE,
  #              linetype = "twodash",
  #              arrow = arrow(length = unit(0.03, "npc"), 
  #                            ends = "both")) +
  # geom_text(data = data.anotate, aes(x = x, y = y, label = label),
  #           inherit.aes = FALSE) + xlab(TeX("$R$"))
  p3 <- p3 + theme(legend.position = c(0.35, 0.75))
} else {
  p3 <- p3 + theme(legend.position="none")
}

p3


pdf(file = file.path(getwd(),
                     paste0("figures/entropy_bound_", n_sub, ".pdf")), 
    width = 6,
    height = 4)

p3
dev.off()



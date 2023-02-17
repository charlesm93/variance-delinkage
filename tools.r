
###############################################################################
## Bounds on the variance shrinkage

eigen_bounds <- function(C) {
  lambda <- eigen(C, symmetric = TRUE, only.values = TRUE)$values
  
  lambda_1 <- min(lambda)
  lambda_n <- max(lambda)
  
  kappa = lambda + 1 / lambda
  kappa_1 = min(kappa)
  kappa_n = max(kappa)
  
  return(list(
    lambda_bound = c(C[1, 1] * lambda_1, C[1, 1] * lambda_n),
    kappa_bound = c(C[1, 1] / (kappa_n - 1), C[1, 1] / (kappa_1 - 1))
  ))
}


###############################################################################
## Bounds on the trace of the shrinkage matrix

Obj_lb <- function(n, R, lambda) {
  - (n - 2)^2 / (n - (1 + R) * lambda) - 1 / lambda - 1 / (R * lambda)
}

Obj_ub <- function(n, R, lambda, k) {
  (n - k) / lambda + (k - 1) / (R * lambda) + 1 / (n - (R * (k - 1) + (n - k)) * lambda)
}

trace_bound <- function(R, n) {
  ## compute lower bound
  a = R * (n - 2)^2 - (1 + R)^2
  b = 2 * n * (1 + R)
  c = - n^2
  lambda_a <- (-b + sqrt(b^2 - 4 * a * c)) / (2 * a)
  lambda_b <- (-b - sqrt(b^2 - 4 * a * c)) / (2 * a)
  
  # Check if the solutions verify the boundary conditions.
  boundary_l <- n / (R * (n - 1) + 1)
  boundary_u <- n / (R + n - 1)
  bc_a <- boundary_l <= lambda_a & lambda_a <= boundary_u
  bc_b <- boundary_l <= lambda_b & lambda_b <= boundary_u
  
  if (bc_a & !bc_b) LbFun <- Obj_lb(n, R, lambda_a)
  if (!bc_a & bc_b) LbFun <- Obj_lb(n, R, lambda_b)
  if (bc_a & bc_b) LbFun <- max(Obj_lb(n, R, lambda_a), Obj_lb(n, R, lambda_b))
  if (!bc_a & !bc_b) LbFun <- max(Obj_lb(n, R, boundary_l),
                                  Obj_lb(n, R, boundary_u))
  
  ## Combute upper bound
  for (k in 1:n) {
    lambda_a <- n / (R * k + (n - k))
    lambda_b <- n / (R * (k - 1) + (n - k + 1))
    Ub_k <- max(Obj_ub(n, R, lambda_a, k), Obj_ub(n, R, lambda_b, k))
    if (k == 1) {
      UbFun <- Ub_k
    } else {
      UbFun <- max(UbFun, Ub_k)
    }
  }
  
  return(list(LbFun = - LbFun, UbFun = UbFun))
}

#####################################################################
## Upper bound for entropy loss (Algorithm 1 + computation of the
## joint bound on entropy)

ObjF <- function(lambda, R, n, k) {
  return ((n - k + (k - 1) / R) / lambda + 
    1 / (n - (R * (k - 1) + n - k) * lambda))
}

ObjG <- function(lambda, R, n) {
  return (log(lambda) + log(R * lambda) + 
            (n - 2) * log((n - (1 + R) * lambda) / (n - 2)))
}

ObjE <- function(p, R, n, k) {
  return (- (k - 1) * log(p) - (n - k) * log(R * p) -
            log(1 - (k - 1 + R * (n - k)) * p))
}

entropy_bound <- function(R, n) {
  # Bound the shrinkage
  for (k in 1:n) {
    lambda_a <- n / (R * k + n - k)
    lambda_b <- n / (R * (k - 1) + n - k + 1)
    Fk <- max(ObjF(lambda_a, R, n, k), ObjF(lambda_b, R, n, k))
    if (k == 1) {
      Fun <- Fk
    } else {
      Fun <- max(Fun, Fk)
    }
  }

  # Bound the delinkage
  lambda <- 2 / (1 + R)
  G <- ObjG(lambda, R, n)
  
  # Jointly bound the entropy
  for (k in 1:n) {
    p_a <- 1 / (k - 1 + R * (n - k + 1))
    p_b <- 1 / (k + R * (n - k))
    Ek <- max(ObjE(p_a, R, n, k), ObjE(p_b, R, n, k))
    if (k == 1) {
      Efun <- Ek
    } else {
      Efun <- max(Efun, Ek)
    }
  }

  return(
    list(shrinkage_ub = 0.5 * n * log(Fun / n),
         delinkage_lb = 0.5 * G,
         joint_entropy_ub = n * log(1 / n) + Efun,
         entropy_ub = 0.5 * n * log(Fun / n) + 0.5 * G)
  )
}

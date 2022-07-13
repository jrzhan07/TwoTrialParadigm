## ------------------------------------------------------------------------- ## 
##
## Script name: Simulations.R
## Purpose of script: Reproduce simulation results from Section 3.2.3 and 
## Figures 6 and 7 in the manuscript
##
##    Zhan SJ, Kunz CU, Stallard N. Should the two-trial paradigm still be the 
##    gold standard in drug assessment? (submitted)
##
## Author: Stella Jinran Zhan
##
## ------------------------------------------------------------------------- ## 


# (0) Initialisation ---------------------------------------------------------- 
# |- Load R packages ----
library(parallel)
library(metafor)

# |- Functions ----
# f.datagen: to obtain the effect estimate, between-trial heterogeneity estimate,
#            and test statistic for a single iteration 
# Parameters: 
#   Delta_RE = true common effect 
#   std = standard deviation
#   n.arm = number of patients per arm 
#   f = fraction of patients in the first trial
#   tau.2 = between-trial heterogeneity 
#   md = estimator for tau.2
f.datagen <- function(Delta_RE = 0, std = 1, n.arm = 650/2, f = 0.5, 
                      alpha = 0.05, tau.2, md = "DL"){
  simul <- matrix(NA, nrow = 1, ncol = 4)
  colnames(simul) <- c("mu", "tau", "Z", "flag")
  
  n1 <- ceiling(n.arm*f)
  n2 <- floor(n.arm*(1-f))
  nmat <- c(n1, n2)
  sk.2 <- 2*(std^2) / nmat
  sk <- sqrt(sk.2)
  Delta1 <- rnorm(n = 1, mean = Delta_RE, sd = sqrt(tau.2 + sk.2[1]))  
  Delta2 <- rnorm(n = 1, mean = Delta_RE, sd = sqrt(tau.2 + sk.2[2]))  
  Deltak <- c(Delta1, Delta2)
  
  rma1 <- try(rma.uni(yi = Deltak, sei = sk, method = md))
  simul[1, "mu"]      <- rma1$b
  simul[1, "tau"]     <- sqrt(rma1$tau2)
  simul[1, "Z"]       <- rma1$zval
  simul[1, "flag"]    <- ifelse(rma1$zval > qnorm(1-(alpha/2)^2), 1, 0)
  return(simul)
}

# f.simul: to get the estimated power based on n.sim simulations from f.datagen
# Parameters: 
#   tau.2_vec = vector of different tau^2 values for simulations
#   n.sim = number of simulations for each tau.2/f combination
#   Delta_RE = true common effect 
#   f = fraction of patients in the first trial
f.simul <- function(tau.2_vec, n.sim, Delta_RE, f){
  sim.row <- unlist(mclapply(tau.2_vec, function(tau.2_value){
    # out: unlisted list of n.sim results from f.datagen
    out <- unlist(mclapply(1:n.sim, function(row){
      f.datagen(Delta_RE = Delta_RE, f = f, tau.2 = tau.2_value)
    }))
    out.mat <- matrix(out, ncol = 4, byrow = T) 
    flag.mean <- mean(out.mat[, 4])             
  }))
}


# (1) Simulation --------------------------------------------------------------
n.sim <- 1000000 # number of simulations
tau.2_vec <- seq(0, 1, by = 0.1)^2 # different tau^2 values for simulations
Delta_RE_vec <- c(0, 0.2, 0.4, 0.8)

simul.out <- matrix(NA, nrow = 1, ncol = 2+length(tau.2_vec))
colnames(simul.out) <- c("Delta_RE", "f", seq(0, 1, by = 0.1)^2)
simul.out[1, 1:2] <- c(Delta_RE_vec[1], 0.1)
simul.out[1, 3:ncol(simul.out)] <- f.simul(tau.2_vec, n.sim = n.sim, 
                                           Delta_RE = Delta_RE_vec[1], f = 0.1)
simul.out <- rbind(simul.out, c(Delta_RE_vec[1], 0.5, 
                                f.simul(tau.2_vec, n.sim = n.sim, 
                                        Delta_RE = Delta_RE_vec[1], f = 0.5)))

for (D in 2:length(Delta_RE_vec)){
  simul.out <- rbind(simul.out, c(Delta_RE_vec[D], 0.1, 
                     f.simul(tau.2_vec, n.sim = n.sim, 
                             Delta_RE = Delta_RE_vec[D], f = 0.1)))
  
  simul.out <- rbind(simul.out, c(Delta_RE_vec[D], 0.5, 
                       f.simul(tau.2_vec, n.sim = n.sim,
                               Delta_RE = Delta_RE_vec[D], f = 0.5)))
}


# (2) Figures ------------------------------------------------------------------
collty = c("black", "gray60", "solid", "solid")

# |- Figure 6 ----
plot(tau.2_vec, simul.out[2, 3:ncol(simul.out)], 
     type = "l", xlab = expression(tau^2), ylab = "Type I error",
     cex.axis = 1.5, cex.lab = 1.5, cex = 0.2, col = collty[1],
     ylim = c(min(simul.out[1:2, 3:ncol(simul.out)]),
              max(simul.out[1:2, 3:ncol(simul.out)])))
points(tau.2_vec, simul.out[2, 3:ncol(simul.out)], col = collty[1], pch = 19)

lines(tau.2_vec, simul.out[1, 3:ncol(simul.out)], col = collty[1], lty = 2)
points(tau.2_vec, simul.out[1, 3:ncol(simul.out)], pch = 18, col = collty[1])

text(0.4, 0.07, expression(Delta == 0), col = 1)
legend('bottomright', title = expression(italic(f)), c("0.1", "0.5"), 
       lty = c(2, 1), col = collty[1], 
       pch = c(18, 19), horiz = F, cex = 1.5)

# |- Figure 7 ----
# Delta_RE = 0.2 (f = 0.1, 0.5)
plot(tau.2_vec, simul.out[4, 3:ncol(simul.out)], 
     type = "l", xlab = expression(tau^2), ylab = "Power", 
     ylim = c(0,1), cex.axis = 1.5, cex.lab = 1.5, cex = 0.2, col = collty[1])
points(tau.2_vec, simul.out[4, 3:ncol(simul.out)], pch = 19, col = collty[1])
lines(tau.2_vec, simul.out[3, 3:ncol(simul.out)], col = collty[1], lty = 2)
points(tau.2_vec, simul.out[3, 3:ncol(simul.out)], pch = 18, col = collty[1])

# Delta_RE = 0.4 (f = 0.1, 0.5)
lines(tau.2_vec, simul.out[6, 3:ncol(simul.out)], col = collty[2], lty = 1)
points(tau.2_vec, simul.out[6, 3:ncol(simul.out)], pch = 19, col = collty[2])
lines(tau.2_vec, simul.out[5, 3:ncol(simul.out)], col = collty[2], lty = 2)
points(tau.2_vec, simul.out[5, 3:ncol(simul.out)], pch = 18, col = collty[2])

# Delta_RE = 0.8 (f = 0.1, 0.5)
lines(tau.2_vec, simul.out[8, 3:ncol(simul.out)], col = 1, lty = 1, lwd = 2)
points(tau.2_vec, simul.out[8, 3:ncol(simul.out)], pch = 19, col = 1)
lines(tau.2_vec, simul.out[7, 3:ncol(simul.out)], col = 1, lty = 2, lwd = 2)
points(tau.2_vec, simul.out[7, 3:ncol(simul.out)], pch = 18, col = 1)

text(0.2, 0.8, expression(Delta == 0.8), col = 1,  font = 2)
text(0.2, 0.4, expression(Delta == 0.4), col = collty[2])
text(0.2, 0.1, expression(Delta == 0.2), col = 1)
legend('topright', title = expression(italic(f)), c("0.1","0.5"), 
       lty = c(2, 1), col = 1, pch = c(18, 19), horiz = F, cex = 1.5)
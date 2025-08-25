### --- Libraries ---
install.packages(c("copula", "energy", "mvtnorm"), dependencies = TRUE)
library(copula)
library(energy)
library(mvtnorm)

### --- Frechet Generator Functions ---
rfrechet <- function(n, shape = 1) {
  u <- runif(n)
  return((1 / (-log(u)))^(1 / shape))
}

qfrechet <- function(p, shape = 1) {
  return((1 / (-log(p)))^(1 / shape))
}

### --- Marginal Spacing Calculator ---
get_spacings <- function(x) {
  sort(x)[-1] - sort(x)[-length(x)]
}

### --- Spacing-Based Kendall's Tau ---
tau_spacing <- function(x, y) {
  sx <- get_spacings(x)
  sy <- get_spacings(y)
  n <- length(sx)
  sum_val <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      sum_val <- sum_val + sign(sx[i] - sx[j]) * sign(sy[i] - sy[j])
    }
  }
  return(2 * sum_val / (n * (n - 1)))
}

### --- Spacing-Based Bergsma-Dassios Tau* ---
tau_star_spacing <- function(x, y) {
  sx <- get_spacings(x)
  sy <- get_spacings(y)
  n <- length(sx)
  if (n < 4) return(NA)
  combs <- combn(n, 4)
  tau_val <- 0
  for (k in 1:ncol(combs)) {
    idx <- combs[, k]
    h <- function(s) {
      a <- s[1]; b <- s[2]; c <- s[3]; d <- s[4]
      sign((a - b) * (c - d)) * sign((a - c) * (b - d))
    }
    tau_val <- tau_val + h(sx[idx]) * h(sy[idx])
  }
  return(tau_val / choose(n, 4))
}

### --- Spacing-Based Distance Covariance ---
dcov_spacing <- function(x, y) {
  sx <- get_spacings(x)
  sy <- get_spacings(y)
  A <- abs(outer(sx, sx, "-"))
  B <- abs(outer(sy, sy, "-"))
  A_centered <- A - rowMeans(A) - colMeans(A) + mean(A)
  B_centered <- B - rowMeans(B) - colMeans(B) + mean(B)
  return(mean(A_centered * B_centered))
}

### --- Power Simulation Function ---
sim_power <- function(n = 100, R = 500, alpha = 1.5, copula_type = "gumbel",

                      alt = c("noncontig", "contig"), 
                      stat = c("tau", "tau_star", "dcov"), delta = 1) 
                      {
  alt <- match.arg(alt)
  stat <- match.arg(stat)
  powers <- numeric(R)
  if (copula_type == "gumbel") {
    cop <- gumbelCopula(param = 2, dim = 2)
  }
  for (r in 1:R) {
    if (alt == "noncontig") {
      u <- rCopula(n, cop)
      X <- qfrechet(u[,1], shape = alpha)
      Y <- qfrechet(u[,2], shape = alpha)
    } else {
      X <- rfrechet(n, shape = alpha)
      Y <- rfrechet(n, shape = alpha)
      h_xy <- delta * sin(X * Y) / sqrt(n)
      Y <- Y + h_xy
    }
    powers[r] <- switch(stat,
                        "tau" = tau_spacing(X, Y),
                        "tau_star" = tau_star_spacing(X, Y),
                        "dcov" = dcov_spacing(X, Y))
  }
  null_vals <- replicate(R, {
    X0 <- rfrechet(n, shape = alpha)
    Y0 <- rfrechet(n, shape = alpha)
    switch(stat,
           "tau" = tau_spacing(X0, Y0),
           "tau_star" = tau_star_spacing(X0, Y0),
           "dcov" = dcov_spacing(X0, Y0))
  })
  crit_val <- quantile(null_vals, 0.95)
  power_est <- mean(powers > crit_val)
  return(power_est)
}

### --- Main Power Analysis Over Alpha ---
alphas <- c(0.5, 0.8, 1.0, 1.2, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0)
results <- data.frame(
  alpha = alphas,
  tau_power_noncontig = NA,
  tau_star_power_noncontig = NA,
  dcov_power_noncontig = NA,
  tau_power_contig = NA,
  tau_star_power_contig = NA,
  dcov_power_contig = NA
)

for (i in seq_along(alphas)) {
  a <- alphas[i]
  cat("Simulating for alpha =", a, "\n")
  results$tau_power_noncontig[i] 
  <- sim_power(alpha = a, alt = "noncontig", stat = "tau")
  results$tau_star_power_noncontig[i] 
  <- sim_power(alpha = a, alt = "noncontig", stat = "tau_star")
  results$dcov_power_noncontig[i] 
  <- sim_power(alpha = a, alt = "noncontig", stat = "dcov")
  results$tau_power_contig[i] 
  <- sim_power(alpha = a, alt = "contig", stat = "tau")
  results$tau_star_power_contig[i] 
  <- sim_power(alpha = a, alt = "contig", stat = "tau_star")
  results$dcov_power_contig[i] 
  <- sim_power(alpha = a, alt = "contig", stat = "dcov")
}

print(results)
library(ggplot2)
library(reshape2)

df_long <- melt(results, id.vars = "alpha",
                variable.name = "Measure",
                value.name = "Power")

df_long$Type 
<- ifelse(grepl("noncontig", df_long$Measure), "Non-Contiguous", "Contiguous")
df_long$Statistic 
<- gsub("_power_(noncontig|contig)", "", df_long$Measure)

ggplot(df_long, aes(x = alpha, y = Power, color = Statistic, linetype = Type)) +
  geom_line(size = 1) +
  labs(title = "Empirical Power vs Tail Parameter (alpha)",
       x = expression(alpha),
       y = "Empirical Power") +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Dark2") +
  ylim(0, 1)

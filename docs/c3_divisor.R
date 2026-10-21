# Why C3 is not a conversion (differences.md D3).
#
# Maximises the two-step polyserial likelihood of Olsson, Drasgow & Dorans
# (1982), eqs. 19-20 (p.342), in rho alone, with the thresholds at the
# inverse-normal marginal proportions of y (as polycor's first step sets
# them), by optimize() at tol 1e-14.  The likelihood is flat in rho at
# about 1e-8, so the fits are good to about 1e-8, not to tol.  It fits once
# with z = scale(x) (n - 1, as polycor does) and once with the 1/N standard
# deviation (as pg-phenotype and eq. 36 do).  The log p(x) term of eq. 20 is constant in
# rho and is left out.  Run after the harness, which writes
# .work/polycor.in.json:
#
#   pixi run Rscript docs/c3_divisor.R

library(jsonlite)
inp <- fromJSON(".work/polycor.in.json", simplifyDataFrame = FALSE)
for (case in inp$cases) {
  if (!(case$case %in% c("polyserial", "biserial"))) next
  x <- as.double(case$x); y <- as.double(case$y); n <- length(x)
  props <- cumsum(table(y)) / n
  thr <- c(-Inf, qnorm(props[-length(props)]), Inf)
  cat_i <- as.integer(as.factor(y))
  lower <- thr[cat_i]; upper <- thr[cat_i + 1]
  # -log L(rho) from eq. 19: P(Y = j | x) = Phi(tau*_j) - Phi(tau*_(j-1)),
  # tau*_j = (tau_j - rho z) / sqrt(1 - rho^2).
  nll <- function(rho, z) {
    sd_c <- sqrt(1 - rho^2)
    -sum(log(pnorm((upper - rho * z) / sd_c) - pnorm((lower - rho * z) / sd_c)))
  }
  z_r <- as.vector(scale(x))
  z_n <- (x - mean(x)) / sqrt(mean((x - mean(x))^2))
  tight <- function(z) optimize(nll, c(-0.9999, 0.9999), z = z, tol = 1e-14)$minimum
  r_r <- tight(z_r); r_n <- tight(z_n)
  poly <- polycor::polyserial(x, y, ML = FALSE, std.err = TRUE)$rho
  cat(sprintf("%s n=%d  polycor %.12f  exact(n-1) %.12f  exact(1/N) %.12f\n", case$case, n, poly, r_r, r_n))
  cat(sprintf("  polycor optimiser error %.3g; divisor effect %.3g; ratio exact(1/N)/exact(n-1) %.12f; sqrt((n-1)/n) %.12f\n",
              poly - r_r, r_n - r_r, r_n / r_r, sqrt((n - 1) / n)))
}

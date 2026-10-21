# Why C3 is not a conversion (differences.md D3).
#
# Refits the two-step polyserial likelihood of polycor::polyserial (copied
# from its own f(), with its own first-step thresholds) by optimize() at
# tol 1e-14, once with z = scale(x) (n - 1, as polycor does) and once with
# the 1/N standard deviation (as pg-phenotype does).  Run after the harness,
# which writes .work/polycor.in.json:
#
#   pixi run Rscript docs/c3_divisor.R

library(jsonlite)
inp <- fromJSON(".work/polycor.in.json", simplifyDataFrame = FALSE)
for (case in inp$cases) {
  if (!(case$case %in% c("polyserial", "biserial"))) next
  x <- as.double(case$x); y <- as.double(case$y); n <- length(x)
  tab <- table(y); s <- length(tab); cuts <- qnorm(cumsum(tab) / n)[-s]
  yi <- as.numeric(as.factor(y))
  nll <- function(rho, z) {
    cts <- c(-Inf, cuts, Inf)
    tau <- (matrix(cts, n, s + 1, byrow = TRUE) - matrix(rho * z, n, s + 1)) / sqrt(1 - rho^2)
    -sum(log(dnorm(z) * (pnorm(tau[cbind(1:n, yi + 1)]) - pnorm(tau[cbind(1:n, yi)]))))
  }
  z_r <- as.vector(scale(x))
  z_n <- (x - mean(x)) / sqrt(mean((x - mean(x))^2))
  tight <- function(z) optimize(nll, c(-0.9999, 0.9999), z = z, tol = 1e-14)$minimum
  poly <- polycor::polyserial(x, y, ML = FALSE, std.err = TRUE)$rho
  cat(sprintf("%s n=%d  polycor %.12f  exact(n-1) %.12f  exact(1/N) %.12f\n", case$case, n, poly, tight(z_r), tight(z_n)))
  cat(sprintf("  polycor optimiser error %.3g; divisor effect %.3g; ratio exact(1/N)/exact(n-1) %.12f; sqrt((n-1)/n) %.12f\n",
              poly - tight(z_r), tight(z_n) - tight(z_r), tight(z_n) / tight(z_r), sqrt((n - 1) / n)))
}

# sandwich (Zeileis): the OLS line and its cluster-robust SE over Mate Networks.
#
# For each fit: slope and se_ols from lm, and se_cl from
# vcovCL(type = "HC0", cadjust = TRUE), whose only adjustment is G/(G-1).
args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "_io.R"))

parent_offspring <- function(case) {
  out <- list()
  for (key in names(case$fits)) {
    f <- case$fits[[key]]
    d <- data.frame(x = nums(f$x), y = nums(f$y), network = unlist(f$cluster))
    fit <- lm(y ~ x, data = d)
    v <- sandwich::vcovCL(fit, cluster = ~network, type = "HC0", cadjust = TRUE)
    out[[paste0(key, ".slope")]] <- unname(coef(fit)[2])
    out[[paste0(key, ".se_ols")]] <- summary(fit)$coefficients[2, 2]
    out[[paste0(key, ".se_cl")]] <- sqrt(v[2, 2])
    out[[paste0(key, ".n")]] <- nrow(d)
  }
  out
}

run_cases(list(parent_offspring = parent_offspring), "sandwich")

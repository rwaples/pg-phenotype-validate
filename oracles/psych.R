# psych (Revelle): polychoric with no continuity correction and per-pair
# thresholds, rho by optimize() at its default tolerance.
args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "_io.R"))

mate <- function(case) {
  d <- data.frame(x = nums(case$x), y = nums(case$y))
  fit <- suppressMessages(suppressWarnings(psych::polychoric(d, correct = 0, global = FALSE)))
  list(rho = fit$rho[1, 2])
}

run_cases(list(mate = mate), "psych")

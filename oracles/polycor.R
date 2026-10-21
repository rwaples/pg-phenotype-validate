# polycor (Fox): polychor on a pair table, polyserial on paired values.
#
#   rho_2step  polychor(tab, ML = FALSE): thresholds from the margins, rho by optimise()
#   rho_ml     polychor(tab, ML = TRUE, std.err = TRUE): joint ML of rho and thresholds
#   se_ml      sqrt(var[1, 1]) of that joint fit (numerical Hessian)
#   se_fixed   polychor(tab, ML = FALSE, std.err = TRUE): rho's SE with thresholds held fixed
#   serial_rho, serial_se
#              polyserial(x, y, ML = FALSE, std.err = TRUE): the two-step fit
#              (without std.err polyserial returns the ad hoc estimate)
args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "_io.R"))

quietly <- function(expr) withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning"))

table_fits <- function(tab) {
  two <- quietly(polycor::polychor(tab, ML = FALSE))
  ml <- quietly(polycor::polychor(tab, ML = TRUE, std.err = TRUE))
  fixed <- quietly(polycor::polychor(tab, ML = FALSE, std.err = TRUE))
  list(rho_2step = two, rho_ml = ml$rho, se_ml = sqrt(ml$var[1, 1]), se_fixed = sqrt(fixed$var[1, 1]))
}

tetrachoric <- function(case) table_fits(matrix(unlist(case$table), 2, byrow = TRUE))

mate <- function(case) {
  x <- nums(case$x)
  y <- nums(case$y)
  if (case$kinds_xy[[1]] == "continuous") {
    fit <- quietly(polycor::polyserial(x, y, ML = FALSE, std.err = TRUE))
    return(list(serial_rho = fit$rho, serial_se = sqrt(fit$var[1, 1])))
  }
  table_fits(table(x, y))
}

run_cases(list(tetrachoric = tetrachoric, mate = mate), "polycor")

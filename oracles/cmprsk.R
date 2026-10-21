# cmprsk (Gray): cuminc's cumulative incidence of cause 1, read at each case
# age by timepoints() (right-continuous).  cuminc has no late-entry argument,
# so only cases without entry come here.
args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "_io.R"))

aalen_johansen <- function(case) {
  if (!is.null(case$entry)) stop("cuminc takes no late entry")
  age <- nums(case$age)
  status <- as.integer(nums(case$status))
  fit <- cuminc(age, status, cencode = 0)
  grid <- sort(unique(age[status == 1]))
  list(ages = grid, cip = unname(timepoints(fit, grid)$est["1 1", ]))
}

run_cases(list(aalen_johansen = aalen_johansen), "cmprsk")

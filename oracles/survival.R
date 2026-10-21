# survival (Therneau): survfit's Aalen-Johansen CIF of cause 1 at each case
# age, called as pg-phenotype's own oracle (tools/survival_aj.R) calls it.
args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "_io.R"))

aalen_johansen <- function(case) {
  age <- nums(case$age)
  code <- as.integer(nums(case$status))
  status <- factor(code, levels = 0:2, labels = c("censored", "case", "death"))
  if (is.null(case$entry)) {
    fit <- survfit(Surv(age, status) ~ 1, timefix = FALSE)
  } else {
    entry <- nums(case$entry)
    keep <- age > entry  # a row censored at entry is at risk nowhere; Surv needs entry < age
    fit <- survfit(Surv(entry[keep], age[keep], status[keep]) ~ 1, id = seq_len(sum(keep)), timefix = FALSE)
  }
  grid <- sort(unique(age[code == 1]))
  list(ages = grid, cip = fit$pstate[match(grid, fit$time), match("case", fit$states)])
}

run_cases(list(aalen_johansen = aalen_johansen), "survival")

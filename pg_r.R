# The pgphenotype R binding on each case, as run_oracles.py runs the Python
# binding (same families, same result names).
#
#   Rscript pg_r.R <input.json> <output.json>
args0 <- commandArgs(trailingOnly = FALSE)
here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE))))
source(file.path(here, "oracles", "_io.R"))

ped_of <- function(case) {
  data.frame(id = nums(case$id), mother = nums(case$mother), father = nums(case$father), sex = nums(case$sex))
}

tetrachoric <- function(case) {
  r <- correlation_tetrachoric(table = matrix(unlist(case$table), 2, byrow = TRUE))
  list(rho = or_null(r$value), se = or_null(r$se), boundary = or_null(r$boundary), n = r$n)
}

mate <- function(case) {
  traits <- lapply(seq_along(case$kinds), function(i) trait(nums(case$traits[[i]]), kind = case$kinds[[i]]))
  r <- assortative_mate_correlation(ped_of(case), traits, permutations = 0)
  cell <- Filter(function(c) c$mother_trait == case$cell[[1]] && c$father_trait == case$cell[[2]], r$cells)[[1]]
  p <- cell$crude[[1]]
  list(estimator = p$estimator, rho = or_null(p$value), se = or_null(p$se), boundary = or_null(p$boundary),
       n = cell$n, n_mate_networks = cell$n_mate_networks)
}

parent_offspring <- function(case) {
  stratum <- if (isTRUE(case$stratified)) nums(case$stratum) else NULL
  r <- regression_parent_offspring(ped_of(case), trait(nums(case$y), kind = "continuous"), stratum = stratum)
  out <- list()
  put <- function(subset, fits) {
    for (line in c("midparent", "mother", "father")) {
      f <- fits[[line]]
      for (field in c("slope", "se", "se_ols", "n", "n_mate_networks")) {
        out[[paste(subset, line, field, sep = ".")]] <<- or_null(f[[field]])
      }
    }
  }
  put("pooled", r$pooled)
  for (s in r$strata) put(paste0("s", s$stratum), s)
  out
}

krebs <- function(case) {
  ped <- data.frame(id = c(1, 2), mother = c(NA, NA), father = c(NA, NA), twin = c(2, 1))
  prep <- pafgrs_prepare(ped, ndegree = 1)
  lapply(case$scenarios, function(s) {
    cip <- pafgrs_cip(c(0, s$cip_age), c(0, s$K))
    pafgrs_score_univariate(prep, trait(c(0, s$status), kind = "binary"), age = c(NA, s$age), cip = cip,
                            h2 = 1)$est[1]
  })
}

aalen_johansen <- function(case) {
  entry <- if (is.null(case$entry)) NULL else nums(case$entry)
  cip <- pafgrs_cip_aalen_johansen(nums(case$age), nums(case$status), entry = entry)
  list(ages = cip$ages, cip = cip$cip)
}

run_cases(
  list(tetrachoric = tetrachoric, mate = mate, parent_offspring = parent_offspring, krebs = krebs,
       aalen_johansen = aalen_johansen),
  "pgphenotype"
)

# The one interface every R script here follows:
#
#   Rscript <script> <input.json> <output.json>
#
# The input holds `cases`, each with a `case` id, a `family` and its data.  The
# script writes `values` (case id -> named results) and `versions` (R and every
# package it loaded).  Doubles are written with 17 significant digits so they
# round-trip; jsonlite's `digits = NA` keeps only 15.
suppressPackageStartupMessages(library(jsonlite))

run_cases <- function(handlers, packages) {
  args <- commandArgs(trailingOnly = TRUE)
  input <- fromJSON(args[1], simplifyVector = TRUE, simplifyDataFrame = FALSE, simplifyMatrix = FALSE)
  for (p in packages) suppressPackageStartupMessages(library(p, character.only = TRUE))
  values <- list()
  for (case in input$cases) {
    handler <- handlers[[case$family]]
    if (is.null(handler)) stop(sprintf("case %s: no handler for family %s", case$case, case$family))
    values[[case$case]] <- tryCatch(handler(case), error = function(e) list(error = conditionMessage(e)))
  }
  versions <- c(list(R = R.version.string), lapply(setNames(packages, packages), function(p) {
    as.character(packageVersion(p))
  }))
  out <- list(values = values, versions = versions)
  write_json(out, args[2], auto_unbox = TRUE, digits = I(17), null = "null", na = "null")
}

# NULL (absent in the JSON) or NA becomes NA, so a missing value stays one.
num <- function(x) if (is.null(x)) NA_real_ else as.double(x)
nums <- function(x) vapply(x, num, numeric(1))
or_null <- function(x) if (is.null(x)) NA else x

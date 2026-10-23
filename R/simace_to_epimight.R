#!/usr/bin/env Rscript
# simace_to_epimight.R
#
# Converts one simACE replicate (cohort.parquet + pedigree.parquet) into an
# EPIMIGHT pool: one row per proband x trait x kind of relative, written as a
# gzipped TSV.
#
# Usage:
#   Rscript simace_to_epimight.R --cohort=rep1/cohort.parquet --pedigree=rep1/pedigree.parquet
#   Rscript simace_to_epimight.R --cohort=... --pedigree=... --out=OUT.tsv.gz \
#       --rels=FS,HS,PO --first-birth-year=1980 --generation-years=25
#   (--key value works as well as --key=value)
#
# Options (all but --cohort and --pedigree have defaults):
#   --out                 output file                              [OUT.tsv.gz]
#   --rels                kinds of relatives, comma separated      [PO,FS,HS,mHS,pHS,Av,1C]
#                           PO  parents (mother, father)
#                           FS  full siblings (MZ co-twins left out)
#                           HS  half siblings, maternal + paternal
#                           mHS maternal half siblings (same mother, other father)
#                           pHS paternal half siblings (same father, other mother)
#                           Av  aunts and uncles (full siblings or MZ co-twin of a parent)
#                           1G  grandparents (not in the default: with 2 phenotyped
#                               generations no grandparent is phenotyped)
#                           1C  first cousins (children of an aunt or uncle, or
#                               a pair sharing two grandparents through half siblings)
#   --traits              simACE traits to export                  [1,2]
#   --trait-names         names written in the trait column        [trait1,trait2]
#   --first-birth-year    birth year of the first phenotyped
#                         generation                               [1960]
#   --generation-years    years between generations                [30]
#   --death-status        trait_status for a death before onset:
#                         2 = competing risk, 0 = censored         [2]
#
# Output columns:
#   person_id          simACE id, as text
#   trait              trait name
#   trait_status       0 = censored at end of follow-up, 1 = affected,
#                      2 = died before onset (see --death-status)
#   trait_age          completed years of age at onset, death or censoring
#   relatives_kind     code from --rels
#   relatives_n        number of relatives of that kind with a phenotype
#   relatives_n_trait  how many of them are affected
#   birth_year         first-birth-year + generation-years x (generation - first phenotyped generation)
#   sex                0 = female, 1 = male
#
# Probands are the phenotyped people (rows of cohort.parquet with affected1
# not missing). Relatives are counted only if they are phenotyped too, so in
# the first phenotyped generation parents and grandparents give relatives_n = 0.
# PO, Av and 1G count only the older relative (the proband is the child,
# niece/nephew or grandchild). Everyone is kept, including people with
# relatives_n = 0, because EPIMIGHT's population CIF uses every row.
# Pairs are classified as pedigree-graph classifies them, each under its
# closest relationship (a parent is never also an aunt or uncle);
# check/relative_counts.py compares the counts with pedigree-graph.
#
# Needs the R packages arrow and data.table.

suppressPackageStartupMessages({
  library(arrow)
  library(data.table)
})

# ---- options -----------------------------------------------------------------
defaults <- list(
  cohort = NULL,
  pedigree = NULL,
  out = "OUT.tsv.gz",
  rels = "PO,FS,HS,mHS,pHS,Av,1C",
  traits = "1,2",
  `trait-names` = "trait1,trait2",
  `first-birth-year` = "1960",
  `generation-years` = "30",
  `death-status` = "2"
)
relatedness <- c(PO = 0.5, FS = 0.5, HS = 0.25, mHS = 0.25, pHS = 0.25, Av = 0.25, `1G` = 0.25, `1C` = 0.125)

usage <- function() {
  file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  if (length(file) != 1 || !file.exists(file)) {
    cat("Usage: Rscript simace_to_epimight.R --cohort=PATH --pedigree=PATH [options]\n",
        "Options: ", paste0("--", setdiff(names(defaults), c("cohort", "pedigree")), collapse = " "), "\n",
        "The comment block at the top of the script describes each option.\n", sep = "")
    return(invisible())
  }
  lines <- readLines(file)
  cat(sub("^# ?", "", lines[2:(which(lines == "")[1] - 1)]), sep = "\n")
}

# Accepts --key=value and --key value
parse_args <- function(args, defaults) {
  opt <- defaults
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a %in% c("-h", "--help")) {
      usage()
      quit(status = 0)
    }
    m <- regmatches(a, regexec("^--([^=]+)(=(.*))?$", a))[[1]]
    if (length(m) == 0 || !(m[2] %in% names(defaults))) stop("unknown argument: ", a, " (see --help)")
    if (nzchar(m[3])) {
      opt[[m[2]]] <- m[4]
    } else {
      if (i == length(args) || startsWith(args[i + 1], "--")) stop("--", m[2], " needs a value")
      i <- i + 1
      opt[[m[2]]] <- args[i]
    }
    i <- i + 1
  }
  if (is.null(opt$cohort) || is.null(opt$pedigree)) stop("--cohort and --pedigree are required (see --help)")
  opt
}
split_opt <- function(x) trimws(strsplit(x, ",", fixed = TRUE)[[1]])

opt <- parse_args(commandArgs(trailingOnly = TRUE), defaults)
rels <- split_opt(opt$rels)
unknown <- setdiff(rels, names(relatedness))
if (length(unknown)) stop("unknown --rels: ", paste(unknown, collapse = ", "), "; choose from ", paste(names(relatedness), collapse = ", "))
traits <- as.integer(split_opt(opt$traits))
trait_names <- split_opt(opt$`trait-names`)
if (length(trait_names) != length(traits)) stop("--trait-names needs one name per --traits entry")
first_year <- as.integer(opt$`first-birth-year`)
gen_years <- as.numeric(opt$`generation-years`)
death_status <- as.integer(opt$`death-status`)
if (!death_status %in% c(0L, 2L)) stop("--death-status must be 0 or 2")

# ---- read --------------------------------------------------------------------
ped <- as.data.table(read_parquet(opt$pedigree, col_select = all_of(c("id", "generation", "sex", "mother", "father", "twin"))))
outcome_cols <- unlist(lapply(traits, function(t) paste0(c("affected", "t_observed", "death_censored"), t)))
coh <- as.data.table(read_parquet(opt$cohort, col_select = all_of(c("id", outcome_cols))))
coh <- coh[!is.na(get(paste0("affected", traits[1])))]
coh <- ped[coh, on = "id", nomatch = NULL]
cat("Probands (phenotyped people):", nrow(coh), "\n")
first_gen <- min(coh$generation)

# ---- relative pairs: proband id -> relative id, every pair once ----------------
kids <- ped[mother >= 0, .(id, mother, father, twin)]
# all children of the same mother or the same father, then classify
by_mother <- kids[kids, on = "mother", allow.cartesian = TRUE, nomatch = NULL][id != i.id, .(id, rel = i.id, same_father = father == i.father, twin, gp = mother)]
by_father <- kids[kids, on = "father", allow.cartesian = TRUE, nomatch = NULL][id != i.id & mother != i.mother, .(id, rel = i.id, gp = father)]
fs_all <- by_mother[same_father == TRUE & rel != twin, .(id, rel)]
parents <- rbind(kids[, .(id, rel = mother)], kids[, .(id, rel = father)])
# a parent's MZ co-twin counts as an aunt or uncle, as in pedigree-graph
parent_sibs <- rbind(fs_all, ped[twin >= 0, .(id, rel = twin)])
# half siblings, with the parent they share
half_sibs <- rbind(by_mother[same_father == FALSE, .(id, rel, gp)], by_father)

pairs_for <- function(kind) {
  switch(kind,
    PO  = parents,
    FS  = fs_all,
    mHS = by_mother[same_father == FALSE, .(id, rel)],
    pHS = by_father[, .(id, rel)],
    HS  = half_sibs[, .(id, rel)],
    `1G` = unique(parents[parents, on = c(id = "rel"), nomatch = NULL, allow.cartesian = TRUE][, .(id = i.id, rel)]),
    # a parent is never also counted as an aunt or uncle (parents who are full siblings);
    # pedigree-graph reports each pair under its closest relationship
    Av  = unique(parent_sibs[parents, on = c(id = "rel"), nomatch = NULL, allow.cartesian = TRUE][, .(id = i.id, rel)])[!parents, on = c("id", "rel")],
    `1C` = {
      pp <- parents[id %in% coh$id]
      av <- unique(parent_sibs[pp, on = c(id = "rel"), nomatch = NULL, allow.cartesian = TRUE][, .(id = i.id, aunt = rel)])
      cousins <- unique(parents[av, on = c(rel = "aunt"), nomatch = NULL, allow.cartesian = TRUE][, .(id = i.id, rel = id)])
      # pairs sharing two grandparents through half-sibling parents also count as
      # first cousins, as in pedigree-graph (e.g. double half cousins)
      half_av <- half_sibs[pp, on = c(id = "rel"), nomatch = NULL, allow.cartesian = TRUE][, .(id = i.id, aunt = rel, gp)]
      half_c <- parents[half_av, on = c(rel = "aunt"), nomatch = NULL, allow.cartesian = TRUE][, .(id = i.id, rel = id, gp)]
      double <- unique(half_c)[, .N, by = .(id, rel)][N >= 2, .(id, rel)]
      cousins <- unique(rbind(cousins, double))
      sibs <- rbind(by_mother[, .(id, rel)], by_father[, .(id, rel)], kids[twin >= 0, .(id, rel = twin)])
      cousins[id != rel][!sibs, on = c("id", "rel")]
    }
  )
}

pairs <- lapply(setNames(rels, rels), function(kind) {
  pr <- pairs_for(kind)
  pr[id %in% coh$id]
})

# ---- one block per trait x kind ------------------------------------------------
blocks <- list()
for (k in seq_along(traits)) {
  t <- traits[k]
  aff <- coh[[paste0("affected", t)]]
  status <- fifelse(aff, 1L, fifelse(coh[[paste0("death_censored", t)]], death_status, 0L))
  probands <- data.table(
    person_id = as.character(coh$id),
    trait = trait_names[k],
    trait_status = status,
    trait_age = as.integer(floor(coh[[paste0("t_observed", t)]])),
    birth_year = as.integer(first_year + gen_years * (coh$generation - first_gen)),
    sex = as.integer(coh$sex)
  )
  rel_status <- data.table(rel = coh$id, rel_aff = aff)
  for (kind in rels) {
    pr <- pairs[[kind]]
    pr <- rel_status[pr, on = "rel", nomatch = NULL]  # relatives must be phenotyped
    counts <- pr[, .(relatives_n = .N, relatives_n_trait = sum(rel_aff)), by = id]
    b <- copy(probands)
    b[, relatives_kind := kind]
    b[, c("relatives_n", "relatives_n_trait") := counts[.(coh$id), on = "id", .(relatives_n, relatives_n_trait)]]
    b[is.na(relatives_n), `:=`(relatives_n = 0L, relatives_n_trait = 0L)]
    b[, `:=`(relatives_n = as.integer(relatives_n), relatives_n_trait = as.integer(relatives_n_trait))]
    cat(sprintf("  %-7s %-4s %8d probands with relatives, %8d with an affected relative\n",
                trait_names[k], kind, sum(b$relatives_n > 0), sum(b$relatives_n_trait > 0)))
    blocks[[length(blocks) + 1]] <- b
  }
}

pool <- rbindlist(blocks)
setcolorder(pool, c("person_id", "trait", "trait_status", "trait_age", "relatives_kind",
                    "relatives_n", "relatives_n_trait", "birth_year", "sex"))
dir.create(dirname(opt$out), showWarnings = FALSE, recursive = TRUE)
fwrite(pool, opt$out, sep = "\t", compress = "gzip")
cat("Wrote", opt$out, ":", nrow(pool), "rows\n")
cat("relatedness to pass to EPIMIGHT run_h2:",
    paste(sprintf("%s=%s", rels, relatedness[rels]), collapse = ", "), "\n")

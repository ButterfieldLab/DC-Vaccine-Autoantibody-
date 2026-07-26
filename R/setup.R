# ============================================================================
# setup.R
# ----------------------------------------------------------------------------
# Single entry point sourced at the top of every analysis script. It replaces
# the original machine-specific paths and the two lost source() files.
#
#   source("R/setup.R")   # from the project root
#
# It defines `project_root`, `data_path` and `results_path`, then loads the
# package stack (R/libraries.R) and the figure colours/themes (R/color_tools.R).
#
# Raw data are NOT shipped with this repository (see data/README.md). Point
# `data_path` at wherever you placed the input files, or set the environment
# variable AUTOANTIBODY_DATA before sourcing this file.
# ============================================================================

# ---- Locate the project root -----------------------------------------------
.find_root <- function() {
  # 1. Honour an explicit override.
  env <- Sys.getenv("AUTOANTIBODY_ROOT", unset = NA)
  if (!is.na(env) && dir.exists(env)) return(normalizePath(env))
  # 2. Walk up from the working directory looking for the R/ + README marker.
  d <- normalizePath(getwd())
  for (i in 1:6) {
    if (file.exists(file.path(d, "R", "setup.R"))) return(d)
    parent <- dirname(d)
    if (parent == d) break
    d <- parent
  }
  # 3. Fall back to the working directory.
  normalizePath(getwd())
}

project_root <- .find_root()

# ---- Paths -----------------------------------------------------------------
data_path <- Sys.getenv("AUTOANTIBODY_DATA",
                        unset = file.path(project_root, "data"))
if (!endsWith(data_path, .Platform$file.sep)) {
  data_path <- paste0(data_path, .Platform$file.sep)
}
# `file_path` is kept for backward compatibility with the original scripts.
file_path <- paste0(project_root, .Platform$file.sep)
results_path <- file.path(project_root, "results")
fig_path   <- file.path(results_path, "figures")
table_path <- file.path(results_path, "tables")

# ---- Load dependencies + figure aesthetics ---------------------------------
source(file.path(project_root, "R", "libraries.R"))
source(file.path(project_root, "R", "color_tools.R"))

message("Autoantibody pipeline setup complete.\n  project_root: ", project_root,
        "\n  data_path:    ", data_path)

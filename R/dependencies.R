# Thin wrapper around requirements.txt for local use
pkgs <- readLines("requirements.txt")
pkgs <- sub("^r::", "", pkgs[startsWith(pkgs, "r::")])
missing <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
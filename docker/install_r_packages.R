options(repos = c(CRAN = "https://cloud.r-project.org"))

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
  }
}

is_available <- function(pkgs) {
  any(vapply(pkgs, requireNamespace, logical(1), quietly = TRUE))
}

has_blink_api <- function() {
  for (ns in c("BLINK", "blink")) {
    if (requireNamespace(ns, quietly = TRUE)) {
      ex <- tryCatch(getNamespaceExports(ns), error = function(e) character(0))
      if (any(c("BLINK", "Blink", "blink") %in% ex)) {
        return(TRUE)
      }
    }
  }
  FALSE
}

install_github_if_missing <- function(repo, aliases) {
  if (is_available(aliases)) {
    return(invisible(TRUE))
  }
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", dependencies = TRUE)
  }
  try(remotes::install_github(repo, upgrade = "never", dependencies = TRUE), silent = TRUE)
  is_available(aliases)
}

has_farmcpu_api <- function() {
  if (!requireNamespace("FarmCPU", quietly = TRUE)) {
    return(FALSE)
  }
  ns <- asNamespace("FarmCPU")
  exists("FarmCPU.Prior", envir = ns, mode = "function") &&
    exists("FarmCPU.LM", envir = ns, mode = "function")
}

must_have <- c("data.table", "optparse", "bigmemory", "mrMLM")
for (pkg in must_have) {
  install_if_missing(pkg)
}

# BLINK requires callable BLINK/Blink API, not just a namespace named blink.
if (!has_blink_api()) {
  try(install.packages("BLINK", dependencies = TRUE), silent = TRUE)
}
if (!has_blink_api()) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", dependencies = TRUE)
  }
  try(remotes::install_github("YaoZhou89/BLINK", upgrade = "never", dependencies = TRUE), silent = TRUE)
}

# FarmCPU implementation can come from FarmCPUpp or FarmCPU. Try both.
if (!requireNamespace("FarmCPUpp", quietly = TRUE)) {
  try(install.packages("FarmCPUpp", dependencies = TRUE), silent = TRUE)
}
if (!requireNamespace("FarmCPU", quietly = TRUE)) {
  try(install.packages("FarmCPU", dependencies = TRUE), silent = TRUE)
}

# CRAN may not provide FarmCPU/FarmCPUpp for all R builds; use GitHub fallbacks.
if (!requireNamespace("FarmCPUpp", quietly = TRUE) && !requireNamespace("FarmCPU", quietly = TRUE)) {
  install_github_if_missing("amkusmec/FarmCPUpp", c("FarmCPUpp"))
}

# BLINK native implementation relies on FarmCPU.Prior/FarmCPU.LM from FarmCPU.
if (!has_farmcpu_api()) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", dependencies = TRUE)
  }
  # Try common FarmCPU sources in order.
  try(remotes::install_github("jiabowang/FarmCPU", upgrade = "never", dependencies = TRUE), silent = TRUE)
}
if (!has_farmcpu_api()) {
  try(remotes::install_github("YaoZhou89/FarmCPU", upgrade = "never", dependencies = TRUE), silent = TRUE)
}
if (!has_farmcpu_api()) {
  try(remotes::install_github("amkusmec/FarmCPU", upgrade = "never", dependencies = TRUE), silent = TRUE)
}

required_checks <- c("data.table", "optparse", "bigmemory", "mrMLM")
missing_required <- required_checks[!vapply(required_checks, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_required) > 0) {
  stop(sprintf("Missing required R packages after installation: %s", paste(missing_required, collapse = ", ")))
}

if (!has_blink_api()) {
  ex_BLINK <- if (requireNamespace("BLINK", quietly = TRUE)) paste(getNamespaceExports("BLINK"), collapse = ", ") else "not installed"
  ex_blink <- if (requireNamespace("blink", quietly = TRUE)) paste(getNamespaceExports("blink"), collapse = ", ") else "not installed"
  stop(sprintf("BLINK API unavailable after CRAN/GitHub attempts. BLINK exports: [%s]; blink exports: [%s]", ex_BLINK, ex_blink))
}

if (!requireNamespace("FarmCPUpp", quietly = TRUE) && !requireNamespace("FarmCPU", quietly = TRUE)) {
  stop("Neither FarmCPUpp nor FarmCPU is available after CRAN/GitHub installation attempts; full-model workflow cannot proceed.")
}

if (!has_farmcpu_api()) {
  message("WARNING: FarmCPU.Prior/FarmCPU.LM not available; BLINK will use runtime compatibility shim.")
}

cat("R package installation complete.\n")

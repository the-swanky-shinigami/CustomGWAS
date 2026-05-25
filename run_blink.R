#!/usr/bin/env Rscript
# Lightweight runner for BLINK/FarmCPU/mrMLM with a low-memory fallback.
# Usage: Rscript run_blink.R --geno_csv path --pheno path --covar path --out prefix --trait TRAIT --chunk_size 2000 --method blink

args <- commandArgs(trailingOnly=TRUE)
parse_args <- function(args){
  out <- list()
  i <- 1
  while(i <= length(args)){
    a <- args[i]
    if(grepl('^--', a)){
      key <- sub('^--', '', a)
      if(i < length(args) && !grepl('^--', args[i+1])){
        out[[key]] <- args[i+1]
        i <- i + 2
      } else {
        out[[key]] <- TRUE
        i <- i + 1
      }
    } else { i <- i + 1 }
  }
  return(out)
}

opts <- parse_args(args)
geno_csv <- opts$geno_csv
pheno <- opts$pheno
covar <- opts$covar
out_pref <- ifelse(!is.null(opts$out), opts$out, 'blink_out')
trait <- ifelse(!is.null(opts$trait), opts$trait, '')
chunk_size <- as.integer(ifelse(!is.null(opts$chunk_size), opts$chunk_size, 2000))
method <- tolower(ifelse(!is.null(opts$method), opts$method, 'blink'))
strict_native <- !is.null(opts$strict_native) && as.character(opts$strict_native) %in% c('1','TRUE','True','true')

detect_blink_namespace <- function(){
  for (ns in c('BLINK', 'blink')) {
    if (requireNamespace(ns, quietly=TRUE)) {
      ex <- tryCatch(getNamespaceExports(ns), error=function(e) character(0))
      if (any(c('BLINK', 'Blink', 'blink') %in% ex)) {
        return(ns)
      }
    }
  }
  return(NULL)
}

pick_blink_fun <- function(ns){
  ex <- tryCatch(getNamespaceExports(ns), error=function(e) character(0))
  for (fn in c('BLINK', 'Blink', 'blink')) {
    if (fn %in% ex && exists(fn, where=asNamespace(ns), mode='function')) {
      return(fn)
    }
  }
  return(NULL)
}

make_blink_with_farmcpu_shim <- function(blink_fun, blink_ns) {
  shim_env <- new.env(parent = asNamespace(blink_ns))

  # BLINK v0.1.0 can call FarmCPU.Prior/FarmCPU.LM. Provide compatible fallbacks
  # when FarmCPU is unavailable in this runtime.
  shim_env$FarmCPU.Prior <- function(GM, P, Prior = NULL, kinship.algorithm = 'FaST-LMM') {
    n <- nrow(GM)
    if (!is.null(Prior)) {
      prior_num <- suppressWarnings(as.numeric(Prior))
      if (length(prior_num) == n) {
        prior_num[!is.finite(prior_num)] <- 1
        return(prior_num)
      }
    }
    p_num <- suppressWarnings(as.numeric(P))
    if (length(p_num) == 1 && is.finite(p_num)) {
      return(rep(p_num, n))
    }
    if (length(p_num) >= n) {
      p_num <- p_num[seq_len(n)]
      p_num[!is.finite(p_num)] <- 1
      return(p_num)
    }
    rep(1, n)
  }

  shim_env$FarmCPU.LM <- function(y, GDP, w = NULL, orientation = 'col') {
    if (exists('Blink.LM', where = asNamespace(blink_ns), mode = 'function')) {
      y_num <- suppressWarnings(as.numeric(y))
      GDP_mat <- as.matrix(GDP)
      suppressWarnings(storage.mode(GDP_mat) <- 'numeric')
      w_mat <- NULL
      if (!is.null(w)) {
        w_mat <- as.matrix(w)
        suppressWarnings(storage.mode(w_mat) <- 'numeric')
      }
      return(get('Blink.LM', envir = asNamespace(blink_ns))(y = y_num, GDP = GDP_mat, w = w_mat, orientation = orientation))
    }
    stop('Blink.LM is unavailable for FarmCPU.LM compatibility shim')
  }

  # BLINK may call GAPIT helpers even when file output is disabled.
  shim_env$GAPIT.Power <- function(...) {
    NULL
  }
  shim_env$GAPIT.Report <- function(...) {
    NULL
  }

  environment(blink_fun) <- shim_env
  blink_fun
}

if(is.null(geno_csv) || is.null(pheno)){
  stop('Required arguments missing: --geno_csv and --pheno')
}

cat(sprintf('run_blink.R starting. method=%s geno_csv=%s pheno=%s chunk_size=%d strict_native=%s\n', method, geno_csv, pheno, chunk_size, strict_native))

suppressMessages({
  require(data.table)
  require(MASS)
})

# read phenotype
ph <- tryCatch({
  data.table::fread(pheno, data.table=FALSE)
}, error=function(e){ stop('Could not read phenotype file: ', e$message) })

if(ncol(ph) < 2){
  stop('Phenotype file must contain sample ID and trait column')
}

sample_id_col <- names(ph)[1]
trait_col <- if(trait != '') trait else names(ph)[2]
if(!(trait_col %in% names(ph))) stop('Trait column not found in phenotype file')
Y <- ph[[trait_col]]
PHENO_SAMPLES <- as.character(ph[[sample_id_col]])

# read geno table (PLINK --recode A raw or CSV). detect PLINK raw style (FID IID ...) and use columns after 6
geno_dt <- tryCatch({
  data.table::fread(geno_csv, data.table=FALSE)
}, error=function(e){ stop('Could not read genotype file: ', e$message) })

col_names <- names(geno_dt)
start_idx <- 1
if(all(c('FID','IID') %in% col_names[1:2])){
  # PLINK --recode A: first 6 columns are FID IID PAT MAT SEX PHEN
  start_idx <- 7
}

sample_ids_geno <- as.character(geno_dt[[1]])
if(start_idx == 7){
  # ensure sample order matches phenotype by IID (second column)
  sample_ids_geno <- as.character(geno_dt[[2]])
}

if(!all(PHENO_SAMPLES %in% sample_ids_geno)){
  warning('Not all phenotype samples found in genotype file; performing inner join by sample ID')
}

# create genotype matrix with rows matching phenotype order
geno_map <- data.frame()
geno_mat <- NULL
try({
  geno_df <- geno_dt
  if(start_idx == 7){
    ids <- as.character(geno_df[[2]])
  } else {
    ids <- as.character(geno_df[[1]])
  }
  # subset and reorder genotype rows to phenotype order by sample ID
  keep_idx <- match(PHENO_SAMPLES, ids)          # phenotype order -> genotype row index
  keep_rows <- which(!is.na(keep_idx))           # phenotype rows that exist in genotype
  if(length(keep_rows) == 0) stop('No overlapping sample IDs between phenotype and genotype')

  # Align phenotype vectors to the overlapping subset so Y and GD have identical nrow.
  PHENO_SAMPLES <- PHENO_SAMPLES[keep_rows]
  Y <- Y[keep_rows]

  geno_sub <- geno_df[keep_idx[keep_rows], start_idx:ncol(geno_df), drop=FALSE]
  rownames(geno_sub) <- PHENO_SAMPLES
  geno_mat <- as.matrix(geno_sub)
  # Build marker map with required 3 columns: SNP, chromosome, position.
  # If coordinates are not encoded in SNP IDs, use stable placeholders.
  snp_ids <- colnames(geno_mat)
  chr_guess <- rep(1L, length(snp_ids))
  pos_guess <- seq_along(snp_ids)
  m <- regexec('^([^_]+)_([0-9]+)$', snp_ids)
  mm <- regmatches(snp_ids, m)
  for(i in seq_along(mm)){
    if(length(mm[[i]]) == 3){
      chr_guess[i] <- suppressWarnings(as.integer(mm[[i]][2]))
      pos_guess[i] <- suppressWarnings(as.integer(mm[[i]][3]))
    }
  }
  chr_guess[!is.finite(chr_guess)] <- 1L
  pos_guess[!is.finite(pos_guess)] <- seq_len(sum(!is.finite(pos_guess)))
  geno_map <- data.frame(SNP=snp_ids, Chromosome=chr_guess, Position=pos_guess, stringsAsFactors=FALSE)
}, silent=TRUE)
if(is.null(geno_mat)) stop('Failed to construct genotype matrix')

# Try to use requested package if available; otherwise fallback to chunked linear regression
use_pkg <- FALSE
blink_ns <- NULL
if(method == 'blink'){
  blink_ns <- detect_blink_namespace()
}
if(method == 'blink' && !is.null(blink_ns)){
  use_pkg <- TRUE
}
if(method == 'farmcpu' && (requireNamespace('FarmCPUpp', quietly=TRUE) || requireNamespace('FarmCPU', quietly=TRUE))){
  use_pkg <- TRUE
}
if(method == 'mrmlm' && requireNamespace('mrMLM', quietly=TRUE)){
  use_pkg <- TRUE
  suppressMessages(library(mrMLM))
}

results <- data.frame(SNP=character(0), beta=numeric(0), se=numeric(0), p=numeric(0), stringsAsFactors=FALSE)

pick_col <- function(df, candidates){
  nm <- names(df)
  nml <- tolower(nm)
  for(cn in candidates){
    hit <- which(nml == tolower(cn))
    if(length(hit) > 0) return(nm[hit[1]])
  }
  return(NULL)
}

if(use_pkg){
  cat(sprintf('%s package found — attempting to run native method.\n', toupper(method)))
  Ynum <- suppressWarnings(as.numeric(Y))
  keep_trait <- which(!is.na(Ynum))
  if(length(keep_trait) == 0) stop('No non-missing trait values after alignment')

  PHENO_SAMPLES <- PHENO_SAMPLES[keep_trait]
  Ynum <- Ynum[keep_trait]
  geno_mat <- geno_mat[keep_trait, , drop=FALSE]

  suppressWarnings(storage.mode(geno_mat) <- 'numeric')
  if(anyNA(geno_mat)){
    cm <- colMeans(geno_mat, na.rm=TRUE)
    idx <- which(is.na(geno_mat), arr.ind=TRUE)
    geno_mat[idx] <- cm[idx[, 2]]
  }

  # Fast invariant marker filter to prevent matrix inversion crashes in mrMLM / FarmCPU
  cm <- colMeans(geno_mat, na.rm=TRUE)
  cm2 <- colMeans(geno_mat^2, na.rm=TRUE)
  vars <- cm2 - cm^2
  keep_snps <- which(vars > 1e-8)
  if(length(keep_snps) < ncol(geno_mat)){
    cat(sprintf('Filtering %d invariant markers for native engine stability.\n', ncol(geno_mat) - length(keep_snps)))
    geno_mat <- geno_mat[, keep_snps, drop=FALSE]
    geno_map <- geno_map[keep_snps, , drop=FALSE]
  }

  snp_ids <- colnames(geno_mat)
  if(is.null(snp_ids) || length(snp_ids) != ncol(geno_mat)){
    snp_ids <- geno_map$SNP
  }

  chr_guess <- suppressWarnings(as.integer(geno_map$Chromosome))
  pos_guess <- suppressWarnings(as.integer(geno_map$Position))
  chr_guess[!is.finite(chr_guess)] <- 1L
  pos_guess[!is.finite(pos_guess)] <- seq_len(sum(!is.finite(pos_guess)))

  out_file <- paste0(out_pref, '_', method, '_results.csv')

  native_ok <- FALSE
  native_err <- NULL

  tryCatch({
    native_out <- NULL

    if(method == 'blink'){
      Ydf <- data.frame(taxon=PHENO_SAMPLES, trait=Ynum)
      GD <- as.matrix(geno_mat)
      suppressWarnings(storage.mode(GD) <- 'numeric')
      GM <- data.frame(SNP=snp_ids, Chromosome=chr_guess, Position=pos_guess, stringsAsFactors=FALSE)
      CV0 <- matrix(numeric(0), nrow=nrow(Ydf), ncol=0)

      if(requireNamespace('FarmCPU', quietly=TRUE)) suppressMessages(library(FarmCPU))
      if (is.null(blink_ns)) {
        # strict guard against unrelated CRAN package named 'blink'
        ex_BLINK <- if (requireNamespace('BLINK', quietly=TRUE)) paste(getNamespaceExports('BLINK'), collapse=', ') else 'not installed'
        ex_blink <- if (requireNamespace('blink', quietly=TRUE)) paste(getNamespaceExports('blink'), collapse=', ') else 'not installed'
        stop(sprintf("BLINK native API unavailable. BLINK exports: [%s] ; blink exports: [%s]", ex_BLINK, ex_blink))
      }
      blink_fun_name <- pick_blink_fun(blink_ns)
      if (is.null(blink_fun_name)) {
        stop(sprintf('BLINK package namespace %s is present but no BLINK/Blink entrypoint was found', blink_ns))
      }

      blink_fun <- get(blink_fun_name, envir=asNamespace(blink_ns))
      need_farmcpu_shim <- FALSE
      if (identical(blink_ns, 'BLINK')) {
        if (!requireNamespace('FarmCPU', quietly=TRUE)) {
          need_farmcpu_shim <- TRUE
        } else {
          farmcpu_ns <- asNamespace('FarmCPU')
          has_prior <- exists('FarmCPU.Prior', envir=farmcpu_ns, mode='function')
          has_lm <- exists('FarmCPU.LM', envir=farmcpu_ns, mode='function')
          need_farmcpu_shim <- (!has_prior || !has_lm)
        }
      }
      if (need_farmcpu_shim) {
        cat('BLINK compatibility mode: FarmCPU API missing; using runtime shim for FarmCPU.Prior/FarmCPU.LM.\n')
        blink_fun <- make_blink_with_farmcpu_shim(blink_fun, blink_ns)
      }
      native_out <- blink_fun(Y=Ydf, GD=GD, GM=GM, CV=CV0, ncpus=1, file.output=FALSE, maxLoop=2)

      gw <- NULL
      if(is.list(native_out) && 'GWAS' %in% names(native_out)) gw <- as.data.frame(native_out$GWAS, stringsAsFactors=FALSE)
      if(is.data.frame(native_out)) gw <- native_out
      if(is.null(gw) || nrow(gw) == 0) stop('BLINK returned no GWAS table')

      c_snp <- pick_col(gw, c('SNP','rs#','RS#'))
      c_p <- pick_col(gw, c('P.value','p.value','p','P'))
      c_beta <- pick_col(gw, c('effect','beta','estimate'))
      c_se <- pick_col(gw, c('se','stderr','std.error'))
      if(is.null(c_snp) || is.null(c_p)) stop('BLINK GWAS table missing SNP or P-value columns')

      results <- data.frame(
        SNP = as.character(gw[[c_snp]]),
        beta = if(!is.null(c_beta)) suppressWarnings(as.numeric(gw[[c_beta]])) else NA_real_,
        se = if(!is.null(c_se)) suppressWarnings(as.numeric(gw[[c_se]])) else NA_real_,
        p = suppressWarnings(as.numeric(gw[[c_p]])),
        stringsAsFactors=FALSE
      )
    }

    if(method == 'farmcpu'){
      if(requireNamespace('FarmCPUpp', quietly=TRUE)){
        suppressMessages(library(bigmemory))
        suppressMessages(library(FarmCPUpp))
        Ypp <- data.frame(taxon=PHENO_SAMPLES, trait=Ynum)
        GMpp <- data.frame(marker=snp_ids, chromosome=chr_guess, position=pos_guess, stringsAsFactors=FALSE)
        GD_big <- bigmemory::as.big.matrix(geno_mat, backingfile='farmcpupp.bin', descriptorfile='farmcpupp.desc', backingpath=tempdir())

        native_out <- FarmCPUpp::farmcpu(
          Y=Ypp, GD=GD_big, GM=GMpp,
          ncores.glm=1, ncores.reml=1,
          iteration.output=FALSE
        )
        trait_key <- names(native_out)[1]
        gw <- native_out[[trait_key]]$GWAS
        if(is.null(gw) || nrow(gw) == 0) stop('FarmCPUpp returned no GWAS table')

        c_snp <- pick_col(gw, c('marker','SNP','rs#','RS#'))
        c_p <- pick_col(gw, c('p.value','P.value','p','P'))
        c_beta <- pick_col(gw, c('estimate','effect','beta'))
        c_se <- pick_col(gw, c('stderr','se','std.error'))
        if(is.null(c_snp) || is.null(c_p)) stop('FarmCPUpp GWAS table missing marker or P-value columns')

        results <- data.frame(
          SNP = as.character(gw[[c_snp]]),
          beta = if(!is.null(c_beta)) suppressWarnings(as.numeric(gw[[c_beta]])) else NA_real_,
          se = if(!is.null(c_se)) suppressWarnings(as.numeric(gw[[c_se]])) else NA_real_,
          p = suppressWarnings(as.numeric(gw[[c_p]])),
          stringsAsFactors=FALSE
        )
      } else {
        stop('FarmCPUpp package is required for robust native FARMCPU execution in this environment')
      }
    }

    if(method == 'mrmlm'){
      gen_mr <- cbind(chr_guess, pos_guess, t(geno_mat))
      phe_mr <- matrix(Ynum, ncol=1)

      keep_n <- min(8, nrow(geno_mat))
      hdr <- c('rs#', 'chrom', 'pos', 'genotype for code 1', PHENO_SAMPLES[1:keep_n])
      body <- cbind(snp_ids, chr_guess, pos_guess, rep(1, length(snp_ids)), t(geno_mat[1:keep_n, , drop=FALSE]))
      genRaw_light <- rbind(hdr, body)

      native_out <- mrMLM::FASTmrMLM(
        gen = gen_mr,
        phe = phe_mr,
        outATCG = NULL,
        genRaw = genRaw_light,
        kk = NULL,
        psmatrix = NULL,
        svpal = 0.01,
        svrad = 20,
        svmlod = 3,
        Genformat = 1,
        CLO = 1
      )

      if(!is.list(native_out) || !('result1' %in% names(native_out))) stop('mrMLM returned no result1 table')
      r1 <- as.data.frame(native_out$result1, stringsAsFactors=FALSE)
      c_snp <- pick_col(r1, c('RS#','rs#','SNP'))
      c_logp <- pick_col(r1, c("'-log10(P) (FASTmrMLM)'", "-log10(P) (FASTmrMLM)", "-log10(P)", "LOD score"))
      c_beta <- pick_col(r1, c('SNP effect (FASTmrMLM)','QTN effect','effect','beta'))
      if(is.null(c_snp) || is.null(c_logp)) stop('mrMLM result1 missing SNP or -log10(P) columns')

      logp <- suppressWarnings(as.numeric(r1[[c_logp]]))
      pvals <- suppressWarnings(10^(-logp))
      results <- data.frame(
        SNP = as.character(r1[[c_snp]]),
        beta = if(!is.null(c_beta)) suppressWarnings(as.numeric(r1[[c_beta]])) else NA_real_,
        se = NA_real_,
        p = pvals,
        stringsAsFactors=FALSE
      )
    }

    if(is.null(results) || nrow(results) == 0) stop(sprintf('%s native execution produced no rows', toupper(method)))
    results <- results[is.finite(results$p) & !is.na(results$SNP), , drop=FALSE]
    if(nrow(results) == 0) stop(sprintf('%s native execution produced no finite p-values', toupper(method)))
    write.csv(results, file=out_file, row.names=FALSE)
    cat(sprintf('%s finished; results saved to %s\n', toupper(method), out_file))
    native_ok <- TRUE
  }, error=function(e){
    native_err <<- e$message
  })

  if(!native_ok){
    msg <- sprintf('%s native execution failed: %s', toupper(method), native_err)
    if(strict_native){
      stop(msg)
    }
    cat(sprintf('%s\n', msg))
    use_pkg <- FALSE
  }
}

if(!use_pkg){
  if(strict_native){
    stop(sprintf('%s package is not available in this R environment. strict_native=TRUE prevents fallback.', toupper(method)))
  }
  cat(sprintf('%s package not available or failed — running chunked linear regression fallback.\n', toupper(method)))
  n_snps <- ncol(geno_mat)
  chunk_starts <- seq(1, n_snps, by=chunk_size)
  for(s in chunk_starts){
    e <- min(n_snps, s + chunk_size - 1)
    mat_chunk <- geno_mat[, s:e, drop=FALSE]
    # for each SNP, fit linear model trait ~ SNP + (optionally covariates if provided)
    for(j in seq_len(ncol(mat_chunk))){
      snpvec <- mat_chunk[, j]
      if(all(is.na(snpvec))) next
      df <- data.frame(trait=Y, snp=snpvec)
      if(!is.null(covar) && nzchar(covar)){
        covdf <- tryCatch({ data.table::fread(covar, data.table=FALSE) }, error=function(e) NULL)
        if(!is.null(covdf)){
          # align covariates to phenotype samples by first column
          cov_ids <- as.character(covdf[[1]])
          cov_keep <- match(PHENO_SAMPLES, cov_ids)
          cov_sub <- covdf[which(!is.na(cov_keep)), , drop=FALSE]
          if(nrow(cov_sub) == nrow(df)){
            df <- cbind(df, cov_sub[,-1, drop=FALSE])
          }
        }
      }
      # quick try; wrap in try
      fit <- tryCatch({ lm(trait ~ snp, data=df) }, error=function(e) NULL)
      if(is.null(fit)) next
      ss <- summary(fit)
      coefs <- ss$coefficients
      if('snp' %in% rownames(coefs)){
        snp_name <- colnames(mat_chunk)[j]
        # PLINK --recode A uses IDs like 1178_T; strip allele suffix to match notebook SNP_ID.
        snp_name <- sub('_[ACGTN]+$', '', snp_name)
        beta <- coefs['snp', 'Estimate']
        se <- coefs['snp', 'Std. Error']
        p <- coefs['snp', 'Pr(>|t|)']
        results <- rbind(results, data.frame(SNP=snp_name, beta=beta, se=se, p=p, stringsAsFactors=FALSE))
      }
    }
    # write intermediate progress
    cat(sprintf('Processed SNPs %d-%d (%d total)\n', s, e, n_snps))
  }
  out_file <- paste0(out_pref, '_', method, '_fallback_results.csv')
  write.csv(results, out_file, row.names=FALSE)
  cat('Fallback results written to', out_file, '\n')
}

cat('run_blink.R completed.\n')

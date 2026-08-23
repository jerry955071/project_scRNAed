# ============================================================
# 1. Extract one locus from the global alt_mtx / ref_mtx
# ============================================================

get_locus_data <- function(locus_id) {
  
  # alt_mtx and ref_mtx are intentionally obtained from
  # the global environment and are NOT passed as arguments.
  
  stopifnot(
    locus_id %in% Features(alt_mtx),
    locus_id %in% Features(ref_mtx)
  )
  
  # Identify trustworthy metacells by cell name
  trustworthy_cells <- rownames(alt_mtx@meta.data)[
    alt_mtx@meta.data$mcRigor_sc == "trustworthy"
  ]
  
  stopifnot(
    all(trustworthy_cells %in% Cells(ref_mtx))
  )
  
  # Retrieve ALT / REF counts in exactly the same cell order
  alt_locus <- LayerData(
    alt_mtx,
    assay = "RNA",
    layer = "counts",
    cells = trustworthy_cells,
    features = locus_id
  )
  
  ref_locus <- LayerData(
    ref_mtx,
    assay = "RNA",
    layer = "counts",
    cells = trustworthy_cells,
    features = locus_id
  )
  
  stopifnot(
    identical(
      colnames(alt_locus),
      colnames(ref_locus)
    )
  )
  
  # Metadata in matching order
  dat <- ref_mtx@meta.data[
    colnames(ref_locus),
    ,
    drop = FALSE
  ] %>%
    mutate(
      alt = as.numeric(alt_locus),
      ref = as.numeric(ref_locus),
      
      treatment = relevel(
        factor(treatment),
        ref = "normal"
      ),
      
      total = alt + ref
    ) %>%
    filter(
      total > 0
    )
  
  dat
}


# ============================================================
# 2. Prepare one lineage dataset
#
# Used by BOTH model fitting and visualization
# ============================================================

prepare_lineage_data <- function(
    dat,
    lineage,
    weight_cutoff = 0.5
) {
  
  pt_col <- paste0(lineage, "_pseudotime")
  wt_col <- paste0(lineage, "_weight")
  
  stopifnot(
    pt_col %in% colnames(dat),
    wt_col %in% colnames(dat)
  )
  
  dat %>%
    filter(
      .data[[wt_col]] > weight_cutoff,
      is.finite(.data[[pt_col]])
    ) %>%
    mutate(
      lineage_pseudotime = .data[[pt_col]],
      editing_level = alt / total
    )
}


# ============================================================
# 3. Helpers
# ============================================================

# ------------------------------------------------------------
# Extract LRT p-value from an ANOVA table
# ------------------------------------------------------------

get_lrt_p <- function(x) {
  
  if (
    is.null(x) ||
    !"Pr(>Chisq)" %in% colnames(x) ||
    nrow(x) < 2
  ) {
    return(NA_real_)
  }
  
  as.numeric(
    x[["Pr(>Chisq)"]][2]
  )
}


# ------------------------------------------------------------
# Check whether a glmmTMB model is suitable for inference
# ------------------------------------------------------------

model_ok <- function(model) {
  
  if (
    is.null(model) ||
    !inherits(model, "glmmTMB")
  ) {
    return(FALSE)
  }
  
  pdhess <- tryCatch(
    isTRUE(model$sdr$pdHess),
    error = function(e) FALSE
  )
  
  ll_ok <- tryCatch(
    is.finite(
      as.numeric(logLik(model))
    ),
    error = function(e) FALSE
  )
  
  pdhess && ll_ok
}


# ------------------------------------------------------------
# Only perform an LRT if BOTH models are valid
# ------------------------------------------------------------

safe_lrt_p <- function(
    model_small,
    model_large
) {
  
  if (
    !model_ok(model_small) ||
    !model_ok(model_large)
  ) {
    return(NA_real_)
  }
  
  test <- tryCatch(
    anova(
      model_small,
      model_large,
      test = "Chisq"
    ),
    error = function(e) NULL
  )
  
  get_lrt_p(test)
}


# ------------------------------------------------------------
# Extract one conditional-model coefficient
# ------------------------------------------------------------

get_coef <- function(
    model,
    term
) {
  
  tab <- tryCatch(
    summary(model)$coefficients$cond,
    error = function(e) NULL
  )
  
  if (
    is.null(tab) ||
    !term %in% rownames(tab)
  ) {
    return(
      c(
        beta = NA_real_,
        se = NA_real_,
        p = NA_real_,
        OR = NA_real_
      )
    )
  }
  
  beta <- tab[term, "Estimate"]
  
  c(
    beta = beta,
    se = tab[term, "Std. Error"],
    p = tab[term, "Pr(>|z|)"],
    OR = exp(beta)
  )
}


# ------------------------------------------------------------
# Extract random-effect SD
# ------------------------------------------------------------

get_random_sd <- function(
    model,
    group
) {
  
  tryCatch(
    {
      vc <- VarCorr(model)$cond[[group]]
      
      as.numeric(
        attr(vc, "stddev")[1]
      )
    },
    error = function(e) NA_real_
  )
}


# ============================================================
# 4. Fit the four models for one lineage
# ============================================================

fit_lineage_models <- function(
    dat,
    lineage = "L1",
    family_use = binomial(link = "logit"),
    weight_cutoff = 0.5,
    spline_df = 3
) {
  
  # ----------------------------------------------------------
  # Prepare lineage-specific observations
  # ----------------------------------------------------------
  
  dat_lineage <- prepare_lineage_data(
    dat = dat,
    lineage = lineage,
    weight_cutoff = weight_cutoff
  )
  
  # ----------------------------------------------------------
  # Basic checks
  # ----------------------------------------------------------
  
  if (nrow(dat_lineage) < spline_df + 5) {
    stop(
      "Insufficient observations for ",
      lineage
    )
  }
  
  if (
    length(
      unique(dat_lineage$lineage_pseudotime)
    ) < spline_df + 1
  ) {
    stop(
      "Insufficient unique pseudotime values for ",
      lineage
    )
  }
  
  # Count treatment groups ACTUALLY represented in the
  # lineage-specific data rather than retained factor levels
  n_treatments <- dplyr::n_distinct(
    as.character(dat_lineage$treatment)
  )
  
  if (n_treatments < 2) {
    stop(
      "Insufficient treatment groups for ",
      lineage
    )
  }
  
  # ----------------------------------------------------------
  # M0: null
  # ----------------------------------------------------------
  
  m0 <- glmmTMB(
    cbind(alt, ref) ~
      (1 | individual) +
      (1 | sample),
    family = family_use,
    data = dat_lineage
  )
  
  # ----------------------------------------------------------
  # MT: treatment only
  # ----------------------------------------------------------
  
  m_treatment <- glmmTMB(
    cbind(alt, ref) ~
      treatment +
      (1 | individual) +
      (1 | sample),
    family = family_use,
    data = dat_lineage
  )
  
  # ----------------------------------------------------------
  # MP: pseudotime only
  # ----------------------------------------------------------
  
  m_pseudotime <- glmmTMB(
    cbind(alt, ref) ~
      splines::ns(
        lineage_pseudotime,
        df = spline_df
      ) +
      (1 | individual) +
      (1 | sample),
    family = family_use,
    data = dat_lineage
  )
  
  # ----------------------------------------------------------
  # MA: treatment + pseudotime
  # ----------------------------------------------------------
  
  m_additive <- glmmTMB(
    cbind(alt, ref) ~
      treatment +
      splines::ns(
        lineage_pseudotime,
        df = spline_df
      ) +
      (1 | individual) +
      (1 | sample),
    family = family_use,
    data = dat_lineage
  )
  
  # ----------------------------------------------------------
  # Check model validity ONCE
  # ----------------------------------------------------------
  
  model_status <- c(
    null =
      model_ok(m0),
    
    treatment =
      model_ok(m_treatment),
    
    pseudotime =
      model_ok(m_pseudotime),
    
    additive =
      model_ok(m_additive)
  )
  
  # ----------------------------------------------------------
  # Likelihood-ratio tests
  #
  # Invalid model pairs automatically return NA.
  # ----------------------------------------------------------
  
  tests <- list(
    
    # Treatment without adjusting for pseudotime
    treatment_raw =
      safe_lrt_p(
        m0,
        m_treatment
      ),
    
    # Pseudotime without adjusting for treatment
    pseudotime_raw =
      safe_lrt_p(
        m0,
        m_pseudotime
      ),
    
    # Pseudotime after controlling for treatment
    pseudotime_adjusted =
      safe_lrt_p(
        m_treatment,
        m_additive
      ),
    
    # Treatment after controlling for pseudotime
    treatment_adjusted =
      safe_lrt_p(
        m_pseudotime,
        m_additive
      )
  )
  
  list(
    lineage = lineage,
    family = family_use$family,
    data = dat_lineage,
    
    models = list(
      null = m0,
      treatment = m_treatment,
      pseudotime = m_pseudotime,
      additive = m_additive
    ),
    
    model_status = model_status,
    
    tests = tests
  )
}


# ============================================================
# 5. Convert one lineage fit into one compact typed result
# ============================================================

extract_lineage_result <- function(
    fit,
    lineage = fit$lineage
) {
  
  dat <- fit$data
  models <- fit$models
  model_status <- fit$model_status
  
  additive <- models$additive
  
  # ----------------------------------------------------------
  # Treatment coefficients
  # ----------------------------------------------------------
  
  opposite <- get_coef(
    additive,
    "treatmentopposite"
  )
  
  tension <- get_coef(
    additive,
    "treatmenttension"
  )
  
  # ----------------------------------------------------------
  # Number of observations per treatment
  # ----------------------------------------------------------
  
  treatment_counts <- table(
    as.character(dat$treatment)
  )
  
  n_normal <- if (
    "normal" %in% names(treatment_counts)
  ) {
    as.integer(
      treatment_counts["normal"]
    )
  } else {
    0L
  }
  
  n_opposite <- if (
    "opposite" %in% names(treatment_counts)
  ) {
    as.integer(
      treatment_counts["opposite"]
    )
  } else {
    0L
  }
  
  n_tension <- if (
    "tension" %in% names(treatment_counts)
  ) {
    as.integer(
      treatment_counts["tension"]
    )
  } else {
    0L
  }
  
  # ----------------------------------------------------------
  # Dispersion
  # ----------------------------------------------------------
  
  dispersion <- tryCatch(
    sigma(additive),
    error = function(e) NA_real_
  )
  
  # ----------------------------------------------------------
  # IMPORTANT:
  # list() preserves numeric/logical/character data types.
  # ----------------------------------------------------------
  
  out <- list(
    
    # Model family
    family =
      fit$family,
    
    # ------------------------------------------
    # Data / QC
    # ------------------------------------------
    
    n_obs =
      as.integer(
        nrow(dat)
      ),
    
    n_samples =
      as.integer(
        dplyr::n_distinct(
          dat$sample
        )
      ),
    
    n_individuals =
      as.integer(
        dplyr::n_distinct(
          dat$individual
        )
      ),
    
    n_normal =
      n_normal,
    
    n_opposite =
      n_opposite,
    
    n_tension =
      n_tension,
    
    mean_depth =
      mean(
        dat$total,
        na.rm = TRUE
      ),
    
    # ------------------------------------------
    # Main likelihood-ratio tests
    # ------------------------------------------
    
    p_treatment_raw =
      fit$tests$treatment_raw,
    
    p_pseudotime_raw =
      fit$tests$pseudotime_raw,
    
    p_treatment_adj =
      fit$tests$treatment_adjusted,
    
    p_pseudotime_adj =
      fit$tests$pseudotime_adjusted,
    
    # ------------------------------------------
    # Treatment effect sizes from additive model
    # ------------------------------------------
    
    beta_opposite =
      unname(
        opposite["beta"]
      ),
    
    se_opposite =
      unname(
        opposite["se"]
      ),
    
    p_opposite =
      unname(
        opposite["p"]
      ),
    
    OR_opposite =
      unname(
        opposite["OR"]
      ),
    
    beta_tension =
      unname(
        tension["beta"]
      ),
    
    se_tension =
      unname(
        tension["se"]
      ),
    
    p_tension =
      unname(
        tension["p"]
      ),
    
    OR_tension =
      unname(
        tension["OR"]
      ),
    
    # ------------------------------------------
    # Random effects
    # ------------------------------------------
    
    sd_individual =
      get_random_sd(
        additive,
        "individual"
      ),
    
    sd_sample =
      get_random_sd(
        additive,
        "sample"
      ),
    
    # ------------------------------------------
    # Model convergence
    # ------------------------------------------
    
    ok_null =
      unname(
        model_status["null"]
      ),
    
    ok_treatment =
      unname(
        model_status["treatment"]
      ),
    
    ok_pseudotime =
      unname(
        model_status["pseudotime"]
      ),
    
    ok_additive =
      unname(
        model_status["additive"]
      ),
    
    dispersion =
      dispersion
  )
  
  names(out) <- paste0(
    lineage,
    "_",
    names(out)
  )
  
  out
}


# ============================================================
# 6. Empty result if one lineage cannot be fitted
# ============================================================

empty_lineage_result <- function(
    lineage,
    status = "fit_failed"
) {
  
  template <- list(
    
    family = NA_character_,
    
    n_obs = NA_integer_,
    n_samples = NA_integer_,
    n_individuals = NA_integer_,
    
    n_normal = NA_integer_,
    n_opposite = NA_integer_,
    n_tension = NA_integer_,
    
    mean_depth = NA_real_,
    
    p_treatment_raw = NA_real_,
    p_pseudotime_raw = NA_real_,
    p_treatment_adj = NA_real_,
    p_pseudotime_adj = NA_real_,
    
    beta_opposite = NA_real_,
    se_opposite = NA_real_,
    p_opposite = NA_real_,
    OR_opposite = NA_real_,
    
    beta_tension = NA_real_,
    se_tension = NA_real_,
    p_tension = NA_real_,
    OR_tension = NA_real_,
    
    sd_individual = NA_real_,
    sd_sample = NA_real_,
    
    ok_null = FALSE,
    ok_treatment = FALSE,
    ok_pseudotime = FALSE,
    ok_additive = FALSE,
    
    dispersion = NA_real_,
    
    status = as.character(status)
  )
  
  names(template) <- paste0(
    lineage,
    "_",
    names(template)
  )
  
  template
}


# ============================================================
# 7. Safely analyze one lineage
# ============================================================

analyze_lineage <- function(
    dat,
    lineage,
    family_use = binomial(link = "logit"),
    weight_cutoff = 0.5,
    spline_df = 3
) {
  
  fit <- tryCatch(
    
    suppressWarnings(
      fit_lineage_models(
        dat = dat,
        lineage = lineage,
        family_use = family_use,
        weight_cutoff = weight_cutoff,
        spline_df = spline_df
      )
    ),
    
    error = function(e) e
  )
  
  # ----------------------------------------------------------
  # Complete fitting failure
  # ----------------------------------------------------------
  
  if (inherits(fit, "error")) {
    
    return(
      empty_lineage_result(
        lineage,
        status = conditionMessage(fit)
      )
    )
  }
  
  # ----------------------------------------------------------
  # Extract compact results
  # ----------------------------------------------------------
  
  result <- extract_lineage_result(
    fit,
    lineage
  )
  
  # model validity was already calculated in
  # fit_lineage_models()
  all_ok <- all(
    fit$model_status
  )
  
  result[
    paste0(
      lineage,
      "_status"
    )
  ] <- if (all_ok) {
    "OK"
  } else {
    "nonconverged"
  }
  
  result
}


# ============================================================
# 8. Analyze one locus across L1, L2, L3
#    Returns exactly ONE ROW
# ============================================================

analyze_locus <- function(
    locus_id,
    lineages = c("L1", "L2", "L3"),
    family_use = binomial(link = "logit"),
    weight_cutoff = 0.5,
    spline_df = 3
) {
  
  dat <- tryCatch(
    get_locus_data(locus_id),
    error = function(e) e
  )
  
  # Failure while extracting locus data
  if (inherits(dat, "error")) {
    
    out <- list(
      locus = locus_id,
      locus_status = conditionMessage(dat)
    )
    
    for (lineage in lineages) {
      
      out <- c(
        out,
        empty_lineage_result(
          lineage,
          status = "locus_data_failed"
        )
      )
    }
    
    return(
      as.data.frame(
        out,
        check.names = FALSE
      )
    )
  }
  
  # Locus-level QC
  out <- list(
    locus = locus_id,
    locus_status = "OK",
    locus_n_obs = nrow(dat),
    locus_mean_depth =
      mean(
        dat$total,
        na.rm = TRUE
      )
  )
  
  # Fit each lineage using the SAME locus dataframe
  for (lineage in lineages) {
    
    lineage_result <- analyze_lineage(
      dat = dat,
      lineage = lineage,
      family_use = family_use,
      weight_cutoff = weight_cutoff,
      spline_df = spline_df
    )
    
    out <- c(
      out,
      lineage_result
    )
  }
  
  as.data.frame(
    out,
    check.names = FALSE
  )
}


# ============================================================
# 9. Run multiple loci
# ============================================================

run_loci <- function(
    loci,
    lineages = c("L1", "L2", "L3"),
    family_use = binomial(link = "logit"),
    weight_cutoff = 0.5,
    spline_df = 3,
    verbose = TRUE
) {
  
  results <- vector(
    "list",
    length(loci)
  )
  
  for (i in seq_along(loci)) {
    
    if (
      verbose &&
      (
        i == 1 ||
        i %% 100 == 0 ||
        i == length(loci)
      )
    ) {
      message(
        "[",
        i,
        "/",
        length(loci),
        "] ",
        loci[i]
      )
    }
    
    results[[i]] <- analyze_locus(
      locus_id = loci[i],
      lineages = lineages,
      family_use = family_use,
      weight_cutoff = weight_cutoff,
      spline_df = spline_df
    )
  }
  
  dplyr::bind_rows(results)
}


# ============================================================
# 10. Run multiple loci in parallel
# ============================================================

run_loci_parallel <- function(
    loci,
    lineages = c("L1", "L2", "L3"),
    family_use = binomial(link = "logit"),
    weight_cutoff = 0.5,
    spline_df = 3,
    n_cores = 40
) {
  
  message(
    "Analyzing ",
    length(loci),
    " loci using ",
    n_cores,
    " parallel workers..."
  )
  
  results <- parallel::mclapply(
    X = seq_along(loci),
    
    FUN = function(i) {
      
      locus <- loci[i]
      
      # Optional progress output.
      # Messages from different workers may appear out of order.
      message(
        "[",
        i,
        "/",
        length(loci),
        "] ",
        locus
      )
      
      analyze_locus(
        locus_id = locus,
        lineages = lineages,
        family_use = family_use,
        weight_cutoff = weight_cutoff,
        spline_df = spline_df
      )
    },
    
    mc.cores = n_cores,
    
    # 40 persistent-ish chunks rather than forking separately
    # for every locus.
    mc.preschedule = TRUE,
    
    # No RNG is used in the fitting pipeline
    mc.set.seed = FALSE
  )
  
  dplyr::bind_rows(results)
}


# ============================================================
# 11. Add FDR and final model classification
# ============================================================

add_fdr <- function(
    results,
    lineages = c("L1", "L2", "L3"),
    alpha = 0.05,
    method = "BH"
) {
  
  for (lineage in lineages) {
    
    p_treatment_col <- paste0(
      lineage,
      "_p_treatment_adj"
    )
    
    p_pseudotime_col <- paste0(
      lineage,
      "_p_pseudotime_adj"
    )
    
    fdr_treatment_col <- paste0(
      lineage,
      "_FDR_treatment"
    )
    
    fdr_pseudotime_col <- paste0(
      lineage,
      "_FDR_pseudotime"
    )
    
    model_col <- paste0(
      lineage,
      "_model"
    )
    
    results[[fdr_treatment_col]] <-
      p.adjust(
        results[[p_treatment_col]],
        method = method
      )
    
    results[[fdr_pseudotime_col]] <-
      p.adjust(
        results[[p_pseudotime_col]],
        method = method
      )
    
    treat_sig <-
      results[[fdr_treatment_col]] < alpha
    
    pt_sig <-
      results[[fdr_pseudotime_col]] < alpha
    
    results[[model_col]] <- dplyr::case_when(
      
      is.na(treat_sig) |
        is.na(pt_sig) ~
        "Not tested",
      
      treat_sig &
        pt_sig ~
        "Additive",
      
      treat_sig ~
        "Treatment",
      
      pt_sig ~
        "Pseudotime",
      
      TRUE ~
        "Null"
    )
  }
  
  results
}


# ============================================================
# 12. Refit one locus × lineage later for plotting
# ============================================================

refit_locus_lineage <- function(
    locus_id,
    lineage,
    family_use = binomial(link = "logit"),
    weight_cutoff = 0.5,
    spline_df = 3
) {
  
  dat <- get_locus_data(
    locus_id
  )
  
  fit_lineage_models(
    dat = dat,
    lineage = lineage,
    family_use = family_use,
    weight_cutoff = weight_cutoff,
    spline_df = spline_df
  )
}


# ============================================================
# VISUALIZATION
# ============================================================


# ============================================================
# Convert stored family name back into a family object
# ============================================================

family_from_result <- function(
    result_row,
    lineage,
    default = binomial(link = "logit")
) {
  
  family_col <- paste0(
    lineage,
    "_family"
  )
  
  if (!family_col %in% colnames(result_row))
    return(default)
  
  family_name <- as.character(
    result_row[[family_col]][1]
  )
  
  if (
    is.na(family_name) ||
    family_name == ""
  ) {
    return(default)
  }
  
  switch(
    family_name,
    
    "binomial" =
      binomial(link = "logit"),
    
    "betabinomial" =
      glmmTMB::betabinomial(
        link = "logit"
      ),
    
    default
  )
}


# ============================================================
# Fit ONLY the model selected in results
# ============================================================

fit_selected_lineage_model <- function(
    dat_lineage,
    model_class,
    family_use = binomial(link = "logit"),
    spline_df = 3
) {
  
  if (
    is.na(model_class) ||
    model_class == "Not tested"
  ) {
    return(NULL)
  }
  
  formula_use <- switch(
    
    model_class,
    
    "Null" =
      cbind(alt, ref) ~
        (1 | individual) +
        (1 | sample),
    
    "Treatment" =
      cbind(alt, ref) ~
        treatment +
        (1 | individual) +
        (1 | sample),
    
    "Pseudotime" =
      cbind(alt, ref) ~
        splines::ns(
          lineage_pseudotime,
          df = spline_df
        ) +
        (1 | individual) +
        (1 | sample),
    
    "Additive" =
      cbind(alt, ref) ~
        treatment +
        splines::ns(
          lineage_pseudotime,
          df = spline_df
        ) +
        (1 | individual) +
        (1 | sample),
    
    stop(
      "Unknown model class: ",
      model_class
    )
  )
  
  glmmTMB(
    formula_use,
    family = family_use,
    data = dat_lineage
  )
}


# ============================================================
# Create population-level fitted values
# ============================================================

get_fitted_curve <- function(
    model,
    dat_lineage,
    model_class,
    n_grid = 200
) {
  
  pt_grid <- seq(
    min(
      dat_lineage$lineage_pseudotime
    ),
    max(
      dat_lineage$lineage_pseudotime
    ),
    length.out = n_grid
  )
  
  # Treatment/Additive models require one trajectory
  # for each treatment
  if (
    model_class %in%
      c(
        "Treatment",
        "Additive"
      )
  ) {
    
    pred_dat <- tidyr::expand_grid(
      lineage_pseudotime = pt_grid,
      treatment =
        levels(
          dat_lineage$treatment
        )
    )
    
    pred_dat$treatment <- factor(
      pred_dat$treatment,
      levels = levels(
        dat_lineage$treatment
      )
    )
    
  } else {
    
    # Null/Pseudotime model:
    # one common fitted curve
    pred_dat <- tibble::tibble(
      lineage_pseudotime = pt_grid
    )
  }
  
  # Required by the random-effect terms in the formula.
  # re.form = NA means these REs are not included in predictions.
  pred_dat$individual <-
    dat_lineage$individual[1]
  
  pred_dat$sample <-
    dat_lineage$sample[1]
  
  pred <- predict(
    model,
    newdata = pred_dat,
    type = "link",
    se.fit = TRUE,
    re.form = NA
  )
  
  pred_dat %>%
    mutate(
      eta = pred$fit,
      eta_se = pred$se.fit,
      
      fit = plogis(
        eta
      ),
      
      lower = plogis(
        eta -
          1.96 * eta_se
      ),
      
      upper = plogis(
        eta +
          1.96 * eta_se
      )
    )
}


# ============================================================
# Helper for getting p/FDR values from one results row
#
# Returns the first FINITE value among candidates.
#
# Example:
#   FDR available -> use FDR
#   FDR is NA     -> use raw adjusted p-value
# ============================================================

get_result_value <- function(
    result_row,
    candidates
) {
  
  for (col in candidates) {
    
    if (!col %in% colnames(result_row)) {
      next
    }
    
    value <- suppressWarnings(
      as.numeric(
        result_row[[col]][1]
      )
    )
    
    if (
      length(value) == 1 &&
      is.finite(value)
    ) {
      return(value)
    }
  }
  
  NA_real_
}


format_p <- function(x) {
  
  if (!is.finite(x))
    return("NA")
  
  if (x < 0.001)
    return(
      format(
        x,
        scientific = TRUE,
        digits = 2
      )
    )
  
  sprintf(
    "%.3f",
    x
  )
}


# ============================================================
# Plot ONE lineage for ONE locus
# ============================================================

plot_lineage_result <- function(
    result_row,
    dat,
    lineage,
    weight_cutoff = 0.5,
    spline_df = 3,
    point_alpha = 0.25
) {
  
  model_col <- paste0(
    lineage,
    "_model"
  )
  
  model_class <- as.character(
    result_row[[model_col]][1]
  )
  
  # Same preparation function used for model fitting
  dat_lineage <- prepare_lineage_data(
    dat = dat,
    lineage = lineage,
    weight_cutoff = weight_cutoff
  )
  
  # --------------------------------------------
  # Statistics displayed in subtitle
  # --------------------------------------------
  
  p_treatment <- get_result_value(
    result_row,
    c(
      paste0(
        lineage,
        "_FDR_treatment"
      ),
      paste0(
        lineage,
        "_p_treatment_adj"
      )
    )
  )
  
  p_pseudotime <- get_result_value(
    result_row,
    c(
      paste0(
        lineage,
        "_FDR_pseudotime"
      ),
      paste0(
        lineage,
        "_p_pseudotime_adj"
      )
    )
  )
  
  subtitle_text <- paste0(
    model_class,
    "\nTreatment = ",
    format_p(
      p_treatment
    ),
    "\nPseudotime = ",
    format_p(
      p_pseudotime
    )
  )
  
  # --------------------------------------------
  # Base scatter plot
  # --------------------------------------------
  
  p <- ggplot(
    dat_lineage,
    aes(
      x = lineage_pseudotime,
      y = editing_level,
      color = treatment
    )
  ) +
    geom_point(
      aes(
        size = total
      ),
      alpha = point_alpha
    ) +
    scale_size_continuous(
      range = c(
        0.3,
        2
      )
    ) +
    theme_classic() +
    labs(
      title = lineage,
      subtitle = subtitle_text,
      x = paste0(
        lineage,
        " pseudotime"
      ),
      y = "Editing level",
      color = "Treatment",
      size = "Read depth"
    )
  
  # Nothing else to fit
  if (
    nrow(dat_lineage) == 0 ||
    is.na(model_class) ||
    model_class == "Not tested"
  ) {
    return(p)
  }
  
  # --------------------------------------------
  # Refit selected model
  # --------------------------------------------
  
  family_use <- family_from_result(
    result_row,
    lineage
  )
  
  fit <- tryCatch(
    fit_selected_lineage_model(
      dat_lineage = dat_lineage,
      model_class = model_class,
      family_use = family_use,
      spline_df = spline_df
    ),
    error = function(e) NULL
  )
  
  if (
    is.null(fit) ||
    !isTRUE(
      fit$sdr$pdHess
    )
  ) {
    
    return(
      p +
        annotate(
          "text",
          x = Inf,
          y = Inf,
          label = "Fit unavailable",
          hjust = 1.1,
          vjust = 1.5
        )
    )
  }
  
  # --------------------------------------------
  # Prediction
  # --------------------------------------------
  
  pred <- tryCatch(
    get_fitted_curve(
      model = fit,
      dat_lineage = dat_lineage,
      model_class = model_class
    ),
    error = function(e) NULL
  )
  
  if (is.null(pred))
    return(p)
  
  # --------------------------------------------
  # Treatment-specific curves
  # --------------------------------------------
  
  if (
    model_class %in%
      c(
        "Treatment",
        "Additive"
      )
  ) {
    
    p <- p +
      geom_ribbon(
        data = pred,
        aes(
          x = lineage_pseudotime,
          ymin = lower,
          ymax = upper,
          fill = treatment,
          group = treatment
        ),
        inherit.aes = FALSE,
        alpha = 0.12,
        color = NA
      ) +
      geom_line(
        data = pred,
        aes(
          x = lineage_pseudotime,
          y = fit,
          color = treatment,
          group = treatment
        ),
        inherit.aes = FALSE,
        linewidth = 1
      )
    
  } else {
    
    # ------------------------------------------
    # Common Null/Pseudotime curve
    # ------------------------------------------
    
    p <- p +
      geom_ribbon(
        data = pred,
        aes(
          x = lineage_pseudotime,
          ymin = lower,
          ymax = upper
        ),
        inherit.aes = FALSE,
        alpha = 0.12
      ) +
      geom_line(
        data = pred,
        aes(
          x = lineage_pseudotime,
          y = fit
        ),
        inherit.aes = FALSE,
        linewidth = 1
      )
  }
  
  p
}


# ============================================================
# Plot all 3 lineages for ONE locus
# ============================================================

plot_locus_result <- function(
    result_row,
    weight_cutoff = 0.5,
    spline_df = 3
) {
  
  locus_id <- as.character(
    result_row$locus[1]
  )
  
  # Important:
  # load locus once, not separately for every lineage
  dat <- get_locus_data(
    locus_id
  )
  
  plots <- lapply(
    c(
      "L1",
      "L2",
      "L3"
    ),
    function(lineage) {
      
      plot_lineage_result(
        result_row = result_row,
        dat = dat,
        lineage = lineage,
        weight_cutoff = weight_cutoff,
        spline_df = spline_df
      )
    }
  )
  
  patchwork::wrap_plots(
    plots,
    nrow = 1,
    guides = "collect"
  ) +
    patchwork::plot_annotation(
      title = locus_id
    ) &
    theme(
      legend.position = "right"
    )
}
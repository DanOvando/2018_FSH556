fit_hwk2 <- function(data, model = "glmm") {

  clam_data <- list(counts = data$counts,
                    sites = data$site - 1)

  clam_params <- list(
    mean_count = 0,
    site_means = rep(0, n_distinct(clam_data$sites)),
    overdisp = rep(0, nrow(data)),
    log_sigma_site = 0,
    log_sigma_overdisp = 0
  )


  fixed_pars <- list()

  if (model == "glm") {
    fixed_pars <- list(
      log_sigma_overdisp = factor(NA),
      overdisp = factor(rep(NA, length(
        clam_params$overdisp
      ))),
      log_sigma_site = factor(NA),
      site_means = factor(rep(NA, length(
        clam_params$site_means
      )))
    )

  } else if (model == "glmm-site") {
    fixed_pars <- list(log_sigma_overdisp = factor(NA),
                       overdisp = factor(rep(NA, length(
                         clam_params$overdisp
                       ))))

    # clam_params$log_sigma_site <- log(0.0001)

  } else if (model == "glmm-overdisp") {
    fixed_pars <- list(log_sigma_site = factor(NA),
                       site_means = factor(rep(NA, length(
                         clam_params$site_means
                       ))))


  }



  obj <-
    MakeADFun(
      data = clam_data,
      parameters = clam_params,
      DLL = "homework_2",
      random = c("overdisp", "site_means"),
      map = fixed_pars,
      silent = TRUE
    )

  opt <-
    nlminb(
      start = obj$par,
      objective = obj$fn,
      gradient = obj$gr
    )

  diagnostics = data.frame(
    "name" = names(obj$par),
    "Est" = opt$par,
    "final_gradient" = as.vector(obj$gr(opt$par))
  )

  clam_report <- obj$report()

  clam_estimates <-
    summary(sdreport(obj)) %>%
    as.data.frame() %>%
    mutate(variable = rownames(.)) %>%
    set_names(tolower) %>%
    rename(std_error = `std. error`) %>%
    mutate(lower = estimate - 1.96 * std_error,
           upper = estimate + 1.96 * std_error)

  return(clam_estimates)



}
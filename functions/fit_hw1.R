fit_hw1 <- function(data,
                    seen_model = "log_normal",
                    variables = "intercept") {
  if (is.null(data$in_sample)) {
    data$in_sample <- 1
  }
  if (is.null(data$intercept)){
    data$intercept <- 1
  }

  data <- list(
    catches = data$catch,
    seeing_catches = as.numeric(data$catch > 0),
    variables = data %>% select(variables) %>% as.matrix(),
    in_sample = data$in_sample,
    seen_dist = case_when(str_detect(seen_model, "log_normal") ~ 1,str_detect(seen_model, "gamma") ~ 2)
  )

  params <-  list(seen_betas = rep(0, ncol(data$variables)), seeing_betas=rep(0, ncol(data$variables)), log_dist_par = 0)


  obj <- MakeADFun(data = data, parameters = params,DLL = "homework_1")

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

  report <- obj$report()

  report$prob_seeing * report$seen_catch_hat

  summary_table <-
    data_frame(
      model = seen_model,
      log_like = obj$fn(),
      n_pars = length(obj$par),
      l_pred_score = ifelse(
        sum(data$in_sample == 0) == 0,
        0,
        report$oob_nll / sum(data$in_sample == 0)
      ),
      report = list(report),
      diagnostics = list(diagnostics)
    )

  return(summary_table)

}
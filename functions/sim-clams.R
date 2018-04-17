sim_clams <- function(mean_counts, sigma_site, sigma_overdisp, n_sites, n_counts){

  log_site_mean <- rnorm(n_sites, mean_counts, sigma_site)

  log_expected_count <- rnorm(n_counts * n_sites, rep(log_site_mean, each = n_counts),sigma_overdisp)

  counts <- rpois(length(log_expected_count), exp(log_expected_count))

  out <-
    data_frame(
      observation = rep(1:n_counts, n_sites),
      site = rep(1:n_sites, each = n_counts),
      counts = counts,
      site_log_mean = rep(log_site_mean, each = n_counts),
      log_expected_count = log_expected_count,
      sigma_site = sigma_site,
      sigma_overdisp = sigma_overdisp
    )


}
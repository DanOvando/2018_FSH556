#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{


// data //

DATA_VECTOR(counts);

DATA_IVECTOR(sites);

// parameters //

PARAMETER(mean_count);

PARAMETER_VECTOR(site_means);

PARAMETER_VECTOR(overdisp);

PARAMETER(log_sigma_site);

PARAMETER(log_sigma_overdisp);

// model //

Type nll;

int n_obs;

int n_sites;

n_obs = counts.size();

n_sites = site_means.size();

nll = 0;

vector<Type> log_expected_counts(n_obs);


for (int i = 0; i < n_obs; i++){

  log_expected_counts(i) = mean_count + site_means(sites(i)) + overdisp(i);

  nll -= dpois(counts(i), exp(log_expected_counts(i)), true);

  nll -= dnorm(overdisp(i), Type(0), exp(log_sigma_overdisp), true);

}

for (int s = 0; s < n_sites; s++){

  nll -= dnorm(site_means(s), Type(0), exp(log_sigma_site), true);

}


// report //

REPORT(log_expected_counts);

REPORT(site_means);

REPORT(overdisp);

REPORT(n_obs);

return nll;
}

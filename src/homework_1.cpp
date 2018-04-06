#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  //// load data ////
  DATA_VECTOR(catches); // vector of pollock catches

  DATA_VECTOR(seeing_catches);

  DATA_MATRIX(variables); // covaraiates

  DATA_IVECTOR(in_sample);

  DATA_INTEGER(seen_dist); //1 = log-normal, 2 = gamma

  //// parameters ////


  PARAMETER_VECTOR(seen_betas);

  PARAMETER_VECTOR(seeing_betas);

  PARAMETER(log_dist_par);


  //// model ////

  int n_t;

  Type nll;

  Type oob_nll;

  Type scale;


  Type dist_par = exp(log_dist_par);

  matrix<Type> seen_catch_hat  = variables * seen_betas;

  matrix<Type> seeing_catch_hat  = variables * seeing_betas;

  vector<Type> prob_seeing =  1 / (1 + exp(-seeing_catch_hat.array()));


  nll = 0;

  oob_nll = 0;

  n_t = catches.size();

  vector<Type> catch_hat(n_t);


  for (int i = 0; i < n_t; i++){

    catch_hat(i) = prob_seeing(i) * seen_catch_hat(i);

    if (in_sample(i) == 1) {

    nll -= dbinom(seeing_catches(i),Type(1),prob_seeing(i), true);

    if (catches(i) > 0){

      if (seen_dist == 1){ //lognormal

        nll -= dnorm(log(catches(i)), seen_catch_hat(i), dist_par, true);
      } // close lognormal if
      if (seen_dist == 2){

        scale = exp(seen_catch_hat(i)) * pow(dist_par,2);

        nll -= dgamma(catches(i), pow(dist_par, -2), scale, true);

      } // close gamma if
    } // close  catch if

    } // close in sample if

    if (in_sample(i) == 0){

      oob_nll -= dbinom(seeing_catches(i),Type(1),prob_seeing(i), true);

      if (catches(i) > 0){


        if (seen_dist == 1){

          oob_nll -= dnorm(log(catches(i)), seen_catch_hat(i), dist_par, true);
        } // close lnorm if
        if (seen_dist == 2){

          scale = exp(seen_catch_hat(i)) * pow(dist_par,2);

          oob_nll -= dgamma(catches(i), pow(dist_par, -2), scale, true);

        } // close gamma if

      } // close catch if

    } // close in sample if


  } // close model loop



  REPORT(dist_par);

  REPORT(catch_hat);

  REPORT(seen_betas);

  REPORT(seen_catch_hat);

  REPORT(prob_seeing);

  REPORT(oob_nll);

  ADREPORT(seen_betas);

  REPORT(seeing_betas);

  ADREPORT(seeing_betas);

  return nll;

}
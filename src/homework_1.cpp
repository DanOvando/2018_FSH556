#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{

  //// load data ////
  DATA_VECTOR(catches); // vector of pollock catches

  DATA_VECTOR(seeing_catches);

  DATA_MATRIX(variables); // covaraiates

  //// parameters ////
  //
  //

  PARAMETER_VECTOR(seen_betas);

  PARAMETER_VECTOR(seeing_betas);

  PARAMETER(log_sigma_seen);


  //// model ////

  int n;

  Type nll;

  Type sigma_seen = exp(log_sigma_seen);

  matrix<Type> seen_catch_hat  = variables * seen_betas;

  matrix<Type> seeing_catch_hat  = variables * seeing_betas;

  vector<Type> prob_seeing =  1 / (1 + exp(-seeing_catch_hat.array()));

  nll = 0;

  n = catches.size();

  for (int i = 0; i < n; i++){

    if (catches(i) > 0){

      nll -= dbinom(seeing_catches(i),Type(1),prob_seeing(i), true);

      nll -= dnorm(log(catches(i)), seen_catch_hat(i), sigma_seen, true);

    } // close if

    if (catches(i) == 0){

      nll -= dbinom(seeing_catches(i),Type(1),prob_seeing(i), true);

    } // close if

  } // close model loop


  REPORT(seen_betas);

  REPORT(seen_catch_hat);

  REPORT(prob_seeing);

  ADREPORT(seen_betas);

  REPORT(seeing_betas);

  ADREPORT(seeing_betas);

  return nll;

}
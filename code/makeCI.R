library(mvrtn)
library(MultinomialCI)

## return confidence intervals(normal based) under independence assumption
ci_indep = function(mat,mu){
  ### input
  # mat : a matrix having independent simluated observation in each row.
  # mu : expected value of generated observations.
  ### output
  # ci : 95 % confidence interval assuming independent observations.

  means = apply(mat,1,mean) 
  vars = apply(mat,1,var)
  nsim = nrow(mat) 
  n = ncol(mat)
  lower = means-(qnorm(.975)*sqrt(vars/n))
  upper = means+(qnorm(.975)*sqrt(vars/n))
  ci = cbind(lower, upper)
  
  return(ci)
}


## coverage rate under independence assumption
cov_indep = function(mat,mu){
  ### input
  # mat : a matrix having independent simluated observation in each row.
  # mu : expected value of generated observations.
  ### output
  # cov : coverage rate of 95 % confidence interval constructed under independence assumption.

  means = apply(mat,1,mean) 
  vars = apply(mat,1,var)
  nsim = nrow(mat) 
  n = ncol(mat)
  lower = means-(qnorm(.975)*sqrt(vars/n))
  upper = means+(qnorm(.975)*sqrt(vars/n))
  cov = sum(lower < mu & upper > mu) / nsim
  
  return(cov)
}

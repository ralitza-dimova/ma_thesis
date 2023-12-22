#simulation functions
library(MBESS)
library(psych)

make_sample_data <- function(pop_d, n) {
  sample1 = rnorm(n = n, mean = 0, sd = 1) #for population 1, setting up the sample d value
  sample2 = rnorm(n = n, mean = pop_d, sd = 1)
  output <- list(sample1 = sample1, sample2 = sample2 )
  return(output)
}
 
#calculation of d value 
calc_d_value <- function(data_samples){
  pop_1_var1 <- var(data_samples$sample1) #for population 1, calculate the standard deviation pooled  
  pop_1_var2 <- var(data_samples$sample2)
  pop_1_variance_pooled <- (pop_1_var1 + pop_1_var2) / 2
  pop_1_sd_pooled <- sqrt(pop_1_variance_pooled)
  pop1_sample_d_value <- (mean(data_samples$sample2) - mean(data_samples$sample1)) / pop_1_sd_pooled #for population 1, calculating the sample d value
  return(pop1_sample_d_value)
}

#calculation of a Hedges g value
hedges_g_value <- function(d_value_to_convert, n1, n2){
  df <- n1 + n2 - 2
  g_value <- d_value_to_convert*(1-(3/((4*df) - 1)))
  return (c(g_value, df))
}

#calculation of Steiger confidence intervals
calc_ci_value <- function(sample_d, n) {
  ci_info <- ci.smd(smd = sample_d, n.1 = n, n.2 = n, conf.level=.95)
  LL <- ci_info$Lower.Conf.Limit.smd
  UL <- ci_info$Upper.Conf.Limit.smd
  return (c(LL, UL, sample_d))
}

#calculation of Hedges confidence intervals
hedges_CI_function <- function(g_value, n, gamma = 0.95){
  lambda <- g_value*sqrt(n/2)
  LL_g_value <- qt(1/2 - gamma/2, df = 2*(n-1), ncp = lambda)
  UL_g_value <- qt(1/2 + gamma/2, df = 2*(n-1), ncp = lambda)
  limits <- c(c(LL_g_value, UL_g_value)/sqrt(n/2), g_value)
}

#calculate zou Lower limit
delta_d_zou <- function(ci1, ci2){
  d1 <- ci1[3]
  LL1 <- ci1[1]
  UL1 <- ci1[2]
  
  d2 <- ci2[3]
  LL2 <- ci2[1]
  UL2 <- ci2[2]
  
  difference_LL <- d1 - d2 - sqrt(((d1 - LL1)^2) + ((UL2 - d2)^2))
  difference_UL <- d1 - d2 + sqrt(((UL1 - d1)^2) + ((d2 - LL2)^2))
  return (c(difference_LL, difference_UL))
}

###PAIRED FUNCTIONS
#function for calculating d value paired
d_paired_value <- function()


#function for making data
make_rep_t_diff_data <- function(pop.d = 1, pop.cor.t1t2 = .5, n = 1000000) {
  error2= 1
  r2 = pop.cor.t1t2^2
  w2 = error2/r2 -error2
  w = sqrt(w2)
  
  diff = rnorm(n = n, mean = pop.d*w, sd = w)
  
  rerror = rnorm(n = n, mean = 0, sd = sqrt(error2))
  t1 = diff + rerror
  t2 = rerror
  
  output <- list(t1 = t1, t2 = t2)
  return(output)
}


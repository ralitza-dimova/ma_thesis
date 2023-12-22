#pull simulation functions

library(tidyverse)

source("simulation_functions.R", echo = FALSE)

set.seed(9)

#samples  
n1set <- c(25, 100)
n2set <- c(25, 100)
population_d_value <- seq(0.1, 0.9, by = 0.2)
number_of_trials <- 40000



#start time 
start_time <- Sys.time()

#logic for results table
results_table <- as.data.frame(expand.grid(population_d_value, 
                                           population_d_value, 
                                           n1set, 
                                           n2set))

names(results_table) <- c("delta_1", "delta_2", "N1", "N2") 
results_table <- results_table %>%
  mutate(delta_diff = delta_1 - delta_2) %>%
  mutate(N1b = N1)%>%
  mutate(N2b = N2)%>%
  select(delta_diff, delta_1, N1, N1b, delta_2, N2, N2b)%>%
  arrange(desc(delta_diff))


results_table$capture_percentage <- "?"
results_table$H_and_O_capture_percentage <- "?"

number_of_simulations <- dim(results_table)[1]

#start here for function  
for (current_sim in 1:number_of_simulations) {
  #print(sprintf("running simulation %f", current_sim))
  
  pop1_current_pop_d_value <- results_table$delta_1[current_sim] #for population 1, determining pop d value, sample size
  pop1_current_n1 <- results_table$N1[current_sim]
  pop1_current_n1b <- results_table$N1b[current_sim]
  
  pop2_current_pop_d_value <- results_table$delta_2[current_sim] #for population 2, determining pop d value, sample size
  pop2_current_n2 <- results_table$N2[current_sim]
  pop2_current_n2b <- results_table$N2b[current_sim]
  
  delta_diff_input <- results_table$delta_diff[current_sim] #bringing in the delta difference of pop 1 and pop 2
  
  did_it_capture_population <- rep(FALSE, number_of_trials) #setting up the capture percentage
  did_it_capture_population_hedges <- rep(FALSE, number_of_trials) #setting up the capture percentage
  
  for (current_trial in 1:number_of_trials) {
    #population 1
    pop1_samples <- make_sample_data(pop_d = pop1_current_pop_d_value, 
                                     n = pop1_current_n1)
    
    pop1_sample_d_value <- calc_d_value(pop1_samples) #calculate d value
    
    pop1_sample_g_value <- hedges_g_value(d_value_to_convert = pop1_sample_d_value, n1 = pop1_current_n1, n2 = pop1_current_n1b)
    
    pop1_ci = calc_ci_value(sample_d = pop1_sample_g_value[1],
                            n = pop1_current_n1) #steiger CI
    
    pop1_ci_info_lower <- pop1_ci[1]#steiger CI LOWER
    pop1_ci_info_upper <- pop1_ci[2]#steiger CI UPPER
    
    
    pop1_ci_H_and_O = hedges_CI_function(g_value = pop1_sample_g_value[1], n = pop1_current_n1, gamma = 0.95) #hedges CI
    
    pop1_ci_info_lower_H_and_O <- pop1_ci_H_and_O[1] #hedges CI lower
    pop1_ci_info_upper_H_and_O <- pop1_ci_H_and_O[2] #hedges CI upper
    
    #population 2 
    
    pop2_samples <- make_sample_data(pop_d = pop2_current_pop_d_value, 
                                     n = pop2_current_n2)
    
    pop2_sample_d_value <- calc_d_value(pop2_samples) #calculate d value
    
    pop2_sample_g_value <- hedges_g_value(d_value_to_convert = pop2_sample_d_value, n1 = pop2_current_n2, n2 = pop2_current_n2b) #calculate g value
    
    pop2_ci = calc_ci_value(sample_d = pop2_sample_g_value[1], 
                            n = pop2_current_n2) #Steiger CI
    
    pop2_ci_info_lower <- pop2_ci[1]#steiger CI LOWER
    pop2_ci_info_upper <- pop2_ci[2]#steiger CI UPPER
    
    
    pop2_ci_H_and_O = hedges_CI_function(g_value = pop2_sample_g_value[1], n = pop2_current_n2, gamma = 0.95) #hedges CI
    
    pop2_ci_info_lower_H_and_O <- pop2_ci_H_and_O[1] #hedges CI lower
    pop2_ci_info_upper_H_and_O <- pop2_ci_H_and_O[2] #hedges CI upper
    
    #formula for LL and UL for steiger
    confidence_interval_created <- delta_d_zou(pop1_ci, pop2_ci)
    deltaLL <- confidence_interval_created[1]
    deltaUL <- confidence_interval_created[2]
    
    #formula for LL and UL for Hedges
    confidence_interval_created_hedges <- delta_d_zou(pop1_ci_H_and_O, pop2_ci_H_and_O)
    delta_hedges_LL <- confidence_interval_created_hedges[1]
    delta_hedges_UL <- confidence_interval_created_hedges[2]
    
    if (delta_diff_input > deltaLL) { 
      if (delta_diff_input < deltaUL) {
        did_it_capture_population[current_trial] <- TRUE
      } 
    }
    if (delta_diff_input > delta_hedges_LL) { 
      if (delta_diff_input < delta_hedges_UL) {
        did_it_capture_population_hedges[current_trial] <- TRUE
      } 
    }
  }
  current_capture_percent <- sum(did_it_capture_population) /number_of_trials * 100
  #print(current_capture_percent)
  results_table$capture_percentage[current_sim] <- current_capture_percent
  
  current_capture_percent_hedges <- sum(did_it_capture_population_hedges) /number_of_trials * 100
  #print(current_capture_percent_hedges)
  results_table$H_and_O_capture_percentage[current_sim] <- current_capture_percent_hedges
}

#end time
end_time <- Sys.time()
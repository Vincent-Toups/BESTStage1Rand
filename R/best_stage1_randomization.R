#' Main Stage 1 Randomization
#'
#' This function performs stage 1 randomization for a clinical trial, assigning treatments to subjects
#' based on eligibility criteria and treatment-specific covariates. It iterates through all subjects in the input data frame,
#' calculates treatment probabilities, and assigns treatments while updating covariate tables.
#'
#' @param rand_input_df A data frame containing subject-specific information, including eligibility criteria
#' and treatment-specific covariates.
#' @param trts A character vector of treatments available for randomization (default: \code{c("act", "dulox", "ebem", "esc")}).
#' @param trt_min_covs A character vector of covariates to be minimized during randomization (default: \code{c("pain", "opioid", "depress", "pheno")}).
#' @param debug_num An integer for debugging purposes. If greater than zero, the function stops processing after this subject ID.
#' Default is -1, meaning no debugging limit is applied.
#'
#' @return A list containing:
#' \item{assigned_stage1_trts}{A data frame with assigned treatments for each subject, including treatment labels.}
#' \item{trt_sample_size}{A summary table showing the sample size and covariate distribution for each treatment.}
#' \item{cumcov}{A cumulative covariate table, showing the distribution of covariates by treatment.}
#' \item{bacpacrand}{A data frame with detailed information for each subject, including probabilities and assigned treatments.}
#'
#' @details
#' This function leverages treatment-specific eligibility and randomization probabilities to ensure
#' balanced covariate distributions across treatments while assigning treatments to subjects.
#' Covariate minimization is achieved through treatment selection optimization, which considers both eligibility and minimizing
#' imbalance in covariates.
#'
#' @examples
#' # Example usage:
#' rand_input_df <- data.frame(
#'   subjectid = 1:100,
#'   act_elig_w0 = sample(c(0, 1), 100, replace = TRUE),
#'   dulox_elig_w0 = sample(c(0, 1), 100, replace = TRUE),
#'   ebem_elig_w0 = sample(c(0, 1), 100, replace = TRUE),
#'   esc_elig_w0 = sample(c(0, 1), 100, replace = TRUE),
#'   rpt1a = rep("1/2", 100),
#'   rpt2a = rep("1/2", 100),
#'   rpt3a = rep("1/2", 100),
#'   rpt4a = rep("1/2", 100),
#'   rpt5a = runif(100),
#'   rpt6a = runif(100),
#'   rpt7a = runif(100),
#'   rpt8a = runif(100),
#'   pain = sample(c("yes", "no"), 100, replace = TRUE),
#'   opioid = sample(c("yes", "no"), 100, replace = TRUE),
#'   depress = sample(c("yes", "no"), 100, replace = TRUE),
#'   pheno = sample(c("yes", "no"), 100, replace = TRUE)
#' )
#'
#' result <- main_stage1_randomization(rand_input_df)
#'
#' @export
main_stage1_randomization <- function(rand_input_df, 
                                      trts = c("act", "dulox", "ebem", "esc"),
                                      trt_min_covs = c("pain", "opioid", "depress", "pheno"), 
                                      debug_num = -1){
  
  trt_min_cov_tables <- init_minimization_objects(minimization_options = trts, 
                                                  trt_min_covs)
  
  
  #Initialize assigned treatment df.
  assigned_stage1_trts <- rand_input_df %>% 
    transmute(subjectid)
  
  #iterate through all rows of input df (subjects)
  for(i in 1:nrow(rand_input_df)){
    
    #Get data for the subject
    subject_df <- rand_input_df %>% 
      slice(i) %>% 
      mutate(trt_elig_vec = list(which(c(act_elig_w0, dulox_elig_w0, ebem_elig_w0, esc_elig_w0) == 1))) 
    
    .subject_id <- subject_df$subjectid
    
    if(debug_num != -1){
      if(.subject_id > debug_num){
        break
      }
    }
    
    
    if(.subject_id == debug_num){
      glimpse(subject_df)
    }
    
    #Get Theta vector for treatment assignment
    theta_vec <- c(subject_df$rpt1a, subject_df$rpt2a,
                   subject_df$rpt3a, subject_df$rpt4a) %>% 
      #Split on /
      str_split(pattern ="/") %>% 
      #Convert to numeric
      lapply(as.numeric) %>% 
      #Compute decimal.
      sapply(FUN = function(a){a[1]/a[2]})
    
    #Get omega vector for treatment assignment
    omega_vec <- c(subject_df$rpt5a, subject_df$rpt6a,
                   subject_df$rpt7a, subject_df$rpt8a)
    
    
    #Assign treatment
    l_subject_i_trt <- assign_trt(.subject_df = subject_df,
                                  .min_cov_tables = trt_min_cov_tables,
                                  .trt_min_cov_list = trt_min_covs,
                                  .assigned_trt_df = assigned_stage1_trts, 
                                  .subject_num = .subject_id, 
                                  .trt_theta_vec = theta_vec,
                                  .trt_omega_vec = omega_vec, 
                                  .debug_num = debug_num)
    
    #Update trt minimization covariate tables
    
    trt_min_cov_tables <- add_covs_to_min_table(.subject_df = subject_df, 
                                                .min_cov_tables = trt_min_cov_tables,
                                                .assignment = l_subject_i_trt[['assigned_option']], 
                                                .cov_list = trt_min_covs)
    
    #Create or append the cumcov table:
    if(i == 1){
      cumcov <- trt_min_cov_tables %>% 
        do.call(cbind, args = .) %>% 
        as_tibble(rownames = "treatment") %>% 
        mutate(subjectid = subject_df$subjectid)
    } else{
      cumcov <- trt_min_cov_tables %>% 
        do.call(cbind, args = .)  %>% 
        as_tibble(rownames = "treatment") %>% 
        mutate(subjectid = subject_df$subjectid) %>% 
        bind_rows(cumcov, .)
    }
    
    #Create or append bacpacrand
    if(i == 1){
      
      bacpacrand <- subject_df %>%
        select(subjectid,
               contains("elig"), 
               -trt_elig_vec,
               u1, u2, u3) %>%
        mutate(aug_prob = NA_character_,
               switch_prob = NA_character_,
               act_prob = l_subject_i_trt$option_probs[1],
               dulox_prob = l_subject_i_trt$option_probs[2],
               ebem_prob = l_subject_i_trt$option_probs[3],
               esc_prob = l_subject_i_trt$option_probs[4],
               tx = l_subject_i_trt$assigned_option,
               action = NA_character_,
               tx_min = tx)
    }else{
      bacpacrand <- subject_df %>%
        select(subjectid,
               contains("elig"), 
               -trt_elig_vec,
               u1, u2, u3) %>%
        mutate(aug_prob = NA_character_,
               switch_prob = NA_character_,
               act_prob = l_subject_i_trt$option_probs[1],
               dulox_prob = l_subject_i_trt$option_probs[2],
               ebem_prob = l_subject_i_trt$option_probs[3],
               esc_prob = l_subject_i_trt$option_probs[4],
               tx = l_subject_i_trt$assigned_option,
               action = NA_character_,
               tx_min = l_subject_i_trt$g) %>%
        bind_rows(bacpacrand, .)
    }
    
    
    assigned_stage1_trts[i,"trt_num"] <- l_subject_i_trt[['assigned_option']]
    
  }
  
  #Format Output
  assigned_stage1_trts <- assigned_stage1_trts %>% mutate(
    trt = case_when(
      trt_num == 1 ~ "Acceptance and Commitment Therapy (ACT)", 
      trt_num == 2 ~ "Duloxetine",
      trt_num == 3 ~ "Evidence-Based Exercise and Manual Therapy (EBEM)", 
      trt_num == 4 ~ "Enhanced Self-Care (ESC)"
    )
  )
  
  #Treatment Sample Size
  trt_sample_size <- trt_min_cov_tables %>% 
    do.call(what = cbind, args = .) %>% 
    as_tibble(rownames = "treatment") %>% 
    mutate(tx_n = pheno_no + pheno_yes) %>% 
    relocate(tx_n, .after = treatment) %>% 
    rename(pain_gte5 = pain_yes,
           pain_lt5 = pain_no)
  
  #Cumcov
  cumcov <- cumcov %>% 
    relocate(subjectid, .before = treatment)%>% 
    rename(pain_gte5 = pain_yes,
           pain_lt5 = pain_no) %>% 
    select(subjectid, treatment, pheno_yes, 
           pheno_no,
           depress_yes,
           depress_no, pain_gte5, 
           pain_lt5, 
           opioid_yes, opioid_no) %>% 
    rowwise() %>% 
    mutate(tx_n = pheno_yes + pheno_no) %>% 
    relocate(tx_n, .after = treatment)
  
  
  #right now, just output assigned stage1 treatments
  list(assigned_stage1_trts, 
       trt_sample_size, 
       cumcov, 
       bacpacrand)
}


#Functions which comprise the main function in stage1_randomization_main.R
assign_trt <- function(.subject_df, 
                       .min_cov_tables, 
                       .assigned_trt_df,
                       .trt_min_cov_list,
                       .trt_theta_vec,
                       .trt_omega_vec,
                       .subject_num = -1, 
                       .debug_num){

  
  l_s1t <- vector("list", length = 3)
  names(l_s1t) <- c("assigned_option", "option_probs", "g")
  
  #Get vector of eligible treatments
  elig_trts <- .subject_df$trt_elig_vec[[1]]
  n_elig_trts <- length(elig_trts)
  
  if(.debug_num == .subject_num){
    cat(sprintf("\nSubject %s has %s eligible treatments.\n\n", .subject_num,
            n_elig_trts))
  }

  
  if(.subject_num == 1){
    
    elig_intervals <- c(0,(1:n_elig_trts)/n_elig_trts)
    picked_trt_index <- findInterval(.subject_df$u1,
                                     elig_intervals)
    if(.debug_num == .subject_num){
      cat(sprintf("Subject %s was the first subject.\nUsing equal probabilities for eligible treatments trt # %s was selected.\n\n", 
              .debug_num, elig_trts[picked_trt_index]))
    }
    
    
    l_s1t[["assigned_option"]] <- elig_trts[picked_trt_index]
    
    .empty_prob_vec <- rep(0,4)
    .empty_prob_vec[elig_trts] <- 1/n_elig_trts
    
    l_s1t[["option_probs"]] <- .empty_prob_vec
    
    l_s1t[['g']] <- NA_real_
    
    return(l_s1t)
  } else if(.subject_num > 1){
    
    #Get minimizing treatment.
    l_s1t[['g']] <- get_minimizer(.subject_df = .subject_df, 
                       .min_cov_tables = .min_cov_tables, 
                       .cov_list = .trt_min_cov_list, 
                       .elig_trts = elig_trts, 
                       .subject_num = .subject_num, 
                       .debug_num = .debug_num)
    
    #Pick from minimizing treatment
    
    #Need to define an argument .rpt_vec containing the correct
    ## RPT probabilities from the minimizer. This should be defined
    ## In such a way that it's flexible for Stage 2. 
    
    #Returns picked treatment for subject_num.
    l_s1t[['assigned_option']] <- randomization_and_minimization(.subject_df = .subject_df,
                                   .minimizing_trt = l_s1t[['g']], 
                                   .theta_vec = .trt_theta_vec, 
                                   .omega_vec = .trt_omega_vec, 
                                   .elig_trts = elig_trts, 
                                   .subject_num = .subject_num, 
                                   .debug_num = .debug_num)
  
    l_s1t[['option_probs']] <- rep(0,4)
    l_s1t[['option_probs']][l_s1t[['g']]] <- .trt_theta_vec[l_s1t[['g']]]
    
    .elig_trts_no_min <- elig_trts[-which(elig_trts == l_s1t[['g']])]
    
    .elig_omegas <- .trt_omega_vec[.elig_trts_no_min]
    l_s1t[['option_probs']][.elig_trts_no_min] <- (.elig_omegas/sum(.elig_omegas))*(1 - .trt_theta_vec[l_s1t[['g']]])
    
    return(l_s1t)
  }
}


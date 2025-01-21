##################################
## R functions for minimization  #
##################################

init_minimization_objects <- function(minimization_options, covs){
  cov_m_list <- vector("list", length = length(covs))
  
  for(i in 1:length(covs)){
    
    colnames <- paste0(covs[[i]], c("_no", "_yes"))
    
    cov_m_list[[i]] <- matrix(0, nrow = length(minimization_options), 
                              ncol = 2, 
                              dimnames =list(minimization_options, 
                                             colnames))
    
  }
  names(cov_m_list) <- covs
  cov_m_list
}

add_covs_to_min_table <- function(.subject_df, 
                                  .min_cov_tables, 
                                  .assignment, 
                                  .cov_list){
  
  for(c in .cov_list){
    cov_col_index <- .subject_df %>% pull(c) + 1
    cov_row_index <- .assignment 
    .min_cov_tables[[c]][cov_row_index, 
                         cov_col_index] = .min_cov_tables[[c]][cov_row_index, 
                                                               cov_col_index] + 1
    
  }
  .min_cov_tables
}

get_minimizer <- function(.subject_df, 
                          .min_cov_tables, 
                          .cov_list,
                          .elig_trts, 
                          .debug_num, 
                          .subject_num){

  #Get n_treatment
  n_treatment <- nrow(.min_cov_tables[[1]])
  
  d_list <- vector("list", length = n_treatment)
  
  #For each treatment, we will compute a vector of n_covariate
  ## d-values
  for(hypothetical_treatment in 1:n_treatment){
    
    #Assign the individual to the hypothetical treatment to see
    ## how it changes sample-size

    
    change_matrices <- add_covs_to_min_table(.subject_df = .subject_df, 
                          .min_cov_tables = .min_cov_tables,
                          .assignment = hypothetical_treatment, 
                          .cov_list = .cov_list)
    
    d_list[[hypothetical_treatment]] <- get_d_vector(
      .change_matrices = change_matrices, 
      .subject_df = .subject_df, 
      .cov_list = .cov_list)
    
  }
  
  if(.debug_num == .subject_num){
    cat("Printing .min_cov_tables before assigning subject.\n\n")
    print(.min_cov_tables)
    cat("Printing list of d-vectors in G computation...\n")
    print(d_list)
    cat("\n\n")
  }
  
  g_values <- sapply(d_list, sum)
  
  if(.debug_num == .subject_num){
    cat("Printing g_values.\n\n")
    print(g_values)
    cat("\n\n")
  }
  
  elig_g_values <- g_values[.elig_trts]

  if(.debug_num == .subject_num){
    cat("Printing elig_g_values. \n\n")
    print(elig_g_values)
    cat("\n\n")
  }
  
  #Get minimizing treatments
  elig_min_idx <- which(elig_g_values == min(elig_g_values))
  elig_min_trts <- .elig_trts[elig_min_idx]
  
  if(.debug_num == .subject_num){
    cat("Printing elig_min_trts \n\n")
    print(elig_min_trts)
    cat("\n\n")
  }

  #If there is only one minimizer, return it as the minimizing treatment
  if(length(elig_min_trts) == 1){
    
    if(.debug_num == .subject_num){
      cat(sprintf("Only one minimizer, %s,...returning it as minimizing treatment.\n\n", 
                  elig_min_trts))
    }
    
    return(elig_min_trts)
  } 
  
  #Otherwise, randomly choose one of the possible minimizers.
  else{
    
    n_elig_minimizers <- length(elig_min_trts)
    
    if(.debug_num == .subject_num){
      cat(sprintf("Picking from %s possible minimizers using U1 value of %s.\n\n", n_elig_minimizers, 
                  .subject_df$u1))
    }
    
    elig_intervals <- c(0,1:n_elig_minimizers/n_elig_minimizers)
    picked_trt_index <- findInterval(.subject_df$u1,
                                     elig_intervals)
    if(.debug_num == .subject_num){
      cat(sprintf("%s picked as minimizer.\n\n", elig_min_trts[picked_trt_index]))
    }
    
    return(elig_min_trts[picked_trt_index])
  }
}

get_d_vector <- function(.change_matrices, 
                         .subject_df, 
                         .cov_list){
  
  d <- rep(-99, length(.cov_list))

  for(i in 1:length(.cov_list)){
    mm <-.change_matrices[[.cov_list[[i]]]][, (.subject_df %>% pull(.cov_list[[i]]) + 1)]
    d[i] <- max(mm) - min(mm)
  }
  d

}

randomization_and_minimization <- function(.subject_df, 
                               .minimizing_trt, 
                               .theta_vec, 
                               .omega_vec, 
                               .elig_trts, 
                               .subject_num, 
                               .debug_num){
  
  #Because of checks in get_minimizer, .minimizing_trt
  ## should always be in .elig_trts
  if(!(.minimizing_trt %in% .elig_trts)){
    cat("Error in randomization_and_minimization function: minimizing treatment not in eligible treatments.")
  }
  
  if(.debug_num == .subject_num){
    cat(sprintf("Using U2 value of %s to pick minimizer with probability %s.\n\n", 
                .subject_df$u2, .theta_vec[.minimizing_trt]))
  }

  
  #Use U2 to see if we pick the minimizer
  if(.subject_df$u2 <= .theta_vec[.minimizing_trt]){
    #Return .minimizing_trt
    .minimizing_trt
  } else{
    
    
    if(.debug_num == .subject_num){
      cat("Minimizer not picked.\n\n")
    }
    
    
    #Use U3 to see which other eligible treatment we pick
    #Get eligible treatments which aren't the minimizer
    elig_trts_not_minimizer <- .elig_trts[-which(.elig_trts == .minimizing_trt)]
    
    if(.debug_num == .subject_num){
      cat("Here are the eligible treatments which aren't the minimizer")
      print(elig_trts_not_minimizer)
    }
    
    #Compute probabilities from relative weights
    omega_probs <- .omega_vec[elig_trts_not_minimizer]/sum(.omega_vec[elig_trts_not_minimizer])
    omega_intervals <- c(0,cumsum(omega_probs))

    if(.debug_num == .subject_num){
      cat(sprintf("U3 value of %s was used to pick from omega intervals of %s", .subject_df$u3,
          str_c(omega_intervals, collapse = ",")))
    }
    
    #Select treatments using U3.
    picked_trt_index <- findInterval(.subject_df$u3,
                                     omega_intervals)
    elig_trts_not_minimizer[picked_trt_index]
  }
}

#' Create Minimizer Data Frame
#'
#' This function generates a data frame (tibble) used for minimization in randomization processes.
#' It initializes uniform random variables for randomization and populates specified probability 
#' and weight variables with default or user-defined values.
#'
#' @param n_rows An integer specifying the number of rows in the output data frame.
#' @param prob_vars A character vector of column names for probability variables (default: \code{paste0('rpt', 1:4, "a")}).
#' @param prob_value A character string representing the default value for the probability variables 
#' (default: \code{"2/3"}).
#' @param weight_vars A character vector of column names for weight variables (default: \code{paste0('rpt', 5:8, "a")}).
#' @param weight_value A numeric value representing the default value for the weight variables 
#' (default: \code{1}).
#'
#' @return A tibble with the following structure:
#' \item{u1, u2, u3}{Uniform random variables sampled from \code{runif(0, 1)}.}
#' \item{prob_vars}{Columns named according to \code{prob_vars}, filled with \code{prob_value}.}
#' \item{weight_vars}{Columns named according to \code{weight_vars}, filled with \code{weight_value}.}
#'
#' @details
#' This function is designed for use in treatment randomization or minimization workflows. It creates a
#' standardized data frame structure with randomization variables and customizable probability and weight
#' columns.
#'
#' @examples
#' # Create a minimizer data frame with default settings
#' minimizer_df <- create_minimizer_df(100)
#'
#' # Create a minimizer data frame with custom probability and weight variables
#' minimizer_df <- create_minimizer_df(
#'   n_rows = 50, 
#'   prob_vars = c("prob1", "prob2"), 
#'   prob_value = "1/2", 
#'   weight_vars = c("weight1", "weight2"), 
#'   weight_value = 0.5
#' )
#'
#' @export
create_minimizer_df <- function(n_rows, 
                                prob_vars = paste0('rpt', 1:4, "a"), 
                                prob_value = "2/3", 
                                weight_vars = paste0('rpt', 5:8, "a"), 
                                weight_value = 1) {
  # Initialize the tibble with uniform random variables
  df <- tibble(
    u1 = runif(n_rows, 0, 1),
    u2 = runif(n_rows, 0, 1),
    u3 = runif(n_rows, 0, 1)
  )
  
  # Add probability variables as character values
  df <- df %>%
    mutate(!!!setNames(nm = prob_vars, 
                       object = rep(prob_value, length(prob_vars))))
  
  # Add weight variables as numeric values
  df <- df %>%
    mutate(!!!setNames(nm = weight_vars, 
                       object = rep(weight_value, length(weight_vars))))
  
  return(df)
}

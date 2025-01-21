library(testthat)
library(tidyverse)
library(janitor)

catf <- function(...){
    cat(sprintf(...));
}

test_that("Stage 1 randomization produces correct results", {
  
  # Set a random seed for reproducibility
  set.seed(12345)
  
  # Load the input test data
  input_path <- "randomization_stage1_test_input.csv"
  catf("wd: %s\n", getwd())
  catf("Input data %s\n", input_path)
  in_df <- read_csv(input_path)
  catf("Names of in_df %s", paste(names(in_df),collapse=", "));
  
  # Load the expected results
  expected_path <- "stage1_randomization_result_sim1.csv"
  expected_df <- read_csv(expected_path)
  
  # Create the minimizer variables using the new function
  minimizer_df <- create_minimizer_df(n_rows = nrow(in_df))
  
  # Combine minimizer variables with the input data
  ready_df <- in_df %>%
    clean_names() %>%
    cbind(minimizer_df)
  
  # Run the main stage 1 randomization function
  results_list <- main_stage1_randomization(ready_df)
  
  # Process results for comparison
  trt_df <- left_join(ready_df, results_list[[1]], 
                      by = 'subjectid') %>%
    mutate(stage1trt = case_when(
      trt == "Acceptance and Commitment Therapy (ACT)" ~ "act",
      trt == "Duloxetine" ~ "dulox",
      trt == "Enhanced Self-Care (ESC)" ~ "esc",
      trt == "Evidence-Based Exercise and Manual Therapy (EBEM)" ~ "ebem"
    )) %>% 
    rename(stage1trt_num = trt_num) %>% 
    select(-trt) %>% as_tibble();
  

  write_csv(trt_df, "fixed.csv");
  expect_equal(trt_df, expected_df)
})

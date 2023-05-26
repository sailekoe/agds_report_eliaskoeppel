get_selection <- function(data){
  data  |>  
    
    # select only the variables we are interested in
    dplyr::select(TIMESTAMP,
                  GPP_NT_VUT_REF,    # the target
                  ends_with("_QC"),  # quality control info
                  ends_with("_F"),   # includes all all meteorological covariates
                  -contains("JSB")   # weird useless variable
    ) |>
    
    # convert to a nice date object
    dplyr::mutate(TIMESTAMP = ymd(TIMESTAMP)) |>
    
    # set all -9999 to NA
    mutate(across(where(is.numeric), ~na_if(., -9999))) |>  
    
    # retain only data based on >=80% good-quality measurements
    # overwrite bad data with NA (not dropping rows)
    dplyr::mutate(GPP_NT_VUT_REF = ifelse(NEE_VUT_REF_QC < 0.8, NA, GPP_NT_VUT_REF),
                  TA_F           = ifelse(TA_F_QC        < 0.8, NA, TA_F),
                  SW_IN_F        = ifelse(SW_IN_F_QC     < 0.8, NA, SW_IN_F),
                  LW_IN_F        = ifelse(LW_IN_F_QC     < 0.8, NA, LW_IN_F),
                  VPD_F          = ifelse(VPD_F_QC       < 0.8, NA, VPD_F),
                  PA_F           = ifelse(PA_F_QC        < 0.8, NA, PA_F),
                  P_F            = ifelse(P_F_QC         < 0.8, NA, P_F),
                  WS_F           = ifelse(WS_F_QC        < 0.8, NA, WS_F)) |> 
    
    # drop QC variables (no longer needed)
    dplyr::select(-ends_with("_QC"))
  
}




get_selection_laegern <- function(data){
  data  |>  
    
    # select only the variables we are interested in
    dplyr::select(GPP_NT_VUT_REF,    # the target
                  ends_with("_QC"),  # quality control info
                  ends_with("_F"),   # includes all all meteorological covariates
                  -contains("JSB")   # weird useless variable
    ) |>
    # set all -9999 to NA
    mutate(across(where(is.numeric), ~na_if(., -9999))) |>  
    
    # retain only data based on >=80% good-quality measurements
    # overwrite bad data with NA (not dropping rows)
    dplyr::mutate(GPP_NT_VUT_REF = ifelse(NEE_VUT_REF_QC < 0.8, NA, GPP_NT_VUT_REF),
                  TA_F           = ifelse(TA_F_QC        < 0.8, NA, TA_F),
                  SW_IN_F        = ifelse(SW_IN_F_QC     < 0.8, NA, SW_IN_F),
                  LW_IN_F        = ifelse(LW_IN_F_QC     < 0.8, NA, LW_IN_F),
                  VPD_F          = ifelse(VPD_F_QC       < 0.8, NA, VPD_F),
                  PA_F           = ifelse(PA_F_QC        < 0.8, NA, PA_F),
                  P_F            = ifelse(P_F_QC         < 0.8, NA, P_F),
                  WS_F           = ifelse(WS_F_QC        < 0.8, NA, WS_F)) |> 
    
    # drop QC variables (no longer needed)
    dplyr::select(-ends_with("_QC"))
  
}






# function for pre-procession (use all variables but LW_IN_F)
get_pp <- function(data){
  recipes::recipe(GPP_NT_VUT_REF ~ SW_IN_F + VPD_F + TA_F, 
                  data = data|> drop_na()) |> 
    recipes::step_BoxCox(all_predictors()) |> 
    recipes::step_center(all_numeric(), -all_outcomes()) |>
    recipes::step_scale(all_numeric(), -all_outcomes())}






# model evaluation 
eval_model <- function(mod, df_test){
  
  # add predictions to the data frames

  
  
  df_test <- df_test |> 
    drop_na()
  df_test$fitted <- stats::predict(mod, newdata = df_test)
  
  # get metrics tables

  
  metrics_test <- df_test |> 
    yardstick::metrics(GPP_NT_VUT_REF, fitted)
  
  # extract values from metrics tables

  rmse_test <- metrics_test |> 
    filter(.metric == "rmse") |> 
    pull(.estimate)
  rsq_test <- metrics_test |> 
    filter(.metric == "rsq") |> 
    pull(.estimate)
  
  # visualise as a scatterplot
  # adding information of metrics as sub-titles
 
  
  plot_2 <- ggplot(data = df_test, aes(GPP_NT_VUT_REF, fitted)) +
    geom_point(alpha = 0.3) +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted") +
    labs(subtitle = bquote( italic(R)^2 == .(format(rsq_test, digits = 2)) ~~
                              RMSE == .(format(rmse_test, digits = 3))),
         title = "Test set") +
    theme_classic()
  

  
  return(plot_2)
}




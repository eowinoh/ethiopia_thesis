#########################################################################################
#########################################################################################
#########################################################################################
##############
##############Assess determinants of child nutritional status (LAZ, WAZ, WLZ)
##############
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

pacman::p_load(
   tidyverse,
   finalfit,
   mice,
   mitools,
   sandwich,
   DHARMa,
   magick,ResourceSelection,
   glmmTMB, gt
)

###==Anthro predictors
anthro_vars_model <- c(
   "cc1_sex",
   "cm3_diarhoea",
   "cf12_exclusive",
   "cf16_solid_food",
   "hc2_age",
   "education_level",
   "occupation",
   "minimum_mdd",
   "muac_mother",
   "weight_mother",
   "hc3_family_size",
   "hc7_private_agric_land",
   "hc9_gro_crops",
   "wash_toilets",
   "wash_drinking_water",
   "wash_handwashing_place",
   "usually_fast",
   "age_cat",
   "urban_rural",
   "hfias_category_insecure"
)


#######Anthro Binary outcomes
fxn_uni_model <- function(outcome, predictors, imputed_df){
   ##loop through all imputed datasets, fit the model, and pool the results using Rubin's rules
   all_df <- mice::complete(imputed_df, action = "all")
   
   output <- lapply(predictors, function(v){
      
      fits <- lapply(all_df, function(dat){
         
         glm(
            reformulate(v, response = outcome),
            family = binomial(),
            data = dat
         )
         
      })
      
      pooled <- mitools::MIcombine(
         results = lapply(fits, coef),
         variances = lapply(fits, function(fit){
            sandwich::vcovHC(fit, type = "HC3")
         })
      )
      s_pooled <- summary(pooled)
      tibble(
         variable  = v,
         term      = rownames(s_pooled),
         estimate  = exp(s_pooled[, "results"]),
         conf.low  = exp(s_pooled[, "results"] - 1.96 * s_pooled[, "se"]),
         conf.high = exp(s_pooled[, "results"] + 1.96 * s_pooled[, "se"]),
         p.value   = 2 * pnorm(-abs(s_pooled[, "results"] / s_pooled[, "se"]))
      ) %>%
         mutate(
            across(
               c(estimate, conf.low, conf.high, p.value),
               ~ round(.x, 4)
            )
         ) %>%
         filter(term != "(Intercept)")
      
   })
   
   bind_rows(output)
}


uni_wasting_binary<- fxn_uni_model("wasting",anthro_vars_model, imp_df_updated)   
uni_underweight_binary<- fxn_uni_model("underweight",anthro_vars_model, imp_df_updated)    
uni_stunting_binary<- fxn_uni_model("stunting",anthro_vars_model, imp_df_updated)  
####===Multivariable

fxn_multi_model <- function(outcome, imputed_df){
   all_df <- mice::complete(imputed_df, action = "all")
   fits <- lapply(all_df, function(dat){
      
      glm(
         as.formula(
            paste(
               outcome,
               "~ cc1_sex +
          cm3_diarhoea +
          cf12_exclusive +
          cf16_solid_food +
          hc2_age +
          education_level +
          occupation +
          minimum_mdd +
          muac_mother +
          weight_mother +
          hc3_family_size +
          hc7_private_agric_land +
          hc9_gro_crops +
          wash_toilets +
          wash_drinking_water +
          wash_handwashing_place +
          usually_fast +
          age_cat +
          urban_rural +
          hfias_category_insecure"
            )
         ),
         family = binomial(),
         data = dat
      )
      
   })
   
   pooled <- mitools::MIcombine(
      results = lapply(fits, coef),
      variances = lapply(fits, function(fit){
         sandwich::vcovHC(fit, type = "HC3")
      })
   )
   
   s_pooled <- summary(pooled)
   
   results_tbl <- tibble(
      term      = rownames(s_pooled),
      estimate  = exp(s_pooled[, "results"]),
      conf.low  = exp(s_pooled[, "results"] - 1.96 * s_pooled[, "se"]),
      conf.high = exp(s_pooled[, "results"] + 1.96 * s_pooled[, "se"]),
      p.value   = 2 * pnorm(-abs(s_pooled[, "results"] / s_pooled[, "se"]))
   ) %>%
      mutate(
         across(c(estimate, conf.low, conf.high, p.value),
                ~ round(.x, 4))
      ) %>%
      filter(term != "(Intercept)")
   
   list(
      pooled_results = results_tbl,
      fits = fits
   )
   
}

multi_wasting_binary<- fxn_multi_model("wasting", imp_df_updated)   
multi_underweight_binary<- fxn_multi_model("underweight", imp_df_updated)    
multi_stunting_binary<- fxn_multi_model("stunting", imp_df_updated) 

######################===MDD model assessemt
#check convergence
sapply(multi_wasting_binary$fits, function(x) x$converged)
sapply(multi_underweight_binary$fits, function(x) x$converged)
sapply(multi_stunting_binary$fits, function(x) x$converged)

#######======DHARMa residual diagnostics plots(Anthro)

set.seed(202506)
models_list <- list(
   wasting     = multi_wasting_binary$fits,
   underweight = multi_underweight_binary$fits,
   stunting    = multi_stunting_binary$fits
)

for (model_name in names(models_list)) {
   
   fits <- models_list[[model_name]]
   
   ## randomly select 4 imputations
   selected <- sample(seq_along(fits), 4)
   
   ## create DHARMa plots for selected imputations
   for (j in seq_along(selected)) {
      
      i <- selected[j]
      
      sim_res <- simulateResiduals(fits[[i]])
      
      png(
         filename = sprintf("8. Temp/%s_bin_%02d.png", model_name, j),
         width = 1600,
         height = 1200,
         res = 150
      )
      
      plot(sim_res, main = "")
      title(main = paste("Imputation", i), line = 1)
      
      dev.off()
   }
   
   ## read the 4 images
   files <- sprintf(
      "8. Temp/%s_bin_%02d.png",
      model_name,
      1:4
   )
   imgs <- image_read(files)
   ## combine into 2 × 2 figure
   combined_plot <- image_montage(
      imgs,
      tile = "2x2",
      geometry = "1600x1200+10+10"
   )
   ## save combined plot
   image_write(
      combined_plot,
      path = sprintf(
         "5. Plots/DHARMa_sample4_%s_bin.png",
         model_name
      ),
      format = "png"
   )
}
######==check hosmer-lemeshow goodness of fit test for each imputed model
hl_results_waste_bin <- lapply(multi_wasting_binary$fits, function(m) {
   hoslem.test(m$y, fitted(m))
})

## check p-values for all the  models

hl_test_all <- function(mi_object) {
   lapply(mi_object$fits, function(m) {
      hoslem.test(m$y, fitted(m))
   })
}

hl_results_all <- list(
   wasting = hl_test_all(multi_wasting_binary),
   underweight = hl_test_all(multi_underweight_binary),
   stunting = hl_test_all(multi_stunting_binary)
)

hl_tables_all <- lapply(names(hl_results_all), function(model_name) {
   hl_results <- hl_results_all[[model_name]]
   
   data.frame(
      Model = model_name,
      Imputation = seq_along(hl_results),
      Chi_square = sapply(hl_results, function(x) unname(x$statistic)),
      DF = sapply(hl_results, function(x) unname(x$parameter)),
      P_value = sapply(hl_results, function(x) x$p.value)
   )
})
names(hl_tables_all) <- names(hl_results_all)

# gtsave(
#    hl_tables_all$wasting %>%
#       gt() %>%
#       fmt_number(
#          columns = c(Chi_square, DF, P_value),
#          decimals = 4
#       ),
#    filename = "4. Tables/hl_wasting.tex"
# )
# 
# ## underweight
# gtsave(
#    hl_tables_all$underweight %>%
#       gt() %>%
#       fmt_number(
#          columns = c(Chi_square, DF, P_value),
#          decimals = 4
#       ),
#    filename = "4. Tables/hl_underweight.tex"
# )
# 
# ## stunting
# gtsave(
#    hl_tables_all$stunting %>%
#       gt() %>%
#       fmt_number(
#          columns = c(Chi_square, DF, P_value),
#          decimals = 3
#       ),
#    filename = "4. Tables/hl_stunting.tex"
# )
# 
# 


####===check VIF
lapply(multi_wasting_binary$fits, vif)
vif_results <- lapply(seq_along(multi_wasting_binary$fits), function(i){
   
   v <- car::vif(multi_wasting_binary$fits[[i]])
   
   data.frame(
      term = rownames(v),
      gvif = v[, "GVIF^(1/(2*Df))"],
      imputation = i
   )
})

vif_results <- do.call(rbind, vif_results)

vif_results

##############################====== OLS

fxn_uni_model_ols <- function(outcome, predictors, imputed_df){
   ## Loop through all imputed datasets, fit the model, and pool the results using Rubin's rules
   all_df <- mice::complete(imputed_df, action = "all")
   
   output <- lapply(predictors, function(v){
      
      fits <- lapply(all_df, function(dat){
         lm(
            reformulate(v, response = outcome),
            data = dat
         )
      })
      
      pooled <- mitools::MIcombine(
         results = lapply(fits, coef),
         variances = lapply(fits, function(fit){
            sandwich::vcovHC(fit, type = "HC3")
         })
      )
      
      s_pooled <- summary(pooled)
      
      tibble::tibble(
         variable  = v,
         term      = rownames(s_pooled),
         estimate  = s_pooled[, "results"],
         conf.low  = s_pooled[, "results"] - 1.96 * s_pooled[, "se"],
         conf.high = s_pooled[, "results"] + 1.96 * s_pooled[, "se"],
         p.value   = 2 * pnorm(-abs(s_pooled[, "results"] / s_pooled[, "se"]))
      ) %>%
         dplyr::mutate(
            across(
               c(estimate, conf.low, conf.high, p.value),
               ~ round(.x, 4)
            )
         ) %>%
         dplyr::filter(term != "(Intercept)")
   })
   
   dplyr::bind_rows(output)
}


uni_WAZ_ols<- fxn_uni_model_ols("zwei",anthro_vars_model, imp_df_updated)   
uni_WLZ_ols<- fxn_uni_model_ols("zwfl",anthro_vars_model, imp_df_updated)    
uni_LAZ_ols<- fxn_uni_model_ols("zlen",anthro_vars_model, imp_df_updated) 

######=========Multivariable OLS
fxn_multi_model_ols <- function(outcome, imputed_df){
   all_df <- mice::complete(imputed_df, action = "all")
   
   fits <- lapply(all_df, function(dat){
      lm(
         as.formula(
            paste(
               outcome,
               "~ cc1_sex +
             cm3_diarhoea +
             cf12_exclusive +
             cf16_solid_food +
             hc2_age +
             education_level +
             occupation +
             minimum_mdd +
             muac_mother +
             weight_mother +
             hc3_family_size +
             hc7_private_agric_land +
             hc9_gro_crops +
             wash_toilets +
             wash_drinking_water +
             wash_handwashing_place +
             usually_fast +
             age_cat +
             urban_rural +
             hfias_category_insecure"
            )
         ),
         data = dat
      )
   })
   
   pooled <- mitools::MIcombine(
      results = lapply(fits, coef),
      variances = lapply(fits, function(fit){
         sandwich::vcovHC(fit, type = "HC3")
      })
   )
   
   s_pooled <- summary(pooled)
   results_tbl <- tibble::tibble(
      term      = rownames(s_pooled),
      estimate  = s_pooled[, "results"],
      conf.low  = s_pooled[, "results"] - 1.96 * s_pooled[, "se"],
      conf.high = s_pooled[, "results"] + 1.96 * s_pooled[, "se"],
      p.value   = 2 * pnorm(-abs(s_pooled[, "results"] / s_pooled[, "se"]))
   ) %>%
      dplyr::mutate(
         dplyr::across(c(estimate, conf.low, conf.high, p.value), ~ round(.x, 4))
      ) %>%
      dplyr::filter(term != "(Intercept)")
   
   list(
      pooled_results = results_tbl,
      fits = fits
   )
}

multi_WAZ_ols<- fxn_multi_model_ols("zwei", imp_df_updated)   
multi_WLZ_ols<- fxn_multi_model_ols("zwfl", imp_df_updated)    
multi_LAZ_ols<- fxn_multi_model_ols("zlen", imp_df_updated)


png("5. Plots/ WAZ_qqplots_all_imputations.png",
    width = 3000,
    height = 4000,
    res = 300)

par(mfrow = c(5, 4))

for(i in seq_along(multi_WAZ_ols$fits)) {
   plot(
      multi_WAZ_ols$fits[[i]],
      which = 2,
      main = paste("Imp", i)
   )
}

dev.off()


car::vif(multi_WAZ_ols$fits[[1]])


png("5. Plots/WAZ_OLS_diagnostics_imp1.png",
    width = 2400,
    height = 2400,
    res = 300)

par(mfrow = c(2, 2))

plot(multi_WAZ_ols$fits[[1]])

dev.off()


png("5. Plots/LAZ_OLS_diagnostics_imp1.png",
    width = 2400,
    height = 2400,
    res = 300)

par(mfrow = c(2, 2))

plot(multi_LAZ_ols$fits[[1]])

dev.off()


png("5. Plots/WLZ_OLS_diagnostics_imp1.png",
    width = 2400,
    height = 2400,
    res = 300)

par(mfrow = c(2, 2))

plot(multi_WLZ_ols$fits[[1]])

dev.off()





#######====combine the tables...
###helper function to clean gt tables for LaTeX outputer
clean_tables1 <- function(tbl_obj) {
   tbl_obj %>%
      gt::as_latex() %>%
      as.character() %>%
      stringr::str_replace_all("<span.*?>", "") %>%
      stringr::str_replace_all("</span>", "") %>%
      stringr::str_replace(
         "\\\\begin\\{table\\}\\[!h\\]",
         "\\\\begin\\{table\\}[htbp]"
      )
}
#####################================================================ Wasting
uni_wasting_tbl <- uni_wasting_binary %>%
   mutate(
      OR_95CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high)
   ) %>%
   select(variable, term, uni_OR_95CI = OR_95CI, uni_p = p.value)

multi_wasting_tbl <- multi_wasting_binary$pooled_results %>%
   mutate(
      OR_95CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high)
   ) %>%
   select(term, multi_OR_95CI = OR_95CI, multi_p = p.value)

# Merge
wasting_final_binary <- full_join(uni_wasting_tbl, multi_wasting_tbl, by = "term")

uni_WLZ_ols_tbl <- uni_WLZ_ols %>% 
   mutate(
      Beta_95CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high)
   ) %>%
   select(variable, term, uni_Beta_95CI = Beta_95CI, uni_ols_p = p.value)


multi_WLZ_ols_tbl <-  multi_WLZ_ols$pooled_results %>% 
   mutate(
      Beta_95CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high)
   ) %>%
   select(term, multi_Beta_95CI = Beta_95CI, multi_ols_p = p.value)
# Merge
wasting_final_ols <- full_join(uni_WLZ_ols_tbl, multi_WLZ_ols_tbl, by = "term")


#####################====================

wasting_final_combine_tbl <- wasting_final_binary %>%
   full_join(
      wasting_final_ols,
      by = c("variable", "term")
   ) %>%
   select(
      variable,
      term,
      #uni_Beta_95CI,
      #uni_ols_p,
      multi_Beta_95CI,
      multi_ols_p,
      #uni_OR_95CI,
      #uni_p,
      multi_OR_95CI,
      multi_p
   ) %>% 
   mutate(
      term = str_remove(term, paste0("^", variable))
   ) %>%
   gt::gt() %>%
   gt::cols_label(
      #uni_Beta_95CI = "Crude Beta (95% CI)", uni_ols_p = "p-value",
      multi_Beta_95CI = "Adjusted Beta (95% CI)", multi_ols_p = "p-value",
      #uni_OR_95CI = "Crude OR (95% CI)", uni_p = "p-value",
      multi_OR_95CI = "Adjusted OR (95% CI)", multi_p = "p-value"
   ) %>%
   gt::tab_spanner(
      label = "OLS regression",
      columns = c(multi_Beta_95CI, multi_ols_p)
   ) %>%
   gt::tab_spanner(
      label = "Logistic regression",
      columns = c(multi_OR_95CI, multi_p)
   ) %>%
   gt::opt_row_striping() %>%
   clean_tables1()


wasting_final_combine_tbl


#####################===========================================Underweight
## Underweight
uni_underweight_tbl <- uni_underweight_binary %>%
   mutate(
      OR_95CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high)
   ) %>%
   select(variable, term, uni_OR_95CI = OR_95CI, uni_p = p.value)

multi_underweight_tbl <- multi_underweight_binary$pooled_results %>%
   mutate(
      OR_95CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high)
   ) %>%
   select(term, multi_OR_95CI = OR_95CI, multi_p = p.value)

underweight_final_binary <- full_join(uni_underweight_tbl, multi_underweight_tbl, by = "term")


uni_WAZ_ols_tbl <- uni_WAZ_ols %>%
   mutate(
      Beta_95CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high)
   ) %>%
   select(variable, term, uni_Beta_95CI = Beta_95CI, uni_ols_p = p.value)

multi_WAZ_ols_tbl <- multi_WAZ_ols$pooled_results %>%
   mutate(
      Beta_95CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high)
   ) %>%
   select(term, multi_Beta_95CI = Beta_95CI, multi_ols_p = p.value)

underweight_final_ols <- full_join(uni_WAZ_ols_tbl, multi_WAZ_ols_tbl, by = "term")

####apend all
underweight_final_combine_tbl <- underweight_final_binary %>%
   full_join(
      underweight_final_ols,
      by = c("variable", "term")
   ) %>%
   select(
      variable,
      term,
      #uni_Beta_95CI,
      #uni_ols_p,
      multi_Beta_95CI,
      multi_ols_p,
      #uni_OR_95CI,
      #uni_p,
      multi_OR_95CI,
      multi_p
   ) %>%
   mutate(
      term = stringr::str_remove(term, paste0("^", variable))
   ) %>%
   gt::gt() %>%
   gt::cols_label(
      #uni_Beta_95CI = "Crude Beta (95% CI)", uni_ols_p = "p-value",
      multi_Beta_95CI = "Adjusted Beta (95% CI)", multi_ols_p = "p-value",
      #uni_OR_95CI = "Crude OR (95% CI)", uni_p = "p-value",
      multi_OR_95CI = "Adjusted OR (95% CI)", multi_p = "p-value"
   ) %>%
   gt::tab_spanner(
      label = "OLS regression",
      columns = c( multi_Beta_95CI, multi_ols_p)
   ) %>%
   gt::tab_spanner(
      label = "Logistic regression",
      columns = c( multi_OR_95CI, multi_p)
   ) %>%
   gt::opt_row_striping() %>%
   clean_tables1()


###########################################====stunting 
## Stunting
uni_stunting_tbl <- uni_stunting_binary %>%
   mutate(
      OR_95CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high)
   ) %>%
   select(variable, term, uni_OR_95CI = OR_95CI, uni_p = p.value)

multi_stunting_tbl <- multi_stunting_binary$pooled_results %>%
   mutate(
      OR_95CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high)
   ) %>%
   select(term, multi_OR_95CI = OR_95CI, multi_p = p.value)

stunting_final_binary <- full_join(uni_stunting_tbl, multi_stunting_tbl, by = "term")


uni_LAZ_ols_tbl <- uni_LAZ_ols %>%
   mutate(
      Beta_95CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high)
   ) %>%
   select(variable, term, uni_Beta_95CI = Beta_95CI, uni_ols_p = p.value)

multi_LAZ_ols_tbl <- multi_LAZ_ols$pooled_results %>%
   mutate(
      Beta_95CI = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high)
   ) %>%
   select(term, multi_Beta_95CI = Beta_95CI, multi_ols_p = p.value)

stunting_final_ols <- full_join(uni_LAZ_ols_tbl, multi_LAZ_ols_tbl, by = "term")

stunting_final_combine_tbl <- stunting_final_binary %>%
   full_join(
      stunting_final_ols,
      by = c("variable", "term")
   ) %>%
   select(
      variable,
      term,
      #uni_Beta_95CI,
      #uni_ols_p,
      multi_Beta_95CI,
      multi_ols_p,
      #uni_OR_95CI,
      #uni_p,
      multi_OR_95CI,
      multi_p
   ) %>%
   mutate(
      term = stringr::str_remove(term, paste0("^", variable))
   ) %>%
   gt::gt() %>%
   gt::cols_label(
      #uni_Beta_95CI = "Crude Beta (95% CI)", uni_ols_p = "p-value",
      multi_Beta_95CI = "Adjusted Beta (95% CI)", multi_ols_p = "p-value",
      #uni_OR_95CI = "Crude OR (95% CI)", uni_p = "p-value",
      multi_OR_95CI = "Adjusted OR (95% CI)", multi_p = "p-value"
   ) %>%
   gt::tab_spanner(
      label = "OLS regression",
      columns = c(multi_Beta_95CI, multi_ols_p)
   ) %>%
   gt::tab_spanner(
      label = "Logistic regression",
      columns = c( multi_OR_95CI, multi_p)
   ) %>%
   gt::opt_row_striping() %>%
   clean_tables1()


 gt::gtsave(
    underweight_final_combine_tbl,
    "4. Tables/underweight_final_combine_tbl_multi.tex"
 )

gt::gtsave(
   stunting_final_combine_tbl,
   "4. Tables/stunting_final_combine_tbl_multi.tex"
)

gt::gtsave(
   wasting_final_combine_tbl,
   "4. Tables/wasting_final_combine_tbl_multi.tex"
)


writeLines(
   as.character(wasting_final_combine_tbl),
   "4. Tables/wasting_final_combine_tbl_multi.tex"
)

writeLines(
   as.character(stunting_final_combine_tbl),
   "4. Tables/stunting_final_combine_tbl_multi.tex"
)

writeLines(
   as.character(underweight_final_combine_tbl),
   "4. Tables/underweight_final_combine_tbl_multi.tex"
)

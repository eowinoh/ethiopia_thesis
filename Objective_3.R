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
   glmmTMB, gt,
   quantreg
)

#######=====write functions to loop through variables
taus <- c(0.1, 0.25, 0.5, 0.8)
fxn_multi_quantile <- function(outcome, imputed_df, tau = 0.5, se = "boot", R = 500){
   
   all_df <- mice::complete(imputed_df, action = "all")
   
   fml <- as.formula(
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
   )
   
   fit_list <- lapply(all_df, function(dat){
      fit <- quantreg::rq(fml, tau = tau, method = "fn", data = dat)
      
      sm <- tryCatch(
         summary(fit, se = se, R = R),
         error = function(e) NULL
      )
      
      if (is.null(sm) || is.null(sm$coefficients)) return(NULL)
      
      coef_vec <- coef(fit)
      se_mat <- sm$coefficients
      
      # check standard error columns
      se_col <- grep("Std", colnames(se_mat), value = TRUE)[1]
      if (is.na(se_col)) se_col <- colnames(se_mat)[2]
      
      se_vec <- se_mat[, se_col]
      
      common <- intersect(names(coef_vec), rownames(se_mat))
      if (length(common) == 0) return(NULL)
      
      coef_vec <- coef_vec[common]
      se_vec   <- se_vec[common]
      
      list(
         coef = coef_vec,
         var  = diag(se_vec^2),
         fit  = fit
      )
   })
   
   fit_list <- Filter(Negate(is.null), fit_list)
   
   if (length(fit_list) == 0) {
      stop("No quantile regression models produced usable standard errors.")
   }
   
   pooled <- mitools::MIcombine(
      results   = lapply(fit_list, `[[`, "coef"),
      variances = lapply(fit_list, `[[`, "var")
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
      fits = lapply(fit_list, `[[`, "fit")
   )
}

#######-====var to clean names downstream
var_names <- c(
   "cc1_sex", "cm3_diarhoea", "cf12_exclusive", "cf16_solid_food", "hc2_age",
   "education_level", "occupation", "minimum_mdd", "muac_mother", "weight_mother",
   "hc3_family_size", "usually_fast", "hc7_private_agric_land", "hc9_gro_crops",
   "wash_toilets", "wash_drinking_water", "wash_handwashing_place",
   "age_cat", "urban_rural", "hfias_category_insecure"
)
###############======WAZ(Underweight)
multi_WAZ_qr_all <- map_dfr(taus, function(t) {
   res <- fxn_multi_quantile("zwei", imp_df_updated, tau = t)
   
   res$pooled_results %>%
      mutate(
         tau = t,
         Beta_95CI = sprintf("%.3f(%.3f, %.3f)",estimate, conf.low, conf.high)
      ) %>%
      select( term, tau, Beta_95CI, p.value)
})

multi_WAZ_qr_wide <- multi_WAZ_qr_all %>%
   mutate(tau = as.character(tau)) %>%
   pivot_wider(
      id_cols = c(term),
      names_from = tau,
      values_from = c(Beta_95CI, p.value),
      names_glue = "{tau}_{.value}"
   )


add_sig_stars <- function(beta_ci, p) {
   star <- case_when(
      is.na(p)    ~ "",
      p < 0.001   ~ "**",
      p < 0.05    ~ "*",
      TRUE        ~ ""
   )
   paste0(beta_ci, star)
}


multi_WAZ_qr_wide2 <- multi_WAZ_qr_wide %>%
   mutate(
      `0.1_Beta_95CI` = add_sig_stars(`0.1_Beta_95CI`, `0.1_p.value`),
      `0.25_Beta_95CI` = add_sig_stars(`0.25_Beta_95CI`, `0.25_p.value`),
      `0.5_Beta_95CI` = add_sig_stars(`0.5_Beta_95CI`, `0.5_p.value`),
      `0.8_Beta_95CI` = add_sig_stars(`0.8_Beta_95CI`, `0.8_p.value`)
   ) %>% 
   select(-ends_with("p.value")) %>% 
   rowwise() %>%
   mutate(
      variable = var_names[which(startsWith(term, var_names))[1]],
      level = sub(paste0("^", variable), "", term)
   ) %>%
   ungroup() %>% 
   select(variable, level, everything(), -term)



#####====gt table
waz_qr_table <- multi_WAZ_qr_wide2 %>%
   gt::gt() %>%
   gt::cols_label(
      `0.1_Beta_95CI` = "Estimate (95% CI)",
      `0.25_Beta_95CI` = "Estimate (95% CI)",
      `0.5_Beta_95CI` = "Estimate (95% CI)",
      `0.8_Beta_95CI` = "Estimate (95% CI)", 
   ) %>%
   gt::tab_spanner(
      label = "Tau = 0.1",
      columns = c(`0.1_Beta_95CI`)
   ) %>%
   gt::tab_spanner(
      label = "Tau = 0.25",
      columns = c(`0.25_Beta_95CI`)
   ) %>%
   gt::tab_spanner(
      label = "Tau = 0.5",
      columns = c(`0.5_Beta_95CI`)
   ) %>%
   gt::tab_spanner(
      label = "Tau = 0.8",
      columns = c(`0.8_Beta_95CI`)
   ) %>%
   gt::opt_row_striping() %>%
   clean_tables1()

writeLines(
   as.character(waz_qr_table),
   "4. Tables/waz_qr_table.tex"
)

##################=====WLZ(Wasting)
multi_wlz_qr_all <- map_dfr(taus, function(t) {
   res <- fxn_multi_quantile("zwfl", imp_df_updated, tau = t)
   
   res$pooled_results %>%
      mutate(
         tau = t,
         Beta_95CI = sprintf("%.4f (%.4f, %.4f)", estimate, conf.low, conf.high)
      ) %>%
      select( term, tau, Beta_95CI, p.value)
})

multi_wlz_qr_wide <- multi_wlz_qr_all %>%
   mutate(tau = as.character(tau)) %>%
   pivot_wider(
      id_cols = c(term),
      names_from = tau,
      values_from = c(Beta_95CI, p.value),
      names_glue = "{tau}_{.value}"
   )

multi_WLZ_qr_wide2 <- multi_wlz_qr_wide %>%
   mutate(
      `0.1_Beta_95CI` = add_sig_stars(`0.1_Beta_95CI`, `0.1_p.value`),
      `0.25_Beta_95CI` = add_sig_stars(`0.25_Beta_95CI`, `0.25_p.value`),
      `0.5_Beta_95CI` = add_sig_stars(`0.5_Beta_95CI`, `0.5_p.value`),
      `0.8_Beta_95CI` = add_sig_stars(`0.8_Beta_95CI`, `0.8_p.value`)
   ) %>% 
   select(-ends_with("p.value")) %>% 
   rowwise() %>%
   mutate(
      variable = var_names[which(startsWith(term, var_names))[1]],
      level = sub(paste0("^", variable), "", term)
   ) %>%
   ungroup() %>% 
   select(variable, level, everything(), -term)
#####====gt table
wlz_qr_table <- multi_WLZ_qr_wide2 %>%
   gt::gt() %>%
   gt::cols_label(
      `0.1_Beta_95CI` = "Estimate (95% CI)",
      `0.25_Beta_95CI` = "Estimate (95% CI)",
      `0.5_Beta_95CI` = "Estimate (95% CI)",
      `0.8_Beta_95CI` = "Estimate (95% CI)", 
   ) %>%
   gt::tab_spanner(
      label = "Tau = 0.1",
      columns = c(`0.1_Beta_95CI`)
   ) %>%
   gt::tab_spanner(
      label = "Tau = 0.25",
      columns = c(`0.25_Beta_95CI`)
   ) %>%
   gt::tab_spanner(
      label = "Tau = 0.5",
      columns = c(`0.5_Beta_95CI`)
   ) %>%
   gt::tab_spanner(
      label = "Tau = 0.8",
      columns = c(`0.8_Beta_95CI`)
   ) %>%
   gt::opt_row_striping() %>%
   clean_tables1()

writeLines(
   as.character(wlz_qr_table),
   "4. Tables/wlz_qr_table.tex"
)


#######=======LAZ(stunting)
multi_laz_qr_all <- map_dfr(taus, function(t) {
   res <- fxn_multi_quantile("zlen", imp_df, tau = t)
   
   res$pooled_results %>%
      mutate(
         tau = t,
         Beta_95CI = sprintf("%.3f(%.3f, %.3f)", estimate, conf.low, conf.high)
      ) %>%
      select( term, tau, Beta_95CI, p.value)
})

multi_laz_qr_wide <- multi_laz_qr_all %>%
   mutate(tau = as.character(tau)) %>%
   pivot_wider(
      id_cols = c(term),
      names_from = tau,
      values_from = c(Beta_95CI, p.value),
      names_glue = "{tau}_{.value}"
   )

multi_LAZ_qr_wide2 <- multi_laz_qr_wide %>%
   mutate(
      `0.1_Beta_95CI` = add_sig_stars(`0.1_Beta_95CI`, `0.1_p.value`),
      `0.25_Beta_95CI` = add_sig_stars(`0.25_Beta_95CI`, `0.25_p.value`),
      `0.5_Beta_95CI` = add_sig_stars(`0.5_Beta_95CI`, `0.5_p.value`),
      `0.8_Beta_95CI` = add_sig_stars(`0.8_Beta_95CI`, `0.8_p.value`)
   ) %>% 
   select(-ends_with("p.value")) %>% 
   rowwise() %>%
   mutate(
      variable = var_names[which(startsWith(term, var_names))[1]],
      level = sub(paste0("^", variable), "", term)
   ) %>%
   ungroup() %>% 
   select(variable, level, everything(), -term)


#####====gt table
laz_qr_table <- multi_LAZ_qr_wide2 %>%
   gt::gt() %>%
   gt::cols_label(
      `0.1_Beta_95CI` = "Estimate (95% CI)",
      `0.25_Beta_95CI` = "Estimate (95% CI)",
      `0.5_Beta_95CI` = "Estimate (95% CI)",
      `0.8_Beta_95CI` = "Estimate (95% CI)", 
   ) %>%
   gt::tab_spanner(
      label = "Tau = 0.1",
      columns = c(`0.1_Beta_95CI`)
   ) %>%
   gt::tab_spanner(
      label = "Tau = 0.25",
      columns = c(`0.25_Beta_95CI`)
   ) %>%
   gt::tab_spanner(
      label = "Tau = 0.5",
      columns = c(`0.5_Beta_95CI`)
   ) %>%
   gt::tab_spanner(
      label = "Tau = 0.8",
      columns = c(`0.8_Beta_95CI`)
   ) %>%
   gt::opt_row_striping() %>%
   clean_tables1()

writeLines(
   as.character(laz_qr_table),
   "4. Tables/laz_qr_table.tex"
)



########Assessing quantile crossing 
taus <- c(0.10, 0.25, 0.50, 0.80)

f <- as.formula(
   "zlen ~ cc1_sex +
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

# Grab just a few completed datasets
subset_idx <- 1:20
imp_list <- lapply(subset_idx, function(i) complete(imp_df_updated, action = i))

crossing_rate <- sapply(imp_list, function(d) {
   fit <- rq(f, tau = taus, data = d)
   p <- predict(fit)
   mean(apply(p, 1, function(r) any(diff(r) < 0)))
})

cat(sprintf("Crossing per imputation: mean %.2f%%, range %.2f–%.2f%%\n",
            mean(crossing_rate) * 100,
            min(crossing_rate) * 100,
            max(crossing_rate) * 100))
crossing_rate * 100   # per-imputation values

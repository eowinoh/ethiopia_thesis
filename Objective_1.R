#########################################################################################
#########################################################################################
#########################################################################################
##############
##############Identify determinants of maternal and child dietary diversity.
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
   glmmTMB
)

#########################################################################################
##############Maternal Diet Diversity
##########################################################################################
#####==Maternal MDD vars
mdd_vars_model <- c(
   "hc2_age",
   "education_level",
   "occupation",
   "hfias_category_insecure",
   "hc3_family_size",
   "usually_fast",
   "md_illness",
   "muac_mother",
   "weight_mother",
   "hc9_gro_crops",
   "hc7_private_agric_land",
   "urban_rural"
)
mdd_desc <- tbl_summary(
   analysis_data_clean1 %>%
      select(minimum_mdd, all_of(mdd_vars_model)),
   by = minimum_mdd,
   missing = "always",
   missing_text = "Missing",
   type = all_dichotomous() ~ "categorical",
   statistic = all_categorical() ~ "{n} ({p}%)"
) %>%
   add_p() %>%    clean_tables()

writeLines(
   as.character(mdd_desc),
   "4. Tables/objective_1_mdd_table.tex"
)


#######Maternal Diet diversity determinants - univariable logistic regression
mdd_uni_results_binomial <- lapply(mdd_vars_model, function(v) {
   #£££££#####=== loop through each imputed dataset, fit the model, and store the result
   mdd_fits <- lapply(mice::complete(imp_df_updated, action = "all"), function(df) {
      df$minimum_mdd <- factor(
         df$minimum_mdd,
         levels = c("LDD", "ADD")
      )
      
      glm(
         as.formula(paste("minimum_mdd ~", v)),
         family = binomial,
         data = df
      )
   })
   ####========== pool using robust variance from each fit
   pool_res <- mitools::MIcombine(
      results = lapply(mdd_fits, coef),
      variances = lapply(mdd_fits, function(fit) {
         sandwich::vcovHC(fit, type = "HC3")
      })
   )
   ######========summarise results and output as a tibble
   mdd_pool_summary <- summary(pool_res)
   
   tibble(
      term = rownames(mdd_pool_summary),
      estimate = exp(mdd_pool_summary[, "results"]),
      std.error = mdd_pool_summary[, "se"],
      conf.low = exp(mdd_pool_summary[, "results"] - 1.96 * mdd_pool_summary[, "se"]),
      conf.high = exp(mdd_pool_summary[, "results"] + 1.96 * mdd_pool_summary[, "se"]),
      p.value = 2 * pnorm(-abs(mdd_pool_summary[, "results"] / mdd_pool_summary[, "se"])),
      variable = v
   )
})

mdd_uni_results_table_bin <- bind_rows(mdd_uni_results_binomial) %>%
   filter(term != "(Intercept)") %>%
   select(variable, term, estimate, std.error, conf.low, conf.high, p.value)


#######Maternal Diet diversity determinants - multivariable logistic regression
mdd_multi_fit_bin <- lapply(mice::complete(imp_df, action = "all"), function(df) {
   df$minimum_mdd <- factor(
      df$minimum_mdd,
      levels = c("LDD", "ADD")
   )
   glm(
      minimum_mdd ~
         hc2_age +
         education_level +
         occupation +
         hfias_category_insecure +
         md_illness +
         hc3_family_size +
         usually_fast +
         muac_mother +
         weight_mother +
         urban_rural +
         hc9_gro_crops,
      family = binomial,
      data = df
   )
})

mdd_multi_fit_pool <- mitools::MIcombine(
   results   = lapply(mdd_multi_fit_bin, coef),
   variances = lapply(mdd_multi_fit_bin, function(fit) vcovHC(fit, type = "HC3"))
)

mdd_multi_fit_summary <- summary(mdd_multi_fit_pool)

mdd_multi_results_table_bin <- tibble(
   term = rownames(mdd_multi_fit_summary),
   estimate = exp(mdd_multi_fit_summary[, "results"]),
   std.error = mdd_multi_fit_summary[, "se"],
   conf.low = exp(mdd_multi_fit_summary[, "results"] - 1.96 * mdd_multi_fit_summary[, "se"]),
   conf.high = exp(mdd_multi_fit_summary[, "results"] + 1.96 * mdd_multi_fit_summary[, "se"]),
   p.value = 2 * pnorm(-abs(mdd_multi_fit_summary[, "results"] / mdd_multi_fit_summary[, "se"]))
) %>%
   filter(term != "(Intercept)")

mdd_multi_results_table_bin
######################===MDD model assessemt
#check convergence
sapply(mdd_multi_fit_bin, function(x) x$converged)
library(car)
vif_list <- lapply(mdd_multi_fit_bin, vif)

diag_results <- lapply(mdd_multi_fit_bin, function(fit) {
   res <- simulateResiduals(fit)
   
   c(
      uniformity_p = testUniformity(res)$p.value,
      dispersion_p = testDispersion(res)$p.value,
      outlier_p = testOutliers(res)$p.value
   )
})

do.call(rbind, diag_results)


#######======DHARMa residual diagnostics plots(MDD- Binomial)

sampled_imps <- sample(seq_along(mdd_multi_fit_bin), 4)

# Create plots
for (i in sampled_imps) {
   
   sim_res <- simulateResiduals(mdd_multi_fit_bin[[i]])
   
   png(
      filename = sprintf("8. Temp/mdd_bin_%02d.png", i),
      width = 1600,
      height = 1200,
      res = 150
   )
   
   plot(sim_res, main = "")
   title(main = paste("Imputation", i), line = 1)
   
   dev.off()
}

# Read only sampled images
files <- sprintf("8. Temp/mdd_bin_%02d.png", sampled_imps)
imgs <- image_read(files)
#####==== 2 × 2 layout
combined_mdd_bin <- image_montage(
   imgs,
   tile = "2x2",
   geometry = "1600x1200+10+10"
)

image_write(
   combined_mdd_bin,
   path = "5. Plots/mdd_bin_dharma_2x2.png",
   format = "png"
)

sampled_imps


library(splines)

fit_linear <- glm(
   minimum_mdd ~ hc2_age + education_level + occupation +
      hfias_category_insecure + wash_toilets + md_illness +
      hc3_family_size + usually_fast + muac_mother +
      weight_mother + urban_rural + hc9_gro_crops,
   family = binomial,
   data = complete(imp_df_updated, 1)
)

#######functional form of linear covariates
library(splines)
fit_spline_age <- glm(
   minimum_mdd ~ ns(hc2_age, df = 3) + education_level + occupation +
      hfias_category_insecure + wash_toilets + md_illness +
      hc3_family_size + usually_fast + muac_mother +
      weight_mother + urban_rural + hc9_gro_crops,
   family = binomial,
   data = complete(imp_df_updated, 1)
)
fit_spline_muac <- glm(
   minimum_mdd ~ hc2_age + education_level + occupation +
      hfias_category_insecure + wash_toilets + md_illness +
      hc3_family_size + usually_fast + ns(muac_mother, df = 3) +
      weight_mother + urban_rural + hc9_gro_crops,
   family = binomial,
   data = complete(imp_df_updated, 1)
)
anova(fit_linear, fit_spline_muac, test = "Chisq")
fit_spline_weight <- glm(
   minimum_mdd ~ hc2_age + education_level + occupation +
      hfias_category_insecure + wash_toilets + md_illness +
      hc3_family_size + usually_fast + ns(weight_mother, df = 3) +
      muac_mother + urban_rural + hc9_gro_crops,
   family = binomial,
   data = complete(imp_df_updated, 1)
)
anova(fit_linear, fit_spline_weight, test = "Chisq")

fit_spline_f_size <- glm(
   minimum_mdd ~ hc2_age + education_level + occupation +
      hfias_category_insecure + wash_toilets + md_illness +
      weight_mother + usually_fast + ns(hc3_family_size, df = 3) +
      muac_mother + urban_rural + hc9_gro_crops,
   family = binomial,
   data = complete(imp_df_updated, 1)
)
anova(fit_linear, fit_spline_f_size, test = "Chisq")

rm(list = ls(pattern = "^fit_"))
######==check hosmer-lemeshow goodness of fit test for each imputed model
hl_results_mdd_bin <- lapply(mdd_multi_fit_bin, function(m) {
   hoslem.test(m$y, fitted(m))
})

## check p-values for all the simulated models
hl_table_mdd_bin <- data.frame(
   Model = names(hl_results_mdd_bin),
   Chi_square = sapply(hl_results_mdd_bin, function(x) unname(x$statistic)),
   DF = sapply(hl_results_mdd_bin, function(x) x$parameter),
   P_value = sapply(hl_results_mdd_bin, function(x) x$p.value)
)

hl_table_mdd_bin
# library(pROC)
# plot(calibrate(mdd_multi_fit_bin))

#####==== cooks distance for each imputed modelan
quartz(width = 14, height = 12)
cooks_list_mdd_bin <- lapply(mdd_multi_fit_bin, cooks.distance)
png("5. Plots/cooks_mdd_binomial.png", width = 3000, height = 2400, res = 300)
par(mfrow = c(5, 4))

# Plot Cook's distance for each imputed model
for(i in seq_along(cooks_list_mdd_bin)) {
   
   plot(
      cooks_list_mdd_bin[[i]],
      type = "h",
      main = paste("Imputation", i),
      xlab = "Observation index",
      ylab = "Cook's distance",
      ylim = c(0, 0.05)
   )
}

# Reset plotting layout
par(mfrow = c(1, 1))
dev.off()

#######Maternal Diet diversity determinants - univariable poisson regression
mdd_uni_results_poisson <- lapply(mdd_vars_model, function(v) {
   #£££££#####=== loop through each imputed dataset, fit the model, and store the result
   mdd_pois_fits <- lapply(mice::complete(imp_df, action = "all"), function(df) {
      glm(as.formula(paste("mdd_score ~", v)),
          family = poisson(link = "log"),
          data = df)
   })
   ####========== pool using robust variance from each fit
   pool_res_pois <- mitools::MIcombine(
      results = lapply(mdd_pois_fits, coef),
      variances = lapply(mdd_pois_fits, function(fitp) {
         sandwich::vcovHC(fitp, type = "HC3")
      })
   )
   ######========summarise results and output as a tibble
   mdd_pool_summary_pois <- summary(pool_res_pois)
   
   tibble(
      term = rownames(mdd_pool_summary_pois),
      estimate = exp(mdd_pool_summary_pois[, "results"]),
      std.error = mdd_pool_summary_pois[, "se"],
      conf.low = exp(mdd_pool_summary_pois[, "results"] - 1.96 * mdd_pool_summary_pois[, "se"]),
      conf.high = exp(mdd_pool_summary_pois[, "results"] + 1.96 * mdd_pool_summary_pois[, "se"]),
      p.value = 2 * pnorm(-abs(mdd_pool_summary_pois[, "results"] / mdd_pool_summary_pois[, "se"])),
      variable = v
   )
})

mdd_uni_results_table_pois <- bind_rows(mdd_uni_results_poisson) %>%
   filter(term != "(Intercept)") %>%
   select(variable, term, estimate, std.error, conf.low, conf.high, p.value)

mdd_uni_results_table_pois


#######==========Maternal Diet diversity determinants - multivariable logistic regression
mdd_multi_fit_pois <- lapply(mice::complete(imp_df_updated, action = "all"), function(df) {
   glm(
      mdd_score ~
         hc2_age +
         education_level +
         occupation +
         hfias_category_insecure +
         md_illness +
         hc3_family_size +
         usually_fast +
         weight_mother +
         muac_mother +
         urban_rural +
         hc9_gro_crops,
      family = poisson(link = "log"),
      data = df
   )
})

mdd_multi_fit_pool_pois <- mitools::MIcombine(
   results   = lapply(mdd_multi_fit_pois, coef),
   variances = lapply(mdd_multi_fit_pois, function(fit) vcovHC(fit, type = "HC3"))
)

mdd_multi_fit_summary_pois <- summary(mdd_multi_fit_pool_pois)

mdd_multi_results_table_pois <- tibble(
   term = rownames(mdd_multi_fit_summary_pois),
   estimate = exp(mdd_multi_fit_summary_pois[, "results"]),
   std.error = mdd_multi_fit_summary_pois[, "se"],
   conf.low = exp(mdd_multi_fit_summary_pois[, "results"] - 1.96 * mdd_multi_fit_summary_pois[, "se"]),
   conf.high = exp(mdd_multi_fit_summary_pois[, "results"] + 1.96 * mdd_multi_fit_summary_pois[, "se"]),
   p.value = 2 * pnorm(-abs(mdd_multi_fit_summary_pois[, "results"] / mdd_multi_fit_summary_pois[, "se"]))
) %>%
   filter(term != "(Intercept)")

mdd_multi_results_table_pois

##########===== Dispersion check for Poisson model

dispersion_checks_pois <- lapply(1:imp_df$m, function(i) {
   dat <- complete(imp_df, i)
   fit_mdd <- glm(
      mdd_score ~
         hc2_age +
         education_level +
         occupation +
         hfias_category_insecure +
         wash_toilets +
         md_illness +
         hc3_family_size +
         usually_fast +
         weight_mother +
         urban_rural +
         hc9_gro_crops,
      family = poisson(link = "log"),
      data = dat
   )
   ######=== Dispersion = residual deviance / df
    data.frame(
      imputation  = i,
      dispersion  = fit_mdd$deviance / fit_mdd$df.residual,
      deviance    = fit_mdd$deviance,
      df          = fit_mdd$df.residual
   )
})
dispersion_summary_pois <- do.call(rbind, dispersion_checks_pois)
print(dispersion_summary_pois)
t.test(dispersion_summary_pois$dispersion)$conf.int


#######======DHARMa residual diagnostics plots(MDD- Poisson)
#######
#######

sampled_imps_pois <- sample(seq_along(mdd_multi_fit_pois), 4)

# Create plots
for (i in sampled_imps_pois) {
   
   sim_res <- simulateResiduals(mdd_multi_fit_pois[[i]])
   
   png(
      filename = sprintf("8. Temp/mdd_pois_%02d.png", i),
      width = 1600,
      height = 1200,
      res = 150
   )
   
   plot(sim_res, main = "")
   title(main = paste("Imputation", i), line = 1)
   
   dev.off()
}

# Read only sampled images
files <- sprintf("8. Temp/mdd_pois_%02d.png", sampled_imps_pois)
imgs <- image_read(files)
#####==== 2 × 2 layout
combined_mdd_bin <- image_montage(
   imgs,
   tile = "2x2",
   geometry = "1600x1200+10+10"
)

image_write(
   combined_mdd_bin,
   path = "5. Plots/mdd_pois_dharma_2x2.png",
   format = "png"
)

#######==Clear under dispersion present: Try CMP model
mdd_multi_fit_cmp <- lapply(1:imp_df_updated$m, function(i) {
   dat <- complete(imp_df_updated, i)
   glmmTMB::glmmTMB(
      mdd_score ~
         hc2_age +
         education_level +
         occupation +
         hfias_category_insecure +
         md_illness +
         hc3_family_size +
         usually_fast +
         weight_mother +
         muac_mother +
         urban_rural +
         hc9_gro_crops,
      family = compois(link = "log"),
      data = dat
   )
})

#####≠=Pool using Rubin's Rules using model based vcovariance as CMP correctly models dispersion issues

pool_mdd_glmmtmb <- function(mods) {
   m <- length(mods)
   coefs <- lapply(mods, function(fit) fixef(fit)$cond)
   vcovs <- lapply(mods, function(fit) vcov(fit)$cond)
   Qbar <- Reduce("+", coefs) / m
   Ubar <- Reduce("+", vcovs) / m
   B    <- var(do.call(rbind, coefs))
   T_var <- Ubar + (1 + 1/m) * B
   se <- sqrt(diag(T_var))
   z  <- Qbar / se
   p  <- 2 * pnorm(abs(z), lower.tail = FALSE)
   data.frame(
      term      = names(Qbar),
      PR        = round(exp(Qbar), 3),
      conf.low  = round(exp(Qbar - 1.96 * se), 3),
      conf.high = round(exp(Qbar + 1.96 * se), 3),
      p.value   = round(p, 3)
   )
}
mdd_multi_results_cmp <- pool_mdd_glmmtmb(mdd_multi_fit_cmp) %>% 
   filter(term != "(Intercept)")
print(mdd_multi_results_cmp)

###=== Assess for equidispersion
###

dispersion_cmp <- sapply(mdd_multi_fit_cmp, function(fit) {
   as.numeric(sigma(fit))
})
nu_cmp <- sapply(mdd_multi_fit_cmp, function(fit) {
   1 / sigma(fit)
})

#######======Diagnostic plots


sampled_imps_cmp <- sample(seq_along(mdd_multi_fit_cmp), 4)

# Create plots
for (i in sampled_imps_cmp) {
   
   sim_res <- simulateResiduals(mdd_multi_fit_cmp[[i]])
   
   png(
      filename = sprintf("8. Temp/mdd_cmp_%02d.png", i),
      width = 1600,
      height = 1200,
      res = 150
   )
   
   plot(sim_res, main = "")
   title(main = paste("Imputation", i), line = 1)
   
   dev.off()
}

# Read only sampled images
files <- sprintf("8. Temp/mdd_cmp_%02d.png", sampled_imps_cmp)
imgs <- image_read(files)
#####==== 2 × 2 layout
combined_mdd_bin <- image_montage(
   imgs,
   tile = "2x2",
   geometry = "1600x1200+10+10"
)

image_write(
   combined_mdd_bin,
   path = "5. Plots/mdd_cmp_dharma_2x2.png",
   format = "png"
)
#########################################################################################
##############Child Diet Diversity
##########################################################################################
#####==Child MDD vars

cdd_vars_model <- c("cc1_sex","cm3_diarhoea","cf12_exclusive",
                    "cf16_solid_food", "usually_fast", "education_level", "occupation", "urban_rural", "hc2_age", "hc3_family_size",
                    "hc9_gro_crops","minimum_mdd",
                    "cm3_diarhoea","age_cat","hfias_category_insecure")



cdd_desc <- tbl_summary(
   analysis_data_clean1 %>%
      filter(cc4_age_months>=6) %>%
      select(minimum_cdd, all_of(cdd_vars_model)),
   by = minimum_cdd,
   missing = "always",
   missing_text = "Missing",
   type = all_dichotomous() ~ "categorical",
   statistic = all_categorical() ~ "{n} ({p}%)"
) %>%
   add_p() %>% 
   clean_tables()

writeLines(
   as.character(cdd_desc),
   "4. Tables/objective_1_cdd_tableR.tex"
)
######====specify methods for imputation(for children only above 6 months)

######=== re-impute for children
analysis_data_clean_children <- analysis_data_clean1 %>%
   filter(cc4_age_months >= 6)


######====specify methods for imputation
ignore_vars <- c(
   "stunting",
   "wasting",
   "underweight",
   "cdd_score",
   "minimum_cdd",
   "age_cat","birth_size_combine",
   "hfias_category_insecure",
   "cdd_score", "minimum_cdd","minimum_mdd","anaemia"
)

meths <- make.method(analysis_data_clean_children)
meths[ignore_vars] <- ""



#####===remove unnecessary vars from predictor matrix
ignore_vars_pred <- c(
   "stunting",
   "wasting",
   "underweight",
   "cdd_score",
   "age_cat","anaemia","birth_size_combine",
   "hfias_category_insecure",
   "minimum_cdd","minimum_mdd"
)
pred <- make.predictorMatrix(analysis_data_clean_children)
pred[ignore_vars_pred, ] <- 0   
pred[, ignore_vars_pred] <- 0   

######===run MI
imp_df_c <- mice(analysis_data_clean_children, 
                 method = meths, seed = 20250520,
                 predictorMatrix = pred,
                 m = 20)

###= Recategorise the variables
completed_df_c <- complete(imp_df_c, "long", include = TRUE) %>% 
   mutate(stunting = if_else(zlen < -2, "Yes", "No"),
          wasting = if_else(zwfl < -2, "Yes", "No"),
          underweight = if_else(zwei < -2, "Yes", "No"),
          stunting = factor(stunting, levels = c("No", "Yes")),
          wasting = factor(wasting, levels = c("No", "Yes")),
          underweight = factor(underweight, levels = c("No", "Yes")),
          minimum_cdd = factor(minimum_cdd, levels = c("LDD", "ADD")),
          minimum_mdd = factor(minimum_mdd, levels = c("LDD", "ADD")),
          hfias_category_insecure = case_when(
             hfias_category %in% c(
                "Mildly food insecure",
                "Moderately food insecure",
                "Severely food insecure"
             ) ~ "Yes",
             
             hfias_category == "Food secure" ~ "No",
             
             TRUE ~ NA_character_
          ),
          age_cat = cut(
             cc4_age_months,
             breaks = c(0, 6, 12, 24),
             right = FALSE,
             labels = c("0-5", "6-11", "12-24")
          ),
          birth_size_combine = case_when(
             cc6_size_born %in% c("Smaller than average", "Very small") ~ "Small",
             cc6_size_born %in% c("Larger than average", "Very large") ~ "Large",
             cc6_size_born == "Average" ~ "Average",
             TRUE ~ NA_character_
          ),
          birth_size_combine = factor(birth_size_combine, levels = c("Small", "Average", "Large")),
          anaemia = case_when(
             hgb_child < 11 ~ "Yes",
             hgb_child >= 11 ~ "No",
             TRUE ~ NA_character_
          )
   )

imp_df_updated_c <- as.mids(completed_df_c)          
 ###########==== Logistics model
#######Child Diet diversity determinants - univariable logistic regression
cdd_uni_results_binomial <- lapply(cdd_vars_model, function(v) {
#£££££#####=== loop through each imputed dataset, fit the model, and store the result
cdd_fits <- lapply(mice::complete(imp_df_updated_c, action = "all"), function(df) {
                df$minimum_cdd <- factor(
                   df$minimum_cdd,
                   levels = c("LDD", "ADD")
                )
                glm(
                   as.formula(paste("minimum_cdd ~", v)),
                   family = binomial,
                   data = df
                )
             })
####========== pool using robust variance from each fit
 pool_res <- mitools::MIcombine(
                results = lapply(cdd_fits, coef),
                variances = lapply(cdd_fits, function(fit) {
                   sandwich::vcovHC(fit, type = "HC3")
                })
             )
 ######========summarise results and output as a tibble
cdd_pool_summary <- summary(pool_res)
             
             tibble(
                term = rownames(cdd_pool_summary),
                estimate = exp(cdd_pool_summary[, "results"]),
                std.error = cdd_pool_summary[, "se"],
                conf.low = exp(cdd_pool_summary[, "results"] - 1.96 * cdd_pool_summary[, "se"]),
                conf.high = exp(cdd_pool_summary[, "results"] + 1.96 * cdd_pool_summary[, "se"]),
                p.value = 2 * pnorm(-abs(cdd_pool_summary[, "results"] / cdd_pool_summary[, "se"])),
                variable = v
             )
          })
          
cdd_uni_results_table_bin <- bind_rows(cdd_uni_results_binomial) %>%
             filter(term != "(Intercept)") %>%
             select(variable, term, estimate, std.error, conf.low, conf.high, p.value)
          
#######Child Diet diversity determinants - multivariable logistic regression
          
cdd_multi_fit_bin <- lapply(mice::complete(imp_df_updated_c, action = "all"), function(df) {
             df$minimum_cdd <- factor(
                df$minimum_cdd,
                levels = c("LDD", "ADD")
             )
             glm(
                minimum_cdd ~
                   cc1_sex +
                   cm3_diarhoea +
                   cf16_solid_food +
                   usually_fast +
                   education_level +
                   occupation +
                   urban_rural +
                   hc2_age +
                   hc3_family_size +
                   hc9_gro_crops +
                   age_cat +
                   minimum_mdd+
                   hfias_category_insecure,
                family = binomial(),
                data = df
             )
          })
          
cdd_multi_fit_pool <- mitools::MIcombine(
             results   = lapply(cdd_multi_fit_bin, coef),
             variances = lapply(cdd_multi_fit_bin, function(fit) vcovHC(fit, type = "HC3"))
          )
          
cdd_multi_fit_summary <- summary(cdd_multi_fit_pool)
          
cdd_multi_results_table_bin <- tibble(
             term = rownames(cdd_multi_fit_summary),
             estimate = exp(cdd_multi_fit_summary[, "results"]),
             std.error = cdd_multi_fit_summary[, "se"],
             conf.low = exp(cdd_multi_fit_summary[, "results"] - 1.96 * cdd_multi_fit_summary[, "se"]),
             conf.high = exp(cdd_multi_fit_summary[, "results"] + 1.96 * cdd_multi_fit_summary[, "se"]),
             p.value = 2 * pnorm(-abs(cdd_multi_fit_summary[, "results"] / cdd_multi_fit_summary[, "se"]))
          ) %>%
             filter(term != "(Intercept)")
          
cdd_multi_results_table_bin

########===============================++++Model assessment
######################===================CDD model assessment
#check convergence

sapply(cdd_multi_fit_bin, function(x) x$converged)
library(car)
vif_list <- lapply(cdd_multi_fit_bin, vif)

diag_results <- lapply(cdd_multi_fit_bin, function(fit) {
   res <- simulateResiduals(fit)
   
   c(
      uniformity_p = testUniformity(res)$p.value,
      dispersion_p = testDispersion(res)$p.value,
      outlier_p = testOutliers(res)$p.value
   )
})

do.call(rbind, diag_results)
hl_results_cdd_bin <- lapply(cdd_multi_fit_bin, function(m) {
   hoslem.test(m$y, fitted(m))
})
hl_results_cdd_bin


#######======DHARMa residual diagnostics plots(CDD- Binomial)

sampled_imps_cdd_bin <- sample(seq_along(cdd_multi_fit_bin), 4)

# Create plots
for (i in sampled_imps_cdd_bin) {
   
   sim_res <- simulateResiduals(cdd_multi_fit_bin[[i]])
   
   png(
      filename = sprintf("8. Temp/cdd_bin_%02d.png", i),
      width = 1600,
      height = 1200,
      res = 150
   )
   
   plot(sim_res, main = "")
   title(main = paste("Imputation", i), line = 1)
   
   dev.off()
}

# Read only sampled images
files <- sprintf("8. Temp/cdd_bin_%02d.png", sampled_imps_cdd_bin)
imgs_cdd_bin <- image_read(files)
#####==== 2 × 2 layout
combined_cdd_bin <- image_montage(
   imgs_cdd_bin,
   tile = "2x2",
   geometry = "1600x1200+10+10"
)

image_write(
   combined_cdd_bin,
   path = "5. Plots/cdd_bin_dharma_2x2.png",
   format = "png"
)




library(splines)

fit_linear <- glm(
   minimum_cdd ~
      cc1_sex +
      cm3_diarhoea +
      cf16_solid_food +
      usually_fast +
      education_level +
      occupation +
      urban_rural +
      hc2_age +
      hc3_family_size +
      hc9_gro_crops +
      age_cat +
      minimum_mdd+
      hfias_category_insecure,
   family = binomial,
   data = complete(imp_df_updated_c, 1)
)


#######functional form of linear covariates
library(splines)
fit_spline_age <- glm(
   minimum_cdd ~
      cc1_sex +
      cm3_diarhoea +
      cf16_solid_food +
      usually_fast +
      education_level +
      occupation +
      urban_rural +
      ns(hc2_age,df = 3) +
      hc3_family_size +
      hc9_gro_crops +
      age_cat +
      minimum_mdd+
      hfias_category_insecure,
   family = binomial,
   data = complete(imp_df_updated_c, 1)
)
anova(fit_linear, fit_spline_age, test = "Chisq")

fit_spline_f_size <- glm(
   minimum_cdd ~
      cc1_sex +
      cm3_diarhoea +
      cf16_solid_food +
      usually_fast +
      education_level +
      occupation +
      urban_rural +
      hc2_age +
      ns(hc3_family_size, df = 3)+
      hc9_gro_crops +
      age_cat +
      minimum_mdd+
      hfias_category_insecure,
   
   family = binomial,
   data = complete(imp_df_updated, 1)
)
anova(fit_linear, fit_spline_f_size, test = "Chisq")

rm(list = ls(pattern = "^fit_"))
#######Maternal Diet diversity determinants - univariable poisson regression
cdd_uni_results_poisson <- lapply(cdd_vars_model, function(v) {
   #£££££#####=== loop through each imputed dataset, fit the model, and store the result
   cdd_pois_fits <- lapply(mice::complete(imp_df_updated_c, action = "all"), function(df) {
      glm(as.formula(paste("cdd_score ~", v)),
          family = poisson(link = "log"),
          data = df)
   })
   ####========== pool using robust variance from each fit
   pool_res_pois <- mitools::MIcombine(
      results = lapply(cdd_pois_fits, coef),
      variances = lapply(cdd_pois_fits, function(fitp) {
         sandwich::vcovHC(fitp, type = "HC3")
      })
   )
   ######========summarise results and output as a tibble
   cdd_pool_summary_pois <- summary(pool_res_pois)
   
   tibble(
      term = rownames(cdd_pool_summary_pois),
      estimate = exp(cdd_pool_summary_pois[, "results"]),
      std.error = cdd_pool_summary_pois[, "se"],
      conf.low = exp(cdd_pool_summary_pois[, "results"] - 1.96 * cdd_pool_summary_pois[, "se"]),
      conf.high = exp(cdd_pool_summary_pois[, "results"] + 1.96 * cdd_pool_summary_pois[, "se"]),
      p.value = 2 * pnorm(-abs(cdd_pool_summary_pois[, "results"] / cdd_pool_summary_pois[, "se"])),
      variable = v
   )
})

cdd_uni_results_table_pois <- bind_rows(cdd_uni_results_poisson) %>%
   filter(term != "(Intercept)") %>%
   select(variable, term, estimate, std.error, conf.low, conf.high, p.value)

cdd_uni_results_table_pois


#######==========Child Diet diversity determinants - multivariable poison regression
cdd_multi_fit_pois <- lapply(mice::complete(imp_df_updated_c, action = "all"), function(df) {
   glm(
      cdd_score ~
         cc1_sex +
         cm3_diarhoea +
         cf16_solid_food +
         usually_fast +
         education_level +
         occupation +
         urban_rural +
         hc2_age +
         hc3_family_size +
         hc9_gro_crops +
         age_cat +
         minimum_mdd+
         hfias_category_insecure,
      family = poisson(link = "log"),
      data = df
   )
})

cdd_multi_fit_pool_pois <- mitools::MIcombine(
   results   = lapply(cdd_multi_fit_pois, coef),
   variances = lapply(cdd_multi_fit_pois, function(fit) vcovHC(fit, type = "HC3"))
)

cdd_multi_fit_summary_pois <- summary(cdd_multi_fit_pool_pois)

cdd_multi_results_table_pois <- tibble(
   term = rownames(cdd_multi_fit_summary_pois),
   estimate = exp(cdd_multi_fit_summary_pois[, "results"]),
   std.error = cdd_multi_fit_summary_pois[, "se"],
   conf.low = exp(cdd_multi_fit_summary_pois[, "results"] - 1.96 * cdd_multi_fit_summary_pois[, "se"]),
   conf.high = exp(cdd_multi_fit_summary_pois[, "results"] + 1.96 * cdd_multi_fit_summary_pois[, "se"]),
   p.value = 2 * pnorm(-abs(cdd_multi_fit_summary_pois[, "results"] / cdd_multi_fit_summary_pois[, "se"]))
) %>%
   filter(term != "(Intercept)")

cdd_multi_results_table_pois

##########===== Dispersion check for Poisson model

dispersion_checks_pois <- lapply(1:imp_df_c$m, function(i) {
   dat <- complete(imp_df_c, i)
   fit_cdd <- glm(
      cdd_score ~
         cc1_sex +
         cm3_diarhoea +
         cf12_exclusive +
         cf16_solid_food +
         usually_fast +
         education_level +
         occupation +
         urban_rural +
         hc2_age +
         hc3_family_size +
         hc9_gro_crops +
         age_cat +
         hfias_category_insecure,
      family = poisson(link = "log"),
      data = dat
   )
   ######=== Dispersion = residual deviance / df
   data.frame(
      imputation  = i,
      dispersion  = fit_cdd$deviance / fit_cdd$df.residual,
      deviance    = fit_cdd$deviance,
      df          = fit_cdd$df.residual
   )
})
dispersion_summary_pois <- do.call(rbind, dispersion_checks_pois)
print(dispersion_summary_pois)

#######======DHARMa residual diagnostics plots(cdd- Poisson)

sampled_imps_cdd_pois <- sample(seq_along(cdd_multi_fit_pois), 4)

# Create plots
for (i in sampled_imps_cdd_pois) {
   
   sim_res <- simulateResiduals(cdd_multi_fit_pois[[i]])
   
   png(
      filename = sprintf("8. Temp/cdd_pois_%02d.png", i),
      width = 1600,
      height = 1200,
      res = 150
   )
   
   plot(sim_res, main = "")
   title(main = paste("Imputation", i), line = 1)
   
   dev.off()
}

# Read only sampled images
files <- sprintf("8. Temp/cdd_pois_%02d.png", sampled_imps_cdd_pois)
imgs_pois <- image_read(files)
#####==== 2 × 2 layout
combined_cdd_pois <- image_montage(
   imgs_pois,
   tile = "2x2",
   geometry = "1600x1200+10+10"
)

image_write(
   combined_cdd_pois,
   path = "5. Plots/cdd_pois_dharma_2x2.png",
   format = "png"
)

#######====Underdispersion use CMP as before
#######==Clear under dispersion present: Try CMP model
cdd_multi_fit_cmp <- lapply(1:imp_df_c$m, function(i) {
   dat <- complete(imp_df, i)
   glmmTMB::glmmTMB(
      cdd_score ~
         cc1_sex +
         cm3_diarhoea +
         cf16_solid_food +
         usually_fast +
         education_level +
         occupation +
         urban_rural +
         hc2_age +
         hc3_family_size +
         hc9_gro_crops +
         age_cat +
         minimum_mdd+
         hfias_category_insecure,
      family = compois(link = "log"),
      data = dat
   )
})

#####≠=Pool using Rubin's Rules using model based vcovariance as CMP correctly models dispersion issues

pool_cdd_glmmtmb <- function(mods) {
   m <- length(mods)
   coefs <- lapply(mods, function(fit) fixef(fit)$cond)
   vcovs <- lapply(mods, function(fit) vcov(fit)$cond)
   Qbar <- Reduce("+", coefs) / m
   Ubar <- Reduce("+", vcovs) / m
   B    <- var(do.call(rbind, coefs))
   T_var <- Ubar + (1 + 1/m) * B
   se <- sqrt(diag(T_var))
   z  <- Qbar / se
   p  <- 2 * pnorm(abs(z), lower.tail = FALSE)
   data.frame(
      term      = names(Qbar),
      PR        = round(exp(Qbar), 3),
      conf.low  = round(exp(Qbar - 1.96 * se), 3),
      conf.high = round(exp(Qbar + 1.96 * se), 3),
      p.value   = round(p, 3)
   )
}
cdd_multi_results_cmp <- pool_cdd_glmmtmb(cdd_multi_fit_cmp)
print(cdd_multi_results_cmp)

###=== Assess for equidispersion
###

dispersion_cmp <- sapply(cdd_multi_fit_cmp, function(fit) {
   as.numeric(sigma(fit))
})
nu_cmp <- sapply(cdd_multi_fit_cmp, function(fit) {
   1 / sigma(fit)
})

#######======Diagnostic plots
set.seed(20250520)
sampled_imps_cdd_cmp <- sample(seq_along(cdd_multi_fit_pois),4)

# Create plots
for (i in sampled_imps_cdd_cmp) {
   
   sim_res <- simulateResiduals(cdd_multi_fit_cmp[[i]])
   
   png(
      filename = sprintf("8. Temp/cdd_pois_cmp_%02d.png", i),
      width = 1600,
      height = 1200,
      res = 150
   )
   
   plot(sim_res, main = "")
   title(main = paste("Imputation", i), line = 1)
   
   dev.off()
}

##read all temporary images
files_cddpois_cmp <- sprintf("8. Temp/cdd_pois_cmp_%02d.png",
                             sampled_imps_cdd_cmp)

imgs_pois_cmp <- image_read(files_cddpois_cmp)

#####==== 2 × 2 layout
combined_cdd_pois_cmp <- image_montage(
   imgs_pois_cmp,
   tile = "2x2",
   geometry = "1600x1200+10+10"
)

image_write(
   combined_cdd_pois_cmp,
   path = "5. Plots/cdd_cmp_dharma_2x2.png",
   format = "png"
)


######====Clean all and combine tables for tables Maternal diet diversity
uni_tbl <- mdd_uni_results_table_bin %>%
   mutate(
      OR_uni = sprintf(
         "%.3f (%.3f, %.3f)",
         estimate, conf.low, conf.high
      ),
      p_uni = sprintf("%.3f", p.value)
   ) %>%
   select(variable,term, OR_uni, p_uni)

multi_tbl <- mdd_multi_results_table_bin %>%
   mutate(
      OR_multi = sprintf(
         "%.3f(%.3f,%.3f)",
         estimate, conf.low, conf.high
      ),
      p_multi = sprintf("%.3f", p.value)
   ) %>%
   select(term, OR_multi, p_multi)


rr_uni_tbl <- mdd_uni_results_table_pois %>%
   mutate(
      RR_uni = sprintf("%.3f(%.3f,%.3f)", estimate, conf.low, conf.high),
      p_rr_uni = sprintf("%.3f", p.value)
   ) %>%
   select(term, RR_uni, p_rr_uni)

rr_multi_tbl <- mdd_multi_results_table_pois %>%
   mutate(
      RR_multi = sprintf("%.3f(%.3f,%.3f)", estimate, conf.low, conf.high),
      p_rr_multi = sprintf("%.3f", p.value)
   ) %>%
   select(term, RR_multi, p_rr_multi)

rr_multi_tbl_cmp <- mdd_multi_results_cmp %>%
   mutate(
      RR_multi_cmp = sprintf("%.3f(%.3f, %.3f)", PR, conf.low, conf.high),
      p_rr_multi_cmp = sprintf("%.3f", p.value)
   ) %>%
   select(term, RR_multi_cmp, p_rr_multi_cmp)

######----merge with binary
final_tbl_mdd <- uni_tbl %>%
   full_join(multi_tbl, by = "term") %>%
   full_join(rr_uni_tbl, by = "term") %>%
   full_join(rr_multi_tbl, by = "term") %>% 
   full_join(rr_multi_tbl_cmp, by = "term")


View(final_tbl_mdd)

######tidy up term and variable names for nice formatting
final_tbl_mdd_main <- final_tbl_mdd %>%
   select(variable, term, dplyr::contains("multi")) %>%
   dplyr::mutate(term = stringr::str_remove(term, paste0("^", variable)))

View(final_tbl_mdd_main)
library(glue)




###helper function to clean gt tables for LaTeX outputer for overleaf editing
clean_tables1 <- function(tbl_obj) {
   tbl_obj %>%
      gt::as_latex() %>%
      as.character() %>%
      stringr::str_remove_all("<span.*?>") %>%
      stringr::str_remove_all("</span>") %>%
      stringr::str_remove("\\\\begin\\{table\\}\\[[^]]*\\]\\n?") %>%
      stringr::str_remove("\\\\end\\{table\\}\\n?$")
}

mdd_mod_main <- final_tbl_mdd_main %>% 
   gt::gt() %>% 
   gt::cols_label(
      term = "Variable",
      OR_multi = "aOR (95% CI)", p_multi = "p-value",
      RR_multi = "aRR (95% CI)", p_rr_multi = "p-value",
      RR_multi_cmp = "aRR (95% CI)*", p_rr_multi_cmp = "p-value*"
   ) %>%gt::tab_spanner(label = "Logistic regression",
                        columns = c(OR_multi, p_multi)) %>%
   gt::tab_spanner(label = "Poisson regression",
                   columns = c(RR_multi, p_rr_multi, RR_multi_cmp, p_rr_multi_cmp)) %>%
   gt::opt_row_striping() %>%
   clean_tables1()

writeLines(
   as.character(mdd_mod_main),
   "4. Tables/objective_1_mdd_multi_models.tex"
)

#####MDD= Unimodels
final_tbl_mdd_uni <- final_tbl_mdd %>%
   select(variable, term, dplyr::contains("uni")) %>%
   dplyr::mutate(term = stringr::str_remove(term, paste0("^", variable)))
mdd_mod_uni <- final_tbl_mdd_uni %>% 
   gt::gt() %>% 
   gt::cols_label(
      term = "Variable",
      OR_uni = "OR (95% CI)", p_uni = "p-value",
      RR_uni = "PR (95% CI)", p_rr_uni = "p-value"
   ) %>%gt::tab_spanner(label = "Logistic regression",
                        columns = c(OR_uni, p_uni)) %>%
   gt::tab_spanner(label = "Poisson regression",
                   columns = c(RR_uni, p_rr_uni)) %>%
   gt::opt_row_striping() %>%
   clean_tables1()

writeLines(
   as.character(mdd_mod_uni),
   "4. Tables/objective_1_mdd_uni_models.tex"
)

#######child diet diversity
######====Clean all and combine tables for tables Maternal diet diversity
uni_tbl <- cdd_uni_results_table_bin %>%
   mutate(
      OR_uni = sprintf(
         "%.3f (%.3f, %.3f)",
         estimate, conf.low, conf.high
      ),
      p_uni = sprintf("%.3f", p.value)
   ) %>%
   select(variable,term, OR_uni, p_uni) %>% 
   filter(term != "cf12_exclusiveYes") %>% 
   filter(term != "hc7_private_agric_landYes")

multi_tbl <- cdd_multi_results_table_bin %>%
   mutate(
      OR_multi = sprintf(
         "%.3f (%.3f, %.3f)",
         estimate, conf.low, conf.high
      ),
      p_multi = sprintf("%.3f", p.value)
   ) %>%
   select(term, OR_multi, p_multi)


rr_uni_tbl <- cdd_uni_results_table_pois %>%
   mutate(
      RR_uni = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high),
      p_rr_uni = sprintf("%.3f", p.value)
   ) %>%
   select(term, RR_uni, p_rr_uni) %>% 
   filter(term != "cf12_exclusiveYes") %>% 
   filter(term != "hc7_private_agric_landYes")  %>% 
   unique()

rr_multi_tbl <- cdd_multi_results_table_pois %>%
   mutate(
      RR_multi = sprintf("%.3f (%.3f, %.3f)", estimate, conf.low, conf.high),
      p_rr_multi = sprintf("%.3f", p.value)
   ) %>%
   select(term, RR_multi, p_rr_multi)

rr_multi_tbl_cmp <- cdd_multi_results_cmp %>%
   mutate(
      RR_multi_cmp = sprintf("%.3f (%.3f, %.3f)", PR, conf.low, conf.high),
      p_rr_multi_cmp = sprintf("%.3f", p.value)
   ) %>%
   select(term, RR_multi_cmp, p_rr_multi_cmp) %>% 
   filter(term != "(Intercept)")

######----merge with binary then tidy up term and variable names for nice formatting
final_tbl_cdd <- uni_tbl %>%
   full_join(multi_tbl, by = "term") %>%
   full_join(rr_uni_tbl, by = "term") %>%
   full_join(rr_multi_tbl, by = "term") %>% 
   full_join(rr_multi_tbl_cmp, by = "term") %>% 
   filter(!is.na(variable)) %>% 
   dplyr::mutate(term = stringr::str_remove(term, paste0("^", variable)))
######Export main and multivariable models separately for better formatting in overleaf
cdd_mod_main <- final_tbl_cdd %>% 
   select(variable, term, dplyr::contains("multi")) %>%
   gt::gt() %>% 
   gt::cols_label(
      term = "Variable",
      OR_multi = "aOR (95% CI)", p_multi = "p-value",
      RR_multi = "aIRR (95% CI)", p_rr_multi = "p-value",
      RR_multi_cmp = "aIRR (95% CI)*", p_rr_multi_cmp = "p-value*"
   ) %>%gt::tab_spanner(label = "Logistic regression",
                        columns = c(OR_multi, p_multi)) %>%
   gt::tab_spanner(label = "Poisson regression",
                   columns = c(RR_multi, p_rr_multi, RR_multi_cmp, p_rr_multi_cmp)) %>%
   gt::opt_row_striping() %>%
   clean_tables1()

writeLines(
   as.character(cdd_mod_main),
   "4. Tables/objective_1_cdd_multi_modelsR.tex"
)


cdd_mod_uni <- final_tbl_cdd %>% 
   select(variable, term, dplyr::contains("uni")) %>% 
   gt::gt() %>% 
   gt::cols_label(
      term = "Variable",
      OR_uni = "OR (95% CI)", p_uni = "p-value",
      RR_uni = "PR (95% CI)", p_rr_uni = "p-value"
   ) %>%gt::tab_spanner(label = "Logistic regression",
                        columns = c(OR_uni, p_uni)) %>%
   gt::tab_spanner(label = "Poisson regression",
                   columns = c(RR_uni, p_rr_uni)) %>%
   gt::opt_row_striping() %>%
   clean_tables1()

writeLines(
   as.character(cdd_mod_uni),
   "4. Tables/objective_1_cdd_uni_modelsR.tex"
)

#######=== End of Objective One

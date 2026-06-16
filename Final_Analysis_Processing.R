###############################################################################
# Project: Final Analysis Preparation
# Purpose: Clean and standardise variables for final analysis dataset.
# Author: Edward
###############################################################################

# Sections:
# 1. Import data and prepare for analysis - convert to factors, recode variables, drop implausible z scores
# 2. Clean and Select Infant Feeding and Child/Maternal Diet Diversity Variables
# 3. Household food insecurity access scale (HFIAS)
# 4. WASH variables Standardisation
# 5. SES variables - education, occupation, land ownership, agricultural production
# 6 . Calculate Z scores.
# 7. Pre Analysis Data Set
# 8. Assessing clustering
###############################################################################
rm(list = ls())
######Packages
pacman::p_load(tidyverse, gtsummary, kableExtra,mice, mitools, sandwich, DHARMa, magick, ResourceSelection,
               glmmTMB, gt, quantreg, broom.mixed)


####===load dataset
load("2. Clean Data/analysis_data_clean.RData")
load("2. Clean Data/main_dictionary.RData")


#===============================================================================
#1. Import data and prepare for analysis - convert to factors, recode variables, drop implausible z scores
#####===convert into factor

factor_vars <- c("minimum_mdd", "marital_status", "education_level",
                 "md_oils", "md_sweets", "md_spices", "usually_fast", "md_illness",
                 "wash_drinking_water","wash_cooking_water","wash_toilets","wash_handwashing_place",
                 "wash_soap_type", "wash_child_stool_disposal", "hfias_category",
                 "hc7_private_agric_land", "hc9_gro_crops", "education_level", "occupation",
                 "stunting", "wasting", "underweight","wa4_treat_water","yesterday_fast","md27_tablet",
                 "cf12_exclusive","cf3_colustrum","cf9_breastfed_yesterday","cc7_early_9months",
                 "cf7_currnt_breastfeeding","usually_fast","cf16_solid_food")

analysis_data_clean1 <- analysis_data_clean %>%
   mutate(across(all_of(factor_vars), as.factor)) %>% 
   mutate(across(where(is.numeric),
                 ~ replace(.x, .x == 98, NA))) %>%
   mutate(across(where(is.character),
                 ~ replace(.x, .x == "DK", NA))) %>% #not working
   mutate(across(where(is.factor),
                 ~ factor(na_if(as.character(.x), "NA")))) %>% 
   mutate(across(
      where(is.factor),
      ~ forcats::fct_recode(.x, NULL = "DK")
   )) %>% 
   mutate(
      occupation = fct_collapse(
         occupation,
         "Unemployed/Small business" = c(
            "No Job",
            "Small business"
         )
      )
   ) %>% filter(
      (zwfl > -5 & zwfl < 5) | is.na(zwfl),
      (zlen > -6 & zlen < 6) | is.na(zlen),
      (zwei > -5 & zwei < 5) | is.na(zwei)
   ) #### drop implausible z scores  leaving final sample of n = 789
   
####==update with rural urban indicator
analysis_data_clean1 <- analysis_data_clean1 %>%
mutate(urban_rural = factor(
   if_else(kebele.x == "K04", "urban", "rural"),
   levels = c("urban", "rural")  # urban as reference
))

####Reorder min dd
analysis_data_clean1$minimum_mdd <- factor(
   analysis_data_clean1$minimum_mdd,
   levels = c("LDD", "ADD")
)

analysis_data_clean1$minimum_cdd <- factor(
   analysis_data_clean1$minimum_cdd,
   levels = c("LDD", "ADD")
)


###== format child vars(Exclusive breasfeeding babies for child diet diversity analysis)
analysis_data_clean1 <- analysis_data_clean1 %>%
   mutate(
      cdd_score = as.numeric(cdd_score),
      cdd_score = if_else(cc4_age_months < 6, NA_real_, cdd_score),
      minimum_cdd = if_else(cc4_age_months < 6, "", minimum_cdd)
   )

###== based on linear relationship
analysis_data_clean1$age_cat <- cut(
   analysis_data_clean1$cc4_age_months,
   breaks = c(0, 6, 12, 24),
   right = FALSE,
   labels = c("0-5", "6-11", "12-24")
)

###=== Combine HFIAS based on discussion with Dr. Kokeb
analysis_data_clean1 <- analysis_data_clean1 %>%
   mutate(
      hfias_category_insecure = case_when(
         hfias_category %in% c(
            "Mildly food insecure",
            "Moderately food insecure",
            "Severely food insecure"
         ) ~ "Yes",
         
         hfias_category == "Food secure" ~ "No",
         
         TRUE ~ NA_character_
      )
   )  %>% 
   mutate(
      birth_size_combine = case_when(
         cc6_size_born %in% c("Smaller than average", "Very small") ~ "Small",
         cc6_size_born %in% c("Larger than average", "Very large") ~ "Large",
         cc6_size_born == "Average" ~ "Average",
         TRUE ~ NA_character_
      ),
      birth_size_combine = factor(birth_size_combine, levels = c("Small", "Average", "Large"))
   ) %>% mutate(
   anaemia = case_when(
      hgb_child < 11 ~ "Yes",
      hgb_child >= 11 ~ "No",
      TRUE ~ NA_character_
   )
   )
#11 patients in LDD category, 0 in ADD

#####Multiple imputations for missing data
# Replace don't know with NA
#imp <- mice(analysis_data_clean1, seed = 123)
#md.pattern(analysis_data_clean1)
#colSums(is.na(analysis_data_clean1))
#md.pairs(analysis_data_clean1)

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

meths <- make.method(analysis_data_clean1)
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
pred <- make.predictorMatrix(analysis_data_clean1)
pred[ignore_vars_pred, ] <- 0   
pred[, ignore_vars_pred] <- 0   

######===run MI
imp_df <- mice(analysis_data_clean1, 
            method = meths, seed = 20250520,
            predictorMatrix = pred,
            m = 20)

###= Recategorise the variables
completed_df <- complete(imp_df, "long", include = TRUE) %>% 
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


imp_df_updated <- as.mids(completed_df)


#########################################################################################
#########################################################################################
#########################################################################################
##############
##############Descriptive Statictisc
##############
##########################################################################################
##########################################################################################
##########################################################################################
##########################################################################################

pacman::p_load("gtsummary","tidyverse")

######==== user defined function to clean tables 
clean_tables <- function(tbl_obj) {
   tbl_obj %>%
      as_kable_extra(
         format = "latex",
         booktabs = TRUE
      ) %>%
      kable_styling(
         latex_options = c("hold_position", "striped")
      ) %>%
      as.character() %>%
      stringr::str_replace_all("<span.*?>", "") %>%
      stringr::str_replace_all("</span>", "") %>%
      stringr::str_replace(
         "\\\\begin\\{table\\}\\[!h\\]",
         "\\\\begin\\{table\\}[htbp]"
      ) %>%
      str_remove_all("\\\\begin\\{table\\}\\[.*?\\]") %>%
      str_remove_all("\\\\end\\{table\\}")   
}

####==== summary statistics
show_one_level <- list(
   wash_child_stool_disposal ~ "Unsafe disposal",
   urban_rural ~ "rural",
   wash_drinking_water ~ "Unimproved",
   wash_cooking_water ~ "Unimproved",
   wash_toilets ~ "Unimproved",
   minimum_mdd ~ "LDD"
)

non_child_desc_table <- analysis_data_clean1 %>%
   mutate(rural_urban = as.factor(urban_rural)) %>%
   select(
      hc2_age,
      marital_status,
      muac_mother,
      mdd_score,
      minimum_mdd,
      weight_mother,
      usually_fast,
      md_illness,
      hfias_category,
      education_level,
      occupation,
      hc9_gro_crops,
      wa4_treat_water,
      wash_toilets,
      wash_child_stool_disposal,
      wash_drinking_water,
      wash_cooking_water,
      urban_rural
   ) %>%
   tbl_summary(
      by = NULL,
      missing = "no",
      value = show_one_level,
      statistic = list(
         all_categorical() ~ "{n} ({p}%)",
         all_continuous() ~ "{median} ({p25}, {p75})"
      ),
      digits = list(
         all_continuous() ~ 1,
         all_categorical() ~ c(0, 1)
      ),
      label = list(
         hc2_age ~ "Maternal age (years)",
         marital_status ~ "Marital status",
         hfias_category ~ "Household food insecurity category",
         wash_drinking_water ~ "Source of drinking water : Unimproved",
         wash_cooking_water ~ "Source of cooking water: Unimproved",
         wash_toilets ~ "Type of toilet facility: Unimproved",
         mdd_score ~ "Maternal dietary diversity score",
         minimum_mdd ~ "MDD-W: No",
         weight_mother ~ "Maternal weight (kg)",
         muac_mother ~ "Maternal MUAC (cm)",
         education_level ~ "Maternal education level",
         usually_fast ~ "Usually fasts: Yes",
         md_illness ~ "Maternal illness: Yes",
         hc9_gro_crops ~ "Household grows crops: Yes",
         wa4_treat_water ~ "Treats drinking water: Yes",
         wash_child_stool_disposal ~ "Child stool disposal practice: Unsafe",
         urban_rural ~ "Residence: Rural",
         occupation ~ "Maternal occupation"
      )
   ) %>%
   modify_header(label ~ "") %>%
   bold_labels()

non_child_desc_table

library(kableExtra)

latex_code_non_child <- non_child_desc_table %>%
   as_kable_extra(
      format = "latex",
      booktabs = TRUE
   ) %>%
   kable_styling(
      latex_options = c("hold_position", "striped")
   ) %>% 
   str_replace_all("<span.*?>", "") %>%
   str_replace_all("</span>", "") %>% 
   stringr::str_replace(
      "\\\\begin\\{table\\}\\[!h\\]",
      "\\\\begin\\{table\\}[htbp]"
   ) %>% 
   str_remove_all("\\\\begin\\{table\\}\\[.*?\\]") %>%
   str_remove_all("\\\\end\\{table\\}") 

# Export as .tex file
writeLines(
   as.character(latex_code_non_child),
   "4. Tables/2. non_child_desc_table.tex"
)

show_one_level_child <- list(
   cf12_exclusive ~ "Yes",
   stunting ~ "Yes",
   wasting ~ "Yes",
   underweight ~ "Yes",
   cm3_diarhoea ~ "Yes",
   cf1_breastfed ~ "Yes",
   anaemia ~ "Yes",
   cc1_sex ~ "Female",
   minimum_cdd ~ "LDD"
)

child_desc_table <- analysis_data_clean1 %>%
   mutate(
      cf1_breastfed = case_when(
         cf1_breastfed == 1 ~ "Yes",
         cf1_breastfed == 2 ~ "No",
         TRUE ~ NA_character_
      ),
      minimum_cdd = as.factor(minimum_cdd),
      cdd_score = as.numeric(cdd_score))%>% 
   select(
      cc1_sex,
      cc4_age_months,
      cc5_place_birth,
      cc6_size_born,
      cf1_breastfed,
      cf7_currnt_breastfeeding,
      cf12_exclusive,
      cf16_solid_food,
      cdd_score,
      minimum_cdd,
      cm3_diarhoea,
      weight_child,
      length_child,
      muac_child,
      hgb_child,
      zlen,
      zwfl,
      zwei,
      stunting,
      wasting,
      underweight,
      anaemia
   ) %>%
   tbl_summary(
      by = NULL,
      type = list(
         cdd_score ~ "continuous2",
         zlen ~ "continuous",
         zwfl ~ "continuous",
         zwei ~ "continuous"
      ),
      statistic = list(
         all_continuous() ~ "{median} ({p25}, {p75})",
         zlen ~ "{mean} ({sd})",
         zwfl ~ "{mean} ({sd})",
         zwei ~ "{mean} ({sd})",
         cdd_score ~ "{median} ({p25}, {p75})",
         all_categorical() ~ "{n} ({p}%)"
      ),
      missing = "no",
      value = show_one_level_child,
      digits = list(
         all_continuous() ~ 1,
         all_categorical() ~ c(0, 1)
      ),
      
      label = list(
         cc1_sex ~ "Child sex: Female",
         cc4_age_months ~ "Child age (months)",
         cc5_place_birth ~ "Place of birth",
         cc6_size_born ~ "Perceived size at birth",
         cf1_breastfed ~ "Ever breastfed: Yes",
         cf7_currnt_breastfeeding ~ "Currently breastfeeding: Yes",
         cf12_exclusive ~ "Exclusive breastfeeding: Yes",
         cf16_solid_food ~ "Introduced solid foods: Yes",
         cdd_score ~ "Child dietary diversity score",
         minimum_cdd ~ "MDD for children 6-23 months: No",
         cm3_diarhoea ~ "Diarrhoea in past 2 weeks: Yes",
         weight_child ~ "Child weight (kg)",
         length_child ~ "Child length (cm)",
         muac_child ~ "Child MUAC (cm)",
         hgb_child ~ "Child haemoglobin (g/dL)",
         zlen ~ "Length-for-age Z-score",
         zwfl ~ "Weight-for-length Z-score",
         zwei ~ "Weight-for-age Z-score",
         stunting ~ "Stunting: Yes",
         wasting ~ "Wasting: Yes",
         underweight ~ "Underweight: Yes",
         anaemia ~ "Anaemia: Yes"
      )
   ) %>%
   modify_header(label ~ "**Child characteristics**") %>%
   bold_labels()

child_desc_table

library(kableExtra)

latex_code <- child_desc_table %>%
   as_kable_extra(
      format = "latex",
      booktabs = TRUE
   ) %>%
   kable_styling(
      latex_options = c("hold_position", "striped")
   ) %>% as.character() %>%
   str_replace_all("<span.*?>", "") %>%
   str_replace_all("</span>", "") %>% 
   stringr::str_replace(
      "\\\\begin\\{table\\}\\[!h\\]",
      "\\\\begin\\{table\\}[htbp]"
   ) %>% 
   str_remove_all("\\\\begin\\{table\\}\\[.*?\\]") %>%
   str_remove_all("\\\\end\\{table\\}")                 

# Export as .tex file
writeLines(
   as.character(latex_code),
   "4. Tables/1.child_desc_table.tex"
)




########=== Check functional form of age 
df_long <- analysis_data_clean1 %>%
   select(cc4_age_months, zwfl, zwei, zlen) %>%
   pivot_longer(cols = c(zwfl, zwei, zlen),
                names_to = "indicator",
                values_to = "zscore")
plotz <- ggplot(df_long, aes(x = cc4_age_months, y = zscore, color = indicator)) +
   geom_smooth(method = "loess", se = TRUE, linewidth = 1.2) +
   scale_color_manual(
      values = c(
         "zlen" = "#6A5ACD",
         "zwei" = "#D95F02",
         "zwfl" = "#1B9E77"
      ),
      labels = c(
         "zlen" = "LAZ",
         "zwei" = "WAZ",
         "zwfl" = "WLZ"
      )
   ) +
   labs(
      x = "Age (months)",
      y = "Z-score",
      color = NULL
   ) +
   theme_classic(base_size = 13) +
   theme(
      plot.title = element_text(face = "bold"),
      legend.position = "bottom"
   )

ggsave(
   filename = "5. Plots/loess_growth_curves.png",
   plot = plotz,
   width = 10,
   height = 6,
   dpi = 300
)

#####===== assessing clustering
############======Assessing clustering
library(lme4)

m_zlen <- lmer(zlen ~ 1 + (1 | kebele.x), data = analysis_data_clean1)

summary(m_zlen)

# ICC
var_kebele <- as.data.frame(VarCorr(m_zlen))$vcov[1]
var_resid  <- attr(VarCorr(m_zlen), "sc")^2

ICC_laz <- var_kebele / (var_kebele + var_resid)
ICC_laz


library(lme4)

outcomes <- c("LAZ", "WAZ", "WLZ")

icc_results <- data.frame(
   outcome = outcomes,
   ICC = NA_real_
)

for(i in seq_along(outcomes)) {
   
   form <- as.formula(
      paste0(outcomes[i], " ~ 1 + (1 | kebele)")
   )
   
   model <- lmer(form, data = analysis_data_clean1)
   
   var_kebele <- as.data.frame(VarCorr(model))$vcov[1]
   var_resid  <- attr(VarCorr(model), "sc")^2
   
   icc_results$ICC[i] <-
      var_kebele / (var_kebele + var_resid)
}

icc_results




library(performance)

for(y in c("zlen","zwfl","zwei")) {
   
   m <- lmer(
      as.formula(paste0(y, " ~ 1 + (1 | kebele.x)")),
      data = analysis_data_clean1
   )
   
   cat("\n", y, "\n")
   print(icc(m))
}

library(lme4)
library(performance)

get_icc <- function(outcome) {
   
   model <- glmer(
      as.formula(
         paste0(outcome, " ~ 1 + (1 | kebele.x)")
      ),
      family = poisson,
      data = analysis_data_clean1
   )
   
   out <- performance::icc(model)
   
   data.frame(
      outcome = outcome,
      ICC = out$ICC_adjusted
   )
}
do.call(
   rbind,
   lapply(c("cdd_score", "mdd_score"), get_icc)
)
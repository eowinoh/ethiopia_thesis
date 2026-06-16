###############################################################################
# Project: Data Preprocessing and Standardisation Pipeline
# Purpose: Clean raw data, filter analytic variables, subset target variables,
#          and compute anthropometric Z-scores for analysis.
# Author: Edward
###############################################################################

# Sections:
# 1. Import Raw data, disctionaries and append raw files
# 2. Clean and Select Infant Feeding and Child/Maternal Diet Diversity Variables
# 3. Household food insecurity access scale (HFIAS)
# 4. WASH variables Standardisation
# 5. SES variables - education, occupation, land ownership, agricultural production
# 6 . Calculate Z scores.
# 7. Pre Analysis Data Set
###############################################################################
rm(list = ls())

library(tidyverse)
library(here)
library(readxl)
library(gt)
library(gtsummary)

#####1. Import Raw data, disctionaries and append raw files
###== Load dictionary
main_dictionary <- read_excel("1. Raw Data/Follow-up_data_cleaning.xlsx", 
                                       sheet = "data_dictionary") %>% 
   filter(!is.na(clean_name))
#####== load clean up and hgb files
readxl::excel_sheets("1. Raw Data/missing_ids_report_v2.xlsx")
manually_matched <- read_excel("1. Raw Data/missing_ids_report_v2.xlsx", 
                                       sheet = "master_ids_lookup")
sup_matched <- read_excel("1. Raw Data/missing_ids_report_v2.xlsx", 
                                       sheet = "ver_sup") %>% 
   filter(!is.na(anthro_id))
########
follow_up_data_cleaning <- read_excel("1. Raw Data/Follow-up_data_cleaning.xlsx", 
                                             sheet = "visit1_all cleaning") %>% 
   # drop figured duplicates using uuid
   filter(!`_uuid` %in% c(
      "21bd71e4-718e-4792-bac0-234db23036b4",
      "f53b04f3-8786-4d92-ab42-0ded3b010cb2",
      "1e8b7afa-1b01-462a-bd55-e07be26ec352",
      "715c5710-507e-49f2-9058-579ab43f81d0",
      "31f58209-7082-4103-bd34-7cb5034e84eb",
      "ba3b26bb-9a8e-4458-8757-8dc47564fd53",
      "a486ba22-e6a3-4f7f-b756-0beb3f461bbf",
      "3eac392b-3939-4206-9a0a-1f7fbe8c98bf",
      "7b647926-716c-4eaa-9284-f31bf1b4819c",
      "e6eac417-d397-43d0-b4a5-527843f28e46"
   ))


### raw anthro hgb file
anthro_and_hgb_data <- read_excel("1. Raw Data/Anthro and hgb data.xlsx", 
                                       sheet = "last checked") %>% 
   rename(household_id = `Household ID`)


######===== clean from hgb data
hgb_ids <- anthro_and_hgb_data %>% 
            select(kebele, household_id) %>% 
   # Clean household_id
   mutate(
      household_id = na_if(trimws(household_id), ""),
      clean_id = if_else(
         !is.na(household_id) & str_detect(household_id, "\\d"),
         paste0(
            format(round(parse_number(household_id), 1), nsmall = 1),
            tolower(coalesce(str_extract(household_id, "[A-Za-z]+$"), ""))
         ),
         household_id
      )
   ) %>% 
   # Fix kebele codes to match with kebeles from main data
   mutate(
      kebele = case_when(
         kebele == "11" ~ "011",
         kebele == "5"  ~ "005",
         kebele == "7"  ~ "007",
         kebele == "8"  ~ "008",
         TRUE ~ kebele
      )
   ) %>% 
   # Remove .0 at end of clean_id and build final key
   mutate(
      clean_id = str_replace(clean_id, "\\.0$", ""),
      ke_hhid  = paste(kebele, clean_id, sep = "/"),
      ke_hhid  = str_replace_all(ke_hhid, " ", "")
   )



all_match <- follow_up_data_cleaning %>% 
   select(kebele, HHID,`Ke/HHID`) %>%
   inner_join(hgb_ids, by = c("Ke/HHID" = "ke_hhid")) %>% 
   select(`Ke/HHID`) %>% 
   mutate(match_status = "Match",
          main_id = `Ke/HHID`,
          anthro_id = `Ke/HHID`) %>% 
   select(main_id, anthro_id, match_status)


########======
cleaned_ids <- manually_matched %>%
   mutate(match_status = "Reviewed",
          main_id = `Ke/HHID`,
          anthro_id = hgb_ids) %>% 
   select(main_id, anthro_id, match_status)

######====== supervisor matched ids

sup_ids <- sup_matched %>% 
   mutate(match_status = "Supervisor Matched",
          main_id = subject_id,
          anthro_id = anthro_id) %>% 
   select(main_id, anthro_id, match_status)

######====combine the matched and reviewed ids
final_id_mapping <- bind_rows(all_match, cleaned_ids,sup_ids) %>% 
   unique()

missing <- hgb_ids %>% 
   anti_join(final_id_mapping, by = c("ke_hhid" = "anthro_id")) %>% 
   select(ke_hhid) %>% 
   mutate(match_status = "Missing",
          main_id = NA_character_,
          anthro_id = ke_hhid) %>% 
   select(main_id, anthro_id, match_status) %>% 
   unique()

write_csv(missing, "missing_anthro_ids.csv")

#####==merge final mapping with main data and anthro data to check for duplicates

anthro_and_hgb_data_with_ids <- anthro_and_hgb_data %>% 
   # Clean household_id
   mutate(
      household_id = na_if(trimws(household_id), ""),
      clean_id = if_else(
         !is.na(household_id) & str_detect(household_id, "\\d"),
         paste0(
            format(round(parse_number(household_id), 1), nsmall = 1),
            tolower(coalesce(str_extract(household_id, "[A-Za-z]+$"), ""))
         ),
         household_id
      )
   ) %>% 
   # Fix kebele codes to match with kebeles from main data
   mutate(
      kebele = case_when(
         kebele == "11" ~ "011",
         kebele == "5"  ~ "005",
         kebele == "7"  ~ "007",
         kebele == "8"  ~ "008",
         TRUE ~ kebele
      )
   ) %>% 
   # Remove .0 at end of clean_id and build final key
   mutate(
      clean_id = str_replace(clean_id, "\\.0$", ""),
      ke_hhid  = paste(kebele, clean_id, sep = "/"),
      ke_hhid  = str_replace_all(ke_hhid, " ", "")
   ) %>% 
   #merge with final id mapping to get main ids
   left_join(final_id_mapping, by = c("ke_hhid" = "anthro_id")) %>% 
   select(-c("...9"))

 dup_anthro <- final_id_mapping %>%
    count(anthro_id) %>%
    filter(n > 1)

######-=====append with main data

follow_up_data_anthro_hgb <- follow_up_data_cleaning %>% 
   left_join(anthro_and_hgb_data_with_ids, by = c("Ke/HHID" = "main_id")) %>% 
   filter(`Ke/HHID` != "NA") %>% 
   mutate(length_child = length,
          weight_child = wt_C.y,
          muac_child = MUAC_C.y,
          hgb_child = HGB_C.y,
          weight_mother = wt_M.y,
          muac_mother = MUAC_M.y) %>% 
   select(-c(length, wt_C.y, MUAC_C.y, HGB_C.y, wt_M.y, MUAC_M.y,
             wt_C.x, MUAC_C.x, HGB_C.x, wt_M.x, MUAC_M.x, length_C)) %>% 
   filter(!is.na(match_status)) %>% 
   distinct()

####=== drop dups due to dups in anthro files
####== Confirmed from dob, weigh, sex and age between the two files
follow_up_data_anthro_hgb <- follow_up_data_anthro_hgb %>%
   mutate(
      clean_id2 = trimws(as.character(clean_id)),
      length_child2 = suppressWarnings(as.numeric(trimws(as.character(length_child)))),
      hgb_child2 = trimws(as.character(hgb_child)),
      drop_row = case_when(
         clean_id2 == "234.4at"  & near(length_child2, 73.1) ~ 1L,
         clean_id2 == "234.4"   & near(length_child2, 69.5) ~ 1L,
         clean_id2 == "129.7"   & near(length_child2, 71.5) ~ 1L,
         clean_id2 == "56.4"    & near(length_child2, 65.1) ~ 1L,
         clean_id2 == "12.1"    & is.na(length_child2) ~ 1L,
         clean_id2 == "612.0at" & near(length_child2, 61.8) ~ 1L,
         clean_id2 == "612.0at" & near(length_child2, 70.0) ~ 1L,
         clean_id2 == "612.0at" & str_detect(hgb_child2, fixed("11.3/9.8")) ~ 1L,
         clean_id2 == "544.1" & near(length_child2,73.3) ~ 1L,
         TRUE ~ 0L
      )
   ) %>% 
   filter(drop_row == 0L)  %>% 
   select(-c(clean_id2, length_child2, hgb_child2, drop_row))


####===818 records.

######==== translate crops into English
#####=== crop look up file

lookup <- c(
   
   # --- NONE / CODES ---
   "0"="none","00"="none","90"="none","98"="none","99"="none",
   "no"="none","o"="none","oo"="none",
   "yelem"="none","y elem"="none",
   "የለም"="none","የለመም"="none","የለሞ"="none","የለምየለም"="none","ዐዐ"="none","የከም"="none","ለ99"="none",
   
   # --- MAIZE ---
   "maize"="maize",
   "በቆ"="maize","በቆሎ"="maize","በቄለ"="maize","በቄላ"="maize",
   "ቦሎቆ"="maize","ቦቆሎ"="maize","ቦቀሎ"="maize",
   
   # --- TEFF ---
   "tef"="teff","teff"="teff","tff"="teff","tafe"="teff","tafa"="teff","taf"="teff","tefe"="teff","tefee"="teff",
   "ቴፈ"="teff","ቴፉ"="teff","ቴፊ"="teff","ቴፍ"="teff",
   "ጠፍ"="teff","ጤ  ፈ"="teff","ጤፉ"="teff","ጤፊ"="teff","ጤፋ"="teff","ጤፍ"="teff","ጦፈ"="teff",
   
   # --- WHEAT ---
   "wheat"="wheat",
   "sede"="wheat","sade"="wheat","send"="wheat","senda"="wheat",
   "sende"="wheat","sendie"="wheat","snde"="wheat","sineda99"="wheat","sinde"="wheat",
   "ሰነደ"="wheat","ሰነዴ"="wheat","ሰንዴ"="wheat","ሲንዴ"="wheat","ስንዴ"="wheat",
   
   # --- SORGHUM ---
   "sorghum"="sorghum","masela"="sorghum","mashella"="sorghum",
   "mashila"="sorghum","mashla"="sorghum","mashl"="sorghum",
   "መሸለ"="sorghum","መሽለ"="sorghum","መሽ"="sorghum",
   "ማሺላ"="sorghum","ማሽለ"="sorghum","ማሽላ"="sorghum",
   "dagusa"="sorghum","ደጉሳ"="sorghum",
   
   # --- BARLEY ---
   "barley"="barley","barly"="barley","barely"="barley",
   "gabs"="barley","gabse"="barley","gbese"="barley","gbs"="barley","gbse"="barley",
   "gebes"="barley","gebs"="barley",
   "ገቢስ"="barley","ገብስ"="barley","ገብሶ"="barley","ገበስ"="barley","ገብሲ"="barley",
   "ጋብስ"="barley","ጋበስ"="barley",
   
   # --- HARICOT BEANS ---
   "haricotbean"="haricot_beans","haricot bean"="haricot_beans",
   "boleka"="haricot_beans","boleke"="haricot_beans","bolke"="haricot_beans",
   "boloke"="haricot_beans","bolokie"="haricot_beans","bolaki"="haricot_beans","bolka"="haricot_beans","bolok"="haricot_beans",
   "fasoliya"="haricot_beans",
   "በለቄ"="haricot_beans","በሎቄ"="haricot_beans",
   "ቦለቂ"="haricot_beans","ቦለቄ"="haricot_beans","ቦሌቄ"="haricot_beans","ቦሎቄ"="haricot_beans","ቦላቄ"="haricot_beans",
   "ባቆላ"="haricot_beans",
   
   # --- FAVA BEANS ---
   "bakela"="fava_beans","bakkela"="fava_beans",
   "broadbean"="fava_beans","broad bean"="fava_beans","broadbeans"="fava_beans","broadbeab"="fava_beans",
   "ባቄላ"="fava_beans",
   
   # --- PEAS ---
   "pea"="peas","ater"="peas","atere"="peas","atre"="peas",
   "አተር"="peas","ኣተር"="peas",
   
   # --- LENTIL ---
   "lentil"="lentil","mesr"="lentil",
   
   # --- CHICKPEA ---
   "chickpea"="chickpea","shimbra"="chickpea",
   
   # --- ENSET ---
   "enset"="enset","ensat"="enset","eset"="enset",
   "እንሰት"="enset","እንሠት"="enset",
   
   # --- NIGER SEED ---
   "zengada"="niger_seed","zangda"="niger_seed",
   "ዘንጋደ"="niger_seed","ዘንጋዳ"="niger_seed","ዘንገደ"="niger_seed","ዘጋደ"="niger_seed",
   
   # --- CHILI ---
   "berbere"="chili","barbare"="chili","barebra"="chili","berebra"="chili","berebrea"="chili",
   "በረበሬ"="chili","በርበሬ"="chili",
   
   # --- KHAT ---
   "chat"="khat","chate"="khat","chaty"="khat","chahte"="khat",
   "ጫት"="khat",
   
   # --- CABBAGE ---
   "gomen"="cabbage","goman"="cabbage","cabbage"="cabbage","ጎመን"="cabbage",
   
   # --- OTHER ---
   "zaf"="other","lemate"="other","semer"="other","ዱባ"="other",
   "shnkora"="sugarcane"
)

####== write fxn to apply look up to columns
clean_crop <- function(x) {
   x <- tolower(trimws(x))
   x <- gsub("\\s+", " ", x)
   out <- lookup[x]
   out[is.na(out)] <- "other"
   return(out)
}

####== cleancols
cols_to_clean <- names(follow_up_data_anthro_hgb)[startsWith(names(follow_up_data_anthro_hgb), "crop/hc10")]

follow_up_data_anthro_hgb[cols_to_clean] <- lapply(follow_up_data_anthro_hgb[cols_to_clean], clean_crop)


##### Drop furtehr varisables that are not needed for analysi

follow_up_data_anthro_hgb <- follow_up_data_anthro_hgb %>% 
   select(-c("iycf/child name","iycf/note_2",`_validation_status`,`_index`))


#####====label variables
rename_map <- setNames(main_dictionary$raw_name, main_dictionary$clean_name)

####== rename variables

follow_up_data_anthro_hgb <- follow_up_data_anthro_hgb %>% 
   rename(any_of(rename_map))


####====Save only matched files
analysis_data <- follow_up_data_anthro_hgb %>% 
   filter(!is.na(match_status)) %>% 
   ### remove dups.....
   distinct(subject_id, "_uuid", .keep_all = TRUE) 

unmatched_rows <- follow_up_data_anthro_hgb %>% 
   filter(is.na(match_status)) %>% 
   select(subject_id,date,hhid, match_status)

#### Export unmatched rows for reviewt u
write_csv(unmatched_rows, "unmatched_rows.csv")

table(follow_up_data_anthro_hgb$match_status, useNA = "ifany")

####===check duplicate files
 dup_main <- analysis_data %>%
    count(subject_id) %>%
    filter(n > 1)
dup_data <- analysis_data %>% 
   filter(subject_id %in% dup_main$subject_id) %>% 
   arrange(subject_id) %>% 
   select(subject_id, match_status,"_uuid","_submission_time" )


####Tidy up evnvironment
rm(list = setdiff(ls(), c("follow_up_data_anthro_hgb", "main_dictionary","analysis_data")))
#####==========2. Clean and Select Infant Feeding and Child/Maternal Diet Diversity Variables
######============================================================
# Clean Child and Diet Diversity Variables
#============================================================
#####==========

#####===user defined fxn to help with labelling columns with yes/no and yes/no/dk

recode_yesno <- function(x) {
   dplyr::case_when(
      x == 1 ~ "Yes",   
      x == 0 ~ "No"
   )
}
recode_yesnodk <- function(x) {
   dplyr::case_when(
      x == 1 ~ "Yes",   
      x == 2 ~ "No",    
      x == 3 ~ "DK"
   )
}
#===============================================================================
######===child diet diversity variables
#===============================================================================
cdd_vars <- analysis_data %>% 
   select(kebele.x,hhid,subject_id,cf9_breastfed_yesterday,
          starts_with("cf12"), starts_with("cf14"), starts_with("cf15"),cf6_animal_milk) %>% 
   mutate( #1. Breast milk
      cdd_breastmilk = if_else(cf9_breastfed_yesterday == "1", 1, 0)
      %>% as.integer(),
      # 2. Grains, roots, tubers
      cdd_grains_roots_tubers = if_any(
         c(cf15_cf15a_injera, cf15_cf15c_tubers),
         ~ . == "1"
      ) %>% as.integer(),
      
      # 3. Vitamin A rich fruits, vegetables
      cdd_vitamin_a = if_any(
         c(cf15_cf15b_pumpkins, cf15_cf15e_fruits),
         ~ . == 1
      ) %>% as.integer(),
      
      # 4. Other fruits, vegetables
      cdd_other_fv = if_any(
         c(cf15_cf15f_other_fruits, cf15_cf15d_veg),
         ~ . == "1"
      ) %>% as.integer(),
      
      # 4. Flesh/Meaty foods (combination organ meat, meat, fish)
      cdd_flesh_foods = if_any(
         c(cf15_cf15g_liver, cf15_cf15h_any_meat, cf15_cf15j_fish),
         ~ . == "1"
      ) %>% as.integer(),
      
      # 6. Eggs
      cdd_eggs = if_else(cf15_cf15i_egg == "1", 1, 0),
      
      # 6. Legumes & nuts
      cdd_legumes_nuts = if_else(cf15_cf15k_legumes == "1", 1, 0),
      
      # 7. Dairy products
      cf6_animal_milk = as.numeric(cf6_animal_milk),
      cf14c_other_milk = as.numeric(cf14c_other_milk),
      cf14h_yogurt = as.numeric(cf14h_yogurt),
      cf14b_infant_formula = as.numeric(cf14b_infant_formula),
      cdd_dairy = if_any(
         c(cf15_cf15l_other_milk,cf14b_infant_formula,cf14c_other_milk,
           cf14h_yogurt,cf6_animal_milk), ~ . == 1
      )%>% as.integer()) %>%  
   mutate(
      cdd_dairy = replace_na(cdd_dairy, 0), # cow milk introducing NAs, but we want to count them as 0s for the score
      cdd_score = rowSums(select(., starts_with("cdd_")), na.rm = TRUE),
      minimum_cdd = if_else(cdd_score < 5, "LDD", "ADD")
   )  %>% 
   select(-c(cf9_breastfed_yesterday, starts_with("cf12"), starts_with("cf14"), 
             starts_with("cf15"),cf6_animal_milk)) %>% 
   mutate(across(c(cdd_grains_roots_tubers, cdd_vitamin_a,cdd_other_fv,cdd_flesh_foods,
                   cdd_eggs,cdd_legumes_nuts,cdd_dairy), recode_yesno)) %>% 
   select(kebele.x,hhid,subject_id, starts_with("cdd"), minimum_cdd)



other_cddvars <- analysis_data %>%
   select(kebele.x,hhid,subject_id, starts_with("cc"),starts_with("cf"),cm3_diarhoea,cf13a) %>% 
   select(-c(starts_with("cf15"), cc2_mother_care,cc3_dob,cf14b_infant_formula,cf14c_other_milk,
             cf14h_yogurt,cf6_animal_milk)) %>%  # cc2_mother_care all mothers apart from 1
   mutate(cc4_age_months = as.numeric(cc4_age_months),
          cc4_age_months = if_else(cc4_age_months == 32, NA_real_, cc4_age_months),
          cc1_sex = ifelse(cc1_sex == 1,"Male","Female"),
          cc5_place_birth = case_when(cc5_place_birth == 1 ~ "Govt facility",
                                      cc5_place_birth == 2 ~ "Not for profit facility",
                                      cc5_place_birth == 3 ~ "Private facilty",
                                      cc5_place_birth == 4 ~ "Health Post",
                                      cc5_place_birth == 5 ~ "Traditional birth attendant",
                                      cc5_place_birth == 6 ~ "Home",
                                      cc5_place_birth == "other" ~ "Other"),
          cc6_size_born = case_when(cc6_size_born == 1  ~ "Very large",
                                    cc6_size_born == 2  ~ "Larger than average",
                                    cc6_size_born == 3  ~ "Average",
                                    cc6_size_born == 4  ~ "Smaller than average",
                                    cc6_size_born == 5  ~ "Very small",
                                    cc6_size_born == 99 ~ " ",
                                    TRUE ~ NA_character_),
          across(c(cc7_early_9months,cc9_weighed_birth,cf3_colustrum,cf9_breastfed_yesterday,cf16_solid_food), recode_yesnodk),
          across(c(cf7_currnt_breastfeeding,cf11_stop_why_childs,cf11_stop_why_preg,
                   cf11_stop_why_mom_sick,cf11_stop_why_separated,cf11_stop_why_mom_no,
                   cf11_stop_why_solid, cf11_stop_why_childsick,cf11_stop_why_enough_milk), recode_yesno),
          cc9_birth_weight =  as.numeric(cc9_birth_weight),
          cc9_birth_weight = if_else(cc9_birth_weight == 99, NA_real_, cc9_birth_weight),
          cf12_exclusive = if_else(cf12_give_6m == "1", "Yes", "No"),
          cf_breastfed = if_else(
             if_any(c(cf9_breastfed_yesterday, cf7_currnt_breastfeeding), ~ . == "Yes"),
             "Yes",
             "No"
          ),
          across(c(starts_with("cf14"),cm3_diarhoea), recode_yesnodk)
   ) %>% 
   select(-c(starts_with("cf12_give_6m"),cf14i_yogurt_times,cf16_solid_food_times,cf14c_infant_formula_times,cf14e_other_milk_times,
             starts_with("cf6_"),starts_with("cf4_"),cf14a_plain_water, cf14f_juice, cf14g_clear_broth,cf14j_thin_porridge,
             cf14k_other_liquids,cf14l_an_other, starts_with("cf11_stop_why"),starts_with("cc9_"),starts_with("cf2_"),
             starts_with("cf13")))

##########====clean up and append ICYF
child_dd <- cdd_vars %>% 
   left_join(other_cddvars, by = c("kebele.x","hhid","subject_id")) %>% 
   # safe to drop individual food groups
   select(-c(cdd_breastmilk, cdd_grains_roots_tubers, cdd_vitamin_a, cdd_other_fv, 
                                 cdd_flesh_foods,cdd_eggs,cdd_legumes_nuts,cdd_dairy))

#===============================================================================
######===Maternal diet diversity variables
#===============================================================================
mdd_vars <- analysis_data %>% 
   select(kebele.x,hhid,subject_id,starts_with("md_")) %>% 
   mutate(
      # 1. Grains, roots, tubers, plantains
      mdd_grains = as.integer(md_cereals == 1 | md_tubers == 1),
      
      # 2. Pulses (beans, peas, lentils)
      mdd_pulses = as.integer(md_legumes == 1),
      
      # 3. Nuts and seeds (included in your legumes var → reuse)
      mdd_nuts_seeds = as.integer(md_legumes == 1),
      
      # 4. Dairy
      mdd_dairy = as.integer(md_milk == 1),
      
      # 5. Meat, poultry, fish (combine all flesh foods)
      mdd_meat_fish = as.integer(
         md_flesh_meat == 1 | md_organ_meat == 1 | md_fish == 1
      ),
      
      # 6. Eggs
      mdd_eggs = as.integer(md_eggs == 1),
      # 7. Dark green leafy vegetables
      mdd_dark_green = as.integer(md_dark_green == 1),
      # 8. Vitamin A-rich fruits and vegetables
      mdd_vitA = as.integer(md_vita == 1 | md_vita_fruits == 1),
      
      # 9. Other vegetables
      mdd_other_veg = as.integer(md_other_vegetables == 1),
      
      # 10. Other fruits
      mdd_other_fruit = as.integer(md_other_fruits == 1)
   ) %>% 
   mutate(
      mdd_score = rowSums(select(., starts_with("mdd_")), na.rm = TRUE),
      minimum_mdd = if_else(mdd_score < 6, "LDD", "ADD") 
   ) %>% 
   select(-c(md_cereals, md_tubers, md_legumes, md_milk, md_flesh_meat, md_organ_meat, md_fish, 
             md_eggs, md_dark_green, md_vita, md_vita_fruits, md_other_vegetables, md_other_fruits,
             md_fast_type)) %>% 
   mutate(across(c(md_oils, md_sweets, md_spices),recode_yesnodk),
          across(c(mdd_grains, mdd_pulses, mdd_nuts_seeds, mdd_dairy,mdd_meat_fish, 
                   mdd_eggs, mdd_dark_green, mdd_vitA, mdd_other_veg,mdd_other_fruit),
                 recode_yesno))

other_mmdvars <- analysis_data %>%
   select(kebele.x,hhid,subject_id,starts_with("md")) %>% 
   select(-c(starts_with("md_"))) %>% 
   
   ###===== clean up foor avoid
   mutate(md_food_avoid = case_when(
      # --- No / None variations ---
      md4_food_avoid %in% c("የለም","የለመ","የለሚ","የለሞ","የለውም","የልም","ዪለሚ",
                            "Yelem","Yelam","Yalame","Yalam","Yelme","Yelme.","Yelm",
                            "Ylme","Ylm","Yleme","Ylam","Ylem","No","የ ለ ም","99","ምንም የለም")  ~      "none",
      # --- Yes variations ---
      md4_food_avoid %in% c("አዎ","አው","አወ","አዉ","አዎ") ~ "yes",
      
      # --- Porridge / Genfo ---
      md4_food_avoid %in% c("ገንፎ","Genfo","Gefo","Gef0","ganfo","Jenefo",
                            "አጀ ገንፎ","አጃ፣ገንፎ") ~ "porridge",
      
      # --- Flatbread / Kita ---
      md4_food_avoid %in% c("ቂጣ","Kita","ትኩስ ቂጣ","Kita.gomen") ~ "flatbread",
      
      # --- Injera ---
      md4_food_avoid %in% c("Hnjera","Engra") ~ "injera",
      
      # --- Bread ---
      md4_food_avoid %in% c("Bread","Berber") ~ "bread",
      
      # --- Pasta ---
      md4_food_avoid %in% c("Pasta","ፓስታ","Past, aja","Paseta") ~ "pasta",
      
      # --- Rice ---
      md4_food_avoid %in% c("ሩዝ","Ruz","Ruze") ~ "rice",
      
      # --- Meat ---
      md4_food_avoid %in% c("ስጋ","Meat","Muk") ~ "meat",
      
      # --- Coffee 
      md4_food_avoid %in% c("Coffee") ~ "coffee",
      
      # --- Juice 
      md4_food_avoid %in% c("juice") ~ "juice",
      
      # --- Honey -
      md4_food_avoid %in% c("Honey") ~ "honey",
      
      # --- Cheese ---
      md4_food_avoid %in% c("Cheese") ~ "cheese",
      
      # --- Yogurt / Milk ---
      md4_food_avoid %in% c("Yogurt","Wetet","Yewetet  zer","የላም") ~ "milk / yogurt",
      
      # --- Kocho ---
      md4_food_avoid %in% c("ቆጮ","Kocho") ~ "kocho",
      
      # --- Kale / Greens ---
      md4_food_avoid %in% c("Kale","ጎመን") ~ "greens",
      
      # --- Butter ---
      md4_food_avoid %in% c("ገንፎ እና ቅቤ") ~ "porridge with butter",
      
      # --- Mixed dishes ---
      md4_food_avoid %in% c("ፓስታ ባቄላ") ~ "pasta with beans",
      
      # --- Bulla ---
      md4_food_avoid %in% c("Bula","Bulla") ~ "bulla",
      
      # --- Other specific items ---
      md4_food_avoid %in% c("Ergowtte") ~ "traditional food",
      md4_food_avoid %in% c("Serotonin.kayserwate") ~ "mixed foods",
      md4_food_avoid %in% c("Genfo.kita Gorman.baketa.almablt.yakmatta.nager.almblt") ~ "mixed foods",
      
      # --- Default ---
      TRUE ~ md4_food_avoid
   )) %>% 
   ####=== maorbidyt
   mutate(
      md_illness = case_when(
         md31_illness == 1 ~ "Yes",
         md31_illness == 2 ~ "No",
         TRUE ~ NA_character_
      )) %>% 
   mutate(usually_fast = case_when(
      md7_fasting == 1 ~ "Yes",
      md7_fasting == 2 ~ "No",
      md7_fasting == 3 ~ "DK"),
      yesterday_fast = case_when(
         md26_fast == 1 ~ "Yes",
         md26_fast == 2 ~ "No",
         md26_fast == 3 ~ "DK")) %>% 
   mutate(across(c(md6_vit_a, md27_tablet, md29_unusual_intake), recode_yesnodk)) %>%
   select(-c(md1_no_meals_before, md2_no_meals_during, md3_non_meals_delivery, 
             md31_illness,md4_food_avoid,md26_fast,md7_fasting, starts_with("md9_"),md5_change_diet,
             md27_specify,md29_reason,md32_appetite ))



#########=====combine necessasry variables for MDD instrument

maternal_dd <- mdd_vars %>% 
   left_join(other_mmdvars, by = c("kebele.x","hhid","subject_id")) %>% 
   select(-c(mdd_grains, mdd_pulses, mdd_nuts_seeds, mdd_dairy,mdd_meat_fish, 
             mdd_eggs, mdd_dark_green, mdd_vitA, mdd_other_veg,mdd_other_fruit))

sum(complete.cases(maternal_dd))
#########===3.  Household food insecurity access scale (HFIAS)
hfas <- analysis_data %>%
   mutate(
      across(c(hifas2_often, hifas4_oftern, hifas6_often, hifas8_often,
               hifas10_often, hifas12_often, hifas14_often,
               hifas16_often, hifas18_often),
             ~ as.numeric(.))
   ) %>%
   mutate(
      # Q1
      hfias_q1 = if_else(hifas1_enough_food == 1, hifas2_often, 0),
      
      # Q2
      hfias_q2 = if_else(hifas3_unable_eat == 1, hifas4_oftern, 0),
      
      # Q3
      hfias_q3 = if_else(hifas5_limited == 1, hifas6_often, 0),
      
      # Q4
      hfias_q4 = if_else(hifas7_not_really == 1, hifas8_often, 0),
      
      # Q5
      hfias_q5 = if_else(hifas9_smaller == 1, hifas10_often, 0),
      
      # Q6
      hfias_q6 = if_else(hifas11_fewer == 1, hifas12_often, 0),
      
      # Q7
      hfias_q7 = if_else(hifas13_no_food == 1, hifas14_often, 0),
      
      # Q8
      hfias_q8 = if_else(hifas15_sleep_hungry == 1, hifas16_often, 0),
      
      # Q9
      hfias_q9 = if_else(hifas17_whole_day == 1, hifas18_often, 0)
   ) %>%
   mutate(
      hfias_score = rowSums(across(hfias_q1:hfias_q9), na.rm = TRUE),
      hfias_category = case_when(
         # Category 1: Food Secure
         (
            hfias_q1 <= 1 &
               hfias_q2 == 0 &
               hfias_q3 == 0 &
               hfias_q4 == 0 &
               hfias_q5 == 0 &
               hfias_q6 == 0 &
               hfias_q7 == 0 &
               hfias_q8 == 0 &
               hfias_q9 == 0
         ) ~ 1,
         # Category 2: Mildly Food Insecure
         (
            (hfias_q1 >= 2 & hfias_q1 <= 3) |
               (hfias_q2 == 1 & hfias_q1 <= 3) |
               (hfias_q3 == 1 & hfias_q1 <= 3) |
               (hfias_q4 == 1 & hfias_q1 <= 3) |
               (hfias_q5 == 1 & hfias_q1 <= 3) |
               (hfias_q6 == 1 & hfias_q1 <= 3) |
               (hfias_q7 == 1 & hfias_q1 <= 3) |
               (hfias_q8 == 1 & hfias_q1 <= 3) |
               (hfias_q9 == 1 & hfias_q1 <= 3)
         ) ~ 2,
         # Category 3: Moderately Food Insecure 
         (
            (
               hfias_q3 %in% c(2,3) |
                  hfias_q4 %in% c(2,3) |
                  hfias_q5 %in% c(1,2) |
                  hfias_q6 %in% c(1,2)
            ) &
               hfias_q7 == 0 &
               hfias_q8 == 0 &
               hfias_q9 == 0
         ) ~ 3,
         # Category 4: Severely Food Insecure
         (
            hfias_q5 == 3 |
               hfias_q6 == 3 |
               hfias_q7 %in% c(1,2,3) |
               hfias_q8 %in% c(1,2,3) |
               hfias_q9 %in% c(1,2,3)
         ) ~ 4,
         TRUE ~ NA_real_
      ),
      hfias_category = factor(
         hfias_category,
         levels = c(1,2,3,4),
         labels = c(
            "Food secure",
            "Mildly food insecure",
            "Moderately food insecure",
            "Severely food insecure"
         )
      )
   )%>% 
   select(subject_id,hfias_category) 

#########4. WASH variables Standardisation

#https://www.who.int/data/nutrition/nlis/info/improved-sanitation-facilities-and-drinking-water-sources
#https://washdata.org/topics/drinking-water

wash_vars <- analysis_data %>%
   select(kebele.x,subject_id,hhid, starts_with("wa")) %>% 
   mutate(
      wash_drinking_water = case_when(
         wa1_source_drink %in% c(1, 2, 3, 4, 5, 6, 8, 11, 14) ~ "Improved",  # improved
         wa1_source_drink %in% c(7, 9, 10, 12, 13, 15,"other") ~ "Unimproved",        # unimproved
         TRUE ~ NA_character_
      ),
      wash_cooking_water = case_when(
         wa2_source_cook %in% c(1, 2, 3, 4, 5, 6, 8, 11, 14) ~ "Improved",  # improved
         wa2_source_cook %in% c(7, 9, 10, 12, 13, 15,"other") ~ "Unimproved",        # unimproved
         TRUE ~ NA_character_
      ),
      across(c(wa4_treat_water,wa13_antimalaria_spray,wa14_mosquito_net,wa7_share_facility),
             recode_yesnodk),
      wash_toilets = case_when(
         str_detect(wa6_toilet, "\\b(1|2|3|6|7|10)\\b") ~ "Improved",
         str_detect(wa6_toilet, "\\b(4|5|8|9|11|12|13)\\b") ~ "Unimproved",
         TRUE ~ NA_character_
      ),
      wash_child_stool_disposal = case_when(
         wa12_dispose_stool %in% c(1, 2) ~ "Safe disposal",
         wa12_dispose_stool %in% c(3, 4, 5, 6, 7, 8, 9, 98) ~ "Unsafe disposal",
         TRUE ~ NA_character_
      ),
      wash_soap_type = case_when(
         wa11_type_soap == 1 ~ "Basic",
         wa11_type_soap == 2 ~ "Limited",
         wa11_type_soap == 3 ~ "No Service",
         wa11_type_soap == 0 ~ "No Service",
         TRUE ~ NA_character_
      ),
      wash_handwashing_place = case_when(
         wa9_wash_hand %in% c(1, 2) ~ "Dedicated place",
         wa9_wash_hand %in% c(3, 4) ~ "No dedicated place",
         wa9_wash_hand %in% c(5, 6) ~ "None/unknown",
         TRUE ~ NA_character_
      ),
   ) %>% 
   select(kebele.x,subject_id,wa4_treat_water,starts_with("wash_"))

############=====5. SES variables - education, occupation, land ownership, agricultural production
SES_variables <- analysis_data %>% 
   select(kebele.x,subject_id,hhid, starts_with("hc")) %>% 
   mutate(
      edu_clean = case_when(
         # Fix obvious typos
         hc5_3_grade_comp %in% c("1o") ~ "10",
         # Degree / diploma
         str_detect(hc5_3_grade_comp, "degree|ድግሪ|digree") ~ "degree",
         str_detect(hc5_3_grade_comp, "diploma") ~ "diploma",
         # Secondary school (Amharic)
         str_detect(hc5_3_grade_comp, "ሁለተኛ") ~ "secondary",
         # Extract numeric grades (e.g. "10 grade", "10ኛ", etc.)
         str_detect(hc5_3_grade_comp, "\\d+") ~ str_extract(hc5_3_grade_comp, "\\d+"),
         TRUE ~ NA_character_
      ),
      across(c(hc7_private_agric_land, hc9_gro_crops), recode_yesno),
      marital_status = case_when(
         hc4_marital_status == 1 ~ "Married",
         TRUE ~ "Not Married"
      ),
      education_level = case_when(
         hc5_highest_level == 1 ~ "Illiterate",
         hc5_highest_level == 2 ~ "Read and Write",
         hc5_highest_level == 3  ~ "Formal education",
         TRUE ~ NA_character_
      ),
      occupation = case_when(
         hc6_income == 1 ~ "Farming",
         hc6_income == 2 ~ "Govt employment",
         hc6_income == 3 ~ "Small business",
         hc6_income == 4 ~ "No Job",
         hc6_income == 5 ~ "Other",
         TRUE ~ NA_character_
      ),
      hc2_age = as.numeric(hc2_age),
      hc2_age = na_if(hc2_age, 99),
      hc3_family_size = as.numeric(hc3_family_size),
      hc3_family_size = na_if(hc3_family_size, -5),
      hc3_family_size = na_if(hc3_family_size, 99)
      
   )  %>% 
   select(-c(starts_with("hc10_"),hc7_land_timad,hc7_land_acre,edu_clean,hc5_3_grade_comp,
             hc4_marital_status,hc5_highest_level,hc6_income ))



########################6 . Calculate Z scores.
library(anthro)
zscore_df <- analysis_data %>% 
   select(subject_id,cc1_sex,cc4_age_months,weight_child,length_child,muac_child,
          weight_mother,muac_mother,hgb_child,cc3_dob) %>% 
   mutate(cc4_age_months = as.numeric(cc4_age_months),
          cc4_age_months = if_else(cc4_age_months == 32, NA_real_, cc4_age_months),
          weight_child = if_else(weight_child < 0.2, NA_real_, weight_child),
          weight_child = if_else(weight_child > 30, NA_real_, weight_child),
          weight_child = if_else((weight_child == 17.4 & cc4_age_months == 17),NA_real_,weight_child)) %>%
   ###== implausible values as data entry errors and missing 
   mutate(length_child = if_else((length_child == 76.1 & cc4_age_months == 1) |
                                    (length_child == 92.5 & cc4_age_months == 2) |
                                    (length_child == 42.3 & cc4_age_months == 12),
                                 NA_real_, length_child))


zscores <- anthro_zscores(
   sex    = zscore_df$cc1_sex,         # 1 = boy, 2 = girl
   age    = zscore_df$cc4_age_months,
   is_age_in_month = TRUE,
   weight = zscore_df$weight_child,   # kg
   lenhei = zscore_df$length_child
)
anthro_df <- cbind(zscore_df, zscores) %>% 
   select(subject_id,weight_child,length_child,muac_child,weight_mother,hgb_child,
          muac_mother,zlen,zwfl,zwei)  %>% 
   mutate(muac_child = as.numeric(muac_child),
          muac_mother = as.numeric(muac_mother),
          hgb_child = as.numeric(hgb_child))


######## 7. Pre Analysis Data Set
analysis_data_clean <- anthro_df %>% 
   left_join(child_dd, by = c( "subject_id")) %>% 
   left_join(maternal_dd, by = c("kebele.x", "subject_id", "hhid")) %>% 
   left_join(SES_variables, by = c("kebele.x", "subject_id", "hhid")) %>% 
   left_join(wash_vars, by = c("kebele.x", "subject_id")) %>% 
   left_join(hfas, by = c("subject_id")) %>%
   #left_join(phq_score, by = c("subject_id"))  %>%  Ignore after chat with Dr. Kokeb
   mutate(stunting = if_else(zlen < -2, "Yes", "No"),
          wasting = if_else(zwfl < -2, "Yes", "No"),
          underweight = if_else(zwei < -2, "Yes", "No"))




######====outsheet out of range outliers based on WHO cutoffs for review and decision on whether to exclude or correct
######== WHO flags outliers as:

zwfl_outliers <- analysis_data_clean %>% 
   filter(zwfl < -5 | zwfl > 5) %>% 
   select(subject_id, zwfl, weight_child,cc1_sex, length_child, cc4_age_months)

zlen_outliers <- analysis_data_clean %>% 
   filter(zlen < -6 | zlen > 6) %>% 
   select(subject_id, zlen, weight_child,cc1_sex, length_child, cc4_age_months)

zwei_outliers <- analysis_data_clean %>% 
   filter(zwei < -5 | zwei > 5) %>% 
   select(subject_id, zwei, weight_child,cc1_sex, length_child, cc4_age_months)


writexl::write_xlsx(
   list(
      WHZ = zwfl_outliers,
      HAZ = zlen_outliers,
      WAZ = zwei_outliers
   ),
   path = "Z_score_outliers.xlsx"
)

######==== Pre Analysis Data Set
save(analysis_data_clean, file = "2. Clean Data/analysis_data_clean.RData")
save(main_dictionary, file = "2. Clean Data/main_dictionary.RData")

#######=== End of R script. Next file is Final analysis Preprocessing====


#https://www.sciencedirect.com/science/article/pii/S0160412021006772?via%3Dihub
##Wood, White & Royston (2008) — Statistics in Medicine — discusses variable selection instability in the context of MI and notes that examining which variables are selected across imputed datasets is a useful diagnostic, though they do not specify a percentage threshold.
#Heymans & Eekhout (2019) — Applied Missing Data Analysis with SPSS and R — describe running models across imputed datasets separately as a way to assess selection consistency, but again without a formal cutoff.
#Sauerbrei & Schumacher (1992) — Statistics in Medicine — introduced the bootstrap inclusion frequency approach, where a variable is retained if selected in ≥ a threshold percentage of bootstrap samples. The 50% and 80% cutoffs originate here, though in a bootstrap context rather than MI.
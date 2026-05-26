## ============================================================
## 02_consolidate_cohorts.R
## Consolida cohortes a 7 categorías biológicamente coherentes
## ============================================================
suppressPackageStartupMessages({ library(data.table) })

dt <- readRDS(here::here("A4_empirical", "data", "master_long.rds"))

classify <- function(x){
  x <- tolower(trimws(as.character(x)))
  if (x %in% c("lean","no obesity-no t2dm cohort","no obesity- no t2dm cohort"))      return("Lean")
  if (grepl("after.*roux", x) || grepl("post-rygb", x))                                return("Post-RYGB")
  if (grepl("after weight loss.*roux", x))                                             return("Post-RYGB")
  if (grepl("before.*rygb", x) || grepl("before.*roux", x))                            return("Pre-RYGB")
  if (grepl("after.*sleeve", x) || grepl("post-sg", x))                                return("Post-SG")
  if (grepl("before sg", x) || grepl("before.*sleeve", x))                             return("Pre-SG")
  if (grepl("after.*calor.*restric", x) || grepl("obesity-cr", x))                     return("Post-CR")
  if (grepl("before.*calor.*restric", x))                                              return("Pre-CR")
  if (grepl("plus.*type 2", x) || grepl("with.*t2dm", x) || grepl("with type 2", x))   return("T2DM-Obesity")
  if (grepl("type 2 diabetes", x) || x == "t2dm")                                      return("T2DM")
  if (grepl("obesity.*middle", x))                                                     return("Obesity")
  if (grepl("obesity\\s+grade\\s*ii", x))                                              return("Obesity-II")
  if (grepl("obesity-iii", x) || grepl("obesity\\s+grade\\s+iii", x))                  return("Obesity-III")
  if (grepl("obesity-i", x) || grepl("obesity\\s+grade\\s+i\\b", x))                   return("Obesity-I")
  if (x == "obesity")                                                                  return("Obesity")
  return("Obesity")
}

dt[, cohort_lvl1 := sapply(cohort, classify)]
dt[, cohort_lvl2 := dplyr::case_when(
  cohort_lvl1 %in% c("Obesity","Obesity-I","Obesity-II","Obesity-III") ~ "Obesity",
  cohort_lvl1 %in% c("Pre-RYGB","Pre-SG","Pre-CR") ~ "Pre-Intervention",
  cohort_lvl1 %in% c("Post-RYGB","Post-SG","Post-CR") ~ "Post-Intervention",
  cohort_lvl1 %in% c("T2DM-Obesity","T2DM") ~ "T2DM",
  TRUE ~ cohort_lvl1
)]

cat("Cohorte nivel 1 (granular):\n")
print(dt[, .N, by=cohort_lvl1][order(-N)])
cat("\nCohorte nivel 2 (consolidada):\n")
print(dt[, .N, by=cohort_lvl2][order(-N)])

cat("\nN por cohort × hormone (nivel 2):\n")
tab <- dcast(dt, cohort_lvl2 ~ hormone, value.var="value", fun.aggregate=length)
print(tab)

saveRDS(dt, here::here("A4_empirical", "data", "master_long.rds"))
fwrite(dt, here::here("A4_empirical", "data", "master_long.csv"))
cat("\nGuardado con cohortes consolidadas.\n")

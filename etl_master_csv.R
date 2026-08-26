# =============================================================================
# ETL: Tabla maestra AUC E and P Hormones → long tidy format
# =============================================================================
# Input:  tabla maestra "AUC E and P Hormones" (fuera del repositorio; ruta en
#         la variable de entorno EPP10_MASTER_CSV — ver epp10_paths.R)
# Output: <raíz del repositorio>/hormones_long_tidy.csv
#
# Operations:
#   1. Read CSV (skip title row), fill-down Author and Cohort
#   2. Drop Excel merged-cell bleed rows (where Time == literal "Time")
#   3. Apply cohort normalization (frozen in preregistration_cohort_map.yaml)
#   4. Unpivot 11 hormone blocks into long format
#   5. Coerce time and value to numeric
#   6. Emit tidy tibble with metadata columns
# =============================================================================

.libPaths(c("~/.R/library", .libPaths()))

# --- Localiza la raíz del repositorio (Rscript, source() o RStudio) ---------
if (!exists("epp10_path")) {
  .epp10_self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  source(file.path(if (length(.epp10_self)) dirname(sub("^--file=", "", .epp10_self[[1]]))
                   else getwd(), "epp10_paths.R"))
}
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr); library(purrr)
})

INPUT_CSV  <- epp10_master_csv()
OUTPUT_CSV <- epp10_path("hormones_long_tidy.csv")
SHA256     <- "ce2e343d29067a79bbefb9c19bf144d0182e977870b43e9a0c6a1a68e06ceefa"

# --- 1. Read + fill-down ----------------------------------------------------
raw <- read_csv(INPUT_CSV, skip = 1, show_col_types = FALSE,
                name_repair = "unique_quiet")
stopifnot(digest::digest(read_file_raw(INPUT_CSV), algo = "sha256") == SHA256)

filled <- raw %>%
  mutate(across(c(Author, Cohort), \(x) if (is.character(x)) str_trim(x) else x)) %>%
  mutate(across(c(Author, Cohort), \(x) na_if(x, ""))) %>%
  fill(Author, Cohort, .direction = "down") %>%
  # Fill down study-level metadata: constant within each (Author, Cohort) block
  group_by(Author, Cohort) %>%
  fill(Reference, `Number of Subjects`, `Mean Body mass index (kg/m2)`,
       `Mean Age`, `Caloric test`, `Kcal of the Caloric test`,
       .direction = "downup") %>%
  ungroup()

# --- 2. Cohort normalization (frozen per YAML) ------------------------------
normalize_cohort <- function(s) {
  x <- str_to_lower(s)
  has_obesity <- str_detect(x, "obesity|overweight")
  has_t2dm    <- str_detect(x, "type 2 diabetes|t2dm|plus type 2|with type 2")
  explicit_no_t2dm <- str_detect(x, "no obesity-no t2dm|no t2dm|non-diabetic|non diabetic|without t2dm")
  has_rygb    <- str_detect(x, "roux|rygb|rox-en-y|gastric bypass|gastric bypas")
  has_sg      <- str_detect(x, "sleeve gastrectomy|sleeve")
  has_calr    <- str_detect(x, "caloric|calóric")
  is_pre      <- str_detect(x, " before ")
  is_post     <- str_detect(x, " after ") | str_detect(x, " at 1[- ]year| weeks|week 13")
  weeks <- case_when(
    str_detect(x, "6 weeks")                                 ~ "6w",
    str_detect(x, "12 weeks|week 13")                        ~ "12w",
    str_detect(x, "1[- ]year|at 1 year|1 year after|after 1 year") ~ "1y",
    TRUE                                                      ~ NA_character_
  )
  is_post_cr <- is_post & has_calr
  cohort_v10_primary <- case_when(
    str_detect(x, "no obesity-no t2dm")  ~ "no_obese_without_T2DM",
    is_post_cr                           ~ "Obesity",
    is_post & has_rygb                   ~ "RYGBP",
    is_post & has_sg                     ~ "SG",
    has_obesity & has_t2dm               ~ "Obesity_T2DM",
    has_obesity                          ~ "Obesity",
    has_t2dm                             ~ "T2DM",
    TRUE                                 ~ "UNCLASSIFIED"
  )
  cohort_v10_sensitivity <- if_else(is_post_cr, "caloric_restriction_post",
                                    cohort_v10_primary)
  modality <- case_when(is_post & has_rygb ~ "RYGB", is_post & has_sg ~ "SG",
                        is_post_cr ~ "caloric_restriction", TRUE ~ "none")
  surg_status <- case_when(is_pre ~ "pre_surgery", is_post ~ "post_surgery",
                           TRUE ~ "not_applicable")
  had_t2dm_pre_surgery <- case_when(has_t2dm ~ "TRUE",
                                    explicit_no_t2dm ~ "FALSE",
                                    TRUE ~ "unknown")
  tibble(cohort_v10_primary     = cohort_v10_primary,
         cohort_v10_sensitivity = cohort_v10_sensitivity,
         surgery_status         = surg_status,
         weeks_post_surgery     = weeks,
         weight_loss_modality   = modality,
         had_t2dm_pre_surgery   = had_t2dm_pre_surgery)
}

meta <- filled %>% distinct(Cohort) %>% filter(!is.na(Cohort)) %>%
  mutate(norm = map(Cohort, normalize_cohort)) %>% unnest(norm)

# --- 3. Hormone block definitions ------------------------------------------
HORMONE_BLOCKS <- tribble(
  ~hormone_name,    ~col_name,  ~col_time, ~col_value, ~units,
  "ghrelin_total", 17L, 18L, 19L, "pmol/L",
  "ghrelin_acyl",  21L, 22L, 23L, "pmol/L",
  "GIP_total",     25L, 26L, 27L, "pmol/L",
  "GIP_active",    29L, 30L, 31L, "pmol/L",
  "GLP1_total",    33L, 34L, 35L, "pmol/L",
  "GLP1_active",   37L, 38L, 39L, "pmol/L",
  "PYY_total",     41L, 42L, 43L, "pmol/L",
  "PYY_3_36",      46L, 47L, 48L, "pmol/L",
  "insulin",       51L, 52L, 53L, "pmol/L",
  "glucagon",      55L, 56L, 57L, "pmol/L",
  "glucose",       59L, 60L, 61L, "mmol/L"
)

# --- 4. Unpivot per hormone block ------------------------------------------
extract_hormone <- function(df, blk) {
  tibble(
    Author     = df$Author,
    Cohort     = df$Cohort,
    Reference  = df$Reference,
    n_subjects = df[["Number of Subjects"]],
    bmi_mean   = df[["Mean Body mass index (kg/m2)"]],
    age_mean   = df[["Mean Age"]],
    caloric_test = df[["Caloric test"]],
    kcal_test    = df[["Kcal of the Caloric test"]],
    hormone_name = blk$hormone_name,
    units        = blk$units,
    time_raw     = df[[blk$col_time]],
    value_raw    = df[[blk$col_value]]
  ) %>%
    # Drop bleed rows (literal "Time" / "Pmol/L" text)
    filter(!(str_detect(time_raw, "^(Time|time)$") |
             str_detect(value_raw, "^(Pmol/L|pmol/L|mmol/L)$"))) %>%
    # Coerce to numeric
    mutate(time_min = suppressWarnings(as.numeric(time_raw)),
           value    = suppressWarnings(as.numeric(value_raw))) %>%
    filter(!is.na(time_min), !is.na(value)) %>%
    select(Author, Reference, source_cohort = Cohort, n_subjects,
           bmi_mean, age_mean, caloric_test, kcal_test,
           hormone_name, units, time_min, value)
}

long <- pmap_dfr(HORMONE_BLOCKS, function(hormone_name, col_name, col_time,
                                          col_value, units) {
  blk <- list(hormone_name = hormone_name, col_time = col_time,
              col_value = col_value, units = units)
  extract_hormone(filled, blk)
})

# --- 5. Join cohort normalization ------------------------------------------
long_tidy <- long %>%
  left_join(meta, by = c("source_cohort" = "Cohort"))

# --- 6. Report + save ------------------------------------------------------
cat(sprintf("=== ETL COMPLETE ===\n"))
cat(sprintf("Input rows:            %d\n", nrow(raw)))
cat(sprintf("After fill-down:       %d (non-NA Author+Cohort)\n",
            sum(!is.na(filled$Author) & !is.na(filled$Cohort))))
cat(sprintf("Long-format rows:      %d (observations)\n", nrow(long_tidy)))
cat(sprintf("Unique (Author,Cohort): %d\n",
            n_distinct(paste(long_tidy$Author, long_tidy$source_cohort))))
cat(sprintf("Hormones observed:     %d / 11\n",
            n_distinct(long_tidy$hormone_name)))

cat("\nRows per hormone:\n")
print(long_tidy %>% count(hormone_name, name = "n_obs") %>% arrange(desc(n_obs)))

cat("\nRows per primary canonical cohort:\n")
print(long_tidy %>% count(cohort_v10_primary, name = "n_obs"))

cat("\nTimepoint distribution:\n")
print(long_tidy %>% count(time_min, name = "n_obs") %>% arrange(time_min) %>% head(20))

cat("\nUnclassified rows (should be 0):\n")
cat(sum(long_tidy$cohort_v10_primary == "UNCLASSIFIED"), "\n")

write_csv(long_tidy, OUTPUT_CSV)
cat(sprintf("\nSaved: %s\n", OUTPUT_CSV))
cat(sprintf("SHA-256 of output: %s\n",
            digest::digest(read_file_raw(OUTPUT_CSV), algo = "sha256")))

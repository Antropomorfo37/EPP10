## ============================================================
## 01_parse_master.R  —  Parser de la tabla maestra AUC E&P
## H.M. Virgen-Ayala, 2026-05-14
## Convierte el CSV ancho (multi-bloque) a long: study|cohort|kcal|test|N|BMI|age|hormone|time|value
## ============================================================
suppressPackageStartupMessages({
  library(tidyverse); library(data.table)
})

## ---------------------------------------------------------------
## DATA SOURCE NOTE (reproducibility):
## This script (01) ingests the primary master table, which is the
## author-curated digitisation of aggregate data from the source studies.
## This raw input file is NOT distributed in the public repository.
## The harmonised output (master_long.csv / .rds) IS included under
## A4_empirical/data/, so the reproducible public pipeline begins at
## 02_consolidate_cohorts.R. To re-run script 01, set CSV_PATH below to
## your local copy of the master table.
## ---------------------------------------------------------------
CSV_PATH <- "/Users/hectormanuelvirgenayala/Library/Mobile Documents/com~apple~CloudDocs/Enteroendocrino/The entero-insular axis/Turings & EPA/***1 ***Tabla maestra AUC E and P Hormones  copia 2 2 2.csv"

raw <- suppressWarnings(read.csv(CSV_PATH, header = FALSE, stringsAsFactors = FALSE,
                                 fileEncoding = "UTF-8", na.strings = c("","NA")))
cat("Dim raw:", dim(raw), "\n")

## Encabezado (fila 1)
hdr <- as.character(raw[1, ])
hdr <- ifelse(is.na(hdr), "", hdr)

## Localiza los bloques de hormona. Cada bloque ocupa 3 columnas reales:
## col_i (etiqueta hormona, "Time"), col_i+1 = Time(min), col_i+2 = Pmol/L (valor)
## En el header la etiqueta de hormona aparece y luego "Time","Pmol/L" (o mmol/L glucosa).
hormone_labels <- c("TOTAL GHRELIN","ACYLATED GHRELIN","TOTAL GIP","ACTIVE GIP",
                    "TOTAL GLP-1","ACTIVE GLP-1","TOTAL PYY","PYY3-36",
                    "Insulin","Glucagon","Glucose")

blocks <- list()
for (h in hormone_labels){
  ix <- which(hdr == h)
  if (length(ix)) {
    # Time col = ix+1; value col = ix+2
    blocks[[h]] <- c(time_col = ix+1, val_col = ix+2)
  }
}
cat("Bloques de hormona localizados en header:\n"); print(blocks)

## Función helper para parsear numérico
num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", as.character(x))))

## Identificar inicios de estudio: filas con Author no-vacío (col1)
study_starts <- which(!is.na(raw[, 1]) & raw[, 1] != "" & raw[, 1] != "Author")
study_starts <- study_starts[study_starts > 1]
# bloque termina antes del siguiente start
study_ends <- c(study_starts[-1] - 1, nrow(raw))

cat("\nN bloques de estudio detectados:", length(study_starts), "\n")

## Para cada bloque, extraer metadatos del primer renglón y datos hormona de los renglones siguientes
parse_study <- function(s, e){
  meta_row <- raw[s, ]
  author <- meta_row[[1]]
  ref    <- meta_row[[2]]
  N      <- num(meta_row[[3]])
  cohort <- trimws(as.character(meta_row[[4]]))
  bmi    <- num(meta_row[[5]])
  bmi_sd <- num(meta_row[[6]])
  Nfem   <- num(meta_row[[7]])
  age    <- num(meta_row[[8]])
  age_sd <- num(meta_row[[9]])
  ## Caloric test info (cols ~14-16)
  ctest  <- trimws(as.character(meta_row[[15]]))
  kcal   <- num(gsub("[^0-9\\.]", "", as.character(meta_row[[16]])))

  ## Para cada bloque hormona, recoger filas s+1:e con time y value definidos
  body <- raw[(s):e, , drop = FALSE]
  out <- list()
  for (h in names(blocks)){
    tc <- blocks[[h]]["time_col"]; vc <- blocks[[h]]["val_col"]
    tt <- num(body[, tc]); vv <- num(body[, vc])
    keep <- which(!is.na(tt) & !is.na(vv))
    if (length(keep)){
      out[[h]] <- data.table(
        study = author, reference = ref, N = N, cohort = cohort,
        BMI = bmi, BMI_sd = bmi_sd, Nfem = Nfem, age = age, age_sd = age_sd,
        caloric_test = ctest, kcal = kcal,
        hormone = h, time_min = tt[keep], value = vv[keep]
      )
    }
  }
  rbindlist(out)
}

dt_long <- rbindlist(lapply(seq_along(study_starts),
                            function(i) parse_study(study_starts[i], study_ends[i])),
                     fill = TRUE)

cat("\nDim long:", dim(dt_long), "\n")
cat("Estudios únicos:", length(unique(dt_long$study)), "\n")
cat("Hormonas:", paste(unique(dt_long$hormone), collapse=", "), "\n")
cat("Cohortes:\n"); print(sort(unique(dt_long$cohort)))

## Normalizar nombres de cohorte
dt_long[, cohort := trimws(cohort)]
dt_long[cohort == "Obesity Cohort", cohort := "Obesity"]
dt_long[cohort == "No obesity-No T2DM Cohort", cohort := "Lean"]
dt_long[cohort == "Type 2 Diabetes Mellitus  Cohort", cohort := "T2DM"]
dt_long[grepl("Obesity grade  I Cohort", cohort), cohort := "Obesity-I"]
dt_long[grepl("Obesity grade II Cohort", cohort), cohort := "Obesity-II"]
dt_long[grepl("Obesity grade III Cohort", cohort), cohort := "Obesity-III"]
dt_long[grepl("Caloric Restriction", cohort), cohort := "Obesity-CR"]
dt_long[grepl("Roux-en-Y", cohort, ignore.case=TRUE), cohort := "Post-RYGB"]
dt_long[grepl("sleeve", cohort, ignore.case=TRUE), cohort := "Post-SG"]
dt_long[grepl("After Roux", cohort, ignore.case=TRUE), cohort := "Post-RYGB"]

cat("\nCohortes consolidadas:\n"); print(sort(unique(dt_long$cohort)))
cat("\nN observaciones por cohorte × hormona:\n")
print(dcast(dt_long, cohort ~ hormone, value.var = "value", fun.aggregate = length))

fwrite(dt_long, here::here("A4_empirical", "data", "master_long.csv"))
saveRDS(dt_long, here::here("A4_empirical", "data", "master_long.rds"))
cat("\nGuardado: data/master_long.csv y .rds\n")

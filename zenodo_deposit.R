# =============================================================================
# zenodo_deposit.R — Automated Zenodo deposit via API
# =============================================================================
# Crea un depósito en Zenodo, sube archivos, y obtiene DOI
# 
# Requisitos:
#   - httr2 package: install.packages("httr2")
#   - Token ZENODO_API_TOKEN en variables de entorno (.env)
#
# Uso:
#   Rscript zenodo_deposit.R --title "My Title" --files file1.pdf file2.csv
# =============================================================================

suppressPackageStartupMessages({
  library(httr2)
  library(jsonlite)
  library(glue)
})

# ============================================================================
# 1. INICIALIZACIÓN
# ============================================================================

ZENODO_TOKEN <- Sys.getenv("ZENODO_API_TOKEN")
ZENODO_API_URL <- Sys.getenv("ZENODO_API_URL", "https://zenodo.org/api")
# ZENODO_API_URL <- Sys.getenv("ZENODO_SANDBOX_URL", "https://sandbox.zenodo.org/api")  # Para testing

if (ZENODO_TOKEN == "") {
  stop("ERROR: ZENODO_API_TOKEN no configurado en .env\n",
       "  Crea un token en: https://zenodo.org/account/settings/applications/")
}

cat("✓ Zenodo token cargado\n")
cat("  API URL:", ZENODO_API_URL, "\n\n")

# ============================================================================
# 2. FUNCIÓN: Crear depósito
# ============================================================================

create_zenodo_deposit <- function(title, description, authors, 
                                  license = "CC-BY-4.0",
                                  upload_type = "software") {
  
  cat("Creando depósito en Zenodo...\n")
  
  # Metadata
  metadata <- list(
    metadata = list(
      title = title,
      description = description,
      upload_type = upload_type,
      creators = authors,  # list(list(name = "Author", affiliation = "Inst"))
      license = license,
      access_right = "open"
    )
  )
  
  # Request
  url <- paste0(ZENODO_API_URL, "/deposit/depositions")
  
  resp <- request(url) %>%
    req_headers(Authorization = paste("Bearer", ZENODO_TOKEN)) %>%
    req_body_json(metadata) %>%
    req_perform()
  
  if (resp$status_code != 201) {
    cat("ERROR: No se pudo crear depósito\n")
    cat("Status:", resp$status_code, "\n")
    cat("Respuesta:", resp_body_string(resp), "\n")
    return(NULL)
  }
  
  result <- resp_body_json(resp)
  deposit_id <- result$id
  
  cat("✓ Depósito creado: ID =", deposit_id, "\n")
  cat("  URL:", result$links$html, "\n\n")
  
  return(result)
}

# ============================================================================
# 3. FUNCIÓN: Subir archivos
# ============================================================================

upload_files_zenodo <- function(deposit_id, file_paths) {
  
  cat("Subiendo archivos a depósito", deposit_id, "...\n")
  
  for (file_path in file_paths) {
    if (!file.exists(file_path)) {
      cat("  ⚠ Archivo no encontrado:", file_path, "\n")
      next
    }
    
    file_name <- basename(file_path)
    file_size <- file.size(file_path)
    
    cat("  Subiendo:", file_name, 
        "(", format(file_size, units = "auto"), ")\n")
    
    url <- paste0(ZENODO_API_URL, "/deposit/depositions/", 
                  deposit_id, "/files")
    
    resp <- request(url) %>%
      req_headers(Authorization = paste("Bearer", ZENODO_TOKEN)) %>%
      req_body_multipart(
        filename = curl::form_file(file_path)
      ) %>%
      req_perform()
    
    if (resp$status_code != 201) {
      cat("    ERROR: No se pudo subir archivo\n")
      cat("    Status:", resp$status_code, "\n")
      next
    }
    
    result <- resp_body_json(resp)
    cat("    ✓ Subido exitosamente\n")
  }
  
  cat("\n")
}

# ============================================================================
# 4. FUNCIÓN: Publicar depósito y obtener DOI
# ============================================================================

publish_zenodo_deposit <- function(deposit_id) {
  
  cat("Publicando depósito", deposit_id, "...\n")
  
  url <- paste0(ZENODO_API_URL, "/deposit/depositions/", 
                deposit_id, "/actions/publish")
  
  resp <- request(url) %>%
    req_headers(Authorization = paste("Bearer", ZENODO_TOKEN)) %>%
    req_method("POST") %>%
    req_perform()
  
  if (resp$status_code != 202) {
    cat("ERROR: No se pudo publicar depósito\n")
    cat("Status:", resp$status_code, "\n")
    cat("Respuesta:", resp_body_string(resp), "\n")
    return(NULL)
  }
  
  result <- resp_body_json(resp)
  
  cat("✓ Depósito publicado\n")
  cat("  DOI:", result$doi, "\n")
  cat("  URL:", result$doi_url, "\n")
  cat("  URL de registro:", result$links$html, "\n\n")
  
  return(result)
}

# ============================================================================
# 5. EJEMPLO DE USO
# ============================================================================

if (interactive()) {
  
  # Configurar metadata
  TITLE <- "Sparse mFACEs + PTP/IEP Classification: EPP10 Pipeline"
  DESCRIPTION <- "Code for sparse multivariate functional principal component analysis with Chiou normalization and per-analyte periprandial trajectory pattern/inter-expectancy pattern classification applied to ecological meta-analysis of six metabolic cohorts."
  
  AUTHORS <- list(
    list(
      name = Sys.getenv("EPP10_AUTHOR_NAME", "Héctor Manuel Virgen Ayala"),
      affiliation = paste0(
        Sys.getenv("EPP10_DEPARTMENT", "Dept. Clínicas Quirúrgicas"), ", ",
        Sys.getenv("EPP10_INSTITUTION", "Universidad de Guadalajara")
      ),
      orcid = Sys.getenv("EPP10_AUTHOR_ORCID", "0009-0006-2081-2286")
    )
  )
  
  # Crear depósito
  deposit <- create_zenodo_deposit(
    title = TITLE,
    description = DESCRIPTION,
    authors = AUTHORS,
    license = "CC-BY-4.0",
    upload_type = "software"
  )
  
  if (!is.null(deposit)) {
    # Subir archivos
    files_to_upload <- c(
      "bootstrap_B2000_results.rds",
      "fit_mfaces_primary_results.rds",
      "ptp_iep_results.rds",
      "fanova_results.rds",
      "stability_classification_stage.rds",
      "preregistration_cohort_map.yaml",
      "CITATION.cff"
    )
    
    upload_files_zenodo(deposit$id, files_to_upload)
    
    # Publicar
    published <- publish_zenodo_deposit(deposit$id)
    
    # Guardar resultado
    if (!is.null(published)) {
      result_file <- file.path(
        Sys.getenv("EPP10_VERIFICATION_DIR", "."),
        "zenodo_deposit_result.json"
      )
      write_json(published, result_file, pretty = TRUE)
      cat("Resultado guardado:", result_file, "\n")
    }
  }
  
} else {
  
  # Script no-interactivo: parsear argumentos
  args <- commandArgs(trailingOnly = TRUE)
  
  # Ejemplo: Rscript zenodo_deposit.R --title "Title" --description "Desc" --files f1.pdf f2.csv
  
  cat("Zenodo deposit automation (non-interactive mode)\n")
  cat("Use dentro de pipelines CI/CD\n")
}

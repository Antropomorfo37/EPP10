# =============================================================================
# epp10_paths.R — resolución de rutas relativa a la raíz del proyecto
# =============================================================================
# Todos los scripts del pipeline direccionaban sus entradas y salidas por ruta
# absoluta al directorio del autor, lo que impedía ejecutar el repositorio en
# cualquier otra máquina. Este helper resuelve la raíz una sola vez y expone
# epp10_path() para construir rutas a partir de ella.
#
# Uso (cabecera estándar en cada script):
#   if (!exists("epp10_path")) {
#     .epp10_self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
#     source(file.path(if (length(.epp10_self)) dirname(sub("^--file=", "", .epp10_self[[1]]))
#                      else getwd(), "epp10_paths.R"))
#   }
#   read_csv(epp10_path("hormones_long_tidy.csv"))
#
# Orden de resolución de la raíz:
#   1. variable de entorno EPP10_ROOT (override explícito)
#   2. dirname() del --file= de Rscript (los scripts viven en la raíz)
#   3. ascenso desde getwd() buscando los ficheros marcadores del repositorio
# Sólo base R: sin dependencias nuevas, no toca renv.lock.
# =============================================================================

# Ficheros versionados y suficientemente distintivos como para identificar la
# raíz del repositorio sin ambigüedad.
EPP10_ROOT_MARKERS <- c("preregistration_cohort_map.yaml", "CITATION.cff")

.epp10_is_root <- function(dir) {
  all(file.exists(file.path(dir, EPP10_ROOT_MARKERS)))
}

.epp10_find_root_upwards <- function(start = getwd()) {
  dir <- normalizePath(start, winslash = "/", mustWork = FALSE)
  repeat {
    if (.epp10_is_root(dir)) return(dir)
    parent <- dirname(dir)
    if (identical(parent, dir)) return(NA_character_)   # llegamos a "/"
    dir <- parent
  }
}

.epp10_resolve_root <- function() {
  # 1. Override explícito por entorno.
  env_root <- Sys.getenv("EPP10_ROOT", unset = "")
  if (nzchar(env_root)) {
    root <- normalizePath(env_root, winslash = "/", mustWork = FALSE)
    if (!.epp10_is_root(root))
      stop(sprintf(paste0("EPP10_ROOT apunta a '%s', que no parece la raíz del ",
                          "repositorio (faltan: %s)."),
                   root,
                   paste(EPP10_ROOT_MARKERS[!file.exists(file.path(root, EPP10_ROOT_MARKERS))],
                         collapse = ", ")),
           call. = FALSE)
    return(root)
  }

  # 2. Bajo Rscript: los scripts del pipeline viven en la raíz del repositorio.
  self <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(self)) {
    root <- dirname(normalizePath(sub("^--file=", "", self[[1]]),
                                  winslash = "/", mustWork = FALSE))
    if (.epp10_is_root(root)) return(root)
  }

  # 3. source() / RStudio / sesión interactiva: ascender desde el directorio
  #    de trabajo hasta encontrar los marcadores.
  root <- .epp10_find_root_upwards(getwd())
  if (!is.na(root)) return(root)

  stop(paste0("No se pudo localizar la raíz del repositorio EPP10 desde '",
              getwd(), "'. Ejecuta los scripts desde el repositorio o exporta ",
              "EPP10_ROOT=/ruta/al/repo."),
       call. = FALSE)
}

EPP10_ROOT <- .epp10_resolve_root()

#' Construye una ruta a partir de la raíz del repositorio.
#'
#' epp10_path("figures", "Figure1.pdf") -> "<raíz>/figures/Figure1.pdf"
epp10_path <- function(...) {
  if (...length() == 0L) return(EPP10_ROOT)
  file.path(EPP10_ROOT, ...)
}

# -----------------------------------------------------------------------------
# Tabla maestra: única entrada del pipeline que vive FUERA del repositorio.
# No se versiona (datos de origen), así que se localiza por entorno.
# -----------------------------------------------------------------------------
#' Ruta a la tabla maestra AUC E and P Hormones (entrada de etl_master_csv.R).
#'
#' El SHA-256 esperado se verifica en etl_master_csv.R, que es su única fuente.
epp10_master_csv <- function() {
  path <- Sys.getenv("EPP10_MASTER_CSV", unset = "")
  if (!nzchar(path))
    stop(paste0(
      "La tabla maestra no está versionada en el repositorio.\n",
      "  Exporta EPP10_MASTER_CSV con la ruta al CSV de origen\n",
      "  ('Tabla maestra AUC E and P Hormones') antes de ejecutar\n",
      "  etl_master_csv.R o generate_compliance_and_appendix.R."),
      call. = FALSE)
  if (!file.exists(path))
    stop(sprintf("EPP10_MASTER_CSV apunta a '%s', que no existe.", path),
         call. = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

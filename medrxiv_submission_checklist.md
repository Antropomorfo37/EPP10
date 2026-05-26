# medRxiv Submission Checklist — Final Package (v1.2, 2026-04-24)

**Manuscrito**: *Dynamic Enteropancreatic Phenotyping via Sparse mFACEs and PTP/IEP Classification: An Ecological Meta-Analysis of Six Metabolic Cohorts*

**Autor**: Héctor Manuel Virgen Ayala, MD — UdeG + IMSS — ORCID 0009-0006-2081-2286 — hectorvirgenmd@gmail.com

**Release actual**: v1.2 (GitHub tag, Zenodo DOI 10.5281/zenodo.19750294)

## Manuscript core

- [x] **Manuscript PDF** — `manuscript/manuscript.pdf` (22 pp, 197 KB) — *physiologic periprandial responsiveness* + dual-track IEP + author identity
- [x] **Abstract estructurado** — en manuscript.pdf pág. 1 (Context/Objective/Design/Results/Conclusions)
- [x] **Introducción** — §1 + §1.1 Conceptual frame + §1.2 Statistical foundations + §1.3 Aim
- [x] **Methods** — §2.1–§2.8 (incluye §2.5 PTP/IEP scope, §2.6 borderline cutoff, §2.7 Validation, §2.8 Computational environment)
- [x] **Results** — §3.1–§3.9 (incluye dual-track §3.6 sin subclass + §3.7 con glycemic subclass)
- [x] **Discusión** — §4 (Principal finding, Conceptual contribution, Post-bariatric, Diabetic/convergent, Glycemic subclass, Methodological, 5 Limitations, Clinical implications)
- [x] **Tablas 1–3** — en manuscript.pdf (págs. 15–17)
- [x] **Tabla 4A** — IEP Types I–V sin glycemic subclass (pág. 18)
- [x] **Tabla 4B** — IEP Types con glycemic subclass a/b/c restringido a III/IV/V (pág. 19)
- [x] **Figuras 1–3** — [figures/Figure{1,2,3}*.pdf + PNG](figures/) (págs. 20–22 del PDF)

## Suplementarios (4 archivos, todos disponibles)

- [x] **S1. STROBE cohort checklist** — [verification/STROBE_checklist_S1.md](verification/STROBE_checklist_S1.md) + [STROBE_checklist_S1.pdf](verification/STROBE_checklist_S1.pdf) (30 pp, 65 KB) — 22 items cross-referenciados a página PDF y sección del manuscrito
- [x] **S2. Compliance tick-sheet** — [verification/compliance_tick_sheet.md](verification/compliance_tick_sheet.md) + [compliance_tick_sheet_S2.pdf](verification/compliance_tick_sheet_S2.pdf) (4 pp, 22 KB) — 32 checks (31 pass, 1 warning documentado: Jensen bias)
- [x] **S3. Verification Appendix** — [verification/Verification_Appendix_S3.md](verification/Verification_Appendix_S3.md) + [Verification_Appendix_S3.pdf](verification/Verification_Appendix_S3.pdf) (9 pp, 55 KB) — pipeline narrativo + DOI manifest con v1.0 + v1.2
- [x] **S4. PTP/IEP Classification Framework v1.0** — [PTP_IEP_Classification_Framework.docx](PTP_IEP_Classification_Framework.docx) (49 KB) — especificación completa §1–§7

## Reproducibility artefacts (Zenodo v1.2 deposit)

### Data artifacts (SHA-256 archivados)
- [x] `hormones_long_tidy.csv` → `93eba357565db8278452f6b16d2ec1e5f0d18f83af9d851a1254123b03c49f70`
- [x] `cohort_normalization_map.csv` → `fd133fb5bab90c1cc2b3e9873a9990cce88f316b0f90261350391c46923e944d`
- [x] `pseudo_ipd_primary_M1000_rho050_cv100.csv` → `c5d398b3a47d5a33430f1889a5b2d8dcba64020ff75b183a7bcf1c2393f1e4ca`
- [x] `pseudo_ipd_subsample_N50_rho050_cv100.csv` → `2e0f066fe51143d7f6ffb3d842cae01831e1c7fa9be8c7c1c368eb3016189cdc`
- [x] `fanova_results.csv` → `a099aa3f6771176325e764a57b2e5645043a6aefeef11d7f48518788d4aa102a`
- [x] `bands_simultaneous.csv` → `9aaa04140bf7c690386807bbf424a1560437a2d7ddde244cf5e2322270b67de6`
- [x] `stability_classification_stage.csv` — B=2000 clasificación
- [x] `bootstrap_stability_results.rds` — B=50 pipeline
- [x] `bootstrap_B2000_results_N300.rds` — frozen at N=300 per empirical convergence (0.29 pp median drift)
- [x] `iep_freq_without_glycemic_subclass.csv` — Tabla 4A provenance
- [x] `iep_freq_with_glycemic_subclass.csv` — Tabla 4B provenance

### Code (R scripts)
- [x] `etl_master_csv.R`, `simulate_pseudo_ipd.R`, `mfaces_dryrun.R`, `fit_mfaces_realdata.R`
- [x] `classify_ptp_iep.R`, `fanova_permutation.R`, `sensitivity_rho.R`, `sensitivity_cv.R`
- [x] `conformal_bands.R`, `bootstrap_stability.R`, `bootstrap_B2000_zenodo.R`, `recover_B2000_results.R`
- [x] `generate_figures.R`, `generate_compliance_and_appendix.R`

### Registro pre-análisis (OSF)
- [x] `preregistration_cohort_map.yaml` frozen 2026-04-22 — ahora incluye v1.0 + v1.2 Zenodo version DOI cross-references
- [x] **OSF DOI**: `10.17605/OSF.IO/3CZRE` (registrado 2026-04-24)

### Environment
- [x] `renv.lock` — generated
- [x] `sessionInfo.txt` — R 4.5.3 arm64 / Apple M4 Pro / macOS Sequoia 15 / BLAS Accelerate

## Licencias y declaraciones (en §6 del manuscript)

- [x] **Licencia CC-BY 4.0** para medRxiv deposit (Plan-S compliant)
- [x] **Declaración de ética** — Not applicable (ecological meta-analysis of publicly available aggregate data)
- [x] **Declaración de consentimiento** — Not applicable
- [x] **Declaración de financiamiento** — "None. The author received no funding for this work."
- [x] **Competing interests** (36-month window) — "The author declares no competing interests within the past 36 months."
- [x] **Contribuciones por autor** (CRediT taxonomy) — Virgen Ayala (sole author): Conceptualization, Methodology, Software, Formal analysis, Investigation, Data curation, Writing — original draft, review and editing, Visualization, Project administration
- [x] **Data availability statement** — §5 Data and Code Availability (concept DOI + v1.2 DOI + v1.0 DOI)
- [x] **Code availability statement** — §5 + GitHub Antropomorfo37/EPP10 tag v1.2
- [x] **Clinical trial registration** — Not applicable (ecological meta-analysis declared in §1.3 Aim)
- [x] **AI/LLM use** — §6: "Analytic pipeline design, R code implementation, and manuscript drafting were supported by Claude (Anthropic, Opus 4.7, April 2026). All statistical outputs, pre-registration decisions, and scientific claims were validated by the author against the v10.0 master prompt and the PTP/IEP framework v1.0 reference documents. AI is disclosed per openRxiv/ICMJE guidance; AI is not listed as an author."
- [ ] **Preprint DOI** — pending medRxiv screening (2–4 business days post-submit)

## DOI manifest (trazabilidad completa)

| Rol | DOI | URL |
|---|---|---|
| Zenodo concept (latest) | `10.5281/zenodo.19743544` | https://zenodo.org/doi/10.5281/zenodo.19743544 |
| **Zenodo v1.2 (current)** | **`10.5281/zenodo.19750294`** | https://zenodo.org/records/19750294 |
| Zenodo v1.0 (historical) | `10.5281/zenodo.19743545` | https://zenodo.org/records/19743545 |
| **OSF pre-registration** | **`10.17605/OSF.IO/3CZRE`** | https://osf.io/3czre |
| GitHub (current) | `Antropomorfo37/EPP10@v1.2` | https://github.com/Antropomorfo37/EPP10/releases/tag/v1.2 |
| medRxiv preprint | pending | — |

## Pipeline operacional para submit (estado actual)

1. [x] **Pre-registration frozen + OSF DOI minted** — 2026-04-22 / 2026-04-24
2. [x] **Code archived at Zenodo v1.0** — 10.5281/zenodo.19743545
3. [x] **Manuscript conceptual revision (v1.2)** — physiologic periprandial responsiveness + dual-track IEP
4. [x] **Zenodo v1.2 minted via GitHub Release webhook** — 10.5281/zenodo.19750294
5. [x] **Supplementaries S1, S2, S3, S4 complete** — 4 archivos listos (S1/S2/S3 with .md + .pdf, S4 .docx)
6. [x] **DOI manifest propagated** — CITATION.cff, preregistration YAML, Verification S3, manuscript §5 all include v1.2
7. [ ] **medRxiv submit** — upload manuscript.pdf + Figures 1-3 PDFs + S1-S4 + metadata form
8. [ ] **Screening** — 2-4 business days → preprint DOI assigned → update OSF + Zenodo + CITATION.cff

## medRxiv submission-form checklist

On https://www.medrxiv.org/submit — uploading the package:

1. **Category**: Endocrinology (primary) + Evidence-Based Medicine (secondary)
2. **Article type**: Research article
3. **License**: CC-BY 4.0
4. **Manuscript file**: `manuscript/manuscript.pdf` (22 pp, 197 KB)
5. **Figures**: `figures/Figure1_cohort_composition.pdf`, `figures/Figure2_classification_inference.pdf`, `figures/Figure3_trajectory_bands.pdf`
6. **Supplementary files**:
   - `verification/STROBE_checklist_S1.pdf`
   - `verification/compliance_tick_sheet_S2.pdf`
   - `verification/Verification_Appendix_S3.pdf`
   - `PTP_IEP_Classification_Framework.docx`
7. **Author**: Héctor Manuel Virgen Ayala — corresponding author (self)
8. **Ethics statement**: Not applicable — ecological meta-analysis of publicly available aggregate data
9. **Funding statement**: None
10. **Competing interests**: None in past 36 months
11. **Data sharing**: "Analysis code and derived cohort-level trajectories archived at Zenodo DOI 10.5281/zenodo.19743544 (concept) / 10.5281/zenodo.19750294 (v1.2) under CC-BY 4.0. Pre-registration at OSF DOI 10.17605/OSF.IO/3CZRE (frozen 2026-04-22)."
12. **AI/LLM use**: Declared per ICMJE — see §6 of manuscript

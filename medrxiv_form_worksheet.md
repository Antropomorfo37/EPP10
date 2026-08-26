# medRxiv submission — ready-to-paste worksheet

**Use este archivo junto con el formulario medRxiv.** Cada sección corresponde a un campo del formulario. Copia literalmente el texto entre los `---` y pégalo en el campo indicado.

**URL de inicio:** <https://www.medrxiv.org/submit-a-manuscript>
**URL directa submission manager:** <https://submit.medrxiv.org>

---

## 0. Antes de empezar — archivos físicos en disco

Ten estos 8 archivos a mano (el formulario te pedirá subir cada uno):

| # | Rol | Path (relativo a la raíz del repositorio) | Tamaño |
|---|---|---|---|
| 1 | Main manuscript (PDF) | `manuscript/manuscript.pdf` | 197 KB, 22 pp |
| 2 | Figure 1 | `figures/Figure1_cohort_composition.pdf` | 25 KB |
| 3 | Figure 2 | `figures/Figure2_classification_inference.pdf` | 31 KB |
| 4 | Figure 3 | `figures/Figure3_trajectory_bands.pdf` | 53 KB |
| 5 | Supplementary S1 — STROBE | `verification/STROBE_checklist_S1.pdf` | 65 KB, 30 pp |
| 6 | Supplementary S2 — Compliance tick-sheet | `verification/compliance_tick_sheet_S2.pdf` | 22 KB |
| 7 | Supplementary S3 — Verification appendix | `verification/Verification_Appendix_S3.pdf` | 55 KB |
| 8 | Supplementary S4 — PTP/IEP framework v1.0 | `PTP_IEP_Classification_Framework.docx` | 49 KB |

---

## 1. Account creation / sign-in

En el primer uso: crear cuenta en medRxiv con email `hectorvirgenmd@gmail.com`.
Validar email de confirmación (puede tardar 1–2 minutos).

Perfil:
- **First name:** Héctor Manuel
- **Last name:** Virgen Ayala
- **Degree:** MD
- **Institution:** Universidad de Guadalajara
- **Department:** Departamento de Clínicas Quirúrgicas
- **ORCID:** 0009-0006-2081-2286 (conectar vía "Link ORCID" — autenticación en orcid.org)

---

## 2. Article type + category

- **Article type:** Research article (Original Article)
- **Primary subject category:** **Endocrinology**
- **Secondary subject category:** **Evidence-Based Medicine** (or "Nutrition" si no aparece EBM)
- **Article format:** Not a clinical trial / systematic review / meta-analysis of intervention trials — this is an *ecological meta-analysis of observational cohort data*. If the form forces a choice, pick "Research article" and clarify in cover letter.

---

## 3. Title (copia literal)

```
Dynamic Enteropancreatic Phenotyping via Sparse mFACEs and PTP/IEP Classification: An Ecological Meta-Analysis of Six Metabolic Cohorts
```

(162 caracteres, sin punto final — medRxiv no pone punto final en el título)

---

## 4. Running title / short title (≤ 75 chars)

```
Dynamic Enteropancreatic Phenotyping via mFACEs + PTP/IEP (6 cohorts)
```

(70 caracteres)

---

## 5. Authors and affiliations

**Único autor (sole author — corresponding author):**

- **Given name:** Héctor Manuel
- **Family name:** Virgen Ayala
- **Degree:** MD
- **Email:** hectorvirgenmd@gmail.com
- **ORCID:** 0009-0006-2081-2286
- **Corresponding author:** YES
- **Presenting author:** YES (único)

**Affiliation 1 (primary):**
```
Departamento de Clínicas Quirúrgicas, Universidad de Guadalajara, Guadalajara, Jalisco, México
```

**Affiliation 2 (secondary):**
```
Departamento de Cirugía General, Hospital General de Zona No. 14, Instituto Mexicano del Seguro Social (IMSS), Guadalajara, Jalisco, México
```

---

## 6. Abstract (339 words — under medRxiv 350-word limit)

Pegar este bloque completo. **Mantener los saltos de línea entre las 5 secciones.**

```
Context. Periprandial trajectories of the enteropancreatic axis (ghrelin, GIP, GLP-1, PYY, glucagon, insulin, glucose — 11 analyte-forms) undergo cohort-specific remodeling under obesity, type 2 diabetes (T2DM), their co-occurrence, and metabolic surgery. Quantitative multi-hormone integration of aggregate published data — capturing the full physiologic periprandial responsiveness across fasting and postprandial phases — remains statistically unresolved.

Objective. To establish six cohort-specific multi-hormone signatures via two complementary routes — joint-multivariate sparse mFACEs with Chiou normalization, and per-analyte Periprandial Transition Profile (PTP) → Integrated Enteropancreatic Pattern (IEP) Type I–V — under a pre-registered protocol. Glycemic subclassification (a/b/c) is applied to dysphysiological/infraphysiological types (III, IV, V) only.

Design. Ecological cohort-time-arm meta-analysis of 71 (Author × Cohort) tuples from 25 publications, harmonized to six canonical cohorts: no_obese_without_T2DM (reference), Obesity, T2DM, Obesity_T2DM, SG, RYGBP. Pseudo-individualization via Gaussian process AR(1) draws (M = 1000; ρ = 0.5 primary; literature-based hormone-specific CV priors). Joint mFACEs fit with Chiou normalization; per-analyte PTP classes aggregated by deterministic precedence to IEP Type, presented in two parallel tracks (with/without glycemic subclassification).

Results. N = 2750 pseudo-subjects; K = 14 mFPC components at FVE 0.90 (N/K = 196); PC1 incretin loading = 0.817. Omnibus Pillai F = 47.6 (p < 0.001; B = 5000 permutations). Pairwise F vs. reference: T2DM 108.3 > RYGBP 86.8 > Obesity_T2DM 59.3 > SG 34.1 > Obesity 22.5. IEP Type IV.II ranking: RYGBP 76 % > Obesity_T2DM 56 % > SG 55 % > Obesity 45 % > T2DM 29 % > reference 19 %. With glycemic subclassification, IV.II.c (dysphysiological + dysglycemic) = 63 % in RYGBP and 56 % in Obesity_T2DM. Dual-route Spearman ρ = 0.50 (complementary, not redundant). Classification stability (B = 2000) = 81.6 % ≥ 0.80 target.

Conclusions. Sparse mFACEs with Chiou normalization plus dual-track PTP/IEP classification is a reproducible framework for extracting multi-hormone signatures from ecological aggregate data. It operationalizes physiologic periprandial responsiveness as a multi-analyte, multi-state alternative to static secretory capacity. The RYGBP signature (IEP IV.II + Enhanced incretin) is operationally distinguishable from T2DM (III.II infra-physiological), supporting functional stratification of bariatric-surgery candidates.
```

Si el formulario rechaza por exceder 350 palabras pese a ser 339, prueba sin la oración final "The RYGBP signature ..." y quedará en ~315 palabras.

---

## 7. Keywords / MeSH terms (6–8 términos)

```
enteropancreatic axis; functional data analysis; sparse mFACEs; PTP/IEP classification; physiologic periprandial responsiveness; bariatric surgery; ecological meta-analysis; type 2 diabetes
```

---

## 8. Declaration statements — texto exacto para pegar

### 8a. Funding statement

```
None. The author received no funding for this work.
```

### 8b. Conflicts of interest / competing interests

```
The author declares no competing interests within the past 36 months.
```

### 8c. Ethics statement / IRB approval

```
Not applicable. This is an ecological meta-analysis of publicly available aggregate cohort-level data from 25 previously published studies; no primary data collection or subject-level contact was performed. Therefore, ethics committee approval and informed consent are not applicable.
```

### 8d. Consent for publication

```
Not applicable (no subject-level data).
```

### 8e. Clinical trial registration

```
Not applicable. This is not a clinical trial; it is a descriptive-inferential ecological meta-analysis.
```

### 8f. Pre-registration / study protocol registration

```
Pre-registered at the Open Science Framework (OSF), registration DOI 10.17605/OSF.IO/3CZRE, frozen 2026-04-22 prior to data lock. Registration URL: https://osf.io/3czre. Full pre-registration YAML and cohort normalization map archived at Zenodo v1.3 DOI 10.5281/zenodo.19758429.
```

### 8g. Data availability statement

```
All derived cohort-level aggregate trajectories, pseudo-IPD draws, classification outputs, and intermediate tables are archived under CC-BY 4.0 at Zenodo: concept DOI 10.5281/zenodo.19743544 (resolves to latest version) and version DOI 10.5281/zenodo.19758429 (v1.3, current). The original master table is retained under the reasonable-controlled-access principle consistent with ICMJE guidance; the full source-study reference list is preserved in Supplementary File S3 (Verification Appendix).
```

### 8h. Code availability statement

```
Complete R analysis code is archived at Zenodo under CC-BY 4.0 (concept DOI 10.5281/zenodo.19743544; version v1.3 DOI 10.5281/zenodo.19758429) and developed openly at https://github.com/sv8wmxnbp8-hash/EPP10 (release tag v1.3). Reproducibility artefacts include renv.lock, sessionInfo.txt, and SHA-256 checksums of all input and output CSVs documented in Supplementary S3 §Post-release DOI assignments.
```

### 8i. Author contributions (CRediT taxonomy)

```
Héctor Manuel Virgen Ayala (sole author): Conceptualization, Methodology, Software, Formal analysis, Investigation, Data curation, Writing — original draft, Writing — review and editing, Visualization, Project administration.
```

### 8j. Acknowledgements

```
The author acknowledges the use of Claude (Anthropic, Opus 4.7, April 2026) for analytic pipeline design, R code implementation support, manuscript drafting assistance, and verification workflows. All statistical outputs, pre-registration decisions, and scientific claims were independently validated by the author against the v10.0 master analysis plan and the PTP/IEP framework v1.0 specification. AI disclosure follows openRxiv and ICMJE 2024 guidance; AI is not listed as an author.
```

### 8k. AI/LLM use disclosure (medRxiv-specific field si existe)

```
Generative AI assistance was used in pipeline design, code implementation, and manuscript drafting (Claude Opus 4.7, Anthropic, April 2026). All scientific decisions and statistical claims were validated by the human author. AI is not listed as an author. Disclosure follows openRxiv/ICMJE 2024 guidance.
```

---

## 9. License selection

- **License:** CC-BY 4.0 (Creative Commons Attribution 4.0 International)
- **Rationale:** Plan-S compliant; required for reuse in systematic reviews and meta-analyses.

---

## 10. Cover letter (optional — recommended)

Algunos formularios medRxiv piden cover letter. Si se pide:

```
Dear medRxiv Editors,

I am submitting the manuscript "Dynamic Enteropancreatic Phenotyping via Sparse mFACEs and PTP/IEP Classification: An Ecological Meta-Analysis of Six Metabolic Cohorts" for consideration as a preprint.

This work develops and validates a reproducible framework for extracting cohort-specific multi-hormone periprandial signatures from ecological aggregate data. It introduces two methodological contributions: (i) operational application of sparse multivariate functional PCA with Chiou normalization to multi-analyte endocrine trajectories, and (ii) the Periprandial Transition Profile (PTP) → Integrated Enteropancreatic Pattern (IEP) Type I–V classification framework. A novel conceptual construct — physiologic periprandial responsiveness — is introduced to replace the static secretory capacity construct, capturing both fasting and postprandial, secretory and inhibitory phases of the enteropancreatic axis under a single multi-analyte multi-state formalism.

The study is fully pre-registered at OSF (DOI 10.17605/OSF.IO/3CZRE, frozen 2026-04-22); analysis code is archived at Zenodo v1.3 (DOI 10.5281/zenodo.19758429) under CC-BY 4.0; reporting follows the STROBE cohort checklist (Supplementary S1); TRIPOD+AI is declared not applicable because the analysis is descriptive-inferential rather than a subject-level prediction model.

The manuscript has not been submitted elsewhere and is not under consideration at any peer-reviewed journal. AI-assisted writing is fully disclosed.

Sincerely,
Héctor Manuel Virgen Ayala, MD
Universidad de Guadalajara & Instituto Mexicano del Seguro Social
Guadalajara, Jalisco, México
ORCID: 0009-0006-2081-2286
Email: hectorvirgenmd@gmail.com
```

---

## 11. Manuscript file upload — order matters

Upload in this order (so the medRxiv web viewer displays the package coherently):

1. **Manuscript (main):** `manuscript.pdf` — role = "Manuscript" / "Main document"
2. **Figures (if separate):** `Figure1_*.pdf`, `Figure2_*.pdf`, `Figure3_*.pdf` — role = "Figure" each. *Note:* medRxiv sometimes auto-embeds figures if they're already in the main PDF. If so, skip and answer "No, figures are embedded in main manuscript."
3. **Supplementary S1:** `STROBE_checklist_S1.pdf` — role = "Supplementary File" — label = "Supplementary File 1 — STROBE cohort checklist"
4. **Supplementary S2:** `compliance_tick_sheet_S2.pdf` — role = "Supplementary File" — label = "Supplementary File 2 — Compliance tick-sheet"
5. **Supplementary S3:** `Verification_Appendix_S3.pdf` — role = "Supplementary File" — label = "Supplementary File 3 — Verification appendix"
6. **Supplementary S4:** `PTP_IEP_Classification_Framework.docx` — role = "Supplementary File" — label = "Supplementary File 4 — PTP/IEP Classification Framework v1.0"

---

## 12. Final submission checklist (re-confirm before clicking "Submit")

- [ ] All 8 files uploaded and roles assigned
- [ ] Title, authors, affiliations match manuscript.pdf page 1
- [ ] Abstract pasted (339 words)
- [ ] Keywords entered
- [ ] Funding / COI / Ethics / Consent / CT / Pre-reg / Data / Code / CRediT / Ack / AI statements all pasted
- [ ] License = CC-BY 4.0 selected
- [ ] ORCID linked
- [ ] Cover letter (if requested) pasted
- [ ] Preview PDF package rendering correctly in medRxiv viewer
- [ ] Email notification enabled
- [ ] Agree to terms (CC-BY, posting policy, no prior publication)

---

## 13. Post-submission

- Screening typically takes **2–4 business days**.
- Once approved, you get a **medRxiv DOI** of form `10.1101/2026.MM.DD.XXXXXXXX`.
- When assigned:
  - [ ] Update CITATION.cff adding the medRxiv DOI
  - [ ] Update OSF registration with the medRxiv link
  - [ ] Update `medrxiv_submission_checklist.md` in the repo
  - [ ] Tag a new Zenodo version v1.3 with the medRxiv DOI injected into the PDF §5

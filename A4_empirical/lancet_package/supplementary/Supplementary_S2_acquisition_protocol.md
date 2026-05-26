# Supplementary S2 — Periprandial trajectory acquisition protocol (FDEP-TP)

## Purpose

This protocol specifies the prospective acquisition protocol that future TRIPOD+AI validation studies of the FDEP-TP framework should follow. It is published under CC-BY 4·0 as a reusable resource for the periprandial endocrinology research community.

## Subject preparation

1. Overnight fast ≥ 10 h prior to morning testing
2. Abstain from alcohol, vigorous exercise ≥ 24 h
3. Patients on glucose-lowering medication: standard washout per local protocol; record washout interval
4. Bedside DPP-IV inhibitor (vildagliptin or sitagliptin) added to all collection tubes for active-GLP-1 / active-GIP preservation

## Caloric stimulus

Standardised liquid mixed-meal test, 500 kcal:

| Macronutrient | % Energy | Mass (g) |
|--------------|----------|----------|
| Carbohydrate (maltodextrin + sucrose) | 56 % | 70 |
| Protein (whey isolate) | 22 % | 27 |
| Fat (sunflower oil) | 22 % | 12 |
| Volume | — | 250 mL water-based drink |

Alternative: 75 g oral glucose tolerance test (OGTT) for restricted-comparison studies; preserve OGTT vs MMT subgroup as analysis covariate.

## Sampling grid

Peripheral venous sampling at the standardised time grid:

| t (min) | Action |
|---------|--------|
| −15, 0 | Fasting baseline (2 samples averaged) |
| 15, 30 | Early postprandial phase |
| 60, 90 | Peak phase |
| 120, 130, 150 | Decay phase |
| 180 | Late postprandial / return-to-baseline phase |

Total: 9 measurement points × 11 hormones = 99 measurements per subject.

## Analyte panel (full FDEP-TP)

| Analyte | Method | Tube | Inhibitor |
|---------|--------|------|-----------|
| Glucose | Enzymatic (hexokinase) | NaF/oxalate | — |
| Insulin | ELISA (Mercodia or equivalent) | Serum | — |
| GLP-1 active (7-36 amide) | ELISA (Mercodia or equivalent) | EDTA + DPP-IV | DPP-IV inhibitor |
| GLP-1 total (7-36 + 9-36) | ELISA (Millipore EZGLP1T-36K) | EDTA + DPP-IV | DPP-IV inhibitor |
| GIP active | ELISA (Millipore EZHGIP-54K) | EDTA + DPP-IV | DPP-IV inhibitor |
| GIP total | ELISA (Millipore EZHGIP-54K) | EDTA | — |
| PYY 3-36 | Multiplex bead-based (Millipore) | EDTA + protease inhibitor | aprotinin |
| PYY total | RIA (Millipore) | EDTA + aprotinin | aprotinin |
| Ghrelin total | RIA (Millipore EZGRT-89K) | EDTA + AEBSF | AEBSF |
| Ghrelin acyl | ELISA (Bertin Pharma A05106) | EDTA + AEBSF | AEBSF |
| Glucagon | Multiplex or sandwich ELISA (Mercodia) | EDTA + aprotinin | aprotinin |

All samples processed within 30 min at 4 °C, plasma aliquoted and stored at −80 °C. **All samples from one subject must be analysed in the same assay batch** to minimise inter-assay variability.

## Cohort recruitment (prospective validation)

| Stratum | Target n | Inclusion |
|---------|----------|-----------|
| LEAN (reference) | ≥ 40 | BMI 18·5-24·9, no T2DM, no GLP-1 RA, no metabolic surgery |
| Obesity (OB) | ≥ 30 | BMI ≥ 30, no T2DM, no GLP-1 RA, no metabolic surgery |
| Obesity + T2DM | ≥ 30 | BMI ≥ 30, T2DM (HbA1c 6·5-9·0 %), no GLP-1 RA |
| Post-CR | ≥ 25 | ≥ 5 % weight loss maintained ≥ 6 months by caloric restriction |
| POST-RYGBP | ≥ 30 | 12 months post-Roux-en-Y gastric bypass |
| POST-SG | ≥ 25 | 12 months post-sleeve gastrectomy |

**Total n ≥ 180** (strata sized to provide ≥ 25 evaluable subjects per joint block after exclusions).

## Analytic pre-registration

- Eigenfunction estimation pipeline: PACE (Yao 2005) per hormone with FACEs sensitivity (Xiao 2018)
- Multivariate joint operator: Happ-Greven (2018) with Chiou normalisation (Chiou 2014)
- Multivariate-FVE component-selection threshold ≥ 0·90 (Golovkine 2025)
- SRSF phase-amplitude registration (Tucker, Wu & Srivastava 2013)
- sup-t bootstrap simultaneous bands (Goldsmith / Degras)

Falsification criteria:

1. **A different K** — multivariate-FVE ≥ 0·90 retains K ≤ 2 or K ≥ 8 components → reject parsimonious four-axis claim
2. **A different leading-axis loading** — Ψ̂₁ not PYY-dominant with integrated L²-norm > 0·5 → reject distal L-cell reading
3. **Absence of cohort separability** — omnibus Pillai p > 0·05 or all pairwise q > 0·05 → reject cohort-discriminating claim
4. **Sensitivity-to-normalisation breakdown** — switching from Chiou normalisation to integrated-variance normalisation reorders leading axes → reject as Chiou artefact rather than physiological feature

## Ethics and consent

- Approval by local Institutional Review Board (IRB) per Declaration of Helsinki
- Written informed consent (model consent at OSF DOI 10.17605/OSF.IO/3CZRE)
- Data shared in fully anonymised form following the Lancet's data-sharing policy

## Repository deposition

All study outputs to be deposited at OSF and Zenodo under CC-BY 4·0. Use the FDEP-TP reference template (Antropomorfo37/EPP10 v1·4) for analysis pipeline reproducibility.

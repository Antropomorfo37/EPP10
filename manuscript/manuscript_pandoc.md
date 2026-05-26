---
title: "Dynamic Enteropancreatic Phenotyping via Sparse mFACEs and PTP/IEP Classification: An Ecological Meta-Analysis of Six Metabolic Cohorts"
author: "Héctor Manuel Virgen Ayala"
date: "2026-04-24"
---

# Introduction

The enteropancreatic axis orchestrates periprandial energy homeostasis through a temporally coordinated choreography among gastric (ghrelin), intestinal L-cell (GLP-1, PYY) and K-cell (GIP), pancreatic (insulin, glucagon) hormones, and the glycemic regulation signal. Cohort-specific remodeling of this choreography under obesity, type 2 diabetes (T2DM), their co-occurrence, and surgical metabolic modalities (sleeve gastrectomy [SG], Roux-en-Y gastric bypass [RYGBP]) defines metabolic phenotypes with greater fidelity than any static marker individually.^[Madsbad, Dirksen & Holst, *Lancet Diabetes Endocrinol* 2014; Jørgensen, Dirksen et al., *Diabetes* 2013; Ryan, Tremaroli et al., *Nature* 2014; Drucker, *Mol Metab* 2022; Pournaras et al., *Obes Surg* 2010; Nauck & Müller, *Diabetologia* 2023.] Two decades of mechanistic characterization have established the qualitative contours, but quantitative multi-hormone integration of heterogeneous aggregate data across published cohorts remains a statistically open problem.

Sparse functional data analysis (FDA) and functional principal component analysis — specifically the sparse PACE algorithm of Yao, Müller and Wang (2005), its multivariate extension by Happ and Greven (2018), and the *mixed-effects fast covariance estimation* of Li, Wang, Liu and Greven (2021) — provide the natural methodological framework to decompose irregularly sampled periprandial trajectories into orthogonal modes of variation (global amplitude PC1, early-vs-late redistribution PC2, higher-order shape PC3). The Chiou, Chen and Yang (2014) normalization preserves balanced per-hormone contribution in the joint decomposition by equalizing integrated variances, operationally indispensable when hormone scales differ across three orders of magnitude (ghrelin ~200 pmol/L vs GLP-1 active ~5 pmol/L). Complementarily, the PTP/IEP classification framework (v1.0, April 2026) formalizes cohort-time-arm taxonomy through six *Periprandial Transition Profiles* per analyte (Preserved, Borderline Impaired, Impaired, Blunted, Borderline Altered, Altered — plus three post-intervention classes Recovered, Borderline Enhanced, Enhanced) and five *Integrated Enteropancreatic Pattern Types* (I–V) with glucose subtypes a/b/c, constituting a deterministic biological taxonomy complementary to the multivariate inferential analysis.

Here we apply both frameworks in integrated form — joint sparse mFACEs with Chiou normalization (primary inferential route via FANOVA permutation, Pillai multivariate contrasts, and sup-t simultaneous bands) and per-analyte PTP → IEP Type I–V (descriptive taxonomy via deterministic precedence) — to an aggregated master table of 71 cohort-time-arm tuples from 25 published studies covering 11 analyte-forms. The operational objective is to establish and statistically validate six cohort-specific multi-hormone signatures — reference `no_obese_without_T2DM`, `Obesity` without T2DM, `T2DM` without obesity, co-occurring `Obesity_T2DM`, post-`SG`, post-`RYGBP` — under a fully pre-registered protocol (OSF DOI [10.17605/OSF.IO/3CZRE](https://osf.io/3czre)), with companion analysis code archived at Zenodo (concept DOI [10.5281/zenodo.19743544](https://zenodo.org/doi/10.5281/zenodo.19743544); version v1.0 [10.5281/zenodo.19743545](https://zenodo.org/records/19743545)). The study is reported in accordance with the STROBE cohort checklist (von Elm et al. 2007); TRIPOD+AI (Collins et al. 2024) is not applicable because the analysis is descriptive-inferential rather than predictive at the individual level.

# Methods

## Design and inferential unit

Ecological meta-analysis at the cohort-time-arm level (v10.0 §4.1). The master table (input SHA-256 `ce2e343d...`) consolidates 71 unique (Author × Cohort) tuples from 25 independent publications. Forty-seven heterogeneous source labels are mapped to six canonical cohorts via pre-registered normalization rules (decisions D1–D4 in the preregistration YAML): pre-surgery cohorts pooled into Obesity or Obesity_T2DM with `surgery_status` metadata (D1); post-bariatric T2DM status captured as `had_t2dm_pre_surgery ∈ {TRUE, FALSE, unknown}` (D2); six caloric-restriction labels pooled into Obesity with a `caloric_restriction_post` tag in a sensitivity column (D3); the “Healthy Middle-aged” obesity descriptor does not define a distinct cohort and is pooled into Obesity (D4).

## Pseudo-individualization and variance priors

The master reports cohort means without SD/SEM; trajectory variance is estimated via pseudo-IPD following v10.0 §4.2. For each (Author, cohort, hormone), M = 1000 synthetic trajectories are sampled from a Gaussian process $\text{GP}(\mu_{\text{cohort}}(t), K_\rho(s, t))$ with continuous-time AR(1) kernel $K_\rho(s, t) = \rho^{|s-t|/30\,\text{min}}$ and marginal variance $\sigma_h^2(t) = (\text{CV}_h \times \mu_h(t))^2$, where $\text{CV}_h$ is the hormone-specific coefficient of variation drawn from literature (upper-bound typical values: ghrelin 0.40, GLP-1 total 0.35, GIP total 0.30, insulin 0.20, glucose 0.12; full table in the preregistration YAML). The primary run uses $\rho = 0.5$ with CV multiplier 1.0; sensitivity analyses cover $\rho \in \{0.3, 0.7, 0.9\}$ and CV multipliers $\{0.75, 1.25\}$. A deterministic subsample of N = 50 pseudo-subjects per arm (2750 total across 55 productive arms) is used for the mFACEs joint fit; the full M = 1000 archive is deposited on Zenodo.

## Three-stage functional analysis

**Stage 1 — univariate (sensitivity, v10.0 §6.2).** Each of 11 analyte-forms (ghrelin total, ghrelin acyl, GIP total, GIP active, GLP-1 total, GLP-1 active, PYY total, PYY 3-36, glucagon, insulin, glucose) is modeled separately by PACE (Yao, Müller & Wang 2005) on the dynamic component $d_{ih}(t) = x_{ih}(t) - B_{ih}$ (v10.0 §7.1), with baseline $B_{ih}$ defined as the strict $t = 0$ value for the reference cohort and the nearest-to-zero value for the study cohorts. Minimum 3 observations per subject per hormone required; FVE 0.95.

**Stage 2 — mFACEs joint (primary, v10.0 §7.4).** Adaptation of Li, Wang, Liu and Greven (*Stat Med* 2021): per-analyte FACE fits via `face::face.sparse`; pairwise cross-covariance blocks estimated by 2-D tensor-product penalized splines (`mgcv::te`); explicit Chiou normalization applied at block assembly (each off-diagonal block $C_{hh'}$ scaled by $1/\sqrt{v_h v_{h'}}$ where $v_h = \int \text{diag}(C_h)(t) \, dt$), preserving balanced contribution across hormones with disparate absolute variances. BLUP score recovery retains measurement-error variance $\sigma^2$ in the observation model $\Sigma_i = C_i + \text{diag}(\sigma^2_i)$ with adaptive ridge regularization (multipliers 1×, 10×, 100×, 1000× to resolve numerical singularity). Joint eigendecomposition of the $H \times n_{\text{Grid}}$ block covariance produces K retained components satisfying FVE ≥ 0.90 (primary) or 0.95 (sensitivity); identifiability requires $N/K_{\text{eff}} > 10$ (v10.0 §7.3).

**Stage 3 — inference and dual classification.** Formal inference by FANOVA permutation (B = 5000 cohort-label permutations per retained PC; BH-FDR α = 0.05 stratified by test type) plus Pillai multivariate contrasts on the retained score block. Simultaneous 95 % bands for $\hat\mu_{\text{cohort}}(t) - \hat\mu_{\text{ref}}(t)$ via sup-t bootstrap (Montiel Olea & Plagborg-Møller 2019, B = 2000). Dual classification layers:

- **Joint-mFPC 6-class** on the first three retained multivariate scores $(\xi_1, \xi_2, \xi_3)$: pseudo-subjects z-standardized against the `no_obese_without_T2DM` reference distribution, then assigned to one of six classes (Preservado, Impairment limítrofe, Impaired, Blunted, Enhanced, Altered) with pre-registered precedence *Altered > Blunted/Enhanced > Impaired > Impairment limítrofe > Preservado*; the *incretin axis* is identified by argmax over the first three PCs of the integrated squared loading on {GIP, GLP-1, PYY} forms (warning flag if loading < 0.50).

- **Per-analyte PTP → IEP Type I–V** (framework v1.0): classification vector $Z_{ih} = [z(B_{ih}), z(\xi_{ih1}), z(\xi_{ih2}), z(\xi_{ih3})]$ from univariate PACE fits referenced to percentiles of the `no_obese_without_T2DM` distribution, six Primary PTP classes assigned with pathophysiological filter per v1.0 §3.6–§3.7 (Altered restricted to GIP, insulin, glucose, ghrelin; upper-tail GLP-1/PYY/glucagon routed to Enhanced/Borderline Enhanced when glucose subtype = a). Post-bariatric cohorts receive Secondary PTPs (Recovered, Borderline Enhanced, Enhanced). Cohort-time-arm IEP Types assigned by deterministic precedence (Type I–V with glucose subtypes a/b/c) per §6.4.

## Validation

Classification-stage bootstrap (B = 2000): scores from the primary fit held fixed; reference cohort resampled with replacement to perturb the z-score reference distribution; modal class per pseudo-subject + stability = Pr(modal). Pre-registered target: ≥ 80 % of subjects with stability ≥ 0.80. Pipeline-stage bootstrap (B = 50 primary + B = 300 extended Zenodo archive): full re-fit PACE → mFACEs → classification on resampled pseudo-IPD; convergence of prevalence distributions verified vs. B = 50. Sensitivity matrices across $\rho$ and CV multipliers report K, incretin loading, and qualitative cohort ordering.

## Computational environment

R 4.5.3 arm64 on Apple M4 Pro (macOS Sequoia 15), BLAS/LAPACK via Apple Accelerate (vecLib); parallelism via `future::plan(multisession, workers = 8)` to circumvent the January 2026 fork-safety bug with `parallel::mclapply`. Key packages: `fdapace` 0.6.0, `face` 0.1-7, `MFPCA` 1.3-11, `funData` 1.3-9, `conformalInference.fd` 1.1.1, `roahd` 1.4.3, `mgcv` 1.9, `pointblank` 0.12, `irr` 0.84. Seed 20260422 across all stochastic procedures. `renv.lock` and session info deposited on Zenodo.

# Results

## Cohorts and data coverage

After pre-registered normalization, the ETL (output SHA-256 `7ed3b058...`) produced 2917 long-format observations across six canonical cohorts: `no_obese_without_T2DM` (13 arms, 692 obs), `Obesity` (14 arms, 777 obs), `T2DM` (4 arms, 280 obs), `Obesity_T2DM` (7 arms, 528 obs), `SG` (7 arms, 264 obs), `RYGBP` (10 arms, 376 obs). Hormone coverage is asymmetric: glucose 508 obs, insulin 460, GIP total 380, GLP-1 active 315, GLP-1 total 300, PYY total 215, ghrelin acyl 199, ghrelin total 192, PYY 3-36 152, glucagon 100, GIP active 96 (Table 1, Figure 1).

## Pseudo-IPD and primary mFACEs fit

M = 1000 pseudo-subject trajectories per arm were drawn from $\text{GP}(\mu, K_{\rho = 0.5})$; deterministic subsampling to N = 50 per arm yielded 2750 pseudo-subjects in 55 productive arms. The joint mFACEs fit with explicit Chiou normalization retained K = 14 components at FVE 0.90 (N/K = 196, §7.3 criterion cleared with large margin) and K = 21 at FVE 0.95 (N/K = 131). FVE top-5 = {0.16, 0.28, 0.38, 0.48, 0.55}, confirming balanced variance distribution post-Chiou (vs. pre-Chiou in which PC1 alone captured 83.7 % of variance by ghrelin/insulin dominance and the incretin loading was 0.000). Eighty-nine percent of the 110 off-diagonal cross-covariance blocks were estimated via the `mgcv::te` penalized spline with ≥ 50 paired observations; the remaining 11 % were retained at zero. Zero BLUP failures across 2750 pseudo-subjects; ridge regularization triggered informationally but never escalated to missing scores.

## Incretin axis and dual classification

The argmax rule identified **PC1 as the incretin axis with integrated squared loading 0.817** on {GIP total, GIP active, GLP-1 total, GLP-1 active, PYY total, PYY 3-36}; the first three components all exceed 0.60 incretin loading (PC1: 0.817, PC2: 0.699, PC3: 0.608), confirming that incretin variation dominates the first joint subspace after balancing hormone scales. The joint 3-PC classifier (Table 2) yields biologically coherent cohort-specific prevalence patterns (Figure 2A).

## Formal inference

Omnibus Pillai test across 14 retained components produced **F = 47.6** (p < 0.001 permutation; all 14 per-PC F-tests significant under BH-FDR < 0.05). Pairwise Pillai contrasts vs. reference ordered the separation as **T2DM (F = 108.3) > RYGBP (86.8) > Obesity_T2DM (59.3) > SG (34.1) > Obesity (22.5)** (Table 3 Panel A, Figure 2B), all significant at the minimum possible permutation p-value 1/B+1 = 0.0002 (saturation expected given N = 2750). The per-PC F surface (Table 3 Panel B, Figure 2C) reveals **PC1 and PC7 as dominant discriminators**: RYGBP vs. reference attains F = 496 on PC1 (Table 3 Panel B), T2DM vs. reference F = 370 on PC7 — with PC7 emerging as a metabolic–diabetic disruption axis specifically salient in T2DM and in Obesity_T2DM. Coverage of significant PCs: omnibus 14/14, Obesity_T2DM 11/14, T2DM 10/14, RYGBP 9/14, Obesity 8/14, SG 8/14.

## Trajectory differences and simultaneous bands

Of 47 hormone × cohort-vs-reference comparisons, **18 were significant at 95 % simultaneous coverage across the entire 0–180 min grid** (sup-t bootstrap, B = 2000). Maximum absolute differences: ghrelin total SG $\Delta = -143$ pmol/L (fundectomy signature); insulin Obesity $\Delta = +975$ pmol/L (hyperinsulinemic resistance); insulin RYGBP $\Delta = +665$ pmol/L (exaggerated post-surgical response); glucose Obesity_T2DM $\Delta = +8.3$ mmol/L; ghrelin acyl SG $\Delta = -97$ pmol/L; GLP-1 total T2DM $\Delta = -12.8$ pmol/L (diabetic incretin deficiency); PYY 3-36 Obesity_T2DM $\Delta = +32.5$ pmol/L (Figure 3).

## PTP/IEP Type I–V taxonomy and cross-framework concordance

Applying the per-analyte PTP framework with the v1.0 §3.7 pathophysiological filter and the deterministic precedence hierarchy yields the IEP Type distribution in Table 4. **Type IV.II (dysphysiological marked) ranking**: Obesity_T2DM (48.0 %) > RYGBP (38.2 %) > T2DM (28.5 %) > Obesity (22.6 %) > SG (15.7 %) > reference (11.5 %). Classifiability rate ranges from 28.6 % (SG — arms lacking insulin + gut hormone combinations) to 100 % (T2DM — fully covered), reflecting heterogeneous source-study hormone panels.

The Spearman rank correlation between joint Pillai F (cohort vs. reference) and IEP Type IV.II prevalence (N = 5 non-reference cohorts) is **ρ = 0.50**, a moderate concordance rather than near-redundancy (Figure 2D). The two rails identify the same extreme cohorts but with interpretable rank discordance: T2DM ranks #1 on Pillai F but #3 on IV.II because its signal is dominated by **infra-physiological** Type III.II (37 %) rather than dysphysiological Type IV.II — consistent with diabetic incretin deficiency. Obesity_T2DM inverts the pattern, ranking #1 on IV.II prevalence while #3 on Pillai F — a convergent multi-analyte dysphysiology with moderate multivariate separation.

## Stability and robustness

Classification-stage bootstrap (B = 2000) yielded **81.6 % of 2750 pseudo-subjects with modal-class stability ≥ 0.80**, exceeding the pre-registered target. All six cohorts individually achieved 77.5–84.0 %. Pipeline-stage bootstrap (B = 50 pipeline re-fits): 50/50 reps converged; K median = 14 (range 13–15); incretin axis PC1 invariant; loading median = 0.813 (IQR 0.795–0.829); prevalence signatures reproducible with IQR 3–5 pp in modal classes. Extended Zenodo archive at B = 300: 0.29 pp median drift vs. B = 50 across 36 cohort × class pairs and IQR ratio 1.02, empirically confirming that B = 50 already saturated the bootstrap variance and justifying the freeze at N = 300 for the v1.0 deposit. Sensitivity matrices: $\rho \in \{0.3, 0.7, 0.9\}$ preserves K $\in$ [11, 14] and incretin PC1 invariant (loadings 0.79–0.82); CV × $\{0.75, 1.25\}$ preserves K $\in$ [12, 16] and incretin loadings 0.79–0.83.

# Discussion

**Principal finding.** Sparse mFACEs with Chiou normalization jointly applied to eleven periprandial analyte-forms from 71 published cohort-time-arms, paired with the per-analyte PTP → IEP Type I–V classification framework, reveals six cohort-specific multi-hormonal signatures that are statistically robust under two orthogonal classification routes. Both routes converge on the same cohort extremes — reference minimum Type IV.II prevalence (11.5 %) versus Obesity_T2DM/RYGBP/T2DM dominance (29 %–48 %) — while preserving interpretable rank discordance in the rails' internal orderings (Spearman ρ = 0.50). The moderate concordance demonstrates that the two frameworks capture complementary rather than redundant information: joint Pillai F measures multivariate distance from reference; IEP Type IV.II prevalence measures per-analyte dysphysiology. Reporting both strengthens the scientific claim by triangulating cohort signatures from orthogonal methodological angles.

**Post-bariatric signature.** RYGBP combines 37 % Altered + 12 % Enhanced in the joint rail and 38 % Type IV.II in the IEP rail — the only cohort with the concurrent dysphysiology-plus-incretin-amplification signature reproducing the exaggerated GLP-1 response documented mechanistically by Jørgensen, Dirksen et al. (*Diabetes* 2013) via Exendin-9 blockade. SG shows an operationally distinguishable profile: the most pronounced ghrelin suppression across all cohorts (sup-t band $|\Delta| = 143$ pmol/L, significant across 100 % of 0–180 min), consistent with the fundectomy–Piezo1 mechanism of Ryan, Tremaroli et al. (*Nature* 2014) and Zhang et al. (*Sci Signal* 2024). The distinguishability of SG from RYGBP supports the v10.0 §3 stance against pooling surgical modalities — the joint signatures reflect genuinely different remodeling phenotypes.

**Diabetic and convergent signatures.** T2DM without obesity shows 28 % Impaired + 9 % Blunted in the joint rail and 37 % Type III.II (infra-physiological moderate-severe) in the IEP rail — a signature of dominant incretin secretion deficiency compatible with the GIP paradox update in Nauck & Müller (*Diabetologia* 2023) and Gutgesell et al. (*Nat Metab* 2025). Obesity_T2DM presents a **convergent** signature with IEP Type IV.II prevalence at 48 % (the highest in the dataset), indicating partially synergistic disruption beyond the sum of individual cohort patterns; this justifies the v10.0 §2.14.3 discordance-escalation rule applied specifically to this arm.

**Methodological contributions.** Three contributions. (i) Chiou normalization within joint mFACEs is operationally indispensable: without it, PC1–PC2 are dominated by ghrelin and insulin absolute variance and the incretin axis vanishes (loading 0.000); after explicit Chiou scaling, K = 14 components are retained at FVE 0.90 and PC1 concentrates 81.7 % of its squared loading on incretins. (ii) Dual classification (joint-mFPC + per-analyte PTP/IEP) triangulates signatures by orthogonal methodologies, producing a defensible primary–secondary reporting structure. (iii) Practical adaptation of mFACEs via per-analyte FACE + 2-D tensor-product cross-covariance splines + BLUP with σ² retention preserves the key property that distinguishes mFACEs from two-stage PACE → Happ-Greven (smoothing away of measurement-error variance in the reconstruction step).

**Limitations.** Four material limitations bound interpretation. *First*, the master reports only cohort means; trajectory variance is inferred from literature-based CV priors, with the CV × {0.75, 1.25} sensitivity bounding the propagated effect on the primary conclusions. *Second*, inference is ecological pseudo-IPD (v10.0 §4.2); N = 2750 saturates permutation p-values at the minimum B⁻¹ ≈ 0.0002, so we report F magnitudes as the comparative metric of cohort separation. *Third*, the per-analyte PTP framework applied across 10 non-glucose analytes with ±P5/P95 thresholds assigns ~40 % of reference subjects to Type IV by multiple-testing arithmetic (Pr(any of 10 analytes > P95) ≈ 40 %); we apply the v1.0 §3.7 pathophysiological filter (Altered restricted to GIP, insulin, glucose, ghrelin; upper-tail GLP-1/PYY/glucagon routed to Enhanced in post-intervention context), which reduces but does not eliminate this arithmetic floor (residual reference Type IV prevalence: 40.9 %). *Fourth*, `had_t2dm_pre_surgery` is operationally binary {TRUE, unknown} in this dataset (0 rows of explicit FALSE), limiting sensitivity analysis on pre-surgical T2DM stratification within the post-bariatric cohorts — a limitation of source reporting that would benefit from IPD access.

**Clinical implications.** The RYGBP signature (IEP Type IV.II with Enhanced incretin component when glucose is preserved) is operationally distinguishable from the T2DM pure signature (Type III.II infra-physiological) and from the Obesity_T2DM convergent signature (Type IV.II dominant). This suggests immediate clinical application: functional baseline stratification of bariatric surgery candidates by PTP/IEP pre-operative signature could identify candidates more likely to achieve Type II (supra-physiological coherent) post-surgery versus persisting in Type IV. External validation in prospective longitudinal IPD from surgical cohorts is the natural next step.

# Data and Code Availability

Analysis code under CC-BY 4.0 at:

- Zenodo concept DOI: [10.5281/zenodo.19743544](https://zenodo.org/doi/10.5281/zenodo.19743544) (resolves to latest version)
- Zenodo v1.0 DOI: [10.5281/zenodo.19743545](https://zenodo.org/records/19743545)
- GitHub repository: [Antropomorfo37/EPP10 v1.0](https://github.com/Antropomorfo37/EPP10/releases/tag/v1.0)

Pre-registration with frozen analysis plan (2026-04-22) at:

- OSF registration DOI: [10.17605/OSF.IO/3CZRE](https://osf.io/3czre)
- OSF project: [https://osf.io/tr469](https://osf.io/tr469)

Aggregate cohort-level data derived via the ETL are archived alongside the code (`hormones_long_tidy.csv`). The original master table (Tabla maestra AUC E and P Hormones, input SHA-256 `ce2e343d...`) is retained under controlled access consistent with ICMJE data-sharing expectations; source-study reference list is preserved in the Verification Appendix S3 bibliography.

# Reporting Guidelines and Declarations

**Reporting guideline.** Reported in accordance with STROBE cohort (von Elm et al., *Lancet* 2007) — completed checklist in Supplementary File S1. TRIPOD+AI (Collins et al., *BMJ* 2024) declared not applicable (this is a descriptive-inferential mFPCA analysis, not a multivariable clinical prediction model).

**Ethics and consent.** Not applicable (ecological meta-analysis of publicly available aggregate data from peer-reviewed sources; no primary data collection).

**Funding.** [To be completed by author.]

**Competing interests (36-month window).** [To be completed by author.]

**Author contributions (CRediT).** [To be completed; sole author: Héctor Manuel Virgen Ayala — Conceptualization, Methodology, Software, Formal analysis, Investigation, Data curation, Writing (original draft, review, editing), Visualization, Project administration.]

**AI/LLM disclosure.** Analytic pipeline design, R code implementation, and manuscript drafting were supported by Claude (Anthropic, Opus 4.7, April 2026). All statistical outputs, pre-registration decisions, and scientific claims were validated by the author against the v10.0 master prompt and the PTP/IEP framework v1.0 reference documents. AI is disclosed per openRxiv/ICMJE guidance; AI is not listed as an author.

**License.** CC-BY 4.0 (Plan-S compliant) for manuscript, code, and derived aggregate data.

# References

*(Bibliography compiled in `references.bib`. Key references cited throughout the manuscript include: Yao, Müller & Wang JASA 2005; Happ & Greven JASA 2018; Li, Wang, Liu & Greven Stat Med 2021; Chiou, Chen & Yang Stat Sin 2014; Goldsmith, Greven & Crainiceanu Biometrics 2013; Madsbad, Dirksen & Holst Lancet Diabetes Endocrinol 2014; Jørgensen et al. Diabetes 2013; Ryan et al. Nature 2014; Drucker Mol Metab 2022; Pournaras et al. Obes Surg 2010; Nauck & Müller Diabetologia 2023; Laferrère Int J Obes 2011/2016; Bojsen-Møller et al. Int J Obes 2023; Collins et al. BMJ 2024; von Elm et al. Lancet 2007; Montiel Olea & Plagborg-Møller QE 2019.)*

# Tables

**Table 1.** Cohort and hormone coverage. Top panel: source-label counts, productive arms, long-format observations, and pseudo-subject sample per canonical cohort. Bottom panel: long-format observations per analyte-form with hormone-specific CV priors (literature-anchored) used for pseudo-IPD variance.

| Cohort v10.0 | Source labels | Arms | Obs (long) | Pseudo-subj. N |
|---|---|---|---|---|
| `no_obese_without_T2DM` | 1 | 13 | 692 | 650 |
| Obesity | 17 | 14 | 777 | 700 |
| T2DM | 2 | 4 | 280 | 200 |
| Obesity_T2DM | 10 | 7 | 528 | 350 |
| SG | 7 | 7 | 264 | 350 |
| RYGBP | 10 | 10 | 376 | 500 |
| **Total** | **47** | **55** | **2917** | **2750** |

| Hormone | n obs | CV prior | Unit |
|---|---|---|---|
| glucose | 508 | 0.12 | mmol/L |
| insulin | 460 | 0.20 | pmol/L |
| GIP_total | 380 | 0.30 | pmol/L |
| GLP1_active | 315 | 0.40 | pmol/L |
| GLP1_total | 300 | 0.35 | pmol/L |
| PYY_total | 215 | 0.35 | pmol/L |
| ghrelin_acyl | 199 | 0.45 | pmol/L |
| ghrelin_total | 192 | 0.40 | pmol/L |
| PYY_3_36 | 152 | 0.40 | pmol/L |
| glucagon | 100 | 0.30 | pmol/L |
| GIP_active | 96 | 0.35 | pmol/L |

\newpage

**Table 2.** Joint-mFPC 6-class classification prevalence per cohort (% median [Q1–Q3] across B = 50 pipeline bootstrap re-fits). Modal class per cohort in bold. Classifier applied to the first three retained joint mFPC components; Chiou-normalized; FVE 0.90 retention threshold.

| Cohort | Preservado | Imp. lim. | Impaired | Blunted | Enhanced | Altered |
|---|---|---|---|---|---|---|
| `no_obese_without_T2DM` | **41.7** [39.5–42.8] | 29.5 [27.5–31.3] | 16.2 [14.8–17.9] | 1.6 [0.9–2.1] | 1.7 [1.1–2.8] | 8.8 [7.3–11.2] |
| Obesity | **33.9** [31.9–36.6] | 29.8 [28.2–31.2] | 22.2 [20.5–23.8] | 3.7 [1.3–6.5] | 2.1 [1.2–5.0] | 6.0 [3.9–9.5] |
| T2DM | 20.2 [18.5–24.0] | 24.2 [21.1–27.0] | **27.8** [25.1–30.0] | 9.0 [7.0–15.8] | 7.5 [1.9–12.9] | 13.8 [10.5–19.2] |
| Obesity_T2DM | 19.4 [16.6–20.8] | 23.9 [22.0–26.2] | **27.6** [25.2–29.6] | 10.9 [0.9–14.9] | 4.9 [0.6–11.1] | 15.0 [11.7–19.2] |
| SG | **25.4** [22.1–29.5] | 24.7 [23.2–27.1] | 22.6 [21.7–25.0] | 2.3 [0.6–8.9] | 4.9 [0.6–13.0] | 14.7 [10.4–20.6] |
| RYGBP | 9.2 [7.6–11.4] | 15.4 [13.6–16.4] | 20.3 [18.8–22.0] | 2.4 [0.4–16.4] | 12.4 [0.5–21.0] | **37.1** [34.4–40.2] |

\newpage

**Table 3.** Formal inference from FANOVA permutation and Pillai multivariate. Panel A: omnibus and pairwise contrasts. Panel B: top-three PC-wise F statistics per contrast. Panel C: coverage of significant PCs under BH-FDR 0.05. All permutation p-values saturate at 1/B+1 ≈ 0.0002 given N = 2750 pseudo-subjects; F magnitudes are the interpretable comparative metric.

**Panel A — Pillai multivariate contrasts**

| Contrast | Pillai F | p_perm | p_adj_BH |
|---|---|---|---|
| Omnibus (6 cohorts) | 47.6 | < 0.001 | < 0.001 |
| T2DM vs. ref | **108.3** | < 0.001 | < 0.001 |
| RYGBP vs. ref | 86.8 | < 0.001 | < 0.001 |
| Obesity_T2DM vs. ref | 59.3 | < 0.001 | < 0.001 |
| SG vs. ref | 34.1 | < 0.001 | < 0.001 |
| Obesity vs. ref | 22.5 | < 0.001 | < 0.001 |

**Panel B — Dominant per-PC F statistics**

| Contrast | Dominant PC | Secondary PC | Tertiary PC |
|---|---|---|---|
| RYGBP vs. ref | PC1 (F = **496**) | PC7 (F = 275) | PC4 (F = 244) |
| T2DM vs. ref | PC7 (F = **370**) | PC1 (F = 153) | PC3 (F = 104) |
| Obesity_T2DM vs. ref | PC1 (F = 273) | PC7 (F = 223) | PC2 (F = 147) |
| SG vs. ref | PC1 (F = 139) | PC2 (F = 88) | PC11 (F = 83) |
| Obesity vs. ref | PC1 (F = 102) | PC2 (F = 97) | PC7 (F = 36) |

**Panel C — Coverage of significant PCs (FDR < 0.05 of 14 retained)**

| Contrast | PCs sig. / 14 |
|---|---|
| Omnibus | **14 / 14** |
| Obesity_T2DM vs. ref | 11 / 14 |
| T2DM vs. ref | 10 / 14 |
| RYGBP vs. ref | 9 / 14 |
| Obesity vs. ref | 8 / 14 |
| SG vs. ref | 8 / 14 |

\newpage

**Table 4.** PTP/IEP Type I–V distribution per cohort (%) under the deterministic precedence hierarchy of framework v1.0 §6.4, with pathophysiological filter §3.6–§3.7 applied. Type IV.II (dysphysiological marked) rankings replicate the joint Pillai F cohort ordering with moderate Spearman concordance (ρ = 0.50 across the 5 non-reference cohorts).

| Cohort | I.I | III.I | III.II | IV.I | **IV.II** | V.I | not_integrable |
|---|---|---|---|---|---|---|---|
| `no_obese_without_T2DM` (ref) | 4.2 | 7.8 | 7.7 | 29.4 | **11.5** | 0.8 | 38.5 |
| Obesity | 1.4 | 1.6 | 2.1 | 22.3 | **22.6** | 0.0 | 50.0 |
| T2DM | 0.0 | 0.0 | **37.0** | 34.5 | 28.5 | 0.0 | 0.0 |
| Obesity_T2DM | 0.0 | 3.1 | 5.7 | 28.9 | **48.0** | 0.0 | 14.3 |
| SG | 0.0 | 0.6 | 0.0 | 12.3 | 15.7 | 0.0 | 71.4 |
| RYGBP | 0.0 | 0.2 | 0.4 | 11.2 | **38.2** | 0.0 | 50.0 |

\newpage

# Figures

**Figure 1.** Cohort composition and data-flow characteristics.

![Figure 1 — Cohort composition and data-flow characteristics](../figures/Figure1_cohort_composition.pdf){#fig-1 width=100%}

\newpage

**Figure 2.** Classification prevalence, Pillai multivariate inference, per-PC F structure, and cross-framework concordance.

![Figure 2 — Classification prevalence and multivariate inference](../figures/Figure2_classification_inference.pdf){#fig-2 width=100%}

\newpage

**Figure 3.** Trajectory differences with 95 % simultaneous sup-t bands for six key analyte-forms.

![Figure 3 — Trajectory differences with sup-t simultaneous bands](../figures/Figure3_trajectory_bands.pdf){#fig-3 width=100%}

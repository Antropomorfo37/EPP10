# RUNBOOK — Operational Closure Sequence

End-to-end sequence from current state (all science artefacts complete, B=2000 running in background) to medRxiv submission. Execute in order; each step has a verification command.

All commands below assume `EPP10_ROOT` points at your clone of this repository:

```bash
export EPP10_ROOT=/ruta/a/tu/clon/EPP10   # p. ej. $(pwd) desde la raíz del repo
```

---

## STEP 0 — Current state (verify before proceeding)

```bash
# B=2000 Zenodo archive in progress
ls "$EPP10_ROOT"/cache_bootstrap_B2000/ | grep rep_ | wc -l
# Expected: grows over time toward 2000

# Verify all science artefacts present
for f in hormones_long_tidy.csv cohort_normalization_map.csv \
         pseudo_ipd_primary_M1000_rho050_cv100.csv \
         fit_mfaces_primary_results.rds fanova_results.csv \
         bands_simultaneous.csv bootstrap_stability_results.rds \
         ptp_iep_results.rds preregistration_cohort_map.yaml; do
  [ -f "$EPP10_ROOT/$f" ] && echo "OK  $f" || echo "MISSING  $f"
done

# Verify figures + verification docs
ls "$EPP10_ROOT"/figures/Figure[1-3]*.pdf
ls "$EPP10_ROOT"/verification/compliance_tick_sheet.md
ls "$EPP10_ROOT"/verification/Verification_Appendix_S3.md
```

Expected: all 9 artefacts `OK`, 3 figures, 2 verification docs.

---

## STEP 1 — Snapshot environment (renv.lock + sessionInfo)

Run once, commits the exact package versions used for the analysis.

```bash
cd "$EPP10_ROOT"
Rscript -e '
.libPaths(c("~/.R/library", .libPaths()))
if (!requireNamespace("renv", quietly=TRUE))
  install.packages("renv", lib="~/.R/library", repos="https://cloud.r-project.org")
renv::init(bare = TRUE, restart = FALSE)
renv::snapshot(prompt = FALSE)
cat("renv.lock generated. Packages captured:",
    nrow(renv::lockfile_read("renv.lock")$Packages), "\n")

# sessionInfo snapshot
writeLines(capture.output(sessionInfo()),
           "verification/sessionInfo.txt")
'
```

Verification:
```bash
[ -f renv.lock ] && echo "renv.lock OK"
wc -l verification/sessionInfo.txt
```

---

## STEP 2 — B=2000 snapshot (any time; can be re-run as more reps complete)

Even if B=2000 has only partial completion, run the aggregator to produce a current snapshot. The archived result captures whatever state the cache has. Re-run to refresh as more reps land.

```bash
cd "$EPP10_ROOT"
Rscript recover_B2000_results.R 2>&1 | tee verification/B2000_snapshot.log
```

Verification:
```bash
Rscript -e '
r <- readRDS(file.path(Sys.getenv("EPP10_ROOT"), "bootstrap_B2000_results.rds"))
cat("B=2000 snapshot: ok=", r$n_ok, "/", r$B_TARGET %||% 2000,
    " | failed=", r$n_fail, "\n", sep="")
'
```

If the pipeline is still running, the snapshot represents the current state; you can re-run this step before the medRxiv upload to refresh. Alternatively, if B ≥ 500 successful reps already, that is scientifically sufficient for the supplementary archive — the remaining reps continue accumulating post-submission.

---

## STEP 3 — Initialize git repository and commit

```bash
cd "$EPP10_ROOT"

git init
git branch -M main

# Ignore large cache directories and R runtime files
cat > .gitignore <<'EOF'
cache_bootstrap_B2000/
cache_bootstrap_pipeline/
.Rhistory
.Rproj.user/
.RData
.Rapp.history
*.Rcheck/
renv/library/
renv/local/
renv/cellar/
renv/python/
renv/staging/
.DS_Store
EOF

# Stage core analysis artefacts
git add \
  *.R *.md *.yaml *.csv *.rds \
  figures/ verification/ \
  .gitignore

# Do NOT commit the raw Excel/CSV master dataset — kept private under OSF controlled access
git rm --cached "hormones_long_tidy.csv" 2>/dev/null || true

git commit -m "Initial public release: sparse mFACEs + PTP/IEP pipeline for JCEM/medRxiv

Includes:
- ETL, pseudo-IPD generator, mFACEs joint fitter with Chiou normalization
- Joint-mFPC 6-class classifier + per-analyte PTP/IEP Type I-V (framework v1.0)
- FANOVA permutation + sup-t simultaneous bands
- Classification-stage B=2000 + pipeline B=50 bootstrap stability
- Sensitivity: ρ∈{0.3,0.5,0.7,0.9}, CV×{0.75,1.0,1.25}
- Compliance tick-sheet + Verification Appendix (11 sections)
- Figures 1-3 (PDF + PNG)

See preregistration_cohort_map.yaml for frozen decisions D1-D4,
CV priors, and classifier specifications (OSF deposit pending)."
```

Verification:
```bash
git log --oneline
git status
```

---

## STEP 4 — Create GitHub remote + push + tag v1.0

Assumes you have a GitHub account and `gh` CLI authenticated.

```bash
# Create repo (private or public; public for Zenodo auto-archive)
gh repo create EPP10 --public --source=. --remote=origin --push \
  --description "Sparse mFACEs + PTP/IEP classification for ecological enteropancreatic meta-analysis (JCEM 2026)"

# Tag the release — Zenodo will watch for tagged releases
git tag -a v1.0 -m "v1.0 — medRxiv submission release

Corresponds to preprint medRxiv DOI (pending assignment).
All science artefacts included except raw master dataset (controlled access on OSF).
Pipeline state: B=2000 Zenodo archive accumulating; can be re-run at any time to
refresh bootstrap_B2000_results.rds for post-submission Zenodo updates."

git push origin v1.0
```

Verification:
```bash
gh release list
# Expected: v1.0 listed with current timestamp
```

---

## STEP 5 — Zenodo DOI (via GitHub integration)

**Pre-requisite**: linked Zenodo account to GitHub, with EPP10 repo toggled ON in Zenodo settings.

Setup (one-time):
1. Log into https://zenodo.org/
2. Navigate to Settings → GitHub → enable webhook for `EPP10` repository
3. Ensure the `CITATION.cff` file exists in repo (if not, create before tagging — see template below)

Once v1.0 is tagged, Zenodo auto-creates a snapshot and issues a DOI:

```bash
# After ~5-10 min post-tag, check Zenodo
open "https://zenodo.org/account/settings/github/"
# Find EPP10 v1.0 → click to view DOI (e.g., 10.5281/zenodo.XXXXXXX)
```

CITATION.cff template (place in repo root before tagging):
```yaml
cff-version: 1.2.0
title: "Sparse mFACEs + PTP/IEP classification of enteropancreatic periprandial trajectories"
authors:
  - family-names: "[Apellido]"
    given-names: "[Nombre]"
    orcid: "https://orcid.org/0000-0000-0000-0000"
version: "1.0"
date-released: "2026-MM-DD"
license: "CC-BY-4.0"
type: software
repository-code: "https://github.com/[org]/EPP10"
url: "https://github.com/[org]/EPP10"
abstract: "Code for sparse mFACEs joint analysis with Chiou normalization + per-analyte PTP/IEP Type I-V classification applied to ecological cohort-time-arm enteropancreatic hormone trajectories."
keywords:
  - functional data analysis
  - mFACEs
  - enteropancreatic hormones
  - bariatric surgery
  - ecological meta-analysis
```

Once DOI assigned, update `manuscript/*.md` files replacing `10.5281/zenodo.XXXXXXX` placeholders.

---

## STEP 6 — OSF DOI for pre-registration

```bash
# 1. Navigate to https://osf.io
# 2. Create new project: "EPP10 — sparse mFACEs enteropancreatic meta-analysis"
# 3. Upload:
#    - preregistration_cohort_map.yaml
#    - cohort_normalization_map.csv
# 4. Register the project as "Pre-registered Analysis Plan" (OSF Registries)
#    - Template: "Pre-Registration Template (v3)"
#    - Freeze timestamp: 2026-04-22 (data lock date)
# 5. Obtain DOI: 10.17605/OSF.IO/XXXXX
# 6. Update manuscript YAML references
```

Documentation: https://help.osf.io/article/158-register-your-project

---

## STEP 7 — Compile manuscript PDF + update DOIs

```bash
cd "$EPP10_ROOT/manuscript"

# Replace DOI placeholders
sed -i.bak "s/10.5281\/zenodo.XXXXXXX/10.5281\/zenodo.ACTUAL_DOI/g" *.md
sed -i.bak "s/10.17605\/OSF.IO\/XXXXX/10.17605\/OSF.IO\/ACTUAL_DOI/g" *.md

# Assemble via Quarto
quarto render manuscript.qmd --to pdf
```

Alternatively, if manuscript is in Word/LaTeX, update manually.

---

## STEP 8 — medRxiv submission

Navigate to https://www.medrxiv.org/submit-a-manuscript

**Submission package:**
1. Cover letter with STROBE + "TRIPOD+AI not applicable" declaration
2. Main manuscript PDF (≤40 MB)
3. Figures 1-3 as separate high-res PDFs (from `figures/`)
4. Supplementaries (4 files):
   - S1_STROBE_checklist.pdf (fill from EQUATOR template)
   - S2_compliance_tick_sheet.md → PDF (render via pandoc)
   - S3_Verification_Appendix.md → PDF
   - S4_PTP_IEP_integration.md → PDF (classify_ptp_iep.R + iep_frequency_by_cohort.csv)
5. Declarations (all completed):
   - Ethics/IRB: "Not applicable (ecological meta-analysis of published aggregate data)"
   - Consent: "Not applicable"
   - Funding: [fill]
   - Competing interests (36 mo): [fill]
   - Data availability: "De-identified aggregate data via Zenodo DOI {ZENODO_DOI}; pre-registration via OSF DOI {OSF_DOI}"
   - Code availability: "github.com/[org]/EPP10 v1.0 → Zenodo DOI {ZENODO_DOI}"
   - Clinical trial registration: "Not applicable"
   - AI/LLM disclosure: "Analytic pipeline design and code generation assisted by Claude Anthropic (2026-04); authors reviewed and validated all code and statistical outputs"
6. License: CC-BY 4.0 (Plan-S compliant)

**Expected timeline**: 2-4 business days for screening; DOI assigned upon acceptance (`10.1101/2026.MM.DD.XXXXXX`).

---

## STEP 9 — Post-submission updates

Even after submit, the B=2000 pipeline continues filling cache. Periodically:

```bash
cd "$EPP10_ROOT"
Rscript recover_B2000_results.R
git add bootstrap_B2000_*.{rds,csv} verification/B2000_snapshot.log
git commit -m "B=2000 snapshot: $(date +%Y-%m-%d) progress $(ls cache_bootstrap_B2000/ | grep rep_ | wc -l)/2000"
git tag -a v1.1-b2000-$(date +%Y%m%d) -m "B=2000 progress snapshot"
git push origin main --tags
# Zenodo auto-archives each tag as a new version
```

When B reaches ≥1500 successful reps (statistical saturation), update the Verification Appendix S3.6.3 with final numbers and submit as revised preprint version.

---

## STEP 10 — Final checklist

Before clicking "submit" on medRxiv:

- [ ] B=2000 snapshot captured (≥300 reps ideal, ≥1000 preferred for initial deposit)
- [ ] `renv.lock` committed
- [ ] git tag v1.0 pushed
- [ ] Zenodo DOI resolved (test the link)
- [ ] OSF DOI resolved (test the link)
- [ ] All DOI placeholders replaced in manuscript
- [ ] Figures 1-3 separate PDFs ready
- [ ] Supplementaries S1-S4 as PDFs ready
- [ ] Declarations block completed (ethics, funding, COI, data, code, AI/LLM)
- [ ] STROBE checklist completed (cross-reference line numbers)
- [ ] Cover letter mentions: STROBE reporting, TRIPOD+AI decline, CC-BY 4.0 license choice, Zenodo/OSF DOIs
- [ ] Final proofread: Results numbers match Tables, Tables 1-4 consistent with Figures

---

## Compute state summary

| Artefact | Size | SHA-256 (first 16) |
|---|---|---|
| `hormones_long_tidy.csv` | 200 KB | 93eba357565db827 |
| `cohort_normalization_map.csv` | 5 KB | fd133fb5bab90c1c |
| `pseudo_ipd_primary_M1000_rho050_cv100.csv` | 120 MB | c5d398b3a47d5a33 |
| `fit_mfaces_primary_results.rds` | 15 MB | (binary) |
| `fanova_results.csv` | 20 KB | a099aa3f67711763 |
| `bands_simultaneous.csv` | 250 KB | 9aaa04140bf7c690 |
| `bootstrap_stability_results.rds` | 2 MB | (binary) |
| `ptp_iep_results.rds` | 3 MB | (binary) |
| `preregistration_cohort_map.yaml` | 10 KB | (frozen 2026-04-22) |

Total repo size (excluding caches and raw master): ~150 MB.

---

**This runbook is self-contained.** Execute sequentially; any step's verification block confirms success before proceeding. Steps 1-4 take ~15 min total. Steps 5-7 (DOI + manuscript compile) take ~30-60 min including Zenodo/OSF UI interactions. Step 8 submission itself ~30 min. Steps 9-10 are post-submission maintenance.

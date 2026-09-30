# QC Dashboard of analyte data

## Quick Start



### Review QC results

1. Select a matrix.
2. Select a biomarker.
3. Review:
   - Imputation eligibility
   - QC status
   - Summary statistics
   - LOD/LOQ information
   - IMP and MEB integrity checks
   - Distribution plots
4. Investigate biomarkers with **WARNING** or **FAIL** status.

---

# Overview

The Biomonitoring QC Dashboard provides an interactive quality control review of biomonitoring datasets after imputation and derivation of secondary variables.

The dashboard is intended as a final validation step before downstream statistical analyses.

All calculations are performed during QC object generation. The dashboard itself only reads pre-computed QC results and visualizes them.

The application runs locally in the user's R session. Uploaded files are not transferred to external servers.

---

# Input Requirements

The dashboard requires a single input file:

```text
qc_objects.rds
```

This file must have been generated using the accompanying `build_qc_objects.R` script.

The uploaded RDS file must contain the following objects:

```text
biomarker_lookup
biomarker_overview
variable_overview
summary_statistics
matrix_summary
lod_loq_combinations
imp_integrity
meb_integrity
plot_data
distribution_data
qq_data
qq_reference
flags
```

If one or more required objects are missing, the dashboard will reject the file.

---

# Data Assumptions

The QC build script assumes:

- Processed datasets are available in `ds.rds`
- Codebooks are available in `cb.rds`
- Only datasets and codebooks with suffix `_lab` are included
- Biomarkers are represented by their raw variable (`BaseVarname`)
- Censored observations are coded as:

```text
-1
-2
-3
```

These values are treated as non-detect observations.

---

# Dashboard Layout

The dashboard contains three main sections:

## QC Summary

Provides detailed QC information for the selected biomarker.

## Plots

Provides visual assessment of imputation results and distributions.

## All Biomarkers

Provides a searchable overview of all biomarkers within the selected matrix.

---

# Matrix Selection

The left panel allows selection of a matrix.

For each matrix, the dashboard reports:

- Number of records
- Number of biomarkers
- Number of PASS biomarkers
- Number of WARNING biomarkers
- Number of FAIL biomarkers

---

# QC Status Categories

## PASS

No QC issues were detected.

## WARNING

A potential issue was detected.

Examples:

Warnings should be reviewed but do not necessarily invalidate the biomarker.
Currently, all issues are registered as FAIL


## FAIL

A FAIL indicates that one or more QC requirements were not met.

Examples include:

- An expected `_imp` variable is missing, empty or only partially populated.
- A `_meb` variable is missing, empty or only partially populated.
- A `_bin` variable is missing, empty or only partially populated.
- A measured observation has no corresponding non-missing value in `_imp` or `_meb`.
- A censored observation has not been imputed in `_imp` when imputation is expected.
- A censored observation has not been substituted in `_meb`.
- A raw biomarker value is present while both the corresponding LOD and LOQ values are missing.

Biomarkers with a FAIL status should be investigated and corrected before further analysis.

---

# Imputation Eligibility

The dashboard reports whether imputation is expected for a biomarker.

Imputation is expected when:

```text
Percent detected ≥ 30%
AND
Number of unique detected values ≥ 10
```

Both criteria must be satisfied.

When imputation is not expected:

- Absence of the `_imp` variable is not considered a QC failure
- IMP integrity checks are not evaluated

---

# Summary Statistics

The dashboard displays summary statistics for:

## Raw variable

Original biomarker values.

## IMP variable

Imputed biomarker values.

## MEB variable

Biomarker values after substitution of censored observations.

For each variable:

- N
- Number of missing observations
- Minimum
- 25th percentile
- Median
- Mean
- 75th percentile
- Maximum

---

# Detection Variable Summary

The dashboard also summarizes the corresponding:

```text
_bin
```

variable.

Reported statistics include:

- Number of 0 values
- Number of 1 values
- Percent detected

This provides a rapid check of detection frequencies.

---

# QC Summary Panel

The QC Summary panel provides information about:

## Censored observations

Number and percentage of observations coded as:

```text
-1
-2
-3
```

## Imputation variables

Expected imputation-related variables:

```text
_imp
_meb
_bin
```

and derived correction variables:

```text
_imp_sg
_meb_sg

_imp_crt
_meb_crt

_imp_lip
_meb_lip
```

For each variable the dashboard reports:

- Completion percentage
- QC status
- Presence in the dataset

---

# LOD / LOQ Review

The dashboard reports:

- Number of unique LOD/LOQ combinations
- Frequency of each combination
- Percentage of observations represented by each combination

This allows verification of analytical method consistency.

### Warning

A warning is generated when:

- A raw biomarker value is present
- Both corresponding LOD and LOQ values are missing

These cases should be investigated because imputation and substitution procedures depend on valid LOD/LOQ information.

---

# IMP Integrity Check

The IMP integrity panel evaluates completeness of imputed variables.

The dashboard compares:

```text
Raw biomarker
↓
IMP variable
```

For each observation it verifies that:

### Measured observations

All measured values have a non-missing value in `_imp`.

### Censored observations

All censored observations have a non-missing value in `_imp` when imputation is expected.

The dashboard reports:

- Number of measured values retained
- Number of censored values imputed
- Number of missing observations

### Important

The integrity check only verifies completeness.

It does **not** verify whether an imputed value is statistically correct.

It only verifies that expected values are present.

---

# MEB Integrity Check

The MEB integrity panel evaluates the substitution-based variable.

The dashboard compares:

```text
Raw biomarker
↓
MEB variable
```

For each observation it verifies:

### Measured observations

Measured observations have a non-missing value in `_meb`.

### Censored observations

Censored observations have a non-missing substituted value in `_meb`.

The dashboard reports:

- Number of measured values retained
- Number of censored values substituted
- Number of missing observations

As for IMP, this check evaluates completeness, not numerical correctness.

---

# Distribution Plots

## Histogram

Displays:

```text
Measured values
vs
Imputed values
```

on a log scale.

Purpose:

- Assess whether imputed values follow the expected range of observed values.
- Detect unusual imputation behaviour.

---

## Density Plot

Displays the complete distribution of the final `_imp` variable.

Purpose:

- Assess distribution shape.
- Identify potential distortions caused by imputation.

---

## Q-Q Plot

Displays:

```text
log(_imp)
```

against theoretical normal quantiles.

Purpose:

- Assess approximate log-normality.
- Detect extreme deviations from expected distributional assumptions.

---

# All Biomarkers Table

The All Biomarkers tab provides a searchable overview of all biomarkers within the selected matrix.

The table includes:

- Detection metrics
- Imputation eligibility
- Numbers of censored and imputed observations
- LOD/LOQ information
- IMP integrity status
- MEB integrity status
- BIN status
- Overall QC status

Biomarkers are automatically ordered:

```text
FAIL
WARNING
PASS
```

to prioritize review of problematic biomarkers.

---

# Recommended QC Workflow

1. Generate `qc_objects.rds`.
2. Open the dashboard.
3. Select a matrix.
4. Filter to WARNING and FAIL biomarkers.
5. Review:
   - Missing derived variables
   - LOD/LOQ warnings
   - IMP integrity
   - MEB integrity
6. Inspect visualizations for biomarkers requiring imputation.
7. Resolve issues in source data.
8. Rebuild QC objects.
9. Repeat until no unexpected QC failures remain.

---

# Notes

- The dashboard is intended for QC review, not statistical analysis.
- All calculations are pre-computed before dashboard loading.
- The dashboard does not modify the uploaded file.
- Closing the browser window or stopping the R session removes the loaded data from memory.
``
# JUMPt Version 1.2.0 (Python Edition)

---

### ✅ Latest Release
**JUMPt v1.2.0** — Deterministic dual–time–unit implementation (supports both **hours** and **days**) with complete Python-based reproducibility.  
*(Previous: JUMPt v1.1.1 — days only)*

> Looking for the previous release? **JUMPt Version 1.1.1** → **[Go to the older repository](https://github.com/abhijitju06/JUMPt-Version-1.1.1-)**.

---

## Table of Contents
- [Introduction](#introduction)
- [Release Notes (v1.2.0)](#release-notes-v120)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Input Data Preparation](#input-data-preparation)
- [Parameter File Configuration](#parameter-file-configuration)
- [Running JUMPt](#running-jumpt)
- [Output Files](#output-files)
- [Authors](#authors)
- [Maintainers](#maintainers)
- [Acknowledgments](#acknowledgments)
- [References](#references)

---

## Introduction

**JUMPt (Jumbo Mass-Spectrometry–based Proteomics Turnover)** is a Python framework for modeling **protein turnover kinetics** from pulse SILAC (Stable Isotope Labeling by Amino Acids in Cell Culture) experiments.  
It implements a **recycling-aware differential equation (ODE)** system that simultaneously fits unlabeled **free lysine (Lys)** and **protein-bound Lys** trajectories to compute *apparent* and *corrected* half-lives.

Version **1.2.0** extends JUMPt’s deterministic solver to handle datasets expressed in **hours (in vitro)** and **days (in vivo)** automatically, ensuring identical algorithmic behavior across both scales.

---

## Release Notes (v1.2.0)

### New Features

1. **Dual time-unit support**  
   Automatically detects and processes input time points in *hours* or *days*, using:
   - `jumpt_hours.py` — optimized solver for in vitro datasets  
   - `jumpt_days.py` — standard solver for in vivo datasets  

2. **Deterministic global optimization**  
   Fixed random seeds and bounded γ-space optimization ensure reproducible half-life estimation.

3. **Adaptive matrix-exponential solver**  
   Dynamic step control prevents numerical overflow or instability in proteins with extremely short (<0.5 h) or very long (>100 d) turnover rates.

4. **Display thresholds only**  
   The code enforces stability internally but reports display limits like `<0.5 h` or `>2400 h` purely for readability.

5. **Performance tuning**  
   Users can adjust runtime and precision with parameters such as:
   - `n_starts`  
   - `max_nfev_global`  
   - `max_nfev_refine`  
   - `expm_safe_magnitude`  

6. **Unified command-line interface**  
   A single main file (`jumpt_main.py`) automatically dispatches to either *hours* or *days* mode based on the input.

7. **Reproducibility**  
   Identical inputs always yield identical outputs across Python versions and operating systems.

---

## System Requirements

- **Python:** 3.9 or newer (tested up to 3.12)
- **Dependencies:**
  ```bash
  pip install numpy pandas scipy openpyxl

  ## Hardware

**Multi-core CPU and ≥16 GB RAM** are recommended for large proteome datasets.

---

## Installation

### Option A – Run in Place
```bash
git clone https://github.com/<your-org>/JUMPt_v1.2.0.git
cd JUMPt_v1.2.0
python jumpt_main.py --params JUMPt_main.params
```
## Input Data Preparation

Prepare an Excel file (`.xlsx`) with the following structure:

| Row | Description |
|-----|--------------|
| 1 | Pulse time points (in hours or days) |
| 2 | Total Lys (concentration) |
| 3 | Free Lys fraction |
| 4 → N | Protein labeling fractions |

**Example files**
- `7951_Rep1_K_Hours.xlsx` → *in vitro* dataset  
- `Small_Cell_Cere_Inputs.xlsx` → *in vivo* dataset  

The program automatically detects the time unit from the parameter file and numeric range.

---

## Parameter File Configuration

Edit `JUMPt_main.params` before running.

| Category | Parameter | Description |
|-----------|------------|-------------|
| Input | `input_file` | Full path to Excel dataset |
| Units | Auto-detected (*hours* or *days*) |
| Bounds | `display_min_hours`, `display_max_hours`, `display_min_days`, `display_max_days` | Display thresholds |
| Physics | `disable_physics_filter`, `physics_filter_mode` | Controls for enforcing physical constraints |
| Runtime | `n_starts`, `max_nfev_global`, `max_nfev_refine`, `use_refine`, `expm_safe_magnitude` | Runtime tuning |
| Reproducibility | `random_seed` | Ensures deterministic behavior |

Lines beginning with `;`, `#`, or `%` are treated as comments.

---

## Running JUMPt

Execute:

```bash
python jumpt_main.py --params JUMPt_main.params
```
## Console Output

The console automatically reports the processing mode:

```bash
Welcome to JUMPt Python (hours mode)...
```
or

```bash
Welcome to JUMPt Python (days mode)...
```
Batch progress messages show optimization status, for example:

```bash
Optimizing proteins 1–30 (of 247)
```
## Output Files

Two Excel files are generated in the same directory as the input file:

| File | Description |
|------|--------------|
| **results_Corrected_T50_<input>.xlsx** | Corrected half-lives (hours/days), confidence intervals, residuals, extremes, and parameter provenance |
| **results_Apparent_T50_<input>.xlsx** | Apparent half-lives from exponential fits |

---

### Display Limits

- `<0.5 h` for sub-hour turnover  
- `>2400 h` or `>100 d` for extremely stable proteins  

---

## Authors

**Dr. Abhijit Dasgupta** — SRM University AP, India  

---

## Maintainers

For bug reports or feature suggestions, please open a GitHub Issue or contact:  
📧 **abhijitju06@gmail.com**

---

## Acknowledgments

We acknowledge the support of **St. Jude Children’s Research Hospital**, **ALSAC**, and the **U.S. National Institutes of Health** for enabling the JUMP software suite.

---

## References

1. Chepyala et al., *Analytical Chemistry*, 93(40): 13495–13504 (2021)  
2. Wang X. et al., *Molecular & Cellular Proteomics*, 13(12): 3663–3673 (2014)  
3. Li W. et al., *Cell*, 188(8): 2267–2287 (2025)  
4. Yarbro J.M. et al., *Nature Communications*, 16(1): 1533 (2025)


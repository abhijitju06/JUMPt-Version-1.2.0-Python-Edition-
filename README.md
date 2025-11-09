# JUMPt Version 1.2.0 (Python Edition)

---

### ✅ Latest Release
**JUMPt v1.2.0** — Deterministic dual–time–unit implementation (supports both **hours** and **days**) with complete Python-based reproducibility.  
*(Previous: JUMPt v1.1.1 — days only)*

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


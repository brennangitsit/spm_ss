# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

SPM_SS (Subject-Specific Analysis Toolbox) is a MATLAB neuroimaging toolbox for ROI-based and voxel-based multi-subject statistical analyses using functional localizers. It accounts for inter-subject variability in brain activation loci by using supra-threshold voxels from a functional localizer contrast to define analysis units within each subject. Version 18.a, developed at the McGovern Institute, MIT.

**Dependencies:** MATLAB 6.5+, SPM5/SPM8/SPM12, CONN toolbox (optional, for QA plots)

## Running the Toolbox

```matlab
% GUI mode
spm_ss

% Batch mode (programmatic) - the standard 3-step workflow
ss = spm_ss_design(ss);    % Step 1: specify analysis
ss = spm_ss_estimate(ss);  % Step 2: estimate parameters
spm_ss_results(ss, Ic);    % Step 3: view results

% QA plots (requires CONN toolbox)
spm_ss qa <foldername>
```

There is no formal build system, test runner, or linter. Batch examples in `batchexamples/` serve as integration tests (e.g., `GcSS_test.m`).

## Architecture

The toolbox follows a 3-step pipeline: **Design -> Estimate -> Results**, controlled by an `ss` struct that flows through each step.

### Core Pipeline

- `spm_ss.m` - Main entry point and GUI launcher. Dispatches to sub-functions via `spm_ss('<command>')` syntax.
- `spm_ss_design.m` - Specifies the analysis: subject SPM.mat files, contrast selection, localizer thresholds, between-subjects model. All `ss` struct fields are documented in this file's header comments (lines 1-39).
- `spm_ss_estimate.m` - Routes to either `spm_ss_estimate_ROI.m` or `spm_ss_estimate_voxel.m` based on `ss.typen`.
- `spm_ss_results.m` - Interactive visualization, thresholding, and CSV report generation.
- `spm_ss_contrast.m` - Routes to `spm_ss_contrast_ROI.m` or `spm_ss_contrast_voxel.m`.

### Three Analysis Types

| `ss.type` | `ss.typen` | Description |
|-----------|-----------|-------------|
| `'voxel'` | 1 | Voxel-based analysis |
| `'GcSS'` | 2 | Group-constrained Subject-Specific (auto-defined ROIs via watershed parcellation) |
| `'mROI'` | 3 | Manually-defined ROI-based analysis |

### Key Supporting Functions

- `spm_ss_glm.m` - Core GLM estimation and hypothesis testing (T-stats, F-stats, Wilks' Lambda)
- `spm_ss_threshold.m` / `spm_ss_createlocalizermask.m` - Localizer mask creation with multiple correction methods (FDR, FWE, none, percentile, Nvoxels, automatic)
- `spm_ss_watershed.m` - Watershed algorithm for ROI parcellation in GcSS analyses
- `spm_ss_overlap.m` - Computes pair-wise subject overlap measures
- `spm_ss_fdr.m` - False discovery rate calculation
- `spm_ss_crossvalidate_sessions.m` - Partitions contrasts across sessions for cross-validation when localizer and effect-of-interest are non-orthogonal
- `spm_ss_roilabels.m` - Parses ROI labels from TXT, CSV, or XLS files (multiple format variants)
- `spm_ss_display.m` - Surface rendering of subject-specific activations
- `spm_ss_importspm.m` - SPM structure import utilities
- `spm_bcc_*.m` - Between-subjects Cross-validation Contrast analysis (separate pipeline)

### The `ss` Struct

The central data structure. Key fields (fully documented in `spm_ss_design.m` header):

- `ss.swd` - Output directory
- `ss.type` - Analysis type (`'voxel'`, `'GcSS'`, `'mROI'`)
- `ss.files_spm` - Cell array of subject SPM.mat paths
- `ss.EffectOfInterest_contrasts` - Effect contrast names
- `ss.Localizer_contrasts` - Localizer contrast names
- `ss.Localizer_thr_type` / `ss.Localizer_thr_p` - Thresholding method and value
- `ss.overlap_thr_vox` / `ss.overlap_thr_roi` - Overlap thresholds for GcSS
- `ss.homologous` - (GcSS only) 1/0 creates bilaterally symmetric parcels by sagittal flip+average of overlap map before watershed
- `ss.model` - Between-subjects model (1=one-sample t, 2=two-sample t, 3=regression)
- `ss.estimation` - `'OLS'` or `'ReML'`
- `ss.ask` - GUI interaction level (`'none'`, `'missing'`, `'all'`)

### Output Files

Estimation produces in `ss.swd`: `SPM_ss*.mat`, brain volumes (`*beta.img`, `*overlap.img`, `*T_####.img`, `*P_####.img`, `*Pfdr_####.img`), and CSV reports (`*data.csv`, `*results_####.Stats.csv`).

## Conventions

- Cross-validation is handled automatically when localizer and effect-of-interest contrasts are non-orthogonal (same session data).
- Functions use `spm_ss init silent` at the top to check for the evlab17 package (non-fatal if missing).
- GUI interaction uses SPM's `spm_input` and `spm_select` functions.
- Some functions use MATLAB persistent variables for state across interactive sessions.

# Selection Breakdown Beyond the Effectively-Neutral Zone

**Empirically Supported DFE Parameters and Deleterious Fixation Rates in Finite Populations**

Christopher Rupe — Independent Researcher, Cleveland, Ohio, USA
ORCID: [0009-0007-8497-4257](https://orcid.org/0009-0007-8497-4257)

This repository contains the SLiM simulation scripts, raw simulation output, and analysis/plotting scripts supporting the manuscript submitted to the *Journal of Mathematical Biology*.

---

## Overview

The study evaluates the distribution-of-fitness-effects (DFE) parameters used by Hancock & Cardinale (2024, "H&C") in their forward-time response to Basener & Sanford (2018). Using H&C's own SLiM script, one parameter at a time was changed to values drawn from the mainstream population-genetics literature. Under empirically supported DFE parameters — even the most conservative published estimates — fitness trajectories decline and populations move toward extinction, with 93–99% of deleterious fixations falling within the effectively-neutral zone (ENZ) and selection-breakdown zone (SBZ).

All simulations use a finite diploid population (K = 1,000) under a non-Wright–Fisher model with density-dependent regulation.

---

## Repository structure

| Folder | Contents |
|---|---|
| `H&C Original SLiM Script/` | `gen_entropy_test.slim` — the original script from Hancock & Cardinale (included unmodified for reference). |
| `Modified SLiM Script/` | `gen_entropy_test_intrinsic_fitness_and_extinction_output_added.slim` — the script used in this study (H&C's script plus intrinsic-fitness logging and extinction output; see below). |
| `H&C Published Output/` | `GE_run1.txt`–`GE_run5.txt` and `GE_run_full.txt` — H&C's own published simulation outputs, used as the baseline comparison. |
| `H&C Parameters Intrinsic Tracked/` | Reproductions of H&C's baseline (mean −0.5, shape 0.5, 1:100 beneficials) with intrinsic-fitness tracking. Source data for the baseline figures and the §3.1 baseline replicate (227 deleterious / 236 beneficial fixations). |
| `Beneficial Ratio 1 to 1000/` | Runs reducing the beneficial rate to 1:1,000, including the 50k/100k/200k/300k progression. |
| `Mean Fitness Isolated/` | "Batch 1" runs isolating the mean deleterious effect (shape held at 0.5). |
| `Shape Isolated/` | "Batch 2" runs isolating the DFE shape parameter (mean held at −0.5). |
| `Mean s Shape Combined/` | "Batch 3" runs varying shape and mean fitness effect together. |
| `Dominance/` | Runs varying the dominance coefficient (h). |
| `Combined Runs/` | Combined-parameter runs (e.g., shape + mean + beneficial rate together). |
| `R Scripts/` | Plotting and DFE-analysis scripts (see below). |

> **Note on file naming.** Each experiment folder contains run logs (`*.txt` with a `cycle, Population_size, …` header) and, where allele-frequency spectra were needed, full-output dumps (`*_Full_Output*.txt`). Run logs and full-output files that share a base name come from the same run. Files with identical names duplicated across folders are byte-identical copies of the same run.

---

## SLiM scripts

Both scripts target **SLiM 4.x** (non-Wright–Fisher model; `initializeSLiMModelType("nonWF")`).

The modified script differs from H&C's original only by **adding output**, not by changing the model:

1. **Intrinsic-fitness logging** — an extra log column,
   `mean(p1.cachedFitness(NULL)) * p1.individualCount / K`, which separates the
   genomic (intrinsic) fitness from the density-dependent scaling.
2. **Extinction handling** — the run prints a full population dump and stops if the
   population reaches zero, and emits a near-extinction full dump (for allele-frequency
   spectra) when the population drops below five individuals.

All DFE parameters (deleterious gamma mean and shape, beneficial exponential mean, the deleterious:beneficial fraction in the genomic element, the dominance coefficient, and the run duration) are set inline in the script and were edited per run as described in the manuscript. Deleterious (`m2`) and beneficial (`m3`) mutation types use `convertToSubstitution = F`, so fixed mutations remain listed in the output (see "Counting fixations" below).

### Running a simulation

```
slim gen_entropy_test_intrinsic_fitness_and_extinction_output_added.slim
```

Edit the parameter lines at the top of the script (and the terminal-generation block) to reproduce a specific run.

---

## R scripts

| Script | Purpose |
|---|---|
| `slim_plot_generator_v15 Thinner Line.R`, `slim_plot_generator_v14.R` | Fitness / mutation-count trajectory plots from run logs. |
| `slim_multipanel_generator_v2_Fixed_Axes_Length.R`, `slim_multipanel_independent_x_v1.R` | Multi-panel figure assembly (fixed vs. independent axes). |
| `slim_popsize_plotter_v2.R` | Population-size trajectory plots. |
| `DFE_fixed_deleterious_JoMB.R` | Allele-frequency spectrum of fixed deleterious mutations. |
| `DFE_comparison_ENZ.R`, `DFE_comparison_ENZ_SBZ.R` | ENZ / SBZ classification and comparison plots. |

---

## Counting fixations

Because `m2` and `m3` are set with `convertToSubstitution = F`, fixed mutations are **not**
moved to a separate substitutions list — they remain in the `Mutations:` block of a
`sim.outputFull()` dump with a derived-allele count equal to `2N` (twice the diploid
population size). A mutation is therefore counted as **fixed** when its count is ≥ 0.999 × 2N.
Run-log columns (`Deleterious_muts`, `Beneficial_muts`) report *segregating* mutations and
should not be confused with fixation counts.

---

## Provenance and attribution

The files in `H&C Original SLiM Script/` and `H&C Published Output/` originate from
**Hancock & Cardinale (2024)** and are included here, unmodified, solely for reference and
reproducibility. They remain the work of their original authors under their original terms;
no claim of authorship or relicensing is made over them. The modified script in this
repository is a derivative of H&C's original script. Users intending to redistribute the
H&C-origin files should consult the licensing terms of the original source.

---

## Citation

If you use this code or data, please cite the manuscript (full citation to be added on
publication) and the original SLiM software:

> Haller, B.C. & Messer, P.W. (2023). SLiM 4: Multispecies eco-evolutionary modeling.
> *The American Naturalist*, 201(5), E127–E139.

---

## License

This repository is dual-licensed so each part carries the license designed for it:

- **Code** — the modified SLiM script and the R analysis/plotting scripts — is released
  under the **MIT License** (see `LICENSE`).
- **Data** — the simulation output files — are released under the **Creative Commons
  Attribution 4.0 International License (CC BY 4.0)** (see `LICENSE-CC-BY-4.0.txt`).
  Reuse is free, including commercially, with attribution.

Neither license extends to the H&C-origin files described under *Provenance and
attribution* above, which remain the work of their original authors.

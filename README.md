# When well-balanced is not well-identified

Minimal reproducibility assets for an empirical methodological stress test of
crystalloid target trial emulation in critical care electronic health record
data. The repository supports inspection of the outcome rules, reuse of the
diagnostic calculations, and reproduction of the three main figures without
redistributing participant-level records.

## Repository contents

- `outcome_definitions.py`: pure functions for measured-creatinine baseline,
  SMART-style MAKE-30 components, new renal replacement therapy, persistent
  renal dysfunction, and creatinine-based KDIGO acute kidney injury.
- `tte_diagnostic_functions.R`: reusable overlap-weight, weighted-risk,
  observation-weight, missing-outcome-bound, and E-value helpers.
- `aggregate_figure_data.csv`: non-disclosive aggregate estimates and labels
  required by the main figures.
- `reproduce_main_figures.R`: regenerates Figures 1–3 from the aggregate file.
- `environment.yml`: compact Python/R environment specification.
- `LICENSE`: MIT terms for code and CC BY 4.0 terms for aggregate data.

## Data availability

The repository contains aggregate results only. Participant-level data and
participant-level derived data are not included. The source critical-care EHR
databases are restricted-access resources available to credentialed users from
PhysioNet after completion of the applicable training and data-use agreements.
They cannot be redistributed through this repository.

## Reproduce the main figures

Create the environment and run the figure script from the repository root:

```bash
conda env create -f environment.yml
conda activate crystalloid-tte-stress-test
Rscript reproduce_main_figures.R
```

Six files (PDF and PNG for Figures 1–3) will be written to
`reproduced_figures/`. No restricted data are required.

The Python outcome functions accept caller-supplied pandas data frames with
measurement timestamps and values. They do not connect to a database or read
local files.

## Interpretation caveats

- The observational estimates are an identification stress test, not estimates
  recommended for clinical treatment decisions.
- Eligibility based on delivered volume after time zero changes the target
  population.
- Forty-eight-hour dominance definitions incorporate post-baseline care and
  should not be interpreted automatically as per-protocol effects.
- Similarity under inverse probability of observation weighting does not remove
  uncertainty from outcome data that may be missing not at random.
- Randomized estimates are used as outcome-specific calibration signals; they
  are displayed rather than pooled with observational estimates.
- The cross-database feasibility gate failed, so no external treatment-effect
  replication is claimed.

## License

Source code is released under the MIT License. `aggregate_figure_data.csv` and
other non-code parameter or aggregate-result content are released under the
[Creative Commons Attribution 4.0 International License](https://creativecommons.org/licenses/by/4.0/).

## Citation

Please cite the accompanying article once its bibliographic record is
available. Until then, cite this repository by title, repository URL, version,
and access date.

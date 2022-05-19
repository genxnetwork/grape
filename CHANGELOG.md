# Changelog

## [v1.6] - 2022-05-20

### Added

- New workflow for simulation of a big relatives dataset (~500k samples) was added. It's available via `simbig` command of the pipeline launcher.
- Support multiple cores for the preprocessing (`preprocess`) workflow.
- IBD segments weighting feature was added, see `compute-weight-mask` workflow and `--weight-mask` parameter of the pipeline launcher.
- Several options for better control of the samples filtering were added: `--missing-samples`, `--alt-hom-samples`, `--het-samples`.
- Random seed parameter was added for the Ped-sim simulation.

### Changed

- GRAPE flows were renamed in the pipeline launcnher: `ibis_king` -> `ibis-king`, `germline -> germline-king`.
- `readme.md` and the GRAPE scheme were updated and actualized.
- Singularity support was removed in favour of conda environments.
- Code refactoring and clean up.

### Fixed

- Fixed `germline-king` simulation flow.
- Fixed `java command not found` during the `reference` workflow evaluation.

## [v1.5.2] - 2021-12-24

### Added

With `--flow ibis_king` grape now calculates IBD1 and IBD2 shares from KING data for the 0-3 degrees.

### Fixed

- Fixed a bug with parsing ERSA output for large datasets.
- Fixed a bug with setting every values for some rows in relatives.tsv to 2.
- Fixed `total_seg_len` and `total_seg_len_ibd2` calculation. Now `total_seg_len` corresponds to only ibd1 segments.

## [v1.5.1] - 2021-12-13

### Fixed

- Bundle downloading hotfix.
- File verification hotfix for reference workflow.

## [v1.5] - 2021-12-01

### Changed

- Removed singularity from all workflows.
- Many intermediate files are now temporary. This significantly reduces working folder size.

### Fixed

- Fixed removal of duplicate SNPs.
- ERSA-only workflow now correctly detects duplicates or monozygotic twins.

## Dockstore Release - 2021-09-28

### Added

- Dockstore support
- Small and full bundle reference downloading from azure

### Changed

- ERSA can handle 100k samples.
- Preprocessing saves phasing information in vcf input.
- MAF filter is now consistent across different inputs.

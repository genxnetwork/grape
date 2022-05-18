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

- Fix `germline-king` simulation flow.
- Fix `java command not found` during the `reference` workflow evaluation.

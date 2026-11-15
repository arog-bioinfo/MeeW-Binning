# Snakemake pipeline: MeeW-Binning

[![Snakemake](https://img.shields.io/badge/Snakemake-%E2%89%A59.12.0-039475?logo=snakemake&logoColor=white)](https://snakemake.readthedocs.io/)
[![run with conda](https://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/projects/conda/en/latest/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

<!-- TEMPORARY CATALOG OPT-OUT: later replace every standalone "pipeline" with the standard term formed by "work" + "flow", then remove this comment. -->

A Snakemake pipeline for metagenomic binning, domain routing, optional prokaryotic refinement, and coverage estimation.

## Usage

The pipeline starts from user-provided contigs and mapping inputs. It runs selected COMEBin, MetaBinner, SemiBin2, or MetaBAT2 binners, assigns source bins with Tiara, and can refine prokaryotic mappings with DAS Tool. CoverM coverage and abundance is mandatory for the configured prokaryotic and/or eukaryotic domains. Selected binners use CheckM internally for scoring, while final downstream genome QA and filtering belong to Annotation.

Install Snakemake `>=9.12.0` and Conda, then edit [`config/config.yaml`](config/config.yaml) and its sample sheet. Supply the mapping or read inputs required by the selected CoverM mode. COMEBin and MetaBinner also require a complete CheckM 1 database; configure automatic setup or a persistent pre-provisioned database as appropriate for the environment.

Results below `results/{sample}/` include binner outputs, Tiara domain assignments, optional DAS Tool refinement, final `prok_bins`, `euk_bins`, and `unknown_bins`, and CoverM tables. Unknown bins are retained with provenance but are excluded from DAS Tool, CoverM, and downstream annotation.

See the [configuration guide](config/README.md) for binner selection, database handling, refinement, input modes, domain policy, and exact output contracts. Run the commands below from the repository root.

## Deployment options

Conda is the supported software deployment method. After configuring `config/config.yaml`, inspect the planned jobs and then execute them with:

```bash
snakemake --cores 2 --dry-run
snakemake \
  --cores 2 \
  --software-deployment-method conda
```

Adjust cores and resource overrides for the selected binners and execution environment. This pipeline does not supply a container deployment.

## Execution profiles

The repository supplies [profile guidance](profiles/README.md), but no configured pipeline profile. Add or select a user-maintained profile for executor and site-specific resource settings.

## License

The pipeline is distributed under the [MIT License](LICENSE).

## References

No module-specific DOI is declared. Cite the repository version used and the selected binning and analysis tools: COMEBin, MetaBinner, SemiBin2, MetaBAT2, CheckM, Tiara, DAS Tool, CoverM, and Snakemake:

> Köster, J., Mölder, F., Jablonski, K. P., Letcher, B., Hall, M. B., Tomkins-Tinch, C. H., Sochat, V., Forster, J., Lee, S., Twardziok, S. O., Kanitz, A., Wilm, A., Holtgrewe, M., Rahmann, S., & Nahnsen, S. _Sustainable data analysis with Snakemake_. F1000Research, 10:33, **2021**. https://doi.org/10.12688/f1000research.29032.2.

## Workflow overview

This workflow starts from precomputed metagenome inputs and runs:

1. **Binning** with one selected tool in direct mode or at least two distinct selected tools in DAS Tool mode (`comebin`, `metabinner`, `semibin2`, `metabat2`)
2. **Source bin assignment with Tiara** using one assembly prediction for every selected binner
3. Optional **DAS Tool 1.1.7 refinement** of prokaryotic mappings, while preserving original eukaryotic and unknown FASTAs from one explicit source binner
4. Mandatory **CoverM coverage and relative abundance at genome and bin level** for configured prokaryotic and/or eukaryotic bins

MeeW-Binning handles metagenomic samples only. Genome quality assessment is performed by MeeW-Annotation when configured separately.

## Input data

The user-authored sample sheet (`sample_sheet` in `config.yaml`) has the following columns:

| sample | contigs                | bam            |
| ------ | ---------------------- | -------------- |
| s1     | path/to/s1.contigs.fna | path/to/s1.bam |

- `sample`, `contigs`, and `bam` are explicit paths required from the user.
- For ComeBin, `bam` follows the current wrapper/script behavior and may be a BAM directory or BAM path accepted by `run_comebin.sh`.
- For MetaBinner, `bam` may be one BAM file or a directory. Sorted direct child `*.bam` files are passed to JGI depth generation with `--noIntraDepthVariance`.
- For SemiBin2, `bam` may be one BAM file or a directory. A directory resolves to the same sorted direct child `*.bam` list used by MetaBAT2, and every resolved BAM is passed explicitly to `single_easy_bin`.
- For MetaBAT2, `bam` may be one BAM file or a directory. A directory resolves to sorted direct child `*.bam` files and fails with an actionable error when none are present. MetaBAT2 generates its own depth file before binning.
- The validated input contract supports ComeBin, MetaBinner, SemiBin2, and MetaBAT2 with `contigs` and `bam`.
- CoverM raw read modes require additional columns with values that are not empty:
  - `coverm.input_mode: coupled`: `reads_1` and `reads_2`
  - `coverm.input_mode: single`: `single`
  - `coverm.input_mode: interleaved`: `interleaved`
  - `coverm.input_mode: bam`: existing `bam` column; if it points to a directory, direct child `*.bam` files are used.
- The workflow does not perform assembly, mapping, or read QC.

`<results>/domain/annotation_samples.tsv` remains an internal per-bin operational inventory. Standalone Annotation continues to use its own user-authored input sheet.

## Key configuration fields

The example `config.yaml` is organized with block comments for the same major sections used below.

### General Settings

- `sample_sheet`: sample metadata TSV path.
- `resources_root`: optional module-wide persistent resource root. A nonblank
  value is used before ambient roots; blank or omitted preserves
  `MEEW_RESOURCES`, `$XDG_DATA_HOME/meew`, then `~/resources`.
- `checkm1.database_path`: optional **exact** CheckM 1 database directory. A
  nonblank value wins over every resource root. Otherwise the selected path is
  `<effective-resource-root>/checkm1_db`.
- `checkm1.download_database`: defaults to `true`, including when the key is
  omitted from a legacy config. When the selected exact directory is absent or
  empty, COMEBin or MetaBinner job startup invokes that environment's installed
  CheckM Bioconda post-link setup in protected staging. Set it to `false` for a
  strictly pre-provisioned database.
- Outputs owned by the workflow use Snakemake's `<results>` pathvar. Standalone runs resolve `<results>` to `results`.
- Configured paths and sample sheet path fields expand `~` and environment variables. Relative values remain relative to the Snakemake working directory.

### Binning

- `binning.tools`: list containing exactly one binner when `refinement.enabled` is `false`; when enabled, it must contain at least two distinct binners.
- `binning.metabat2.min_contig_length`: minimum contig length passed to MetaBAT2; it must be at least `1500`.
- `binning.metabat2.seed`: deterministic MetaBAT2 seed (default `1`). MetaBAT2's native `*.fa` bins and the membership provenance requested with `--saveCls` are retained under `<results>/binning/metabat2/bins`.
- `binning.metabinner.install_path`: optional source/tool root or Bioconda `bin` directory containing `run_metabinner.sh`, `scripts/`, and `auxiliary/`. A normal Conda prefix containing those entries under `bin/` is also accepted. The shell rule checks this field, then `METABINNER_HOME`, then derives the adjacent tool root from `command -v run_metabinner.sh`; a Bioconda environment on `PATH` therefore needs no explicit path.
- `binning.metabinner.min_contig_length`: strict contig filter threshold (default `1000`, so only lengths greater than 1000 are retained). Composition is fixed at tetranucleotide size four and has 136 features collapsed by reverse complement.
- `binning.metabinner.scale`: upstream MetaBinner dataset scale (`small`, `large`, or `huge`).
- `binning.semibin2.environment`: pretrained SemiBin2 habitat for one BAM (for example `global`). Set it to `null` for self supervised training. Multiple BAMs always force explicit self supervised mode because a pretrained model cannot consume their combined abundance features.
- `binning.semibin2.sequencing_type`: canonical SemiBin2 algorithm, either `short_read` or `long_read`. Hybrid assemblies use `long_read`.
- `binning.semibin2.extra`: unrelated advanced `single_easy_bin` options. Arguments managed by the wrapper for input, coverage input, output, threads, environment or model mode, sequencing type, and compression are rejected here with an actionable error.
- SemiBin2 clears its output root before execution, forces uncompressed FASTA with `--compression none`, and normalizes final bins to `<results>/binning/semibin2/output_bins` and metadata to `<results>/binning/semibin2/bins_info.tsv` for downstream domain separation.

MetaBinner preprocessing uses restartable shell rules in the MetaBinner 1.4.4/Python 3.7.6 environment. Bundled JGI generates the raw table with `--noIntraDepthVariance`; installed `Filter_tooshort.py` applies the strict length threshold; a small inline `awk` command preserves the JGI header and selects contig name plus abundance columns; installed `gen_kmer.py` creates the fixed k=4 composition table. The raw JGI table, filtered FASTA, coverage table, and tetranucleotide table are written under `<results>/binning/metabinner/profiles`.

Before COMEBin augmentation/training or MetaBinner component binning/scoring and
output cleanup, the Binning-owned runtime locks the selected cache and concrete
Conda environment. It validates the installed setup/scoring executables
(including `hmmfetch` and `hmmalign`), checks each tool's own bundled bacterial
and archaeal markers (MetaBinner uses `<install_root>/scripts/checkm_ms`), and
performs bounded format, manifest size/hash, and native HMM-index probes for
`.dmanifest`, `hmms/checkm.hmm`, `hmms/checkm.hmm.ssi`, and
`pfam/Pfam-A.hmm.dat`. A complete selected cache is reused without a downloader.
A nonempty incomplete or corrupt cache is never merged, removed, or replaced.

When automatic setup is needed, the wrapper reuses the **installed CheckM
Bioconda post-link hook**, including its pinned official archive SHA-256. The
legacy hook is confined to unique staging through a space-free alias, its
temporary package registration is confined to a private CheckM package copy,
the complete staged database is validated and atomically published, and native
CheckM registration plus a fresh `DefaultValues`/marker-parser subprocess must
resolve the final exact directory. Shared `DATA_CONFIG` hard links are detached
before final registration so another Conda prefix/package cache is not mutated.
The environment lock is stored inside the canonical environment and remains
held until COMEBin or MetaBinner completes, so jobs sharing one environment are
intentionally serialized even when they select different database paths. This
design gives one shared lock inode to path aliases on a shared filesystem; it
does not claim that local tests prove distributed lock semantics.

Bioconda also runs this hook automatically while it installs CheckM. Setting
`CHECKM_DATA_DIR` during creation of a new environment selects that install-time
destination, while `CHECKM_DATA_PATH` means a local archive to that hook. Merely
setting either variable after an environment has already been cached does not
rerun the installation hook. MeeW's job-time wrapper is the separate mechanism
that prepares the configured persistent destination before COMEBin or
MetaBinner. It deliberately removes `CHECKM_DATA_PATH` so an ambient local
archive cannot redirect configured runtime setup.

The upstream execution root is cleaned only after CheckM preflight and before
each direct `run_metabinner.sh` invocation. The final mapping is
`<results>/binning/metabinner/upstream/metabinner_res/metabinner_result.tsv`;
installed `gen_bins_from_tsv.py` materializes its FASTA bins under
`<results>/binning/metabinner/bins`. Native CheckM HMM commands and MetaBinner's
upstream shell scripts do not safely quote every whitespace-containing path, so
the runtime rejects a whitespace-containing CheckM database path and, for
MetaBinner, input, output, and installation paths before expensive processing.

### Refinement

- The `refinement` block is required and must explicitly set `enabled` and `method`.
- `refinement.enabled: false` requires `refinement.method: none` and exactly one selected binner. This uses direct Tiara domain separation.
- `refinement.enabled: true` requires `refinement.method: das_tool`, at least two distinct selected binners, and `refinement.euk_source_tool` contained in `binning.tools`. Despite its established name, this source-tool field owns preservation of both eukaryotic and unknown original FASTAs.
- `refinement.das_tool.score_threshold` is a number from `0` through `1`; `refinement.das_tool.threads` is a positive integer.
- DAS Tool is pinned to Bioconda `das_tool=1.1.7`, always uses DIAMOND, and receives only input, label, assembly, output, threshold, and thread options managed by the workflow.

Tiara runs once for each sample. One preparation job for the sample processes selected binners in configured order and each binner's direct child FASTAs in sorted order with the same majority base pair and fallback semantics. It writes only metadata: deterministic normalized bin identities use `<sample>__<tool>__bin_<NNNNNN>`, and every contig in a whole bin assigned prokaryotic is included in that binner's mapping. Preparation rejects contigs absent from the assembly, duplicate FASTA identifiers, and contigs assigned to more than one bin within a binner.

Only prokaryotic mappings enter DAS Tool. Final prokaryotic FASTA bins come from DAS Tool; final eukaryotic and unknown FASTA bins are direct copies that preserve every byte from `refinement.euk_source_tool` original FASTAs selected through consolidated assignment metadata. Assignments from every binner remain in the refinement provenance table, but unknown bins from other binners are not materialized into the final domain directory. The final outputs include `domain/prok_bins`, `domain/euk_bins`, and `domain/unknown_bins`; the assignment table and rename map include all three domains, while `domain/annotation_samples.tsv` includes only prokaryotic and eukaryotic bins.

Refinement provenance is retained under `<results>/refinement/`: the consolidated `bin_domain_assignments.tsv`, one headerless `das_tool/inputs/{tool}_prok_bins.tsv` for each selected binner, DAS Tool input provenance, `all_bins.eval.tsv`, `DASTool_summary.tsv`, selected contig mapping, run status, tool log, benchmark, and selected FASTA bins. No prokaryotic inputs and an outcome where DAS Tool selects no bins are valid deterministic outcomes. Declared directories and mappings remain empty, score and summary files contain only headers where upstream did not produce them, and no synthetic bin is created.

The preparation job creates no staged FASTA artifacts.

### Domain Separation

- `domain_separation.fallback_domain`: domain used when prokaryotic and eukaryotic evidence is tied or absent, and for unknown Tiara labels. It defaults to `unknown`; explicit `prok` and `euk` remain accepted for compatibility.

### Coverage and Abundance

- CoverM is always active and its outputs for every configured domain are final targets; the removed `coverm.enabled` key is rejected.
- `coverm.domains`: non-empty unique selection containing `prok`, `euk`, or both; defaults to `[prok, euk]`. `unknown` is never accepted or processed by CoverM.
- `coverm.input_mode`: one of `bam`, `coupled`, `single`, or `interleaved`.
- `coverm.mapper`: raw read mapper choice; one of `minimap2-sr`, `minimap2-ont`, `minimap2-pb`, `minimap2-hifi`, `minimap2-no-preset`, or `bwa-mem2`. This setting applies only to raw read modes, not BAM mode.
- `coverm.methods`: CoverM genome methods. The default recommended list is `relative_abundance`, `mean`, `covered_fraction`, `length`, and `count`.
- `coverm.extra`: additional arguments passed to `coverm genome`.

## CoverM coverage and abundance

Coverage targets are produced after `separate_bins_by_domain` for only the configured prokaryotic and/or eukaryotic domains:

- `<results>/coverm/prok_genome_coverage.tsv`
- `<results>/coverm/euk_genome_coverage.tsv`

Raw reads are recommended for relative abundance estimation because CoverM can map reads consistently with the selected mapper. BAM mode is provided as a fallback/reuse mode for existing mapped reads and does not use `coverm.mapper`.

CoverM is invoked only when the domain directory contains a direct-child regular file matching the configured FASTA extension. Otherwise, the target is a deterministic zero-byte TSV with no genome rows and the log records the skip reason. Native output and failure behavior are unchanged when matching bins exist.

## Tiara separation policy

Bin assignment is deterministic and based on **majority base pairs** per bin:

- prok if bp labeled prok is a strict majority
- euk if bp labeled euk is a strict majority
- otherwise tied or without prokaryotic/eukaryotic evidence and assigned to `domain_separation.fallback_domain` (`unknown` by default; `prok` and `euk` remain compatible)

Unknown Tiara labels use the same fallback, while contigs absent from the Tiara table do not vote. A positive prokaryotic/eukaryotic imbalance remains assigned to that positive domain. In direct mode, unknown bins are deterministically named `<sample>_unknown_bin_<NNN>.fa` and retained under `<results>/domain/unknown_bins` with summary and rename provenance.

In DAS Tool mode this policy is applied before refinement to every selected binner. Unknown assignments remain in `<results>/refinement/bin_domain_assignments.tsv` but are excluded from every DAS Tool mapping. DAS Tool output is routed to the final prokaryotic contract, and only `refinement.euk_source_tool` contributes final eukaryotic and unknown bins. Unknown bins are excluded from DAS Tool, CoverM, `domain/annotation_samples.tsv`, and Annotation, so the workflow performs no unknown-bin QA, taxonomy, coverage, or annotation.

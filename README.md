# CP113 HDOCK screening code

A compact code release for the historical KPNA4/MC1R peptide screening analysis.
Includes docking command preparation, score extraction, score distributions and
ranked score tables. This is a selected code subset, not an end-to-end reproduction
of the manuscript, structure generation, clustering, or subsequent simulations.

## Quick start

Tested with Python 3.12 and the pinned dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
python scripts/analyze_scores.py examples/scores_sample.csv --output results/example --top 10
```

For the full analysis, replace the sample path with your local
`merged_scores_sp_id.csv`. The preserved table has 152,867 peptide rows, with
one score column for each target; this is not 152,867 independent experiments.
Outputs: score summary, source counts, a distribution figure, and top-ranked
tables for KPNA4, MC1R and the sum of both scores. Input rows are not altered.
The two-target sum is a descriptive ranking, not calibrated affinity.

## Prepare docking jobs (optional)

Install HDOCK and its `createpl` utility separately and place both on `PATH`.
Provide prepared receptor PDBs and ligand PDBs under `ligands/part_001/`, etc.
Run from the repository root, using paths without spaces or shell metacharacters:

```bash
python scripts/generate_dock_commands.py -r inputs/KPNA4.pdb -l ligands -o docking/KPNA4 -cf kpna4_jobs.sh -df kpna4_done.log
bash kpna4_jobs.sh
python scripts/extract_scores.py docking/KPNA4 --receptor KPNA4 --output results/kpna4_raw.csv
```

Repeat with MC1R and separate output, command and done files. `createpl` exports
the top pose (`-nmax 1 -complex`). Ligand filenames must be unique within a run:
the historical generator tracks completion by filename only. Inspect the generated
commands before execution. GNU Parallel may run the lines with a user-chosen job
count. The release validation does not execute docking.

Extraction retains relative source paths and the numeric text in PDB remarks.
Raw extraction does not reconstruct the historical merged ID assignments; do not
join different libraries on ligand basename alone.

## Data and interpretation

`examples/scores_sample.csv` contains up to 10 existing rows per library source,
in original order and with original numeric text. It is a demonstration subset.
The full table, receptor/ligand structures, model weights and raw docking archives
are not included. See [data schema and limits](docs/DATA.md),
[script provenance](docs/PROVENANCE.md), and [validation](docs/VALIDATION.md).

HDOCK scores are docking scores, not physical binding free energies. `Ligand`
contains identifiers, not amino-acid sequences. The original analysis module is
preserved for traceability; use the supplied entry point, which only calls its
score-related functions. Its legacy sequence/heatmap routines are not validated
for this schema and are not part of the release workflow.

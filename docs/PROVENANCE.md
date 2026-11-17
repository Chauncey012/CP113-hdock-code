# Script provenance

Paths below are relative to the source CYC113 project, not this standalone package.
Hashes and sample selection details are recorded in `../source_manifest.json`.

| Release file | Origin and changes |
| --- | --- |
| scripts/generate_dock_commands.py | Exact copy of `05_provenance/HDOCK/generate_dock_commands.py` |
| scripts/score_analysis_original.py | Exact copy of `03_results/screening/scripts.py`; its standalone main is not the supported entry point |
| scripts/analyze_scores.py | New portable wrapper calling the original score distribution, ranking and CSV export functions; validates scores and IDs, adds summary/source counts |
| scripts/extract_scores.py | Small replacement informed by `05_provenance/HDOCK/extract_hdock.py`; recursively reads receptor-prefixed filenames, keeps original numeric text and source paths, fails on missing scores |

The historical extractor assumed ligand-first filenames, while the command
generator emits receptor-first filenames. The release extractor explicitly accepts
`--receptor` to resolve this naming mismatch. It reads only first-model header
remarks and leaves absent RMSD values blank.

This package excludes internal scheduling, credentials, host-specific deployment,
manuscript drafts and large scientific outputs. No third-party HDOCK executable or
model weights are redistributed. No new software license is assigned on behalf of
the authors; the authors can add their chosen license before broader redistribution.

# Input data

The merged CSV schema is:

| Column | Meaning |
| --- | --- |
| Ligand_ID | Identifier in the full screening-table namespace |
| Ligand | Original ligand identifier, not a sequence |
| Score_KPNA4 | KPNA4 HDOCK score |
| Score_MC1R | MC1R HDOCK score |
| Root | Library/run source |
| Part | Input partition |

Historical library counts: DB9 11,979; DRAMP 12,788; UniProt 128,097;
cosmetics 3. Total: 152,867 rows, each with two target scores.
The table's saved decimal text is retained, without claiming additional precision
beyond the source file. The full source SHA-256 is in `source_manifest.json`.

Selected-candidate tables use a different SP identifier namespace. The project
audit also records two unresolved selected-versus-screening score conflicts
(SP_0024833 and SP_0034834 in the selected namespace). Do not merge those tables
by identical SP identifiers or silently correct their scores. This release uses
the full-screening schema and does not reconstruct historical candidate selection.

The complete dataset is not bundled or deposited by this code release. A separate
data deposit and its verified DOI would be needed to support a public-data claim.

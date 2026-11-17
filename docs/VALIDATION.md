# Validation

Executed locally on 2026-11-16 (America/Los_Angeles), using Python 3.12 and the
versions in `requirements.txt`, in a separate temporary virtual environment.

- Ran the documented analysis entry point on the 33-row example (`--top 10`).
- Ran it on all 152,867 rows of the historical table (`--top 100`).
- Both runs exited successfully and produced the six documented outputs.
- Independently sorted the full table with Python's CSV reader and confirmed
  every identifier in all three top-100 tables, including original-order tie handling.
- Generated a docking command from dummy receptor/ligand paths; checked its
  receptor-prefixed filename and top-pose arguments, and ran `bash -n`.
- Extracted a synthetic first-model score and RMSD; confirmed receptor-prefix
  parsing and exact scientific-notation text preservation.
- Confirmed a PDB with no score causes extraction to fail without creating its CSV.
- Compared the two preserved scripts byte-for-byte against their source files.

No HDOCK calculation, structure prediction, historical selection reconstruction,
wet-lab validation or data-repository deposit was performed. Analysis outputs and
temporary environments remain outside the release archive. The copied historical
module's unused functions are not covered by these checks. Matplotlib is pinned
to 3.10.7 to retain the historical boxplot `labels` argument.

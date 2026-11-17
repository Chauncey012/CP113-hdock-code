#!/usr/bin/env python3
"""Extract first-model REMARK scores from receptor-prefixed HDOCK PDB files."""
import argparse
import csv
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--receptor", required=True, help="Receptor PDB basename, e.g. KPNA4")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.input.is_dir():
        parser.error("Input must be a directory")
    number = r"([-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?)"
    rows, missing = [], []
    prefix = args.receptor + "_"
    for path in sorted(args.input.rglob("*.pdb")):
        if not path.stem.startswith(prefix):
            continue
        values = {}
        with path.open() as handle:
            for line in handle:
                if line.startswith(("ATOM", "HETATM", "ENDMDL")):
                    break
                for field in ("Score", "RMSD"):
                    match = re.search(r"REMARK\s+" + field + r":\s*" + number, line)
                    if match:
                        values.setdefault(field, match.group(1))
        if "Score" not in values:
            missing.append(str(path.relative_to(args.input)))
            continue
        rows.append({"Ligand": path.stem[len(prefix):], "Receptor": args.receptor,
                     "hdock_score": values["Score"], "hdock_RMSD": values.get("RMSD", ""),
                     "source_file": str(path.relative_to(args.input))})
    if missing or not rows:
        parser.error(f"Found {len(rows)} scored files; missing scores in {len(missing)} files: {missing[:5]}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"Extracted {len(rows)} rows to {args.output}")


if __name__ == "__main__":
    main()

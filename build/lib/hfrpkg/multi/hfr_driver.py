#!/usr/bin/env python3

import os
from hfrpkg.runner import run_reaction
import argparse

reaction_map = {
    "1": "isogyric",
    "2": "isodesmic",
    "3": "hypohomodesmotic",
    "4": "homodesmotic"
}

def main(input_file, method, basis):
    with open(input_file, "r") as f:
        for i, line in enumerate(f, start=1):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) != 2:
                print(f"[SKIP] Line {i} malformed: {line}")
                continue

            level, smiles = parts
            reaction_type = reaction_map.get(level)
            if not reaction_type:
                print(f"[SKIP] Invalid level {level} on line {i}")
                continue

            outfolder = f"{i}.mhfr"
            #print(f"[RUN] Line {i}: level={level}, type={reaction_type}, SMILES={smiles}, folder={outfolder}")

            try:
                run_reaction(
                    action_type="write",
                    reaction_type=reaction_type,
                    input_smiles=smiles,
                    lhs=None,
                    rhs=None,
                    substruct=None,
                    replacement=None,
                    outfolder=outfolder,
                    method = method,
                    basis = basis
                )
                #print(f"[DONE] Folder written: {outfolder}")
            except Exception as e:
                print(f"[ERROR] Line {i} failed: {e}")
def main_cli():
    parser = argparse.ArgumentParser(description="Run multiple HFR jobs from a .txt list")
    parser.add_argument("input_file", help="Path to input .txt file")
    parser.add_argument("--m", help="METHOD")
    parser.add_argument("--b", help="BASIS")
    args = parser.parse_args()
    main(args.input_file, method=args.m, basis=args.b)

if __name__ == "__main__":
    main_cli()

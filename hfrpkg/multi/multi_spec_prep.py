import sys
import os
import glob
import argparse
import shutil
from hfrpkg.multi.compute_multiple import combined_postprocess as compute_multiple
from hfrpkg.multi.unique_spec import make_spec

def multi_spprep(method, basis, extension):
#  compute_multiple()
  
  make_spec(method, basis, extension)


  





def main_cli():
    parser = argparse.ArgumentParser(description="Generate and submit SP .com files from optimized .log files")
    parser.add_argument("--m", "--method", dest="method", default="B3LYP", help="DFT method (default: B3LYP)")
    parser.add_argument("--b", "--basis", dest="basis", default="6-31G", help="Basis set (default: 6-31G)")
    parser.add_argument(
        "--s", choices=["g", "o", "p"], default="g",
        help="Software to use: 'g' for Gaussian, 'o' for Orca, 'p' for Psi4 (default: 'g')"
    )
    args = parser.parse_args()

    multi_spprep(method=args.method, basis=args.basis, extension = args.s)

if __name__ == "__main__":
    main_cli()

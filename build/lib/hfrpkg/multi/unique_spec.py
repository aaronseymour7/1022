import sys
import os
import glob
import argparse
import shutil
from AaronTools.input import Theory, FileWriter
from AaronTools.geometry import Geometry
from AaronTools.job_control import SubmitProcess
from AaronTools.theory.job_types import SinglePointJob
from hfrpkg.run_unique import run_jobs
from hfrpkg.utils import get_extensions


def make_spec(method, basis, extension):
    extension_map = {
        "g": ".com",
        "o": ".inp",
        "p": ".in"
    }
    out_map = {
        "g": ".log",
        "o": ".out",
        "p": ".dat"
    }
    let_map = {
        "g": "Gaussian",
        "o": "ORCA",
        "p": "Psi4"
    }
    ext = extension_map.get(extension.lower())
    out = out_map.get(extension.lower())
    
    if ext is None or out is None:
        print(f"Unknown software code '{extension}'. Please use 'g', 'o', or 'p'.")
        sys.exit(1)
    os.makedirs("unique_files/spec", exist_ok=True)
    folder = "unique_files/spec/"
    index_path = os.path.join(folder, "index.txt")
    shutil.copyfile("unique_files/index.txt", index_path)
    optin, optout = get_extensions("1.mhfr/index.txt")
    log_files = glob.glob(os.path.join(folder, "*"+ optout))

    if not log_files:
        print(f"No {optout} files found in folder.")
        
        sys.exit(1)

    level = Theory(
        method=method,
        basis=basis,
        job_type=SinglePointJob()
    )
    
    
    for log_path in log_files:
        name = os.path.splitext(os.path.basename(log_path))[0]
        com_path = os.path.join(folder, name + ext)

        try:
            geom = Geometry(log_path)
            geom.write(outfile=com_path, theory=level)
        except Exception as e:
            print(f"Error processing {log_path}: {e}")
    
    with open("unique_files/index.txt", "r") as fin, open(index_path, "w") as fout:
        for line in fin:
            if line.startswith("Filename") or line.startswith("Level"):
                fout.write(line)
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                parts[3] = let_map.get(extension.lower())  
            fout.write("\t".join(parts) + "\n")
    cwd = os.getcwd()
    os.chdir(folder)
    try:
        run_jobs()
    finally:
        os.chdir(cwd)
        



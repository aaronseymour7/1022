from rdkit import Chem
from rdkit.Chem import AllChem
import os
import sys
import shutil
import glob
from AaronTools.job_control import SubmitProcess
from AaronTools.theory import Theory, OptimizationJob, FrequencyJob
from AaronTools.fileIO import FileReader

def ext_from_soft(header_line: str):

    parts = header_line.strip().split("\t")
    try:
        software = parts[1].lower()
        method = parts[3].lower()
        basis = parts[5].lower()
    except IndexError:
        raise ValueError("Header line does not have enough columns for software/method/basis")

    if software == "gaussian":
        return ".log", method, basis
    elif software == "orca":
        return ".dat", method, basis
    elif software == "psi4":
        return ".out", method, basis
    else:
        raise ValueError(f"Unknown software: {software}")

def get_enthalpy(logfile):
    try:
        reader = FileReader(logfile, just_geom=False)
        if 'E_ZPVE' in reader.keys():
            return reader['E_ZPVE']
        else:
            raise ValueError(f"Not finished job {logfile}")
    except Exception:
        c = os.getcwd()
        raise ValueError(f"Error in{c} ,{logfile}")
        
        
def compConf():
    cwd = os.getcwd()
    unique_dir = os.path.join(cwd, "unique_files")

    if not os.path.isdir(unique_dir):
        raise FileNotFoundError("unique_files is not found")

    os.chdir(unique_dir)

    try:
        with open("index.txt", "r") as f:
            lines = [line.strip() for line in f if line.strip()]
    except Exception:
        raise FileNotFoundError("index.txt not found in unique_files")

    extension, method, basis = ext_from_soft(lines[0])
    conf_dirs = glob.glob("*.conf")
    for cd in conf_dirs:
        E= 100000000
        dest = cd.split(".")[0]
        name = "x"
        os.chdir(cd)
        outs = glob.glob(f"*{extension}")
        for file in outs:
            enthalpy = get_enthalpy(file)
            if enthalpy<E:
                E= enthalpy
                name = file
        if name == "x":
            raise ValueError(f"No valid enthalpy found in {cd}")
        shutil.copy(name, os.path.join(unique_dir, f"{dest}{extension}"))
        os.chdir(unique_dir)
    print("Lowest energy conformers successfully extracted.")


def main_cli():
    compConf()

if __name__ == "__main__":
    compConf()


from rdkit import Chem
from rdkit.Chem import AllChem
import os
import sys
import shutil
import glob
from AaronTools.job_control import SubmitProcess
from AaronTools.theory import Theory, OptimizationJob, FrequencyJob
from AaronTools.fileIO import FileReader


def checkFile(logfile,cd):
    try:
        reader = FileReader(logfile, just_geom=False)
        if 'E_ZPVE' in reader.keys():
            pass
        else:
            print(f"[NOT FINISHED]{cd} {logfile}")
    except Exception as e:
        print(f"[ERROR]{cd} {logfile} ({e})")
def main():

    cwd = os.getcwd()
    unique_dir = os.path.join(cwd, "unique_files")

    if not os.path.isdir(unique_dir):
        raise FileNotFoundError("unique_files is not found")

    os.chdir(unique_dir)
    conf_dirs = glob.glob("*.conf")
    for cd in conf_dirs:
        os.chdir(cd)
        outs = glob.glob(f"*.log")
        for file in outs:
            checkFile(file,cd)
        os.chdir(unique_dir)

def main_cli():
    main()

if __name__ == "__main__":
    main()

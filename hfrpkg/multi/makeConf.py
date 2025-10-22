from rdkit import Chem
from rdkit.Chem import AllChem
import os
import glob
from AaronTools.job_control import SubmitProcess
from AaronTools.theory import Theory, OptimizationJob, FrequencyJob
from AaronTools.atoms import Atom
from rdkit.Chem import rdMolAlign, rdMolDescriptors
from AaronTools.geometry import Geometry
import argparse

def geom_from_rdkit(rdkitmol, confId=0):

    atom_list = []
    conf = rdkitmol.GetConformer(confId)
    for i, atom in enumerate(rdkitmol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        atom_list.append(Atom(element=atom.GetSymbol(),
                              coords=(pos.x, pos.y, pos.z)))
    return Geometry(atom_list)


def prune_conformers_by_rmsd_keep_order(mol, candidate_conf_ids, rmsd_cutoff=0.15):

    kept = []
    for cid in candidate_conf_ids:
        keep = True
        for kcid in kept:
            rmsd = rdMolAlign.GetBestRMS(mol, mol, prbId=cid, refId=kcid)
            if rmsd < rmsd_cutoff:
                keep = False
                break
        if keep:
            kept.append(cid)
    return kept


def generate_conformers_from_inchi(inchi, outfolder, ext, rmsd_cutoff=0.15):

    mol = Chem.MolFromInchi(inchi)
    if mol is None:
        raise ValueError(f"Could not parse InChI: {inchi!r}")
    mol = Chem.AddHs(mol)

    level = Theory(
        method=ext[1],
        basis=ext[2],
        job_type=[OptimizationJob(), FrequencyJob()]
    )

    n_rot = rdMolDescriptors.CalcNumRotatableBonds(mol)
    n_confs = max(5, n_rot * 5)

    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_confs, randomSeed=42))
    if not conf_ids:
        raise RuntimeError("Embedding produced no conformers")

    conf_ids = prune_conformers_by_rmsd_keep_order(mol, conf_ids, rmsd_cutoff=rmsd_cutoff)


    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant="MMFF94s")
    energies = []
    for cid in conf_ids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
        if ff is not None:    
            ff.Minimize()
            energies.append((cid, ff.CalcEnergy()))

    if energies:
        energies.sort(key=lambda x: x[1])
        ordered_conf_ids = [cid for cid, e in energies]
        energy_dict = {cid: e for cid, e in energies}
    else:
        ordered_conf_ids = conf_ids


    final_conf_ids = prune_conformers_by_rmsd_keep_order(mol, ordered_conf_ids, rmsd_cutoff=rmsd_cutoff)

    os.makedirs(outfolder, exist_ok=True)

    for i, conf_id in enumerate(final_conf_ids, start=1):
        geom = geom_from_rdkit(mol, confId=conf_id)
        outfile = os.path.join(outfolder, f"{i}{ext[0]}")
        geom.write(outfile=outfile, theory=level)
        #print(f"Wrote {outfile} (conf {conf_id})")


def run_jobs():
    com_files = glob.glob("*.com") + glob.glob("*.inp") + glob.glob("*.in")
    if not com_files:
        print("No runnable files found in current directory.")
        return

    for f in com_files:
        submit_process = SubmitProcess(f, 24, 12, 24)
        try:
            submit_process.submit()
            #print(f"Submitted {f}")
        except Exception as e:
            print(f"Failed to submit {f}: {e}")
def count():
    com_files = glob.glob("*.com") + glob.glob("*.inp") + glob.glob("*.in")
    return len(com_files)

def ext_from_soft(header_line: str):

    parts = header_line.strip().split("\t")
    try:
        software = parts[1].strip().lower()
        method = parts[3].strip()
        basis = parts[5].strip()
    except IndexError:
        raise ValueError("Header line does not have enough columns for software/method/basis")

    if not method:
        method = "B3LYP"
    if not basis:
        basis = "6-31G"

    method = method.lower()
    basis = basis.lower()

    if software == "gaussian":
        return ".com", method, basis
    elif software == "orca":
        return ".inp", method, basis
    elif software == "psi4":
        return ".in", method, basis
    else:
        raise ValueError(f"Unknown software: {software!r}")


def main(rc):
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

    ext = ext_from_soft(lines[0])


    for line in lines[2:]:
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        idx, inchi = parts[0], parts[1]
        outfolder = f"{idx}.conf"
        try:
            generate_conformers_from_inchi(inchi, outfolder, ext, rmsd_cutoff=rc)
        except Exception as e:
            print(f"Failed for {idx}: {e}")

    conf_dirs = glob.glob("*.conf")
    num_jobs = 0
    for cd in conf_dirs:
        os.chdir(cd)
        num_jobs = num_jobs +count()
        run_jobs()
        os.chdir(unique_dir)

    os.chdir(cwd)
    print(f"{num_jobs} jobs")

def main_cli():
    parser = argparse.ArgumentParser(description="Make unique conformers")
    parser.add_argument("--rmsd", default=".15", help="RMSD threshold, default = .15Å")
    args = parser.parse_args()
    main(args.rmsd)

if __name__ == "__main__":
    main(.15)


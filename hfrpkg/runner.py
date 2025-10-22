import os
from rdkit import Chem
from hfrpkg.utils import (
    geom_from_rdkit,
    isogyric_count,
    isodesmic_count,
    hypohomodesmotic_count,
    homodesmotic_count,
    display_reaction_counts,
)
from AaronTools.theory import Theory, OptimizationJob, FrequencyJob
from AaronTools.fileIO import FileWriter
import importlib.resources
#from hfrpkg.core import Isogyric, Isodesmic, Hypohomodesmotic, Homodesmotic
from hfrpkg.sandbox.Hfcore import Isogyric, Isodesmic, Hypohomodesmotic, Homodesmotic
reaction_map = {
    "homodesmotic": Homodesmotic,
    "isodesmic": Isodesmic,
    "isogyric": Isogyric,
    "hypohomodesmotic": Hypohomodesmotic,
}


extension_map = {
    "g": ".com",
    "o": ".inp",
    "p": ".in"
}
software_map = {
    "g": "Gaussian",
    "o": "ORCA",
    "p": "Psi4"
}

def get_Hf(inchi):
        try:
            with importlib.resources.open_text("hfrpkg.data", "ATcT_lib.txt", encoding="utf-8") as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) >= 6 and parts[3] == inchi:
                        return float(parts[5])
        except Exception:
            pass
        return None
def run_reaction(action_type, reaction_type, input_smiles, lhs=None, rhs=None, substruct=None, replacement=None, outfolder=None,method=None, basis=None, extension=None):

    mol = (Chem.MolFromInchi(input_smiles))
    if mol is not None:
        mol = Chem.AddHs(mol)
    if mol is None:

        mol = Chem.MolFromSmiles(input_smiles)
        if mol is not None:
            mol = Chem.AddHs(mol)
        else:
            raise ValueError(f"Invalid SMILES: {input_smiles}")
    if method is None:
        method = "B3LYP"
    if basis is None:
        basis = "6-31G"
    if extension is None:
        extension = "g"
    software = software_map[extension.lower()]
    ext = extension_map[extension.lower()]
    level = Theory(
            method=method,
            basis=basis,
            job_type=[OptimizationJob(), FrequencyJob()]
        )
    reaction_fn = reaction_map[reaction_type.lower()]
    rhs_mols, lhs_mols, status = reaction_fn(mol, lhs, rhs, substruct, replacement)

    if status != 'Optimal':
        print("Infeasible")
        return
    
    if action_type == "count":
        display_reaction_counts(mol, reaction_fn)

    elif action_type == "view":
        print("-----------Reactants-----------")
        for mol, coeff in lhs_mols:
            print(f"({coeff})*{Chem.MolToSmiles(mol)}")
        print("-----------Products-----------")
        for mol, coeff in rhs_mols:
            print(f"({coeff})*{Chem.MolToSmiles(mol)}")

    elif action_type == "write":
        if not outfolder:
            raise ValueError("Outfolder must be provided for write action")
        os.makedirs(outfolder, exist_ok=True)
        
        index_file = os.path.join(outfolder, "index.txt")
        with open(index_file, "w") as idx:
            idx.write(f"Level:\t{reaction_fn.__name__}\tSoftware:\t{software}\tMethod:\t{method}\tBasis:\t{basis}\n")
            idx.write(f"Input SMILES:\t{input_smiles}\n")
            idx.write("Filename\tInChI\tSMILES\n")
            Ri = Li = 1
            for mol, coeff in lhs_mols:
                geom = geom_from_rdkit(mol)
                smiles = Chem.MolToSmiles(mol)
                inchi = Chem.MolToInchi(mol)
                Hf = get_Hf(inchi)
                name = f"R{Li}_{coeff}"
                outfile=os.path.join(outfolder, name + ext)
                geom.write(outfile=outfile, theory=level)    
                idx.write(f"{name}\t{inchi}\t{smiles}\t{Hf}\n")
                Li += 1
            for mol, coeff in rhs_mols:
                geom = geom_from_rdkit(mol)
                smiles = Chem.MolToSmiles(mol)
                inchi = Chem.MolToInchi(mol)
                Hf = get_Hf(inchi)
                name = f"P{Ri}_{coeff}"
                outfile=os.path.join(outfolder, name + ext)
                geom.write(outfile=outfile, theory=level)
                idx.write(f"{name}\t{inchi}\t{smiles}\t{Hf}\n")
                Ri += 1
        print("Reaction written to", outfolder)

    print("Complete")

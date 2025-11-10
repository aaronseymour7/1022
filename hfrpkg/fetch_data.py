import os
import sys
import argparse
from rdkit import Chem
from AaronTools.fileIO import FileReader
from AaronTools.geometry import Geometry

def find_filename_by_inchi(inchi_query, index_file="index.txt"):
    with open(index_file, "r") as f:
        lines = f.readlines()

    data_started = False
    for line in lines:
        if line.strip().startswith("Filename"):
            data_started = True
            continue
        if not data_started:
            continue

        parts = line.strip().split("\t")
        if len(parts) < 2:
            parts = line.strip().split()
        if len(parts) < 2:
            continue

        filename = parts[0]
        inchi = parts[1]

        if inchi.strip() == inchi_query.strip():
            return filename
    return None

def get_extensions(index_path="index.txt"):
    ext_map = {
        "gaussian": ('.com', '.log'),
        "orca": ('.inp', '.out'),
        "psi4": ('.in', '.out')
    }
    try:
        with open(index_path) as f:
            first_line = f.readline().strip().lower()
            if "software:" in first_line:
                parts = first_line.replace(":", "").split()
                if "psi4" in parts:
                    return ext_map["psi4"]
                elif "gaussian" in parts:
                    return ext_map["gaussian"]
                elif "orca" in parts:
                    return ext_map["orca"]
    except FileNotFoundError:
        print(f"[ERROR] Could not find {index_path}")
    except Exception as e:
        print(f"[ERROR] Reading {index_path}: {e}")
    return ('.com', '.log')


def find_filename_by_smiles(smiles_query, index_file="index.txt"):
    with open(index_file, "r") as f:
        lines = f.readlines()

    data_started = False
    for line in lines:
        if line.strip().startswith("Filename"):
            data_started = True
            continue
        if not data_started:
            continue

        parts = line.strip().split("\t")
        if len(parts) < 3:
            parts = line.strip().split()
        if len(parts) < 3:
            continue

        filename = parts[0]
        smiles = parts[2]

        if smiles.strip() == smiles_query.strip():
            return filename
    return None


def main(input_smiles, info):
    cwd = os.getcwd()
    unique_path = os.path.join(cwd, "unique_files")
    if not os.path.exists(unique_path):
        print("Missing unique_files/")
        return

    os.chdir(unique_path)

    inext, outext = get_extensions()
    mol = Chem.MolFromSmiles(input_smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {input_smiles}")

    input_inchi = Chem.MolToInchi(mol)
    file = find_filename_by_inchi(input_inchi)


    if file is None:
        smiles = Chem.MolToSmiles(Chem.AddHs(mol))
        file = find_filename_by_smiles(smiles)

    if file is None:
        print("Molecule not found in unique_files/index.txt")
        return

    outfile = file.strip() + outext 
    if not os.path.exists(outfile):
        raise FileNotFoundError(f"Output file not found: {outfile}")

    if info == 'energy':
        try:
            reader = FileReader(outfile, just_geom=False)
            if 'E_ZPVE' in reader.keys():
                print(reader['E_ZPVE'])
                return reader['E_ZPVE']
            else:
                raise ValueError(f"No E_ZPVE found in {outfile}")
        except Exception as e:
            raise ValueError(f"Error reading energy: {e}")

    elif info == 'geom':
        try:
            geom = Geometry(outfile)
            if geom is None:
                raise ValueError(f"Failed to read geometry from {outfile}")
            print(geom)
            return geom
        except Exception as e:
            raise ValueError(f"Error reading geometry: {e}")


def main_cli():
    parser = argparse.ArgumentParser(description="Fetch output data from unique_files/index.txt")
    parser.add_argument("information", choices=["energy", "geom"], help="Information to retrieve")
    parser.add_argument("query", type=str, help="Input molecule SMILES string")
    args = parser.parse_args()

    main(args.query, args.information) 


if __name__ == "__main__":
    main_cli()

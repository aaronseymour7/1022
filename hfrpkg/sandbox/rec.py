import sys
import argparse
from rdkit import Chem
from hfrpkg.sandbox.pred import recurHomodesmotic as recurs
from hfrpkg.sandbox.pred import get_Hf
from hfrpkg.core import Homodesmotic

def Combo(input_mol):
    missing_inchis = []
    
    reaction = recurs(input_mol)
    if reaction is None:
        reaction = Homodesmotic(input_mol)
    if len(reaction) != 3:
        print('Infeasible')
        return None
    rhs, lhs, status = reaction
    if status != 'Optimal':
        print('Infeasible')
        return None
    for mol, coeff in lhs[:-1]:
        inchi = Chem.MolToInchi(mol)
        Hf = get_Hf(inchi)
        if Hf is None:
            missing_inchis.append(inchi)
    for mol, coeff in rhs:
        inchi = Chem.MolToInchi(mol)
        Hf = get_Hf(inchi)
        if Hf is None:
            missing_inchis.append(inchi)
    
    for i in missing_inchis:
        feas = recurs(Chem.MolFromInchi(i))
        if feas is None:
            print(f"cannot solve for {inchi}")
            print("Uncalculable")
            return None
        print(f"Homodesmotic can solve for {i}")
    return rhs,lhs,status

def main_cli():
    parser = argparse.ArgumentParser(description="Recursively caculable")
    parser.add_argument("input", type=str, help="Input molecule InChI")
    args = parser.parse_args()
    mol = Chem.MolFromInchi(args.input)
    if mol is None:
        print("Error: could not parse input InChI")
        return
    result = Combo(input_mol=mol)
    if result is None:
        print('Infeasible')
        return
   
    rhs_mols, lhs_mols, status = result

    if status != 'Optimal':
        print("Infeasible")
        return


    print("-----------Reactants-----------")
    for mol, coeff in lhs_mols:
        print(f"({coeff})*{Chem.MolToSmiles(mol)}")
    print("-----------Products-----------")
    for mol, coeff in rhs_mols:
        print(f"({coeff})*{Chem.MolToSmiles(mol)}")


if __name__ == "__main__":
    main_cli()

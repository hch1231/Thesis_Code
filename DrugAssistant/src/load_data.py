from rdkit import Chem
import os


def read_molecule_file(file_path):
    """
    Reads a molecular file and returns a list of RDKit molecule objects.

    Supports .smi, .sdf, and .xyz files.

    Args:
        file_path (str): Path to the molecular file.

    Returns:
        list of Chem.Mol: A list of RDKit molecule objects, or None if the file format is unsupported.
    """
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"The file '{file_path}' does not exist.")

    file_extension = os.path.splitext(file_path)[-1].lower()

    if file_extension == ".smi":
        supplier = Chem.SmilesMolSupplier(file_path)
        molecules = [mol for mol in supplier if mol is not None]

    elif file_extension == ".sdf":
        supplier = Chem.SDMolSupplier(file_path)
        molecules = [mol for mol in supplier if mol is not None]

    elif file_extension == ".xyz":
        molecules = read_xyz_file(file_path)

    else:
        raise ValueError(f"Unsupported file format: {file_extension}")

    return molecules


def read_xyz_file(file_path):
    """
    Reads an .xyz file and converts it into RDKit molecule objects.

    Args:
        file_path (str): Path to the .xyz file.

    Returns:
        list of Chem.Mol: A list containing a single RDKit molecule object.
    """
    with open(file_path, "r") as file:
        lines = file.readlines()

    atom_count = int(lines[0].strip())
    atoms = []
    coords = []

    for line in lines[2:2 + atom_count]:  # Skip the first two lines
        parts = line.split()
        atoms.append(parts[0])  # Element symbol
        coords.append([float(x) for x in parts[1:4]])  # Coordinates

    mol = Chem.RWMol()
    conf = Chem.Conformer(len(atoms))

    for i, (atom, coord) in enumerate(zip(atoms, coords)):
        atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(atom)
        mol.AddAtom(Chem.Atom(atomic_num))
        conf.SetAtomPosition(i, coord)

    mol.AddConformer(conf)
    return [mol]


# Example usage
if __name__ == "__main__":
    file_path = "D:\WorkSpace\PycharmFile\chemAgent\data\smi\AAAA.smi"  # Replace with your file path
    sdf_path = "D:\WorkSpace\PycharmFile\chemAgent\data\sdf\AAAARM.xaa.sdf\AAAARM.xaa.sdf"
    try:
        molecules = read_molecule_file(sdf_path)
        print(f"Successfully read {len(molecules)} molecules from {file_path}.")
    except Exception as e:
        print(f"Error: {e}")

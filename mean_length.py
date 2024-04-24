import os

def count_ca_atoms_in_pdb(filepath):
    """ Count the number of CA atoms in a single PDB file. """
    count = 0
    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith("ATOM") and " CA " in line:
                count += 1
    return count

def mean_residues_in_pdbs(directory):
    """ Calculate the mean number of residues in all PDB files in the specified directory. """
    files = [f for f in os.listdir(directory) if f.endswith('.pdb')]
    total_residues = 0
    for file in files:
        filepath = os.path.join(directory, file)
        total_residues += count_ca_atoms_in_pdb(filepath)
    
    if len(files) > 0:
        return total_residues / len(files)
    else:
        return 0

def unique_chains_in_pdb(filepath):
    """ Count the number of unique chains in a single PDB file. """
    chains = set()
    with open(filepath, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                parts = line.split()
                if len(parts) > 5:
                    chain_id = parts[4]  # Fifth column for chain ID
                    chains.add(chain_id)
    return len(chains)

def mean_chains_in_pdbs(directory):
    """ Calculate the mean number of unique chains in all PDB files in the specified directory. """
    files = [f for f in os.listdir(directory) if f.endswith('.pdb')]
    total_chains = 0
    for file in files:
        filepath = os.path.join(directory, file)
        total_chains += unique_chains_in_pdb(filepath)
    
    if len(files) > 0:
        return total_chains / len(files)
    else:
        return 0

# Specify the directory containing your PDB files
directory_path = '/home/simon.herman/Téléchargements/peptides'
mean_chains = mean_chains_in_pdbs(directory_path)
print(f"The mean number of unique chains per protein is: {mean_chains:.2f}")

# Specify the directory containing your PDB files
directory_path = '/home/simon.herman/Téléchargements/peptides'
mean_length = mean_residues_in_pdbs(directory_path)
print(f"The mean length of proteins is: {mean_length:.2f} residues")

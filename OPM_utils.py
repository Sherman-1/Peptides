
__all__ = ["size_picker", "read_pdb"]

from Bio import SeqIO

def size_picker(fasta_file, min_size=20, max_size=70):
    """Reads a FASTA file and returns a randomly picked size from its sequences."""
    import random
    sizes = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        size = len(record.seq)
        if min_size <= size <= max_size:
            sizes.append(size)
    return random.choice(sizes)


from collections import OrderedDict

def read_pdb(input_file):
    res_coord_dic = OrderedDict() # Store all residues along with their z values

    membrane_x_values = []
    membrane_y_values = []
    membrane_z_values = []


    # The two membrane z-values
    inner_membrane_z = min(membrane_z_values) - margin
    outer_membrane_z = max(membrane_z_values) + margin
    # The membrane x and y limit values to avoid peptide segments that are far from the represented membrane (this is indicative of abnormal peptide segments).
    min_membrane_x = min(membrane_x_values) - 10
    max_membrane_x = max(membrane_x_values) + 10
    min_membrane_y = min(membrane_y_values) - 10
    max_membrane_y = max(membrane_y_values) + 10
    with open(input_file, 'r') as f_in:
        for line in f_in:
            if line.startswith('ATOM') and line[13:15].strip() == 'CA':
                chain = line[21]
                x_value = float(line[30:38])
                y_value = float(line[38:46])
                z_value = float(line[46:54])
                resSeq = int(line[22:26])
                res_coord_dic[chain,resSeq] = {'x': x_value, 'y': y_value, 'z': z_value}
            elif line.startswith('HETATM') and 'DUM' in line:
                membrane_x_values.append(float(line[30:38]))
                membrane_y_values.append(float(line[38:46]))
                membrane_z_values.append(float(line[46:54]))

    return res_coord_dic
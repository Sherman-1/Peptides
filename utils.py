import numpy as np
from collections import defaultdict
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random
import warnings

IORF_PATH = "Scer_NCBI_iORF.faa"


class SuperOD(OrderedDict):
    def __getitem__(self, key):
        if isinstance(key, slice):
            
            keys = list(self.keys())[key]
            return SuperOD((k, self[k]) for k in keys)
        else:
            
            return OrderedDict.__getitem__(self, key)
        
class ResidueError(ValueError):
    """Error to handle non canonic residues in files"""
    def __init__(self, res, message="Residue not in the list of canonic residues"):
        self.res = res
        self.message = message
        super().__init__(f"{self.message}: {self.res}")

def size_picker_v3(fasta_file = IORF_PATH, min_length = 0, max_length = 1000, n_samples = 1):

    """
    Same thing as v2 but draw n_samples sizes from the distribution
    all at once, avoiding to parse the fasta file n_samples times
    """

    sizes = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        length = len(record.seq)
        if min_length <= length <= max_length:

            sizes.append(length)

    return random.choices(sizes, k=n_samples)

def pdb_struct_to_fasta(pdb_struct, writing_path = ".", write = False):
    
    records = []
    
    for chain_id in pdb_struct["CA"]:

        sequence = ""

        for res_number, data in pdb_struct["CA"][chain_id].items():

            sequence += data["res_name"]
            
        records = SeqRecord(Seq(sequence), id=f"{pdb_struct['protein_name']}_{chain_id}", description=f"")

    if write == True:

        if writing_path.endswith("/"):
             
             writing_path = writing_path[:-1]

        SeqIO.write(records, f"{writing_path}/{pdb_struct['protein_name']}.fasta", "fasta")

    return records
    

    
def read_pdb(file_path, secondary_structure_path = None):

    """
    Function to process PDB format records.

    PDB Record Format:
    
    Line starts with ATOM
    ----
    1-4    : "ATOM"                            (character)
    7-11   : Atom serial number                (right, integer)
    13-16  : Atom name                         (left*, character)
    17     : Alternate location indicator      (character)
    18-20  : Residue name                      (right, character)
    22     : Chain identifier                  (character)
    23-26  : Residue sequence number           (right, integer)
    27     : Code for insertions of residues   (character)
    31-38  : X orthogonal Å coordinate         (right, real (8.3))
    39-46  : Y orthogonal Å coordinate         (right, real (8.3))
    47-54  : Z orthogonal Å coordinate         (right, real (8.3))
    55-60  : Occupancy                         (right, real (6.2))
    61-66  : Temperature factor                (right, real (6.2))
    73-76  : Segment identifier                (left, character)
    77-78  : Element symbol                    (right, character)
    79-80  : Charge                            (character)
    
    Line starts with HETATM
    ------
    1-6    : "HETATM"                          (character)
    7-80   : Same as ATOM records

    Line starts with TER
    ---
    1-3    : "TER"                             (character)
    7-11   : Serial number                     (right, integer)
    18-20  : Residue name                      (right, character)
    22     : Chain identifier                  (character)
    23-26  : Residue sequence number           (right, integer)
    27     : Code for insertions of residues   (character)  
    """
    

    aa_dict = {
                    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
                    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
                    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
                }

    pdb_struct = {}

    pdb_struct["protein_name"] = file_path.split("/")[-1].split(".")[0]

    pdb_struct["full"] = defaultdict(SuperOD)
    pdb_struct["CA"] = defaultdict(SuperOD)
    pdb_struct["membrane_coord"] = []

    array = []

    with open(file_path,"r") as f:

        while line := f.readline():

            suspected_membrane_line = False

            if line.startswith("ATOM"):

                # Line format is given in the docstring

                atom_number = int(line[6:11].strip())
                atom_name = line[12:16].strip()

                res_number = int(line[22:26].strip())
                res_name = line[17:20].strip()

                chain_id = line[21]

                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                if res_name == "DUM":

                    suspected_membrane_line = True

                elif res_name not in aa_dict:

                    raise ResidueError(res_name)

                if not suspected_membrane_line:

                    # Is the chain_id empty ?
                    if not bool(chain_id.strip()):

                        chain_id = "A"

                    if chain_id not in pdb_struct["full"]:

                        pdb_struct["full"][chain_id] = SuperOD()
                        pdb_struct["CA"][chain_id] = SuperOD()

                        pdb_struct["full"][chain_id][atom_number] = {

                            "coord" : [x,y,z],
                            "atom_name" : atom_name,
                            "res_name" : aa_dict[res_name],
                            "res_number" : res_number
                        }


                    else:

                        pdb_struct["full"][chain_id][atom_number] = {

                            "coord" : [x,y,z],
                            "atom_name" : atom_name,
                            "res_name" : aa_dict[res_name],
                            "res_number" : res_number
                        }

                    
                        if atom_name == "CA":

                            pdb_struct["CA"][chain_id][res_number] = {

                                "coord" : [x,y,z],
                                "res_name" : aa_dict[res_name],
                                "res_number" : res_number,
                            }

            if (line.startswith("HETATM") and "DUM" in line) or suspected_membrane_line:

                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())

                array.append([x,y,z])


    pdb_struct["membrane_coord"] = np.array(array)
    pdb_struct["protein_length"] = { }

    # For each chain, compute the length.
    # The pdbs are never evenly formated, thus we can't 
    # infer the length from the residues numbers
    # We just count the number of CA atoms for each chain
    for chain_id in pdb_struct["CA"]:
       
        pdb_struct["protein_length"][chain_id] = len(pdb_struct["CA"][chain_id])

    if secondary_structure_path:

        with open(secondary_structure_path, "r") as f:

            
            
            # Line format given by extract_SS.cpp 
            # chain_id << '\t' << res_number << '\t' << secondary_structure << '\n';
            while line := f.readline():

                line = line.split()
                chain_id = line[0]
                res_number = int(line[1])
                secondary_structure = line[2]

                pdb_struct["CA"][chain_id][res_number].update({"secondary_structure" : secondary_structure})

                if secondary_structure == "H" or secondary_structure == "S":

                    pdb_struct["CA"][chain_id][res_number].update({"folded" : True})

                else:

                    pdb_struct["CA"][chain_id][res_number].update({"folded" : False})

    return pdb_struct

def binarize_structure(pdb_struct: dict, lower_margin=0, margin=10):
    
    min_z_membrane = np.min(pdb_struct["membrane_coord"][:, 2])
    max_z_membrane = np.max(pdb_struct["membrane_coord"][:, 2])

    min_x_membrane = np.min(pdb_struct["membrane_coord"][:, 0]) - 10
    max_x_membrane = np.max(pdb_struct["membrane_coord"][:, 0]) + 10

    min_y_membrane = np.min(pdb_struct["membrane_coord"][:, 1]) - 10
    max_y_membrane = np.max(pdb_struct["membrane_coord"][:, 1]) + 10

    in_membrane_binaries = { chain_id: "" for chain_id in pdb_struct["CA"].keys() }
    in_margin_binaries = { chain_id: "" for chain_id in pdb_struct["CA"].keys() }
    
    for chain_id, residues in pdb_struct["CA"].items():

        last_res_number = list(residues.keys())[0]
        print(f"First residue number : {last_res_number}")
        for res_number, data in residues.items():

            x = data["coord"][0]
            y = data["coord"][1]
            z = data["coord"][2]

            # Sometimes, residues are missing from the 3D structure and we don't have their coordinates
            # We can monitor this by checking if the current residue number is not the previous residue number + 1

            if int(res_number) != last_res_number + 1:

                # Fill in the gaps with zeros, by default
                for i in range(int(res_number) - last_res_number - 1):
                    in_membrane_binaries[chain_id] += "0"
                    in_margin_binaries[chain_id] += "0"
                
                # Then compute the binary string for the current residue
                in_membrane = "1" if min_z_membrane <= z <= max_z_membrane + lower_margin else "0"
                in_margin = "1" if min_z_membrane - margin <= z <= max_z_membrane + margin else "0"

                # Update the pdb_struct dict on the fly 
                pdb_struct["CA"][chain_id][res_number].update({"in_membrane" : in_membrane, "in_margin" : in_margin})
                
            else:

                # Residues that are out of the x / y membrane bounds are usually 
                # pathologic predictions and should be discarded

                is_x = min_x_membrane <= x <= max_x_membrane
                is_y = min_y_membrane <= y <= max_y_membrane
                is_z = min_z_membrane <= z <= max_z_membrane
                is_z_plus_margin = min_z_membrane - margin <= z <= max_z_membrane + margin

                in_membrane = "1" if is_x and is_y and is_z else "0"
                in_margin = "1" if is_x and is_y and is_z_plus_margin else "0"

                # Update the pdb_struct dict on the fly
                pdb_struct["CA"][chain_id][res_number].update({"in_membrane" : in_membrane, "in_margin" : in_margin})


            in_membrane_binaries[chain_id] += in_membrane
            in_margin_binaries[chain_id] += in_margin

            last_res_number = int(res_number)

    
    return in_membrane_binaries, in_margin_binaries

def define_tm_segments(binary_dict : dict):

    segment_indices = {chain_id : [] for chain_id in binary_dict.keys() }

    start_index = None

    for chain_id in binary_dict:
        for i, bit in enumerate(binary_dict[chain_id]):
            if bit == "1":
                if start_index is None:
                    start_index = i+1
            else:
                if start_index is not None:
                    length = i-start_index+1
                    if length >= 15: # minimum length of a TM segment, although 20 is the length of a typical alpha helical TM segment
                        
                        # python indices are 0-based, so we add 1 to match 
                        # the 1-based residue numbering in the PDB file
                        segment_indices[chain_id].append((start_index, i+1, length))

                    # Wether the segment is long enough or not, we reset the start_index
                    start_index = None

        # If we ended in the middle of a TM segment, we add it to the list
        if start_index is not None:
            length = len(binary_dict[chain_id]) - start_index + 1
            if length >= 15:
                segment_indices[chain_id].append((start_index, len(binary_dict[chain_id]), length))

    return segment_indices

def elongate_tm_segments(tm_indices, pdb_struct, min_length=20, max_length=70, verbose = False):
    """
    This function elongates transmembrane segments to a random size drawn from a given size distribution,
    distributing the elongation both upstream and downstream without considering overlapping boundaries between segments.

    Parameters:
    - tm_indices (dict): Dictionary where keys are chain_ids and values are lists of tuples (start, end, length)
      for transmembrane segments.
    - pdb_struct (dict): Dictionary containing protein-related information such as protein lengths.
    - min_length (int): Minimum length of the elongated segment.
    - max_length (int): Maximum length of the elongated segment.

    Returns:
    - int: Function completion status (0 for successful execution).
    """


    for chain_id in tm_indices:

        # Draw the random elongation lengts once, as it needs to read a fasta file each time
        desired_lengths = size_picker_v3(min_length=min_length, max_length=max_length, n_samples=len(tm_indices[chain_id]))
        chain_length = pdb_struct["protein_length"][chain_id]

        for i, (start_current, end_current, length_current) in enumerate(tm_indices[chain_id]):

            desired_length = desired_lengths[i]
            old_start = start_current
            old_end = end_current

            if desired_length <= length_current:

                tm_indices[chain_id][i] = (start_current, end_current, old_start, old_end)
                
            elongation_needed = desired_length - length_current

            # Randomly determine the amount of elongation downstream and upstream
            downstream = random.randint(0, elongation_needed)
            upstream = elongation_needed - downstream

            if verbose:
                print(f"Chain {chain_id}, segment {i+1}")
                print(f"Desired length: {desired_length}")
                print(f"Current length: {length_current}")
                print(f"Upstream elongation: {upstream}")
                print(f"Downstream elongation: {downstream}")

            # Calculate new coordinates ensuring they stay within the protein bounds
            new_end_coordinates = min(end_current + downstream, chain_length)
            new_start_coordinates = max(start_current - upstream, 1)

            # Update the segment information with new start and end positions
            tm_indices[chain_id][i] = (new_start_coordinates, new_end_coordinates, old_start, old_end)

    return 0

def extract_elongated_sequences(tm_indices : dict, pdb_struct : dict, gaps : int = 2):
    
    records = []
    records_shorts = []
    protein_name = pdb_struct["protein_name"]

    log = open(f"{protein_name}.log","w")
    

    for chain_id in tm_indices:

        for i, (start, end, old_start, old_end) in enumerate(tm_indices[chain_id]):

            sequence = ""
            sequence_short = ""
            buffer = "" 
            have_left_margin = False
            have_left_membrane = False

            # Look in pdb_struct for the current elongated segment
            for res_number in range(start, end+1):

                # Security check
                if res_number in pdb_struct["CA"][chain_id]:

                    # We are before the TM segment that is elongated
                    if res_number < old_start:
                        
                        # We are before the TM and in the margin : 
                        if pdb_struct["CA"][chain_id][res_number]["in_margin"] == "1":
                            
                            if have_left_margin:

                                sequence += pdb_struct["CA"][chain_id][res_number]["res_name"]
                                print(f"Residue {res_number} in margin and in sequence ")
                                print(f"Buffer at this point : {buffer}")

                            else:
                                buffer += pdb_struct["CA"][chain_id][res_number]["res_name"]
                                print(f"Residue {res_number} in margin and in buffer ")
                                print(f"Buffer at this point : {buffer}")

                        # We are still before the segment but not in the margin
                        if pdb_struct["CA"][chain_id][res_number]["in_margin"] == "0":
                            
                            if 1 < res_number:

                                margin_plus_one = pdb_struct["CA"][chain_id][res_number+1]["in_margin"] == "1"
                                margin_minus_one = pdb_struct["CA"][chain_id][res_number-1]["in_margin"] == "1"
                                margin_plus_two = pdb_struct["CA"][chain_id][res_number+2]["in_margin"] == "1"
                                margin_minus_two = pdb_struct["CA"][chain_id][res_number-2]["in_margin"] == "1"

                                # The current residue is not in the margin but the surrounding ones
                                # are => we can accept it as a margin residue
                                if (margin_plus_one or margin_plus_two) and (margin_minus_one or margin_minus_two):
                                    if have_left_margin:
                                        sequence += pdb_struct["CA"][chain_id][res_number]["res_name"]
                                        print(f"Residue {res_number} surrounded and we have left the margin so in sequence")
                                    else:
                                        buffer += pdb_struct["CA"][chain_id][res_number]["res_name"]
                                        print(f"Residue {res_number} surrounded and we have not left the margin so in buffer")

                                # The current residue is not in the margin and the residues around it 
                                # are not either : 
                                else:
                                    have_left_margin = True
                                    buffer = ""
                                    print(f"#############")
                                    print(f"Residue {res_number} not surrounded, we left the margin")
                                    print(f"#############")

                            # We are not in the margin but 
                            else:
                                have_left_margin = True
                                buffer = ""
                                print(f"#############")
                                print(f"Residue {res_number} less than 1, we left the margin")
                                print(f"#############")

                    elif res_number == old_start:

                        # We know arrived at the beginning of the TM segment
                        # before it was elongated by the elongate_tm_segments function
                        # If the residues between this one and the previous one 
                        # were all in the margin, the buffer is not empty
                        print("\n\n")
                        print(f"We arrived at TM, buffer : {buffer}")
                        

                        if buffer and not have_left_margin:
                            sequence += buffer
                            buffer = ""
                            print(f"Buffer added to sequence : {sequence}")
                        
                        sequence += pdb_struct["CA"][chain_id][res_number]["res_name"]
                        print("\n\n")

                    # we are in the TM segment
                    elif old_start < res_number <= old_end:

                        sequence += pdb_struct["CA"][chain_id][res_number]["res_name"]

                    # We know left the TM segment and are in the C-terminal portion of the elongated tm-segment
                    elif res_number > old_end:

                        
                        # If we are in the membrane but we already left it, 
                        # it means we are BACK in it : we are in the next TM segment
                        if pdb_struct["CA"][chain_id][res_number]["in_membrane"] == "1":

                            # old_end + 2 to be sure that we really left the tm segment that is elongated
                            # so far i've only seen inter-tm segments up to 2 amino acids ...
                            if res_number > old_end + 2:
                                
                                if have_left_membrane:
                                    break

                            else:
                                # We are in the membrane so we are in the margin, by definition
                                sequence += pdb_struct["CA"][chain_id][res_number]["res_name"]

                        # Not in the membrane 
                        elif pdb_struct["CA"][chain_id][res_number]["in_membrane"] == "0":
                            
                            # We are passed the coordinates of the original TM segment
                            # => we left the membrane
                            have_left_membrane = True

                            # As usual
                            if pdb_struct["CA"][chain_id][res_number]["in_margin"] == "1":

                                sequence += pdb_struct["CA"][chain_id][res_number]["res_name"]

                            # We give ourselves some margin ( lol ), if a residue is not in the margin
                            # but the 2 surrounding ones are, the current one must be " floating " 
                            # away from the margin limit => we take it 
                            elif pdb_struct["CA"][chain_id][res_number]["in_margin"] == "0":

                                margin_plus_one = pdb_struct["CA"][chain_id][res_number+1]["in_margin"]
                                margin_minus_one = pdb_struct["CA"][chain_id][res_number-1]["in_margin"]

                                if margin_plus_one == "1" and margin_minus_one == "1":

                                    sequence += pdb_struct["CA"][chain_id][res_number]["res_name"]

                                # i_th AA is not in margin, surrounding ones neither, we know we 
                                # already left the principal tm segment 
                                # => we are leaving the margin for good, break the loop
                                else:
                                    break

                else:

                    if pdb_struct["CA"][chain_id][res_number-1]["in_membrane"] == "1" and pdb_struct["CA"][chain_id][res_number+1]["in_membrane"] == "1":
                        
                        sequence = ""
                        sequence_short = ""
                        buffer = ""

                        log.write(f"Residue {res_number} is in the membrane and not defined !")
                        warnings.warn(f'{protein_name} has {res_number} not defined and in the membrane !')
                        continue

                    if pdb_struct["CA"][chain_id][res_number-1]["in_margin"] == "1" and pdb_struct["CA"][chain_id][res_number+1]["in_margin"] == "1":

                        if pdb_struct["CA"][chain_id][res_number + 1]["folded"] and pdb_struct["CA"][chain_id][res_number - 1]["folded"]:

                            log.write(f"Residue {res_number} is missing and surrounded by secondary structure inside the margin\n")
                            sequence += "X"
                        
                    sequence += "X"

            # Just write the amino acids of the tm-segments without elongation on a separate fasta
            for res_number in range(old_start, old_end+1):
                
                if res_number in pdb_struct["CA"][chain_id]:
                    sequence_short += pdb_struct["CA"][chain_id][res_number]["res_name"]
                else:
                    sequence_short += "X"
            
            
            sequence = sequence.strip("X")
            sequence_short = sequence_short.strip("X")
            
            if not "X" * gaps in sequence:
            
                record = SeqRecord(Seq(sequence), id=f"{protein_name}_{chain_id}_{i+1}", description=f"")
                records.append(record)
                
            if not "X" * gaps in sequence_short:
            
                record = SeqRecord(Seq(sequence_short), id=f"{protein_name}_{chain_id}_{i+1}_short", description=f"")
                records_shorts.append(record)

    log.close()

    return records, records_shorts


def test(pdb_path = "all_pdbs/1uaz.pdb"):

    random.seed(random.randint(0,1000))

    pdb = read_pdb(pdb_path)
    print(list(pdb["CA"]["A"].keys())[0])
    in_membrane_binaries, in_margin_binaries = binarize_structure(pdb)
    tm_indices = define_tm_segments(in_membrane_binaries)
    elongate_tm_segments(tm_indices, pdb, verbose = False, min_length = 60, max_length=70)
    records, records_short = extract_elongated_sequences(tm_indices, pdb)
    pdb_struct_to_fasta(pdb, write = True)
    SeqIO.write(records, "test.fasta", "fasta")

    return records









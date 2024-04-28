import numpy as np
from collections import defaultdict
from collections import OrderedDict
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import SeqRecord
import random
from scipy.spatial.distance import pdist

random.seed(25032024)

def compute_pairwise_distances(pdb_struct):
    # Extract coordinates from pdb_struct
    coords = []
    for chain in pdb_struct["full"].values():
        for residue in chain.values():
            coords.append(residue["coord"])

    # Return mean pairwise distance
    return np.mean(pdist(coords))

def size_picker_v2(fasta_file = IORF_PATH, min_length = 0, max_length = 1000):

    sizes = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        length = len(record.seq)
        if min_length <= length <= max_length:

            sizes.append(length)
    
    return random.choice(sizes)

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

def pdb_struct_to_fasta(pdb_struct):

    aa_dict = {
                    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
                    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
                    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
                }
    fasta = ""

    for res_number in pdb_struct["CA"]:

        res_name = pdb_struct["CA"][res_number]["res_name"]

        fasta += aa_dict[res_name]

    record = SeqRecord.SeqRecord(Seq(fasta), id="1uaz", description="1uaz")
    filename = f"{pdb_struct['protein_name']}.fasta"
    SeqIO.write(record, filename, "fasta")

    return 0

class SuperOD(OrderedDict):
    def __getitem__(self, key):
        if isinstance(key, slice):
            
            keys = list(self.keys())[key]
            return SuperOD((k, self[k]) for k in keys)
        else:
            
            return OrderedDict.__getitem__(self, key)
    
def read_pdb(file_path):
    
    aa_dict = {
                    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
                    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
                    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
                }

    pdb_struct = {}

    pdb_struct["protein_name"] = file_path.split(".")[0]

    pdb_struct["full"] = defaultdict(SuperOD)
    pdb_struct["CA"] = defaultdict(SuperOD)
    pdb_struct["membrane_coord"] = []

    array = []

    with open(file_path,"r") as f:

        line = f.readline()

        while line:

            line = line.split()

            if line[0] == "ATOM":

                # Line format : 
                # ATOM      2  CA  MET A   1      24.767  -2.102  13.513  1.00  0.00      A1A9 C  

                x = float(line[6])
                y = float(line[7])
                z = float(line[8])

                atom_name = line[2]
                atom_number = int(line[1])
                res_name = line[3]
                res_number = int(line[5])
                
                chain_id = line[4]

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

                
                    if line[2] == "CA":

                        pdb_struct["CA"][chain_id][res_number] = {

                            "coord" : [x,y,z],
                            "res_name" : aa_dict[res_name],
                            "res_number" : res_number,
                        }

            elif line[0] == "HETATM" and "DUM" in line:

                # Line format :
                # HETATM  643  O   DUM   643     -24.000  -6.000  14.200   
            
                x = float(line[-3])
                y = float(line[-2])    
                z = float(line[-1])

                array.append([x,y,z])

            line = f.readline()

    pdb_struct["membrane_coord"] = np.array(array)
    pdb_struct["protein_length"] = { }

    # For each chain, compute the length.
    # The length is the highest value of the residue numbers
    for chain_id in pdb_struct["CA"]:
       
        pdb_struct["protein_length"][chain_id] = max([int(res_number) for res_number in pdb_struct["CA"][chain_id].keys()])
        print(f"Chain {chain_id} has length {pdb_struct['protein_length'][chain_id]}")

    return pdb_struct

def return_binaries(pdb_struct: dict, lower_margin=0, margin=5):
    
    min_z_membrane = np.min(pdb_struct["membrane_coord"][:, 2])
    max_z_membrane = np.max(pdb_struct["membrane_coord"][:, 2])

    in_membrane_binaries = { chain_id: "" for chain_id in pdb_struct["CA"].keys() }
    in_margin_binaries = { chain_id: "" for chain_id in pdb_struct["CA"].keys() }

    # Assuming that the residues start at 1
    last_res_number = 0
    for chain_id, residues in pdb_struct["CA"].items():
        for res_number, data in residues.items():

            z = data["coord"][2]

            # Sometimes, residues are missing from the 3D structure and we don't have their coordinates
            # We can monitor this by checking if the current residue number is not the previous residue number + 1
            # If some residues are missing, fill the binary string with zeros corresponding to the number of missing residues
            # By default we assume that the missing residues are not of interest and are not in the membrane

            if int(res_number) != last_res_number + 1:

                # Fill in the gaps with zeros
                for i in range(int(res_number) - last_res_number - 1):
                    in_membrane_binaries[chain_id] += "0"
                    in_margin_binaries[chain_id] += "0"
                
                # Then compute the binary string for the current residue
                in_membrane = "1" if min_z_membrane <= z <= max_z_membrane + lower_margin else "0"
                in_margin = "1" if min_z_membrane - margin <= z <= max_z_membrane + margin else "0"

                # Update the pdb_struct dict on the fly 
                pdb_struct["CA"][chain_id][res_number].update({"in_membrane" : in_membrane, "in_margin" : in_margin})
                
            else:

                in_membrane = "1" if min_z_membrane <= z <= max_z_membrane + lower_margin else "0"
                in_margin = "1" if min_z_membrane - margin <= z <= max_z_membrane + margin else "0"

                # Update the pdb_struct dict on the fly
                pdb_struct["CA"][chain_id][res_number].update({"in_membrane" : in_membrane, "in_margin" : in_margin})

            in_membrane_binaries[chain_id] += in_membrane
            in_margin_binaries[chain_id] += in_margin

            last_res_number = int(res_number)

    # Now before returning, we can assume that the sequence 101 
    # corresponds to an error. We could like to change the 101 to 111 
    return in_membrane_binaries, in_margin_binaries

def extract_tm_segments_indices(binary_dict : dict):

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

def elongate_tm_segments(tm_indices : dict, pdb_struct : dict, min_length=20, max_length=70):
    """
    This function takes a list of tuples containing the start and end indices of putative transmembrane (tm) segments
    Extracted for the same multiple-fragments transmembrane protein.
    For example, GPCR proteins have 7 transmembrane segments, they will end up in a list of 7 tuples.

    For each tm segment, the function will elongate the segment to a random size drawn from a given size distribution,
    given by the size_picker_v2 function.

    The function will elongate the segment only if the size of the segment is smaller than the size drawn from the distribution.
    The goal here is to "draw" from the parts of the sequence that are not transmembrane segments, and elongate the tm segments.
    One main goal is to avoid drawing twice from the same region to elongate two tm segments that are adjacent to each other.

    Input:

    tm_indices[chain_id] : list of tuples 
                # [ (12,26,15), (45, 60, 16), (80, 100, 21) ...]
                # [ (start, end, length), ... ]

    min_length : int
                # minimum length of the elongated segment

    max_length : int
                # maximum length of the elongated segment
    """

    for chain_id in tm_indices:

        
        protein_length = pdb_struct["protein_length"][chain_id]

        ##### Treat first TM Segment separately ##### 


        desired_length = size_picker_v2(min_length=min_length, max_length=max_length)

        
        # First TM Segment
        start_current = tm_indices[chain_id][0][0]
        end_current = tm_indices[chain_id][0][1]
        length_current = tm_indices[chain_id][0][2]


        if desired_length > length_current:


            # Second TM Segment
            start_next = tm_indices[chain_id][1][0]

            elongation_left_to_do = desired_length - length_current


            downstream = random.randint(0, elongation_left_to_do)

            
            lefts = None

            # The new end of this tm segment should not exceed the start of the next tm segment
            if downstream + end_current > start_next:

                new_end_coordinates = start_next - 1

                lefts = downstream - ( start_next - end_current )



            else:

                new_end_coordinates = end_current + downstream


            upstream = elongation_left_to_do - downstream

            

            if lefts:

                upstream += lefts



            if start_current - upstream < 1:

                new_start_coordinates = 1



            else:

                new_start_coordinates = start_current - upstream

            tm_indices[chain_id][0] = (new_start_coordinates, new_end_coordinates, new_end_coordinates - new_start_coordinates)


        ##### Treat from the second TM Segment to the penultimate one ( n-1 ) #####

        for i in range(1, len(tm_indices[chain_id]) - 1):

            # Target size that the current tm should reach
            desired_length = size_picker_v2(min_length=min_length, max_length=max_length)

            # ith TM Segment
            start_current = tm_indices[chain_id][i][0]
            end_current = tm_indices[chain_id][i][1]
            length_current = tm_indices[chain_id][i][2]

            # check before anything else to save computation time
            if desired_length <= length_current:

                # If there is no elongation to do, we skip to the next segment
                # and the coordinates of the ith segment are not modified
                continue
            
            # (i+1)th TM Segment
            start_next = tm_indices[chain_id][i+1][0]


            # (i-1)th TM Segment
            end_previous = tm_indices[chain_id][i-1][1]
            
            # Compute the number of residues that are required to elongate the current segment
            elongation_left_to_do = desired_length - length_current


            # Randomly choose the number of residues to elongate downstream ( toward the C-terminal )
            downstream = random.randint(0, elongation_left_to_do)

            lefts = None

            # The new end of this tm segment should not exceed the start of the next tm segment
            if downstream + end_current > start_next:

                # Hence take everyting that is between the end of the current tm segment and the start of the next one
                new_end_coordinates = start_next - 1

                # What is " left " from downstream that could not be taken cause of the next tm ? 
                lefts = downstream - (start_next - end_current)

            else:

                new_end_coordinates = end_current + downstream

            ## If there is elongation that was not taken from downstream, add it to the upstream
            upstream = elongation_left_to_do - downstream
            if lefts:

                upstream += lefts


            # The new start of this tm segment should not be lower than the end of the previous tm segment
            if start_current - upstream < end_previous:

                new_start_coordinates = end_previous + 1 

            else:

                new_start_coordinates = start_current - upstream


            tm_indices[chain_id][i] = (new_start_coordinates, new_end_coordinates, new_end_coordinates - new_start_coordinates)

            


        ##### Treat the last TM Segment #####

        # Target size that the current tm should reach
        desired_length = size_picker_v2(min_length=min_length, max_length=max_length)


        # Last TM Segment
        start_current = tm_indices[chain_id][-1][0]
        end_current = tm_indices[chain_id][-1][1]
        length_current = tm_indices[chain_id][-1][2]

        # check before anything else to save computation time
        if desired_length <= length_current:

            # If there is no elongation to do, we skip to the next segment
            # and the coordinates of the ith segment are not modified
            return 0

        # (i-1)th TM Segment
        end_previous = tm_indices[chain_id][-2][1]

        # Compute the number of residues that are required to elongate the current segment
        elongation_left_to_do = desired_length - length_current



        # Randomly choose the number of residues to elongate downstream ( toward the C-terminal )
        downstream = random.randint(0, elongation_left_to_do)

        lefts = None

        # The new end of this final tm should not exceed the protein length
        if downstream + end_current > protein_length:

            # Hence take everyting that is between the end of the current tm segment and the start of the next tm segment
            new_end_coordinates = protein_length

            # What is " left " from downstream that could not be taken because the protein is too short after the last tm ? 
            lefts = downstream - (protein_length - end_current)


        else:

            new_end_coordinates = end_current + downstream


        upstream = elongation_left_to_do - downstream
        if lefts:

            upstream += lefts


        # The new start of this tm segment should not be lower than the end of the previous tm segment
        if start_current - upstream < end_previous:

            new_start_coordinates = end_previous + 1 

        else:

            new_start_coordinates = start_current - upstream     


        tm_indices[chain_id][-1] =(new_start_coordinates, new_end_coordinates, new_end_coordinates - new_start_coordinates + 1)

    return 0

def elongate_tm_segments_overlap(tm_indices, pdb_struct, min_length=20, max_length=70):
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

            # Calculate new coordinates ensuring they stay within the protein bounds
            new_end_coordinates = min(end_current + downstream, chain_length)
            new_start_coordinates = max(start_current - upstream, 1)

            # Update the segment information with new start and end positions
            new_length = new_end_coordinates - new_start_coordinates + 1
            tm_indices[chain_id][i] = (new_start_coordinates, new_end_coordinates, old_start, old_end)

    return 0


def write_segments_v2(tm_indices : dict, pdb_struct : dict):
    
    records = []
    records_shorts = []
    protein_name = pdb_struct["protein_name"]
    output_file = f"{protein_name}_TM_segments.fasta"
    output_file_short = f"{protein_name}_TM_segments_short.fasta"
    

    for chain_id in tm_indices:

        protein_length = pdb_struct["protein_length"][chain_id]

        for i, (start, end, old_start, old_end) in enumerate(tm_indices[chain_id]):

            sequence = ""
            sequence_short = ""
            buffer = "" 
            have_left_margin = False
            have_left_membrane = False

            # Look in pdb_struct for the current elongated segment
            for res_number in range(start, end+1):

                if res_number in pdb_struct["CA"][chain_id]:

                    # We are in elongated residues, i.e the ones that were not in the TM segment
                    if res_number < old_start:
                        
                        # We are before the TM and in the margin : 
                        if pdb_struct["CA"][chain_id][res_number]["in_margin"] == "1":
                            
                            if have_left_margin:
                                sequence += pdb_struct["CA"][chain_id][res_number]["res_name"]

                            else:
                                buffer += pdb_struct["CA"][chain_id][res_number]["res_name"]

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
                                    else:
                                        buffer += pdb_struct["CA"][chain_id][res_number]["res_name"]

                                # The current residue is not in the margin and the residues around it 
                                # are not either : 
                                else:
                                    have_left_margin = True

                                    

                            else:
                                have_left_margin = True

                    elif res_number == old_start:

                        # We know arrived at the beginning of the TM segment
                        # before it was elongated by the elongate_tm_segments function
                        # If the residues between this one and the previous one 
                        # were all in the margin, the buffer is not empty

                        if buffer:
                            sequence += buffer
                            buffer = ""
                        
                        sequence += pdb_struct["CA"][chain_id][res_number]["res_name"]

                    # we are in the TM segment
                    elif old_start < res_number <= old_end:

                        sequence += pdb_struct["CA"][chain_id][res_number]["res_name"]

                    elif res_number > old_end:

                        # We know left the TM segment  
                        
                        if pdb_struct["CA"][chain_id][res_number]["in_membrane"] == "1":

                            # If we are in the membrane but we already left it, 
                            # it means we are BACK in it : we are in the next TM segment

                            if res_number > old_end + 2:
                                
                                if have_left_membrane:
                                    break

                            else:
                                # We are in the membrane so we are in the margin, by definition
                                sequence += pdb_struct["CA"][chain_id][res_number]["res_name"]

                        elif pdb_struct["CA"][chain_id][res_number]["in_membrane"] == "0":

                            have_left_membrane = True
                            if pdb_struct["CA"][chain_id][res_number]["in_margin"] == "1":

                                sequence += pdb_struct["CA"][chain_id][res_number]["res_name"]

                            elif pdb_struct["CA"][chain_id][res_number]["in_margin"] == "0":

                                margin_plus_one = pdb_struct["CA"][chain_id][res_number+1]["in_margin"]
                                margin_minus_one = pdb_struct["CA"][chain_id][res_number-1]["in_margin"]

                                if margin_plus_one == "1" and margin_minus_one == "1":

                                    sequence += pdb_struct["CA"][chain_id][res_number]["res_name"]

                                else:
                                    break

                        

                else:
                    print(f"Residue {res_number} not found in the PDB file")
                    sequence += "X"

            for res_number in range(old_start, old_end+1):
                
                if res_number in pdb_struct["CA"][chain_id]:
                    sequence_short += pdb_struct["CA"][chain_id][res_number]["res_name"]
                else:
                    sequence_short += "X"
            
            record = SeqRecord.SeqRecord(Seq(sequence), id=f"{protein_name}_{chain_id}_{i+1}", description=f"")
            records.append(record)
            
            record = SeqRecord.SeqRecord(Seq(sequence_short), id=f"{protein_name}_{chain_id}_{i+1}_short", description=f"")
            records_shorts.append(record)

    SeqIO.write(records, output_file, "fasta")
    SeqIO.write(records_shorts, output_file_short, "fasta")

    return 0
                    
                
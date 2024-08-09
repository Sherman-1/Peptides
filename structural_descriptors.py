import numpy as np
from scipy.spatial.distance import pdist, squareform
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
import polars as pl
import glob
import multiprocessing
import warnings
import os
import subprocess
import tempfile



from bin.HullRadV9 import compute_metrics as hullrad
from bin.hydromoment import analyze_sequence

import argparse 

def suppress_warnings(func):
    def wrapper(*args, **kwargs):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return func(*args, **kwargs)
    return wrapper


PDBParser.get_structure = suppress_warnings(PDBParser.get_structure)
parser = PDBParser(QUIET=True)
ppb=PPBuilder()


AA = { "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V" }

TEMP_DIR = "/store/EQUIPES/BIM/MEMBERS/simon.herman/Peptides/tmp"

def extract_matrix(structure):
    """
    Returns distance matrix of a PDB file.
    """

    try:
        full_barycentres = []
        for model in structure:
            for chain in model:     
                coords = []
                for residue in chain:
                    if residue.id[0] == ' ':  # Only consider standard residues (ATOM records)
                        # Compute the barycentre of the residue
                        # Only consider atoms that are part of the side chain
                        barycentre = np.mean([atom.get_coord() for atom in residue.get_unpacked_list()], axis=0)
                        coords.append(barycentre)
                full_barycentres.extend(coords)
            # Only compute the first model
            break

        return squareform(pdist(full_barycentres))
        
    except Exception as e:
        print(f"Error processing {structure.id}: {e}")
        return None





def parse_assignment(lines) -> dict:
    
    """
    Parse an ASG line from STRIDE
    
    STRIDE Format as in official documentation : 
    
    ASG    Detailed secondary structure assignment

    Format:  
    
    6-8     Residue name
    10-10   Protein chain identifier
    12-15   PDB	residue	number
    17-20   Ordinal residue number
    25-25   One	letter secondary structure code	**)
    27-39   Full secondary structure name
    43-49   Phi	angle
    53-59   Psi	angle
    65-69   Residue solvent accessible area

    Returns
    -------
    Assignment namedtuple
    """
    Assignments = {}
    
    for line in lines: 
        
        if line.startswith("ASG"):
            
            resname = line[5:8].strip()
            chain_id = line[9:10].strip()
            res_id = line[11:15].strip()
            res_num = line[16:20].strip()
            structure_code = line[24:25].strip()
            structure_name = line[26:39].strip()
            phi = line[42:49].strip()
            psi = line[52:59].strip()
            area = line[64:69].strip()
            
            Assignments[chain_id] = { 
                                    
                                "resname" : resname,   
                                "res_id" : res_id,
                                "res_num" : res_num,
                                "structure_code" : structure_code,
                                "structure_name" : structure_name,
                                "phi" : phi,
                                "psi" : psi,
                                "area" : area
                                
                            }
            
    if Assignments:
        
        return Assignments
    
    else:
        
        return None

def stride(pdbpath):
  
    """
    Perform STRIDE analysis 
    
    Parameters
    ----------
    pdbpath : str
      path to the pdb file to analyse

    Returns
    -------
    output : decoded standard output of 
    """
    
    stride_path = 'stride'
    p = subprocess.Popen([stride_path, pdbpath],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    
    try:
        stdout, stderr = p.communicate()
        
        lines = stdout.decode('utf-8')
        
        if lines :
            
            return lines.split("\n")
        
        else: 
            
            return None
    
    except Exception as e:
        print(e)
        return None
  


def describe_matrix(distance_matrix, k=5):

    """
    Returns a set of descriptive statistics of a distance matrix.
    """

    n = distance_matrix.shape[0]

    frobenius = np.linalg.norm(distance_matrix, 'fro')
    
    eigvals, _ = np.linalg.eigh(distance_matrix)
    eigens = { f"eig_{i}": eigvals[i] for i in range(k) }

    triu_indices = np.triu_indices_from(distance_matrix, k=1)  # Upper triangular indices
    distances = distance_matrix[triu_indices]

    mean_distance = np.mean(distances)
    median_distance = np.median(distances)
    std_distance = np.std(distances)
    min_distance = np.min(distances)
    max_distance = np.max(distances)
    q25_distance = np.percentile(distances, 25)
    q75_distance = np.percentile(distances, 75)

    return {

        "mean_distance": mean_distance,
        "median_distance": median_distance,
        "std_distance": std_distance,
        "min_distance": min_distance,
        "max_distance": max_distance,
        "q25_distance": q25_distance,
        "q75_distance": q75_distance,
        "frobenius": frobenius,
        **eigens
    }

def compute_structural_metrics(pdb_path):

    basename = os.path.basename(pdb_path).split(".")[0]

    try:
        
        structure = parser.get_structure("protein", pdb_path)
        distance_matrix = extract_matrix(structure)
        
    except Exception as e:
        return None

    seq = []
    for pp in ppb.build_peptides(structure):
        seq.extend(pp.get_sequence())
        
    if len(seq) == 0 or set(seq) <= AA:
        return None
        
    aa_composition = { aa : seq.count(aa) / len(seq) for aa in AA }

    chem = analyze_sequence(sequence = seq, window = -1)

    chem_features = {

        "av_h" : np.mean([d[6] for d in chem]),
        "av_uH" : np.mean([d[7] for d in chem]),
        "d" : np.mean([d[8] for d in chem]),
        "n_tot_pol" : np.mean([d[9] for d in chem]),
        "n_tot_apol" : np.mean([d[10] for d in chem]),
        "n_charged" : np.mean([d[11] for d in chem]),
        "n_aromatic" : np.mean([d[12] for d in chem])

    }

    geometric_features = describe_matrix(distance_matrix)

    phys_features = hullrad(pdb_path)
    
    second_struct_percent = stride(pdb_path)
    
    if phys_features is None or geometric_features is None or chem_features is None or second_struct_percent is None:

        return None

    else: 

        return { **{"id": basename, "second_struct_percent": second_struct_percent}, 
                **geometric_features, **phys_features, **chem_features, **aa_composition}


def pool_process(pdb_paths, num_processes):
    
    with multiprocessing.Pool(processes=num_processes) as pool:
        print(f"Processing {len(pdb_paths)} PDB files using {num_processes} processes")
        results = pool.map(compute_structural_metrics, pdb_paths)

    return [ res for res in results if res is not None ]


def call_pool(pdb_paths, num_processes, category, output):
    
    results = pool_process(pdb_paths, num_processes)

    return pl.DataFrame(results).with_columns(category = pl.lit(category))


    
def main(categories):

    categories = {

        #"poly_shorts": "/store/EQUIPES/BIM/MEMBERS/simon.herman/Peptides/save/pdbs/polytopics/shorts",
        "bi_shorts": "/store/EQUIPES/BIM/MEMBERS/simon.herman/Peptides/save/pdbs/bitopics/shorts",
        "peripherals": "/store/EQUIPES/BIM/MEMBERS/simon.herman/Peptides/save/pdbs/peripherals/longs",
        # "associated": subprocess.check_output(["python3", "horizontals.py"]).decode().strip(),
        "S3": "/store/EQUIPES/BIM/MEMBERS/paul.roginski/OLD/ORFPRED/new_data/SCOPe/classes_a_b_c_d_e_20_to_70_uniq_pdb",
        "small": "/store/EQUIPES/BIM/MEMBERS/paul.roginski/OLD/ORFPRED/new_data/SCOPe/class_g_20_to_70_uniq_pdb",

    }

    dfs = {}
    for category, path in categories.items():
        
        print(f"Processing {category} ...")
        
        pdbs = [file for ext in ('pdb', 'ent') for file in glob.glob(f"{path}/*.{ext}")]   
        
        pdbs = np.random.choice(pdbs, 200, replace = False)
        
        df = call_pool(pdbs, 60, category, f"{category}.tsv")

        dfs[category] = df
        
    return pl.concat(dfs.values())

    
    
if __name__ == "__main__":
    
    df = main()
    
    df.write_csv("structural_descriptors.tsv", separator = "\t", include_header = True)
    













#!/usr/bin/env python3


"""
## type_id correspondance : 

#   1 : Multitopic transmembrane
#   2 : Monotopic / Peripherals
#   3 : Peptides 

## classtype_id correspondance :

###### Multitopic transmembrane ######
#  1  : Alpha-helical polytopic 
#  11 : Bitopic proteins 
#  2  : Beta barrel 

###### Monotopic ######
#  4  : All alpha 
#  3  : All beta
#  5  : Alpha / beta
#  6  : Alpha + beta

###### Peptides #######
#  7  : Alpha-helical peptides
#  8  : Beta-helical peptides
#  9  : Beta hairpins 
#  10 : Non-regular peptides
"""

import polars as pl 
import multiprocessing as mp 
from bin.MMseqs import MMseqs2API
from pathlib import Path, PosixPath
from Bio import SeqIO
from tqdm import tqdm
import os
import shutil
import random


from bin.peripheral import peripheral
from bin.transmembrane import transmembrane

MIN_LENGTH = 20
MAX_LENGTH = 70
MIN_SEGMENT_LENGTH = 15
INNER_MARGIN = 0
CLOSE_MARGIN = 5
MARGIN = 20
GAPS = 1
IORF_CSV = "input/iORFs.csv"


def get_paths():


    metadata = (

        pl.read_csv("input/proteins-2024-05-07.csv", separator = ",", infer_schema_length = 20000)
            .with_columns(
                pl.concat_str([
                        pl.lit("input/OPM"),
                        pl.concat_str([pl.col("pdbid"), pl.lit(".pdb")], separator = "")           
                    ], separator = "/",
                ).alias("pdb_path")
            )
    )

    ## FILTERS ## 


    # Annotated bitopic proteins, put threshold of 20 to avoid some false positives
    bitopic_proteins = ((pl.col("classtype_id") == 11) & (pl.col("thickness") >= 20)) 

    # All peptides that are crossing the membrane ( > 20 ), regardless of their folding type
    bitopic_peptides = ((pl.col("type_id") == 3) & (pl.col("thickness") >= 20)) 
    # Multitopic proteins to be cut 
    polytopic_proteins = (pl.col("classtype_id") == 1) & (pl.col("thickness") >= 20)


    # Some "peripheral" proteins are anchored to the mb through lipidation, empirical threshold of 12 to remove them
    peripheral_proteins = (pl.col("type_id") == 2) & (pl.col("thickness") <= 12) 
        
    # All peptides that are not beta-helical or non-regular and that are not crossing the membrane
    peripheral_peptides = ((pl.col("classtype_id") == 7) & (pl.col("classtype_id") == 9) & (pl.col("thickness") < 20)) 

    horizontal_peripheral_peptides = peripheral_peptides & (pl.col("tilt") > 75)
        
    # Mis-annotated bitopic proteins
    misannotated_proteins = ((pl.col("classtype_id") == 11) & (pl.col("thickness") < 20))

    tm_paths = {
        
        "bitopic_proteins" : metadata.filter(bitopic_proteins)["pdb_path"].to_list(),
        "bitopic_peptides" : metadata.filter(bitopic_peptides)["pdb_path"].to_list(),
        "polytopic_proteins" : metadata.filter(polytopic_proteins)["pdb_path"].to_list()

    }

    peripheral_paths = {

        "peripheral_proteins" : metadata.filter(peripheral_proteins)["pdb_path"].to_list(),
        "peripheral_peptides" : metadata.filter(peripheral_peptides)["pdb_path"].to_list(),
        "horizontal_peripheral_peptides" : metadata.filter(horizontal_peripheral_peptides)["pdb_path"].to_list(),
        "misannotated_proteins" : metadata.filter(misannotated_proteins)["pdb_path"].to_list()

    }

    return tm_paths, peripheral_paths

def write_pdb(structure: dict, id : str, writing_dir : PosixPath) -> int:

    lines = []
    try:
        for res_number, atom_dict in structure.items():
            for atom_number, atom_line in atom_dict.items():
                lines.append(atom_line)

        if lines:
            file_name = f"{id}.pdb"
            file_path = f"{writing_dir}/{file_name}"

            with open(file_path, "w") as output:
                output.write("".join(lines)) # Lines are already \n terminated

        else:
            
            return ValueError("Empty structure")

    except Exception as e:
        
        return e


    return 0

def try_tm(args : list):

    try:
        return transmembrane(*args)
    except Exception as e:
        return None
    
def try_peripheral(args : list):

    try:
        return peripheral(*args)
    except Exception as e:
        return None
    
def process_results(results):

    sequences = [ seq for res in results if res is not None for seq in res["sequences"] ]

    structures = {

        k2: v2
        for res in results if res is not None
        for _dict in res["structures"].values()
        for k2, v2 in _dict.items()
    }

    return sequences, structures
    
def main():

    random.seed(42)

    transmembranes, peripherals = get_paths()

    mmseqs = MMseqs2API(threads = max(1, os.cpu_count() - 2), cleanup=True)

    database = Path("database")
    database.mkdir(exist_ok = True)

    # RM everyting under database, dir or file
    for path in database.iterdir():

        try:
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                shutil.rmtree(path)
        except Exception as e:
            print(f"Error deleting file {path}: {e}")

    pbar = tqdm(transmembranes.items(), leave = False, total = len(transmembranes))

    covs = [0.5, 0.7, 0.9]
    idens = [0.5, 0.7, 0.9]

    for category, paths in pbar:

        pbar.set_description(f"Processing {category}")

        with mp.Pool(max(1, mp.cpu_count() - 4)) as pool:

            results = pool.map(try_tm, [(path, None, MARGIN, INNER_MARGIN, MIN_LENGTH, MAX_LENGTH, GAPS, IORF_CSV, None, False) for path in paths])

        results = [ res for res in results if res is not None ]

        sequences, structures = process_results(results) 

        SeqIO.write(sequences, f"tmp/{category}.fasta", "fasta")

        cat_path = database / Path(category)
        cat_path.mkdir(exist_ok = True)

        for cov in covs:

            for id in idens:

                representatives = mmseqs.fasta2representativeseq(fasta_file = f"tmp/{category}.fasta", writing_dir = "tmp", cov = cov, iden = id, cov_mode = 0)

                param_paths = cat_path / Path(f"cov_{cov}_iden_{id}")
                param_paths.mkdir(exist_ok = True)

                pdb_path = param_paths / Path("pdb")
                pdb_path.mkdir(exist_ok = True)

                fasta_path = param_paths / Path("fasta")
                fasta_path.mkdir(exist_ok = True)

                common_keys = set(representatives.keys()) & set(structures.keys())

                representatives = { k : representatives[k] for k in common_keys }
                structures = { k : structures[k] for k in common_keys }

                print(f"Number of representatives : {len(representatives)}")
                print(f"Number of structures : {len(structures)}")

                seq_to_write = []
                wrong_keys = []
                for key in common_keys:

                    returned = write_pdb(structures[key], key, pdb_path)

                    if returned == 0:

                        seq_to_write.append(representatives[key])
                    
                    else:

                        print(f"Error writing pdb {key} : {returned}")
                        wrong_keys.append(key)

                SeqIO.write(seq_to_write, f"{fasta_path}/{category}.fasta", "fasta")

    pbar.close()

if __name__ == '__main__':

    main()

        
        


    

        
    
    


    
        

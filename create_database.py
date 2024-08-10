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
from bin.peripheral import peripheral
from bin.transmembrane import transmembrane
from pathlib import Path
from Bio import SeqIO
from tqdm import tqdm
import os

import multiprocessing

MARGIN = 50
MIN_LENGTH = 20
MAX_LENGTH = 70
MIN_SEGMENT_LENGTH = 15
INNER_MARGIN = 0
GAPS = 1
iORF_csv = "input/iORFs.csv"


def write_pdb_files(structures: dict, category : str, length : str) -> int:

    curr_dir = os.getcwd()

    writing_dir = Path(os.path.join(curr_dir, "pdbs", category, length))

    writing_dir.mkdir(parents=True, exist_ok=True)

    for protein_name, struct_dict in structures.items():
        for chain_id, segment_dict in struct_dict.items():
            if not segment_dict:
                continue
            for segment_id, residue_dict in segment_dict.items():
                if not residue_dict:
                    continue    
                lines = []
                try:
                    for res_number, atom_dict in residue_dict.items():
                        for atom_number, atom_line in atom_dict.items():
                            lines.append(atom_line)

                except Exception as e:
                    print(f"Error writing pdb {protein_name}_{chain_id}_{segment_id}: {e}")

                if lines:
                    file_name = f"{protein_name}_{chain_id}_{segment_id}_{length}.pdb"
                    file_path = f"{writing_dir}/{file_name}"

                    with open(file_path, "w") as output:
                        output.write("".join(lines)) # Lines are already \n terminated

    return 0


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


# Some "peripheral" proteins are anchored to the mb through lipidation, empirical threshold of 12, we don't want them
# Although we might miss through candidates between 12 and 20 
peripheral_proteins = (pl.col("type_id") == 2) & (pl.col("thickness") <= 12) 
    
# All peptides that are not beta-helical or non-regular and that are not crossing the membrane
peripheral_peptides = ((pl.col("classtype_id") == 7) & (pl.col("classtype_id") == 9) & (pl.col("thickness") < 20)) 

horizontal_peripheral_peptides = peripheral_peptides & (pl.col("tilt") > 75)
    
# Mis-annotated bitopic proteins
misannotated_proteins = ((pl.col("classtype_id") == 11) & (pl.col("thickness") < 20))


## PATHS ##

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


if __name__ == "__main__":  

    n_cpus = multiprocessing.cpu_count() - 1

    pool = multiprocessing.Pool(processes=n_cpus)

    for category, paths in tqdm(tm_paths.items(), desc = "transmembrane proteins", leave = False):
        
        print(f'Processing {category} category')
        
        args = [(path, MARGIN, INNER_MARGIN, MIN_LENGTH, MAX_LENGTH, MIN_SEGMENT_LENGTH, GAPS, iORF_csv, iORF_csv, False) for path in paths]
        
        results = pool.starmap(transmembrane, args)
        
    
    


    
        

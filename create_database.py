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

    pl.read_csv("/store/EQUIPES/BIM/MEMBERS/simon.herman/Peptides/input/proteins-2024-05-07.csv", separator = ",", infer_schema_length = 20000)
    .with_columns(
        pl.concat_str([
                pl.lit("/store/EQUIPES/BIM/MEMBERS/simon.herman/Peptides/input/OPM_database"),
                pl.concat_str([pl.col("pdbid"), pl.lit(".pdb")], separator = "")           
            ], separator = "/",
        ).alias("pdb_path")
    )
)

## FILTERS                                                     

# Bitopics 
bitopics = (

    # Annotated bitopic proteins, put threshold of 20 to avoid some false positives
    ((pl.col("classtype_id") == 11) & (pl.col("thickness") >= 20)) | 
    # All peptides that are crossing the membrane ( > 20 ), regardless of their folding type
    ((pl.col("type_id") == 3) & (pl.col("thickness") >= 20)) |
    # Multitopic proteins to be cut 
    (pl.col("classtype_id") == 1) & (pl.col("thickness") >= 20)

) 

# Peripherals 
peripherals = (

    # Some "peripheral" proteins are anchored to the mb through lipidation, empirical threshold of 12, we don't want them
    # Although we might miss through candidates between 12 and 20 
    (pl.col("type_id") == 2) & (pl.col("thickness") <= 12) |
    # All peptides that are not beta-helical or non-regular and that are not crossing the membrane
    ((pl.col("classtype_id") == 7) & (pl.col("classtype_id") == 9) & (pl.col("thickness") < 20)) |
    # Mis-annotated bitopic proteins
    ((pl.col("classtype_id") == 11) & (pl.col("thickness") < 20))

)

peripherals_path = metadata.filter(peripherals).select("pdb_path").to_series().to_list()[:100]
bitopics_path = metadata.filter(bitopics).select("pdb_path").to_series().to_list()

if __name__ == "__main__":
        
    nb_extracted_segments = []
    for inner_margin in tqdm(range(1,31), desc = "Margin", leave = False):
        
        MARGIN = inner_margin + 10
        
        counter = 0

        periphs_args = [(path, MARGIN, inner_margin, MIN_LENGTH, MAX_LENGTH, MIN_SEGMENT_LENGTH, iORF_csv, iORF_csv, GAPS) for path in peripherals_path]

        with multiprocessing.Pool(processes = 30) as pool:

            results = [ res for res in pool.starmap(peripheral, periphs_args) if type(res) == dict ]
            
            
        counter += sum([len(res) for res in results])
        
        nb_extracted_segments.append(counter)
    
    
    print(nb_extracted_segments)
        
    


    
        

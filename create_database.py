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
from bin.MMseqs import MMseqs2API
from pathlib import Path
from Bio import SeqIO
from tqdm import tqdm
import os
import shutil


MIN_LENGTH = 20
MAX_LENGTH = 70
MIN_SEGMENT_LENGTH = 15
INNER_MARGIN = 0
CLOSE_MARGIN = 5
MARGIN = 20
GAPS = 1
IORF_CSV = "input/iORFs.csv"


def write_pdb_files(structures: dict, length : str, writing_dir : str) -> int:

   
    for segment_id, residue_dict in structures.items():
        
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
            file_name = f"{segment_id}_{length}.pdb"
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


## SETUP RESULTS DIRECTORY ##

results_dir = Path("results")
results_dir.mkdir(parents=True, exist_ok=True)

# If there are subdirectories, remove them

for sub_dir in results_dir.iterdir():
    
    shutil.rmtree(sub_dir)


## TRANSMEMBRANE ##

pbar = tqdm(tm_paths.items(), total = len(tm_paths), leave = False)

structures = {}
sequences = {}
for category, pdb_paths in pbar:
    
    pbar.set_description(f"Processing {category}")
    
    sub_pbar = tqdm(pdb_paths, leave = False, total = len(pdb_paths))
    
    i = 0
    
    structures[category] = {}
    sequences[category] = []

        
    for pdb_path in sub_pbar:
        
        try:
        
            sub_pbar.set_description(f"Processing {pdb_path}")
            
            protein_name = os.path.basename(pdb_path).split(".")[0]
            
            # transmembrane(file_path, secondary_structure_path, margin, inner_margin, min_length, max_length, gaps, iorf_path, csv_path, verbose = False ):
            
            result = transmembrane(pdb_path, None, MARGIN, INNER_MARGIN, MIN_LENGTH, MAX_LENGTH, GAPS, IORF_CSV, False)
            
            for chain_id, segment_dict in result["structures"].items():
                for segment_id, residue_dict in segment_dict.items():
                    structures[category].update({segment_id : residue_dict})
            
            sequences[category].extend(result["records"])
            
            i += 1  
            
            if i == 1:
                break # Debug
                   
        except Exception as e:
            
            print(f"Error processing {pdb_path}: {e}")
            continue
        
    SeqIO.write(sequences[category], f"tmp/{category}.fasta", "fasta")

pbar.close()

## PERIPHERAL ##

pbar = tqdm(peripheral_paths.items(), total = len(peripheral_paths), leave = False)

for category, pdb_paths in pbar:
    
    pbar.set_description(f"Processing {category}")
    
    sub_pbar = tqdm(pdb_paths, leave = False, total = len(pdb_paths))
    
    structures[category] = {}
    sequences[category] = []
    
    i = 0
    
    for pdb_path in sub_pbar:
        
        try:
        
            sub_pbar.set_description(f"Processing {pdb_path}")
            
            protein_name = os.path.basename(pdb_path).split(".")[0]
            
            # peripheral(pdb_path, close_margin, outer_margin, min_length, max_length, min_segment_length, iorf_csv, iorf_fasta, gaps, verbose = False):
            
            result = peripheral(pdb_path, CLOSE_MARGIN, MARGIN, MIN_LENGTH, MAX_LENGTH, MIN_SEGMENT_LENGTH, IORF_CSV, None, GAPS, False)
                        
            for chain_id, segment_dict in result["structures"].items():
                for segment_id, residue_dict in segment_dict.items():
                    structures[category].update({segment_id : residue_dict})
            
            sequences[category].extend(result["records"])
            
            i += 1
            
            if i == 100:
                break # Debug
            
        except Exception as e:
            
            print(f"Error processing {pdb_path}: {e}")
            continue
        
    SeqIO.write(sequences[category], f"tmp/{category}.fasta", "fasta")

pbar.close()

coverage = [0.5, 0.7, 0.9]
identity = [0.5, 0.7, 0.9]

mmseqs = MMseqs2API(threads = max(os.cpu_count() - 2, 1), cleanup = True)

for cov in coverage: 
    
    for id in identity: 
        
        for category in sequences.keys():
            
            try:
                
                mmseqs.fasta2representativeseq(f"tmp/{category}.fasta", f"results/{category}_cov_{cov}_id_{id}", cov, id)
                
            except Exception as e:
                
                print(f"Error creating representatives for {category}: {e}")
                continue
            
            break   
        
        break
    
    break
        
        
        
        


    

        
    
    


    
        

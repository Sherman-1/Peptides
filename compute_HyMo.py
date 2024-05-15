from hydromoment import read_fasta_file, analyze_sequence
from utils import read_pdb, ResidueError, pdb_struct_to_fasta

import polars as pl 

# Column typeid wrongly formated : ="1qjp" instead of 1qjp
metadata = (
    pl.scan_csv("proteins-2024-05-07.csv", has_header = True)
    .with_columns(
        pl.col("pdbid").str.replace(pattern  = r'="(\w+)"', value = "${1}")
    )
)

## type_id correspondance : 

#   1 : Multitopic transmembrane
#   2 : Monotopic / Peripherals 
#   3 : Peptides ( either TM or peripherals )


## classtype_id correspondance : 

#   8 : Beta-helical peptides ( we do not want them ) 


peptides = (
    metadata
    .filter(
        (pl.col("type_id") == 3) & 
        (pl.col("classtype_id") != 8)
    )
    .collect()
)


for pdb_id in peptides.select(pl.col("pdbid")).to_series():

    pdb_path = f"peptides/{pdb_id}.pdb"
    wrong_pdbs = 0

    try:
        pdb = read_pdb(pdb_path)
        pdb_struct_to_fasta(pdb_path)
    except ResidueError as e:
        print(f"Error reading {pdb_id}: {e}")
        wrong_pdbs += 1
        continue
    except FileNotFoundError:
        print(f"File not found: {pdb_id}")
        print(peptides.filter(pl.col("pdbid") == pdb_id))
        wrong_pdbs += 1
        continue
    
print(f"Total wrong PDBs: {wrong_pdbs}")
    


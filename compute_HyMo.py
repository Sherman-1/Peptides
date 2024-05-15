from hydromoment import read_fasta_file, analyze_sequence
from utils import read_pdb, ResidueError, pdb_struct_to_fasta, test

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

test()
    


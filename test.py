import glob
from bin.data.utils import read_pdb
from pathlib import Path

files = glob.glob("database/horizontal/cov_0.3_iden_0.3/pdb/*")

i = 0
for pdb_file in files:


    path = Path(pdb_file)
    pdb = read_pdb(None,path)

    
    i += 1
    print(pdb["protein_name"])
    print(pdb["sequence"])
    
    print(i)



        
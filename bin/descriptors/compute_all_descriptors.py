from pTrans_embeddings import main as pTrans_embeddings
from sequence_descriptors import process_data

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from pathlib import Path
import os
import glob


def search_main_directory():

    """
    Iteratively searches for the peptide directory and returns its absolute path.
    """

    global main_directory
    main_directory = None
    for i in range(3): 
        backward = "../" * i
        main_directory = Path(f"{backward}Peptides").resolve()
        if main_directory.exists():
            break

    if main_directory is None:
        raise FileNotFoundError("Peptide directory not found")
    
    print(f"Working on the main directory : {main_directory}")
    
    return main_directory

def main():

    search_main_directory()

    writting_dir = main_directory / "database" / "descriptors"
    writting_dir.mkdir(exist_ok = True)

    print("Bonjour")

    for fasta in glob.glob(f"{main_directory}/database/full/*.fasta"):

        print(f"Working on {fasta}")
        
        records = list(SeqIO.parse(fasta, "fasta"))

        category = Path(fasta).stem

        if category != "horizontal":
            continue

        seq_based_descriptors = process_data(records, category) 

        embeddings = pTrans_embeddings(records)

        seq_based_descriptors.join(embeddings, how = "inner", on = "id").write_parquet(f"{writting_dir}/{category}.parquet")

        

if __name__ == "__main__":

    main()
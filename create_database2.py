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
from Bio import SeqIO, SeqRecord
from tqdm import tqdm
import os
import shutil
import random
from glob import glob
import logging

from bin.peripheral import peripheral
from bin.transmembrane import transmembrane

# Constants
MIN_LENGTH = 20
MAX_LENGTH = 70
MIN_SEGMENT_LENGTH = 10
INNER_MARGIN = 0
CLOSE_MARGIN = 5
MARGIN = 20
GAPS = 1
IORF_CSV = "input/iORFs.csv"
COVERAGE_THRESHOLDS = [0.3, 0.5, 0.7]
IDENTITY_THRESHOLDS = [0.3, 0.5, 0.7]


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
    peripheral_peptides = ((pl.col("type_id") == 3) & (pl.col("thickness") < 20) & (pl.col("tilt") <= 80)) 

    horizontal_peripheral_peptides = ((pl.col("type_id") == 3) & (pl.col("thickness") < 20) & (pl.col("tilt") > 80))
        
    # Mis-annotated bitopic proteins
    misannotated_proteins = ((pl.col("classtype_id") == 11) & (pl.col("thickness") < 20))


    membranome = list(glob("input/membranome/*.pdb"))
    bitopic_proteins = metadata.filter(bitopic_proteins)["pdb_path"].to_list()
    bitopic_peptides = metadata.filter(bitopic_peptides)["pdb_path"].to_list()


    tm_paths = {
        
        "bitopic" : membranome + bitopic_proteins + bitopic_peptides,
        "polytopic" : metadata.filter(polytopic_proteins)["pdb_path"].to_list()

    }

    peripheral_proteins = metadata.filter(peripheral_proteins)["pdb_path"].to_list()
    peripheral_peptides = metadata.filter(peripheral_peptides)["pdb_path"].to_list()
    misannotated_proteins = metadata.filter(misannotated_proteins)["pdb_path"].to_list()

    peripheral_paths = {

        "peripheral" : peripheral_proteins + peripheral_peptides + misannotated_proteins,
        "horizontal" : metadata.filter(horizontal_peripheral_peptides)["pdb_path"].to_list(),

    }

    return tm_paths, peripheral_paths

def write_pdb(structure: dict, id: str, writing_dir: PosixPath) -> int:
    lines = []
    try:
        for res_number, atom_dict in structure.items():
            if not isinstance(atom_dict, dict):
                return 1 
            for atom_number, atom_line in atom_dict.items():
                lines.append(atom_line)
        if lines:
            file_name = f"{id}.pdb"
            file_path = writing_dir / file_name
            with open(file_path, "w") as output:
                output.write("".join(lines))  # Lines are already \n terminated
        else:
            logging.error(f"Empty structure for ID: {id}")
            return 2  
    except Exception as e:
        logging.error(f"Error writing PDB for ID: {id}, Error: {e}")
        return 3  
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

def setup_environment():
    """Sets up the random seed, logging, and database directory."""
    random.seed(42)
    logging.basicConfig(level=logging.DEBUG, format='%(levelname)s - %(message)s')
    create_database_directory()

def create_database_directory():
    """Creates the database directory and removes any existing contents."""
    database = Path("database")
    database.mkdir(exist_ok=True)
    for path in database.iterdir():
        try:
            if path.is_file():
                path.unlink()
            elif path.is_dir():
                shutil.rmtree(path)
        except Exception as e:
            logging.error(f"Error deleting file {path}: {e}")

def get_paths():
    """Fetches and filters paths for transmembrane and peripheral proteins."""
    metadata = pl.read_csv("input/proteins-2024-05-07.csv", separator=",", infer_schema_length=20000)
    metadata = metadata.with_columns(
        pl.concat_str([
            pl.lit("input/OPM"),
            pl.concat_str([pl.col("pdbid"), pl.lit(".pdb")], separator="")
        ], separator="/").alias("pdb_path")
    )

    bitopic_proteins = (pl.col("classtype_id") == 11) & (pl.col("thickness") >= 20)
    bitopic_peptides = (pl.col("type_id") == 3) & (pl.col("thickness") >= 20)
    polytopic_proteins = (pl.col("classtype_id") == 1) & (pl.col("thickness") >= 20)
    peripheral_proteins = (pl.col("type_id") == 2) & (pl.col("thickness") <= 12)
    peripheral_peptides = (pl.col("type_id") == 3) & (pl.col("thickness") < 20) & (pl.col("tilt") <= 80)
    horizontal_peripheral_peptides = (pl.col("type_id") == 3) & (pl.col("thickness") < 20) & (pl.col("tilt") > 80)
    misannotated_proteins = (pl.col("classtype_id") == 11) & (pl.col("thickness") < 20)

    membranome = list(glob("input/membranome/*.pdb"))
    tm_paths = {
        "bitopic": membranome + metadata.filter(bitopic_proteins)["pdb_path"].to_list() + metadata.filter(bitopic_peptides)["pdb_path"].to_list(),
        "polytopic": metadata.filter(polytopic_proteins)["pdb_path"].to_list(),
    }

    peripheral_paths = {
        "peripheral": metadata.filter(peripheral_proteins)["pdb_path"].to_list() + metadata.filter(peripheral_peptides)["pdb_path"].to_list() + metadata.filter(misannotated_proteins)["pdb_path"].to_list(),
        "horizontal": metadata.filter(horizontal_peripheral_peptides)["pdb_path"].to_list(),
    }

    return tm_paths, peripheral_paths

def process_results(results):
    """Processes the results from the multiprocessing pool."""
    sequences = [seq for res in results if res is not None for seq in res["sequences"]]
    structures = {
        k2: v2
        for res in results if res is not None
        for _dict in res["structures"].values()
        for k2, v2 in _dict.items()
    }
    return sequences, structures

def process_category(category, paths, process_func, mmseqs):
    """Processes a single category of transmembranes or peripherals."""
    results = run_multiprocessing(paths, process_func)
    sequences, structures = process_results(results)
    write_sequences(category, sequences)
    cat_path = setup_category_directory(category)
    for cov in COVERAGE_THRESHOLDS:
        for iden in IDENTITY_THRESHOLDS:
            representatives = mmseqs.fasta2representativeseq(fasta_file=f"tmp/{category}.fasta", writing_dir="tmp", cov=cov, iden=iden, cov_mode=0)
            handle_representatives(category, representatives, structures, cov, iden, cat_path)

def run_multiprocessing(paths, process_func):
    """Runs multiprocessing on the given paths with the provided processing function."""
    with mp.Pool(max(1, mp.cpu_count() - 4)) as pool:
        return pool.map(process_func, paths)

def write_sequences(category, sequences):
    """Writes sequences to a FASTA file."""
    SeqIO.write(sequences, f"tmp/{category}.fasta", "fasta")

def setup_category_directory(category):
    """Creates the necessary directories for the given category."""
    cat_path = Path("database") / Path(category)
    cat_path.mkdir(exist_ok=True)
    return cat_path

def handle_representatives(category, representatives, structures, cov, iden, cat_path):
    """Handles representative sequences and writes PDB files."""
    param_paths, pdb_path, fasta_path = create_param_directories(cat_path, cov, iden)
    common_keys = set(representatives.keys()) & set(structures.keys())
    write_pdbs(common_keys, structures, pdb_path)
    write_fasta_for_common_keys(category, common_keys, representatives, fasta_path, cov, iden)

def create_param_directories(cat_path, cov, iden):
    """Creates directories for specific coverage and identity parameters."""
    param_paths = cat_path / Path(f"cov_{cov}_iden_{iden}")
    param_paths.mkdir(exist_ok=True)
    pdb_path = param_paths / Path("pdb")
    pdb_path.mkdir(exist_ok=True)
    fasta_path = param_paths / Path("fasta")
    fasta_path.mkdir(exist_ok=True)
    return param_paths, pdb_path, fasta_path

def write_pdbs(common_keys, structures, pdb_path):
    """Writes PDB files for the common keys."""
    for key in common_keys:
        returned = write_pdb(structures[key], key, pdb_path)
        if returned != 0:
            logging.error(f"Error writing PDB for key {key}")

def write_fasta_for_common_keys(category, common_keys, representatives, fasta_path, cov, iden):
    """Writes FASTA sequences for the common keys."""
    seq_to_write = []
    for key in common_keys:
        pdb_file_path = f"{fasta_path}/../pdb/{key}.pdb"
        if validate_pdb_file(pdb_file_path):
            if key in representatives:
                seq_to_write.append(representatives[key])
            else:
                logging.error(f"Key {key} not found in representatives")
    SeqIO.write(seq_to_write, f"{fasta_path}/{category}_cov_{cov}_iden_{iden}.fasta", "fasta")

def validate_pdb_file(pdb_file):
    """Validates the existence and size of the PDB file."""
    try:
        if os.path.getsize(pdb_file) != 0:
            return True
        else:
            logging.error(f"Empty PDB file: {pdb_file}")
    except FileNotFoundError:
        logging.error(f"File not found: {pdb_file}")
    except IOError as e:
        logging.error(f"IO error reading {pdb_file}: {e}")
    return False

def main():
    setup_environment()
    transmembranes, peripherals = get_paths()
    mmseqs = MMseqs2API(threads=max(1, os.cpu_count() - 2), cleanup=True)

    process_data("Transmembranes", transmembranes, try_tm, mmseqs)
    process_data("Peripherals", peripherals, try_peripheral, mmseqs)

def process_data(data_type, data_paths, process_func, mmseqs):
    """Processes transmembrane or peripheral data."""
    pbar = tqdm(data_paths.items(), leave=False, total=len(data_paths))
    for category, paths in pbar:
        pbar.set_description(f"Processing {category} ({data_type})")
        process_category(category, paths, process_func, mmseqs)
    pbar.close()

if __name__ == '__main__':
    main()

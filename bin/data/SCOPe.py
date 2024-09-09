import requests
import tarfile
from io import BytesIO
import os
import shutil
from Bio import SeqIO, SeqRecord
from pathlib import Path, PosixPath
from bin.MMseqs import MMseqs2

MIN_LENGTH = 20
MAX_LENGTH = 70
COVS = [0.3]
IDENS = [0.3]

def create_directories(cat_path : PosixPath, cov, id):

    """Creates the directories for the category, coverage, and identity."""

    param_paths = cat_path / Path(f"cov_{cov}_iden_{id}")
    param_paths.mkdir(exist_ok=True)
    pdb_path = param_paths / Path("pdb")
    pdb_path.mkdir(exist_ok=True)
    fasta_path = param_paths / Path("fasta")
    fasta_path.mkdir(exist_ok=True)


    return pdb_path, fasta_path

def download_and_extract_tarball(url):
    
    print(f"Downloading {url} ... ")
    response = requests.get(url)
    print("Done")
    tar = tarfile.open(fileobj=BytesIO(response.content), mode="r:gz")
    extracted_files = {}
    
    print("Extracting data")
    for member in tar.getmembers():
        if member.isfile():
            file_content = tar.extractfile(member).read().decode('utf-8')
            extracted_files[member.name] = file_content
    
    return extracted_files

def write_pdbstyle_data_to_disk(pdbstyle_data, output_dir="tmp"):

    os.makedirs(output_dir, exist_ok=True)

    for name, content in pdbstyle_data.items():

        file_path = os.path.join(output_dir, name)
        os.makedirs(os.path.dirname(file_path), exist_ok=True)
        with open(file_path, 'w') as file:
            file.write(content)
    
    print(f"All files have been written to {output_dir}")

def download_pdbstyle(write = False):

    pdbstyle_tree = []
    
    try:
        for i in range(1, 9):
    
            print(f"Processing pdbstyle-2.08-{i}")
            url = f"https://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-2.08-{i}.tgz"
            extracted_files = download_and_extract_tarball(url)
            pdbstyle_tree.extend([name for name in extracted_files.keys() if name.endswith('.ent')])

            if write:
                print(f"Writing pdbstyle-2.08-{i} to disk")
                write_pdbstyle_data_to_disk(extracted_files, output_dir=f"tmp/pdbstyle-2.08-{i}")
                print("Writing pdbstyle_tree in tmp")
                with open("tmp/pdbstyle_tree.txt", "w") as f:
                    for line in pdbstyle_tree:
                        f.write(f"{line}\n")
    
    except requests.exceptions.HTTPError as e:
        print(f"Error downloading file : {e}")
    except Exception as e:
        print(f"Error: {e}")
    
    return pdbstyle_tree

def download_fasta_file(url, output_filename):

    print("Downloading fasta file ... ")
    response = requests.get(url)
    with open(output_filename, 'w') as fasta_file:
        fasta_file.write(response.text)


def extract_class_sequences(fasta_path, min_length, max_length):

    classes = {}
    AAs = "ACDEFGHIKLMNPQRSTVWY"

    for record in SeqIO.parse(fasta_path, "fasta"):
        # Check if the first character of the second word in the description is 'g'
        # Format : d6iyia_ a.1.1.0 (A:) automated matches {Acipenser stellatus [TaxId: 7903]}
        class_ = record.description.split()[1][0]   
        seq = record.seq.upper()

        if set(seq).issubset(AAs) and min_length <= len(seq) <= max_length:
            if class_ not in classes:
                classes[class_] = {}
            classes[class_][record.id] = record

    return classes

def main():

    #pdbstyle_tree = download_pdbstyle(write = True)

    with open("pdbstyle_tree.txt", "r") as file:
        
        pdbstyle_tree = [line.strip() for line in file.readlines()]

    id_to_path = { 
        os.path.basename(path).split('.')[0] : Path(path)
        for path in pdbstyle_tree
    }

    print(id_to_path)

    fasta_url = "https://scop.berkeley.edu/downloads/scopeseq-2.08/astral-scopedom-seqres-gd-all-2.08-stable.fa"
    fasta_path = "tmp/astral.fasta"

    #download_fasta_file(fasta_url, fasta_path)

    records_per_class = extract_class_sequences(fasta_path, MIN_LENGTH, MAX_LENGTH)

    mmseqs = MMseqs2(threads = max(1, os.cpu_count() - 4), cleanup = True)

    index_to_name = {"g" : "Small", "abcde" : "Structured_solution"}

    for indexes in ("g","abcde"):

        category = index_to_name[indexes]

        db_dir = Path(f"database/{category}")
        db_dir.mkdir(parents = True, exist_ok = True)
        
        records = [ record for index in indexes if index in records_per_class for record in records_per_class[index].values() ]

        SeqIO.write(records, f"tmp/{category}.fasta","fasta")

        for cov in COVS: 

            for id in IDENS:
    
                repr = mmseqs.fasta2representativeseq(f"tmp/{category}.fasta",db_dir, cov, id)

                pdb_dir, fasta_dir = create_directories(db_dir, cov = cov, id = id)

                sequences = []
                for id in repr.keys(): 

                    try:

                        pdb_path = id_to_path[id] 
                        destination_dir = pdb_dir / pdb_path.name
                        seq = repr[id]
                        assert type(seq) == SeqRecord
                        shutil.copy(pdb_path, destination_dir )
                        sequences.append(seq)
                    
                    except Exception as e:
                        print(e)

                SeqIO.write(sequences, fasta_dir, "fasta")


                        










            


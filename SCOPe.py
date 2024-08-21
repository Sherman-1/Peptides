import requests
import tarfile
from io import BytesIO
import re
from collections import Counter
from Bio import SeqIO


MIN_LENGTH = 20
MAX_LENGTH = 70

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

def download_and_process_pdbstyle():


    pdbstyle_data = {}
    pdbstyle_tree = []

    for i in range(1, 2):
        print(f"Processing pdbstyle-2.08-{i}")
        url = f"https://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-2.08-{i}.tgz"
        extracted_files = download_and_extract_tarball(url)
        pdbstyle_data.update(extracted_files)
        
        pdbstyle_tree.extend([name for name in extracted_files.keys() if name.endswith('.ent')])
    
    return pdbstyle_data, pdbstyle_tree

def download_fasta_file(url, output_filename):

    response = requests.get(url)
    with open(output_filename, 'w') as fasta_file:
        fasta_file.write(response.text)


# Unused
def count_classes_in_fasta(path):

    class_counts = Counter()
    for record in SeqIO.parse(path, "fasta"):

        # Assuming the class is represented by the first character of the second word in the description
        first_char = record.description.split()[1][0]
        class_counts[first_char] += 1
    
    return class_counts

def extract_class_sequences(fasta_path, _class):

    sequences = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        # Check if the first character of the second word in the description is 'g'
        first_char = record.description.split()[1][0]   
        if first_char == _class:

            sequences[record.id] = str(record.seq)
    
    return sequences

def calculate_sequence_sizes(sequences):
    sizes = {header: len(sequence) for header, sequence in sequences.items()}
    return sizes

def filter_sequences_by_size(sequences, min_length, max_length):

    filtered_sequences = {
        header: seq for header, seq in sequences.items()
        if min_length <= len(seq) <= max_length
    }
    return filtered_sequences

def main():

    pdbstyle_data, pdbstyle_tree = download_and_process_pdbstyle()

    fasta_url = "https://scop.berkeley.edu/downloads/scopeseq-2.08/astral-scopedom-seqres-gd-all-2.08-stable.fa"
    fasta_path = "tmp/astral.fasta"
    download_fasta_file(fasta_url, fasta_path)

    class_g_sequences = extract_class_sequences(fasta_path, "g")

    filtered_sequences = filter_sequences_by_size(class_g_sequences)

if __name__ == "__main__":
    main()

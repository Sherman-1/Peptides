#!/usr/bin/env python3

import argparse
import logging
import random
from Bio import SeqIO
import sys, os
from tqdm import tqdm
from collections import defaultdict
from pathlib import Path
 

from utils import read_pdb, binarize_transmembrane, define_tm_segments, elongate_tm_segments, extract_elongated_sequences_v2
from utils import setup_logger, exception_catcher, ascii_separator, format_pdb_line, generate_pdb


transmembrane_logger = setup_logger('transmembrane_logger', 'transmembrane.log')

@exception_catcher(logging.getLogger('transmembrane_logger'))
def transmembrane(file_path, secondary_structure_path, margin, inner_margin, min_length, max_length, gaps, iorf_path, csv_path ):

    random.seed(random.randint(0,1000))

    transmembrane_logger.info("Starting transmembrane")
    
    pdb_struct = read_pdb(logger = transmembrane_logger,file_path = file_path, secondary_structure_path = secondary_structure_path)
    protein_name = pdb_struct["protein_name"]

    transmembrane_logger.info(f"{ascii_separator}")

    in_membrane_binaries, _ = binarize_transmembrane(transmembrane_logger, pdb_struct, margin, inner_margin)

    transmembrane_logger.info(f"{ascii_separator}")

    tm_indices = define_tm_segments(transmembrane_logger,in_membrane_binaries, pdb_struct)

    transmembrane_logger.info(f"{ascii_separator}")

    elongate_tm_segments(transmembrane_logger, tm_indices, pdb_struct, iorf_path, csv_path, min_length, max_length)

    transmembrane_logger.info(f"{ascii_separator}")

    res_dict = extract_elongated_sequences_v2(transmembrane_logger,tm_indices = tm_indices, pdb_struct = pdb_struct, gaps = gaps)
       
    return res_dict


def process_transmembrane_file(pdb_path: str, secondary_structure_path: str, margin: int, inner_margin: int, 
                     min_length: int, max_length: int, gaps: int, iorf_path: str, csv_path: str) -> None:

    protein_name = os.path.basename(pdb_path).split(".")[0]
    sequences = []
    sequences_short = []

    try:
        res_dict = transmembrane(file_path=pdb_path, 
                                 secondary_structure_path=secondary_structure_path,
                                 margin=margin, 
                                 inner_margin=inner_margin,
                                 min_length=min_length, 
                                 max_length=max_length, 
                                 gaps=gaps,
                                 iorf_path=iorf_path,
                                 csv_path=csv_path)
        
        sequences.extend(res_dict["records"])
        sequences_short.extend(res_dict["records_shorts"])
        structure = res_dict["structures"]
        structure_short = res_dict["structures_shorts"]

    except Exception as e:
        transmembrane_logger.error(f"Error processing file {pdb_path}: {e}")
        #print(f"Error processing file {pdb_path}: {e}")
        return 1

    return sequences, sequences_short, structure, structure_short, protein_name

def main(): 

    parser = argparse.ArgumentParser()

    i_inputs = parser.add_mutually_exclusive_group(required=True)

    i_inputs.add_argument("--input", help="Input file", default = None)
    i_inputs.add_argument("--input_files", nargs='+', help="Input files", default = None)

    parser.add_argument("--secondary_structures", help="Secondary structure file")
    parser.add_argument("--min_length", type=int, help = "Mininimum target size for transmembrane segments elongation", default = 20)
    parser.add_argument("--max_length", type=int, help = "Maximum target size for transmembrane segments elongation", default = 70)
    parser.add_argument("--margin", type=int, help = "Supplementary distance from the membrane to consider amino acids of interest", default = 10)
    parser.add_argument("--inner_margin", 
                        type=int,
                        help = "Distance to artificially reduce the membrane width to capture more amino acids as out-of-membrane", 
                        default = 0)
    parser.add_argument("--gaps", type=int, help = "Minimum number of missing AA in a segment to discard it", default = 1)
    parser.add_argument("--output", help = "Output file", default = sys.stdout)
    i_opts = parser.add_mutually_exclusive_group(required=True)

    i_opts.add_argument('--fasta', help='File with iORFs in FASTA format', default = None)
    i_opts.add_argument('--csv', help='File with iORFs lengths in CSV format', default = None)

    args = parser.parse_args()

    if args.input:
        
        protein_name = os.path.basename(args.input).split(".")[0]

        sequences, sequences_short, structure, structure_short, protein_name = process_transmembrane_file(pdb_path = args.input,
                                                                                                            secondary_structure_path = None,
                                                                                                            margin = args.margin, 
                                                                                                            inner_margin = args.inner_margin, 
                                                                                                            min_length = args.min_length, 
                                                                                                            max_length = args.max_length, 
                                                                                                            gaps = args.gaps,
                                                                                                            iorf_path  = None, 
                                                                                                            csv_path = args.csv)
                
        print([len(record.seq) for record in sequences])

    elif args.input_files:

        sequences = []
        sequences_short = []
        structures = {}
        structures_short = {}

        for file in args.input_files:
            
            protein_name = os.path.basename(file).split(".")[0]

            try:

                res_dict = transmembrane(file_path = file, 
                        secondary_structure_path = args.secondary_structures,
                        margin = args.margin, 
                        inner_margin = args.inner_margin,
                        min_length = args.min_length, 
                        max_length = args.max_length, 
                        gaps = args.gaps,
                        iorf_path = args.fasta,
                        csv_path = args.csv)
                
                sequences.extend(res_dict["records"])
                sequences_short.extend(res_dict["records_shorts"])
                structures[protein_name] = res_dict["structures"]
                structures_short[protein_name] = res_dict["structures_shorts"]

            except Exception as e:

                transmembrane_logger.error(f"Error processing file {file}: {e}")
                print(f"Error processing file {file}: {e}")
                
                
    
    
                        
                        
if __name__ == "__main__":
    
    main()
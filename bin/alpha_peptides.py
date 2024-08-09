#!/usr/bin/env python3

import argparse
import logging
import sys
import os

from Bio import SeqIO

from utils import read_pdb, extract_horizontal_alpha_helix
from utils import setup_logger, exception_catcher

import argparse
import os

class ValidatePDBSecondaryStructure(argparse.Action):
    def __init__(self, option_strings, dest, nargs=None, **kwargs):
        if nargs is not None and nargs != '+':
            raise ValueError("Expecting number of args to be n+")
        
        super(ValidatePDBSecondaryStructure, self).__init__(option_strings, dest, nargs=nargs, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        if not hasattr(namespace, 'pdb_dict'):
            namespace.pdb_dict = {}
        if not hasattr(namespace, 'sec_struct_dict'):
            namespace.sec_struct_dict = {}
        
        is_pdb = '--pdbs' in option_string
        
        for value in values:
            base_name = os.path.splitext(os.path.basename(value))[0]
            if is_pdb:
                namespace.pdb_dict[base_name] = value
            else:
                namespace.sec_struct_dict[base_name] = value

        setattr(namespace, self.dest, values)


def pair_files(parser, namespace):

    pairs = []
    pdb_dict = namespace.pdb_dict
    sec_struct_dict = namespace.sec_struct_dict

    for base_name in pdb_dict:
        if base_name in sec_struct_dict:
            pairs.append((pdb_dict[base_name], sec_struct_dict[base_name]))
        else:
            parser.warning(f"No corresponding secondary structure file found for PDB file: {pdb_dict[base_name]}")
    
    for base_name in sec_struct_dict:
        if base_name not in pdb_dict:
            parser.warning(f"No corresponding PDB file found for secondary structure file: {sec_struct_dict[base_name]}")
    
    setattr(namespace, 'file_pairs', pairs)


def check_args(parser, args):

    """
    Input files should come as pairs, hence pdbS have to be passed along with structureS, 
    same for pdb and structure ! 
    """

    if args.pdb and args.pdbs:
        parser.error("Cannot use --input and --input_files at the same time.")
    
    if args.pdb and args.secondary_structures:
        parser.error("Unique secondary structure file is required when using --input.")
    
    if (args.pdbs and not args.secondary_structures) or (args.pdbs and args.secondary_structure):
        parser.error("Multiple secondary structure files are required when using --pdbs.")
    
    if not args.pdb and not args.pdbs:
        parser.error("Either --input or --input_files is required.")


alpha_peptides_logger = setup_logger('alpha_peptides_logger', 'alpha_peptides.log')
@exception_catcher(logging.getLogger('alpha_peptides_logger'))
def extract_peptides_alpha_helixes(file_path, secondary_structure_path, checking_depth, min_length, gaps):

    pdb_struct = read_pdb(alpha_peptides_logger, file_path, secondary_structure_path)
    
    alpha_helix = extract_horizontal_alpha_helix(alpha_peptides_logger, pdb_struct, checking_depth, min_length, gaps)

    return alpha_helix


def main():

    parser = argparse.ArgumentParser(description="Process PDB files and their corresponding secondary structure files.")
    
    parser.add_argument(
        '--pdbs',
        metavar='PDB',
        type=str,
        nargs='+',
        action=ValidatePDBSecondaryStructure,
        help='List of PDB files'
    )

    parser.add_argument(
        '--secondary_structures',
        metavar='STRUCT',
        type=str,
        nargs='+',
        action=ValidatePDBSecondaryStructure,
        help='List of secondary structure files'
    )

    parser.add_argument(
        '--pdb',
        metavar='PDB',
        type=str,
        help='Single PDB file'
    )

    parser.add_argument(
        '--secondary_structure',
        metavar='STRUCT',
        type=str,
        help='Single secondary structure file'
    )

    parser.add_argument(
        '--min_length',
        metavar='N',
        type=int,
        default=10,
        help='Minimum length of the alpha helix'
    )

    parser.add_argument(
        '--gaps',
        metavar='N',
        type=int,
        default=1,
        help='Minimum number of gaps in an helix to discard it',
    )

    parser.add_argument(
        '--depth',
        metavar='N',
        type=int,
        default=3,
        help='Number of AA between two predicted alpha helixes to consider them as a single one'
    )

    parser.add_argument(
        '--output',
        metavar='OUTPUT',
        type=str,
        help='Output file'
    )

    args = parser.parse_args()

    logger = logging.getLogger('alpha_peptides_logger')

    check_args(parser, args)

    if args.pdbs and args.secondary_structures:

        pair_files(parser, args)

        records = []

        for pdb, sec_struct in args.file_pairs:

            logger.info(f"Processing pair {pdb} and {sec_struct}")

            try:
            
                helix = extract_peptides_alpha_helixes(pdb, sec_struct, args.depth, args.min_length, args.gaps)

                if helix : records.extend(helix)

            except Exception as e:

                alpha_peptides_logger.error(f"Error processing PDB file: {pdb}. {e}")

    elif args.pdb and args.secondary_structure:

        records = extract_peptides_alpha_helixes(args.pdb, args.secondary_structure, args.depth, args.min_length, args.gaps)

    else:

        raise ValueError("Unexpected combination of arguments.")

    if args.output and records:
        
        SeqIO.write(records, args.output, "fasta")

    else:

        SeqIO.write(records, sys.stdout, "fasta")

if __name__ == '__main__':
    main()




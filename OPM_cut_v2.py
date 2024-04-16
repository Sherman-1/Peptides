#!/usr/bin/env python3

"""
This program takes an OPM pdb file, where a protein is embedded in a membrane, and isolates the transmembrane fragments.
It also reads a FASTA of yeast's iORFs and extracts their size.
Thus it can let a few residues on each fragment to obtain 'natural' sized peptides. 
"""

import argparse
import sys
from Bio import SeqIO
import random
from collections import OrderedDict
import os


def remove_outside_residues(input_file, size_fasta, min_size, max_size, output_format, output_file=None, margin=3.0, max_gap=5, keepmb=False, verbose=False):
    
    
    # Second pass: identify the intra-membrane fragments
    with open(input_file, 'r') as f_in:

        fragments = []
        current_fragment = []
        last_resSeq = None
        last_chain = None

        for (chain, resSeq), coords in res_coord_dic.items():
            x_value = coords['x']
            y_value = coords['y']
            z_value = coords['z']

            #### WARNING ####
            # If the residue is in the membrane,
            if inner_membrane_z <= z_value <= outer_membrane_z and min_membrane_x <= x_value <= max_membrane_x and min_membrane_y <= y_value <= max_membrane_y:

                # If the fragment is empty, start a new one with the current residue
                if current_fragment == []:
                    current_fragment.append((chain, resSeq, z_value))
                    last_resSeq = resSeq
                    last_chain = chain

                # If the current fragment is not empty,
                else :

                    # If the current residue belongs to the current fragment, add it.
                    if abs(resSeq - last_resSeq) <= max_gap+1 and chain == last_chain:
                        current_fragment.append((chain, resSeq, z_value))
                        last_resSeq = resSeq
                        last_chain = chain

                    # If the current residue belongs to a new fragment, close the current one and start a new one.
                    else :
                        fragments.append(current_fragment)
                        current_fragment = []
                        current_fragment.append((chain, resSeq, z_value))
                        last_resSeq = resSeq
                        last_chain = chain

            # If the residue is not in the membrane, close the current fragment if it is not empty.
            else:
                if current_fragment:
                    fragments.append(current_fragment)
                    current_fragment = []
                    last_resSeq = None
                    last_chain = None

        # At the end of the file, close the current fragment if it is not empty.
        if current_fragment:
            fragments.append(current_fragment)

    if verbose:
        print(f'res_coord_dic: from {min(res_coord_dic.keys())} to {max(res_coord_dic.keys())}')
        print(f'Fragments:')
        for fragment in fragments:
            # print([ele[0] for ele in fragment])
            chain = fragment[0][0]
            frag_start = fragment[0][1]
            frag_end = fragment[-1][1]
            frag_z_min = min(z_value for chain, resSeq, z_value in fragment)
            frag_z_max = max(z_value for chain, resSeq, z_value in fragment)
            print(f'Fragment from chain {chain} from {frag_start} to {frag_end}, z_min: {frag_z_min}, z_max: {frag_z_max}')

    # The transmembrane fragments are the ones that cross the center of the membrane (i.e. the z=0 plane)
    tm_fragments = [fragment for fragment in fragments if min(z_value for chain, resSeq, z_value in fragment) < 0 and max(z_value for chain, resSeq, z_value in fragment) > 0]
    if verbose: print(f'Transmembrane fragments:')

    def complete_tm_fragment(tm_fragment, frag_ind):
        """Complete a transmembrane fragment with its flanking residues to reach a target size.
        The target size is picked from a model FASTA file."""

        frag_chain = tm_fragment[0][0]

        if verbose:
            frag_start = tm_fragment[0][1]
            frag_end = tm_fragment[-1][1]
            frag_z_min = min(z_value for chain, resSeq, z_value in tm_fragment)
            frag_z_max = max(z_value for chain, resSeq, z_value in tm_fragment)
            print(f'\n***TM frag {frag_ind+1} : chain {frag_chain} from {frag_start} to {frag_end}, z_min: {frag_z_min}, z_max: {frag_z_max}***')

        # Pick a size for the peptide
        target_size = size_picker(size_fasta, min_size, max_size)
        fragment_size = len(tm_fragment)
        if verbose : print(f'target size: {target_size}, fragment size: {fragment_size}')
        
        # If the fragment is smaller than the target size, complete it
        if fragment_size < target_size: 
            extra_size = target_size - fragment_size
            if verbose : print(f'extra size: {extra_size}')
            
            # Determine how many extra residues will be added to each side of the fragment
            upstream_extra = random.choice(range(0, extra_size+1))
            downstream_extra = extra_size - upstream_extra
            if verbose : print(f'extra upstream: {upstream_extra}, extra downstream: {downstream_extra}')

            # Then if possible, add the extra residues to the fragment
            # Get the size of the fragment and if possible, add extra residues to it.
            def get_extra_fragment(extra_nb, side):
                extra_fragment = []
                gap = -1
                first_residue_index = tm_fragment[0][1]
                last_residue_index = tm_fragment[-1][1]
                for i in range(extra_nb):
                    if gap >= max_gap:
                        break
                    if side == 'upstream':
                        extra_residue_index = first_residue_index - i - 1
                    elif side == 'downstream':
                        extra_residue_index = last_residue_index + i + 1
                    if (frag_chain, extra_residue_index) in res_coord_dic:
                        gap = -1
                        extra_residue = (frag_chain, extra_residue_index, res_coord_dic[(frag_chain, extra_residue_index)])
                        if side == 'upstream':
                            extra_fragment.insert(0, extra_residue)
                        elif side == 'downstream':
                            extra_fragment.append(extra_residue)
                    else:
                        gap += 1

                return extra_fragment
            
            upstream_extra_fragment = get_extra_fragment(upstream_extra, 'upstream')
            if verbose :
                if upstream_extra_fragment:
                    print(f'upstream extra fragment: {upstream_extra_fragment[0][1]} to {upstream_extra_fragment[-1][1]}')
                else:
                    print(f'upstream extra fragment: None')
            downstream_extra_fragment = get_extra_fragment(downstream_extra, 'downstream')
            if verbose : 
                if downstream_extra_fragment:
                    print(f'downstream extra fragment: {downstream_extra_fragment[0][1]} to {downstream_extra_fragment[-1][1]}')
                else:
                    print(f'downstream extra fragment: None')
            tm_fragment = upstream_extra_fragment + tm_fragment + downstream_extra_fragment

            # If some residues could not be added upstream, try to add them downstream
            if len(upstream_extra_fragment) < upstream_extra and len(downstream_extra_fragment) == downstream_extra:
                remaining_upstream_extra = upstream_extra - len(upstream_extra_fragment)
                downstream_extra_fragment = get_extra_fragment(remaining_upstream_extra, 'downstream')
                if verbose : 
                    if downstream_extra_fragment:
                        print(f'completed downstream extra fragment: {downstream_extra_fragment[0][1]} to {downstream_extra_fragment[-1][1]}')
                    else:
                        print(f'completed downstream extra fragment: None')
                tm_fragment = tm_fragment + downstream_extra_fragment
            elif len(downstream_extra_fragment) < downstream_extra and len(upstream_extra_fragment) == upstream_extra:
                remaining_downstream_extra = downstream_extra - len(downstream_extra_fragment)
                upstream_extra_fragment = get_extra_fragment(remaining_downstream_extra, 'upstream')
                if verbose : 
                    if upstream_extra_fragment:
                        print(f'completed upstream extra fragment: {upstream_extra_fragment[0][1]} to {upstream_extra_fragment[-1][1]}')
                    else:
                        print(f'completed upstream extra fragment: None')
                tm_fragment = upstream_extra_fragment + tm_fragment

            # The final fragment
            if verbose : print(f'final fragment: {tm_fragment[0][1]} to {tm_fragment[-1][1]}')

        # Third pass: write the output by only keeping the crossing fragment
        with open(input_file, 'r') as f_in:
            if output_file:
                path_no_ext = os.path.splitext(output_file)[0]
                if output_format == 'pdb':
                    extension = '.pdb'
                elif output_format == 'fasta':
                    extension = '.faa'
                f_out = open(path_no_ext + "_frag_" + str(frag_ind+1) + extension, 'w')
            else:
                f_out = sys.stdout

            if output_format == 'pdb':
                for line in f_in:
                    if line.startswith('ATOM'):
                        chain = line[21]
                        if chain == frag_chain:
                            resSeq = int(line[22:26])
                            if (chain, resSeq) in [(chain, resSeq) for chain, resSeq, z_value in tm_fragment]:
                                f_out.write(line)
                    elif not line.startswith('HETATM') or (keepmb and 'DUM' in line):
                        f_out.write(line)

                if output_file:
                    f_out.close()

            elif output_format == 'fasta':
                basename = os.path.splitext(os.path.basename(input_file))[0]
                # Dictionary mapping three-letter codes to one-letter codes
                aa_dict = {
                    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E',
                    'PHE': 'F', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                    'LYS': 'K', 'LEU': 'L', 'MET': 'M', 'ASN': 'N',
                    'PRO': 'P', 'GLN': 'Q', 'ARG': 'R', 'SER': 'S',
                    'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
                }
                f_out.write(f'>{basename}_{frag_chain}_TM_frag_{frag_ind+1}\n')
                last_resSeq = None
                for line in f_in:
                    if line.startswith('ATOM') :
                        chain = line[21]
                        if chain == frag_chain:
                            resSeq = int(line[22:26])
                            if (chain, resSeq) in [(chain, resSeq) for chain, resSeq, z_value in tm_fragment] and resSeq != last_resSeq:
                                aa = line[17:20]
                                if aa in aa_dict:
                                    f_out.write(aa_dict[aa])
                                    last_resSeq = resSeq
                f_out.write('\n')

            if output_file:
                f_out.close()

    for i,tm_fragment in enumerate(tm_fragments):
        complete_tm_fragment(tm_fragment, i) # complete and write down each transmembrane fragment

def main():
    parser = argparse.ArgumentParser(description='Remove residues outside the membrane from a Membranome PDB file.')
    parser.add_argument('input_file', help='The input PDB file.')
    parser.add_argument('-s','--size_fasta', help='Model FASTA file to pick a size from.')
    parser.add_argument('-l', '--min_size', type=int, default=20, help='The minimal size of the peptide. Default is 20.')
    parser.add_argument('-u', '--max_size', type=int, default=70, help='The maximal size of the peptide. Default is 70.')
    parser.add_argument('-f','--output_format', choices=['pdb', 'fasta'], default='pdb', help='The output format. Default is pdb.')
    parser.add_argument('-o', '--output_file', help='The output PDB file. If not provided, output will be printed to stdout.')
    parser.add_argument('-m', '--margin', type=float, default=3.0, help='The distance threshold for removing residues. Default is 3.0.')
    parser.add_argument('-g', '--gap', type=int, default=5, help='The maximal gap between residues to consider them part of the same fragment. Default is 5')
    parser.add_argument('--keepmb', action='store_true', help='Keep the membrane atoms. Default is to remove them.')
    parser.add_argument('-v','--verbose', action='store_true', help='Print a lot a stuff.')
    args = parser.parse_args()

    if args.verbose:
        print(args)

    # Cut a peptide of the target size that includes the single-crossing membrane fragment of the PDB (write a new PDB).
    remove_outside_residues(args.input_file, args.size_fasta, args.min_size, args.max_size, args.output_format, args.output_file, args.margin, args.gap, args.keepmb, args.verbose)

if __name__ == "__main__":
    main()
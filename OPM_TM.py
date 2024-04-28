from utils import *


def OPM_TM(input_file): 
    
    pdb_struct = read_pdb(input_file)
    in_membrane_binaries, in_margin_binaries = return_binaries(pdb_struct)
    tm_indices = extract_tm_segments_indices(in_membrane_binaries)
    elongate_tm_segments(tm_indices, pdb_struct)
    segments_strings = segments_to_string(tm_indices, pdb_struct)
    combined_strings = combine_strings_v2(segments_strings, in_margin_binaries)
    coordinates = extract_coordinates_v3(combined_strings)
    
    return coordinates



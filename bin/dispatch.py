#!/usr/bin/env python3

import polars as pl 
import os 
import argparse

def main(metadata_path):

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
    
    ###### Peripherals ######
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

    metadata = (
        pl.read_csv(metadata_path, separator = ",", has_header = True)
        
        #.with_columns(
        #    pdbid = pl.col("pdbid").str.replace(pattern  = r'="(\w+)"', value = "${1}")
        #)   
    ) 

    ## FILTERS                                                     
    
    # Peptides or shorts proteins for which an alpha helix is completely crossing the membrane
    bitopics = ((pl.col("classtype_id") == 11) & (pl.col("thickness") >= 20) | ((pl.col("classtype_id") == 7) & (pl.col("thickness") >= 20))) # + membranome
    
    # Multitopic transmembranes, only alpha-helical 
    multitopics = (pl.col("classtype_id") == 1)

    # Peptides or short proteins located near the membrane surface
    associated_proteins = (pl.col("type_id") == 2) | ((pl.col("classtype_id") == 11) & (pl.col("thickness") < 20))

    # Alpha helical horizontal peptides, laying on the membrane surface in the " Interface " zone
    horizontal_alpha_peptides = (pl.col("classtype_id") == 7) & (pl.col("thickness") < 20) & (pl.col("tilt").is_between(70, 90, closed = "none"))

    # Alpha helical peptides not horizontal and just a little bit crossing the membrane or beta hairpins
    associated_peptides = ((pl.col("classtype_id") == 7) & (pl.col("thickness") < 20) & (pl.col("tilt") < 70)) | (pl.col("classtype_id") == 9)

    ## DISPATCHING

    # Bitopic proteins
    (metadata
        .filter(bitopics)
        .write_csv("bitopics.csv")
    )

    # Multitopic transmembranes
    (metadata
        .filter(multitopics)
        .write_csv("multitopics.csv")
    )

    # Peripherals
    (metadata
        .filter(associated_proteins)
        .write_csv("associated_proteins.csv")
    )

    # Horizontal alpha helical peptides
    (metadata
        .filter(horizontal_alpha_peptides)
        .write_csv("horizontal_alpha_peptides.csv")
    )

    # Associated alpha helix
    (metadata
        .filter(associated_peptides)
        .write_csv("associated_peptides.csv")
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("--metadata", type = str, required = True)

    args = parser.parse_args()

    main(args.metadata)

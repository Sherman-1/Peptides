process horizontal_batch { 

    input:
    path pdb_files
    path secondary_structures 

    output:
    path "horizontal.fasta", emit : horizontal_fasta

    script:
    """
    alpha_peptides.py --pdbs ${pdb_files} --secondary_structures ${secondary_structures} --min_length 4 >> horizontal.fasta
    """
}
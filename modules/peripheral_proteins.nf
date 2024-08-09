process cut_peripheral_proteins {

    clusterOptions '-q lowprio -l ncpus=2 -l mem=2gb'

    errorStrategy 'ignore'

    input:
    path pdb_path

    output:
    path "*.fasta", emit : peripheral
    

    script:
    """
    peripheral.py --input ${pdb_path} 
    """
}
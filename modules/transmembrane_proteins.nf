process cut_transmembrane_proteins {

    clusterOptions '-q lowprio -l ncpus=2 -l mem=2gb'

    errorStrategy 'ignore'

    input:
    path pdb_path
    //path secondary_structures
    each path(iORFs)

    output:
    path "*_elongated.fasta", emit : elong_fasta
    path "*_short.fasta", emit : short_fasta

    // --secondary_structures ${secondary_structures}
    script:
    """
    transmembrane.py --input ${pdb_path} --csv ${iORFs} 
    """


}

process cut_transmembrane_proteins_batch {

    clusterOptions '-q common -l ncpus=2 -l mem=2gb -l walltime=1000:00:00'

    
    publishDir "save/pdbs/polytopic/longs", pattern: "*_long.pdb", mode: 'copy'
    publishDir "save/pdbs/polytopic/shorts", pattern: "*_short.pdb", mode: 'copy'


    input:
    path pdb_paths
    //path secondary_structures
    path(iORFs)

    output:
    path "polytopics.fasta", emit : segment_elong
    path "polytopics_short.fasta", emit : segment_short
    path "*_long.pdb", emit : pdbs_elongs // several files
    path "*_short.pdb", emit : pdbs_shorts // several files

    
    script:
    """
    cat ${iORFs} > local_iORFs.csv
    transmembrane.py --input_files ${pdb_paths} --csv local_iORFs.csv 
    """

}

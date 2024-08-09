process check_format {

    clusterOptions '-l ncpus=1 -l mem=2gb -q lowprio'

    input:
    path pdb_file

    output:
    tuple path(pdb_file), stdout

    script:
    """
    count=\$(awk '\$1 == "ATOM" || \$1 == "HETATM"' ${pdb_file} | grep -c "\\*" || true)

    echo "\$count"
    """
}

process dispatch { 

    clusterOptions '-l ncpus=10 -l mem=20gb -q common'

    input:
    path metadata

    output:
    path "bitopics.csv", emit : bitopics
    path "multitopics.csv", emit : multitopics
    path "associated_proteins.csv", emit : associated_proteins
    path "horizontal_alpha_peptides.csv", emit : horizontal_alpha_peptides
    path "associated_peptides.csv", emit : associated_peptides

    script:
    """
    dispatch.py --metadata ${metadata}
    """
}

process remove_duplicate_segment { 

    clusterOptions '-q lowprio -l ncpus=2 -l mem=2gb'

    input:
    path pdbFile

    output:
    path "*_filtered.pdb"

    script:

    id = params.overlap_identity
    """
    pdb2fasta ${pdbFile} > ${pdbFile}.fasta
    cd-hit -i ${pdbFile}.fasta -o ${pdbFile}.cdhit -c ${id} -n 5 -d 0 -M 2000 -T 2 
    grep '^>' ${pdbFile}.cdhit | cut -d ':' -f2 | cut -f1 | sort | uniq > uniq_chain_IDs.txt
    awk '\$1 != "ATOM"' ${pdbFile} > ${pdbFile.baseName}_filtered.pdb
    check_chain_ID uniq_chain_IDs.txt ${pdbFile} >> ${pdbFile.baseName}_filtered.pdb
    """

}

process remove_duplicate_segment_batch { 

    input:
    path pdbFiles

    output:
    path "*_filtered.pdb"

    script:
    id = params.overlap_identity
    """
    for pdb_file in ${pdbFiles}; do

        base=\$(echo "\${pdb_file%.txt}")
        pdb2fasta \$pdb_file > \${base}.fasta
        cd-hit -i \${base}.fasta -o \${base}.cdhit -c ${id} -n 5 -d 0 -M 2000 -T 2 
        grep '^>' \${base}.cdhit | cut -d ':' -f2 | cut -f1 | sort | uniq > uniq_chain_IDs.txt
        awk '\$1 != "ATOM"' \$pdb_file > \${base}_filtered.pdb
        check_chain_ID uniq_chain_IDs.txt \${pdb_file} >> \${base}_filtered.pdb

    done
    """
}
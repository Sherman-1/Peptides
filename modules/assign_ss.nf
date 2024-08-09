process mkdssp { 

    clusterOptions '-q common -l ncpus=2 -l mem=6gb'   

    label 'mkdssp'

    errorStrategy 'ignore'

    input:
    path pdb 

    output:
    tuple path("*.dssp"), path(pdb)

    script:
    """
    pdb_tidy="/home/simon.herman/.local/bin//pdb_tidy"
    pdb_tocif="/home/simon.herman/.local/bin//pdb_tocif"
    mkdssp_container="/store/EQUIPES/BIM/MEMBERS/simon.herman/Peptides/singularity_cache/hermansimon-dssp-4.4.img"

    awk '\$1 == "ATOM"' ${pdb} > temp.pdb 
    \$pdb_tidy temp.pdb > temp_clean.pdb
    singularity exec --bind \$(pwd):/data/ \$mkdssp_container mkdssp --write-other --verbose /data/temp_clean.pdb > temp.dssp 2> log
    extract_ss temp.dssp ${pdb.baseName}.dssp 2>> log
    """

}

process stride_batch {

    clusterOptions '-q common -l ncpus=2 -l mem=4gb'

    input:
    path pdbs 

    output:
    path "*.stride", emit : strides
    path pdbs, emit : pdbs

    script:
    """
    #!/bin/bash

    for pdb in ${pdbs} ; do
        base=\$(echo "\${pdb%.*}")
        stride -f \${pdb} > temp || true 
        if [ -s temp ] ; then
            # Extract chain_id, res_number and secondary structure
            awk '/^ASG/ {print \$3 "\t" \$4 "\t" \$6}' temp > \${base}.stride
        fi
    done
    """
}


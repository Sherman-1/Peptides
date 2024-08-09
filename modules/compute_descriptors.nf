
process hydromoment {

    clusterOptions '-q lowprio -l ncpus=2 -l mem=2gb'

    input:
    path fasta_file

    output:
    path "*.csv", emit : hydro

    script:
    """
    hydromoment.py --seqfile ${fasta_file} --window -1 --scale Fauchere-Pliska
    """
    
}


process HIT { 

    clusterOptions '-q lowprio -l ncpus=6 -l mem=16gb'

    label 'ORFMine'

    input:
    path fasta

    output:
    path "*.csv"


    script:
    """
    orfold -faa ${fasta} -options HIT > ${fasta.baseName}_HIT.csv
    """


}

process topology {

    clusterOptions '-q lowprio -l ncpus=2 -l mem=2gb'

    input:
    path pdbs // multiple pdb files

    output:
    path "hydromoment.csv", emit : hydro
    path "topology.csv", emit : topo

    script:
    """
    for file in $pdbs; do 

        base=\$(basename \$file .ent) 
        echo -n "\$base," > buffer


        # Compute topology metrics
        echo -n "\$(HullRadV9.py \$file)," >> buffer

        # Compute secondary structure content percentage
        stride \$file | awk '\$1 == "ASG" {print \$6}' | sort | uniq -c > temp 
        awk '{if (\$2=="H" || \$2=="E") sum+=\$1; total+=\$1} END {print (sum/total)*100}' temp >> buffer

        # grep command to look for error messages in the current buffer 
        if grep -q -e "missing" -e "found" -e "residue" buffer; then
            echo "Error in \$base"
            continue
        fi

        # If no error, append the buffer to the final file
        cat buffer >> topology.csv

        pdb2fasta \$file >> seq.fasta

    done

    hydromoment.py --window -1 --seqfile seq.fasta --outfile hydromoment.csv
    """

}

process hullrad { 

    clusterOptions '-q lowprio -l ncpus=2 -l mem=2gb'

    input:
    path pdb 


    output:
    path csv 

    script
    """
    HullRadV9.py ${pdb} > ${pdb.baseName}.csv
    """
}
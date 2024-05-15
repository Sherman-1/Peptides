#!/bin/bash

# Define the input PDB file
pdb_file="1uaz.pdb"

# Command pipeline
awk '$1 == "ATOM"' "$pdb_file" | \
    mkdssp --output-format mmcif --write-other --quiet /dev/stdin 2> /dev/null | \
    awk '{
    # Check for the line that contains the specific pattern
    if (/^_struct_conf.end_auth_seq_id/) {
        flag = 1;
        prevLine = $0;
        next;
    }
    
    # Print the previous line if we are past the specified line and not on the first subsequent line
    if (flag && NR > 1) {
        print prevLine;
    }

    # Update the previous line with the current one
    if (flag) {
        prevLine = $0;
    }
}
END {
    # Do not print the last stored line
}' | tail -n +2 | ./extract_SS /dev/stdin > 1uaz.ss




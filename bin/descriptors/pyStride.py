#!/usr/local/bin/python3

import subprocess


def _parse_output(stdout) -> dict:
    
    """
    Parse an ASG line from STRIDE
    
     ASG Detailed secondary structure assignment

    Format:  
    6-8  Residue name
    10-10 Protein chain identifier
    12-15 PDB	residue	number
    17-20 Ordinal residue number
    25-25 One	letter secondary structure code	**)
    27-39 Full secondary structure name
    43-49 Phi	angle
    53-59 Psi	angle
    65-69 Residue solvent accessible area

    Input : 
    
        - stdout of stride
        
    Output : 
    
        - Dict of secondary structure assignment for each AA for each chain.
        
    """
    
    
    
    _ =  [ line for line in stdout.split('\n') if line.startswith('ASG') ]
    
    _dict = {}
    
    for line in _:
        
        chain_id = line[9].strip()
        res_number = int(line[11:15])
        secondary_structure = line[24]
        
        if chain_id not in _dict:
            
            _dict[chain_id] = {}
            
        else:
            
            _dict[chain_id][res_number] = secondary_structure
            
    return _dict
    
    
def stride(pdbpath) -> dict:
    
    """
    Perform STRIDE analysis

    Parameters
    ----------
    pdbpath : str
      path to the pdb file to analyse

    Returns
    -------
    output : str
      The raw stdout from stride as a string
    """
    
    p = subprocess.Popen(['stride', pdbpath],
                         stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE)
    stdout, stderr = p.communicate()
    
    output = stdout.decode(encoding = "utf-8")
    
    return _parse_output(output)
    

if __name__ == "__main__": 
    
    import argparse
    
    
    parser = argparse.ArgumentParser(
                    prog='PyStride'
    )
    
    parser.add_argument("--input", type = str, required = True,
                        help = "PDB path")
    
    
    args = parser.parse_args()
    
    print(stride(args.input))
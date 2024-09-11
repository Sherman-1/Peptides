#!/usr/bin/env python3

def generate_random_protein_sequences(n, min_length, max_length) -> list[str]:
    
    sequences = []
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    
    for _ in range(n):
        length = random.randint(min_length, max_length)
        sequence = ''.join(random.choices(aa, k=length))
        sequences.append(sequence)
    
    return sequences
    
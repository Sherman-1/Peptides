#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 18:33:22 2020

@author: mheinzinger
"""

print("================ Loading Libraries =================")

import warnings
import logging
warnings.filterwarnings("ignore")

import time
from tqdm import tqdm 
import numpy as np
import h5py

import torch
from transformers import T5EncoderModel, T5Tokenizer

from Bio import SeqIO
import polars as pl 

import random


logger = logging.getLogger("memory_profile")

"""
torch.cuda.memory._record_memory_history(
    max_entries=100000
)


logger.info("=====================================")
logger.info("Memory allocated before loading model: {} GB".format(torch.cuda.memory_allocated()/1e9))
logger.info("Mem info before loading model: {}".format(torch.cuda.mem_get_info()))  
logger.info("Free memory before loading model: {} GB".format(torch.cuda.memory_reserved()/1e9))
"""


def generate_random_protein_sequences(n, min_length, max_length) -> list[str]:
    
    sequences = []
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    
    for _ in range(n):
        length = random.randint(min_length, max_length)
        sequence = ''.join(random.choices(aa, k=length))
        sequences.append(sequence)
    
    return sequences


def get_T5_model(device):

    model = T5EncoderModel.from_pretrained("Rostlab/prot_t5_xl_half_uniref50-enc")
    tokenizer = T5Tokenizer.from_pretrained('Rostlab/prot_t5_xl_half_uniref50-enc', do_lower_case=False)

    
    model = model.to(device) 
    model = model.eval()
    
    if device==torch.device("cpu"):
        print("Casting model to full precision for running on CPU ...")
        model.to(torch.float32)

    model = model.to(device)
    model = model.eval()

    return model, tokenizer

def get_sequences(seq_path, max_seq):

    AA = { "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y" }

    _ = { record.id : str(record.seq.upper()) for record in SeqIO.parse(seq_path, "fasta") if set(record.seq.upper()) <= AA and 20 <= len(record.seq) <= 70 }

    if len(_) > max_seq:

        sampled_keys = np.random.choice(list(_.keys()), size=max_seq, replace=False)
        seqs = {key : _[key] for key in sampled_keys}

    else:

        seqs = _

    if not seqs:
        raise ValueError("No valid sequences found in the input file")

    return seqs, max([len(seq) for seq in seqs.values()])

def get_embeddings( model, tokenizer, device, seqs, per_residue, per_protein, 
                   max_residues=4000, max_seq_len=1000, max_batch=100 ):

    """
    # per_residue indicates that embeddings for each residue in a protein should be returned.
    # per_protein indicates that embeddings for a whole protein should be returned (average-pooling)
    # max_residues gives the upper limit of residues within one batch
    # max_seq_len gives the upper sequences length for applying batch-processing
    # max_batch gives the upper number of sequences per batch
    """

    results = {
                "residue_embs" : dict(), 
                "protein_embs" : dict(),
                "sec_structs" : dict() 
               }

    # sort sequences according to length (reduces unnecessary padding --> speeds up embedding)
    seq_dict   = sorted( seqs.items(), key=lambda kv: len( seqs[kv[0]] ), reverse=True )
    start = time.time()
    batch = list()

    for seq_idx, (pdb_id, seq) in tqdm(enumerate(seq_dict,1), desc = "Batching sequences", total=len(seq_dict)):
        seq = seq
        seq_len = len(seq)
        seq = ' '.join(list(seq))
        batch.append((pdb_id,seq,seq_len))

        # count residues in current batch and add the last sequence length to
        # avoid that batches with (n_res_batch > max_residues) get processed 
        n_res_batch = sum([ s_len for  _, _, s_len in batch ]) + seq_len 
        if len(batch) >= max_batch or n_res_batch>=max_residues or seq_idx==len(seq_dict) or seq_len>max_seq_len:

            pdb_ids, seqs, seq_lens = zip(*batch)
            batch = list()

            # add_special_tokens adds extra token at the end of each sequence
            token_encoding = tokenizer.batch_encode_plus(seqs, add_special_tokens=True, padding="longest")
            input_ids      = torch.tensor(token_encoding['input_ids']).to(device)
            attention_mask = torch.tensor(token_encoding['attention_mask']).to(device)

            
            try:

                with torch.no_grad():

                    logger.info("Embedding for {} (L={})".format(pdb_id, seq_len))
                    logger.info("Batch size: {}".format(input_ids.size(0)))

                    # returns: ( batch-size x max_seq_len_in_minibatch x embedding_dim )
                    embedding_repr = model(input_ids, attention_mask=attention_mask)
                    #logger.info("Memory allocated: {:.2f} GB".format(torch.cuda.memory_allocated()/1e9))
                    #logger.info("Mem info : {}".format(torch.cuda.mem_get_info()))


            except RuntimeError:
                print("RuntimeError during embedding for {} (L={})".format(pdb_id, seq_len))
                continue


            for batch_idx, identifier in enumerate(pdb_ids): # for each protein in the current mini-batch
                s_len = seq_lens[batch_idx]
                # slice off padding --> batch-size x seq_len x embedding_dim  
                emb = embedding_repr.last_hidden_state[batch_idx,:s_len]

                if per_residue: # store per-residue embeddings (Lx1024)
                    results["residue_embs"][ identifier ] = emb.detach().cpu().numpy().squeeze()
                if per_protein: # apply average-pooling to derive per-protein embeddings (1024-d)
                    protein_emb = emb.mean(dim=0)
                    results["protein_embs"][identifier] = protein_emb.detach().cpu().numpy().squeeze()


    passed_time=time.time()-start
    avg_time = passed_time/len(results["residue_embs"]) if per_residue else passed_time/len(results["protein_embs"])
    print('\n############# EMBEDDING STATS #############')
    print('Total number of per-residue embeddings: {}'.format(len(results["residue_embs"])))
    print('Total number of per-protein embeddings: {}'.format(len(results["protein_embs"])))
    print("Time for generating embeddings: {:.1f}[m] ({:.3f}[s/protein])".format(
        passed_time/60, avg_time ))
    print('\n############# END #############')

    logger.info("Total number of per-residue embeddings: {}".format(len(results["residue_embs"])))
    logger.info("Total number of per-protein embeddings: {}".format(len(results["protein_embs"])))
    logger.info("Time for generating embeddings: {:.1f}[m] ({:.3f}[s/protein])".format(
        passed_time/60, avg_time ))
    
    #logger.info("Memory state after all embeddings are generated: {} GB".format(torch.cuda.memory_allocated()/1e9))
    


    return results


def main():

    device = torch.device('cuda:0' if torch.cuda.is_available() else 'cpu')

    if not torch.cuda.is_available():
        print("Warning: No GPU found")

    print("Using device: {}".format(device))


    print("Loading model ...")
    model, tokenizer = get_T5_model(device)

    categories = {

        #"Small" : "../christos/out_small.faa",
        #"Small" : "../christos/fastas/Small_christos.fasta",
        #"Folded" : "../christos/fastas/Folded_christos.fasta",
        #"Disordered" : "../christos/fastas/Disordered_christos.fasta",
        #"Bound" : "/store/EQUIPES/BIM/MEMBERS/simon.herman/Peptides/save/fasta/test.fasta",
        "Bitopic_short" : "../save/fasta/bi_shorts.fasta",
        "Bitopic_long" : "../save/fasta/bi_longs.fasta",
        "Polytopic_short" : "../save/fasta/poly_shorts.fasta",
        "Polytopic_long" : "../save/fasta/poly_longs.fasta",
        

    }

    dfs = {}
    
    for category, fasta in categories.items():
        
        print(f"\n############# Processing category : {category} #############")

        sequences, max_protein_length = get_sequences(fasta, 5000)

        per_residue = False
        per_protein = True
        
        print("\n ### Generating embeddings ###")

        embeddings = get_embeddings(model, tokenizer, device, sequences, per_residue, per_protein, max_residues=4000, max_seq_len=max_protein_length, max_batch=100)

        per_prot = embeddings["protein_embs"]

        data = [ [id, *list(emb)] for id, emb in per_prot.items()]

        dfs[category] =  pl.DataFrame(data, schema = ["id", *[f"emb_{i+1}" for i in range(1024)]]).with_columns(category = pl.lit(category))

    pl.concat(dfs.values()).write_parquet("random_and_reversed_christos.parquet", compression = "zstd")    

    return 0


if __name__ == '__main__':

    main()
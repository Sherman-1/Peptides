#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 19/07/2024 

@author : Simon HERMAN

Based on Christos Papadopoulos work
"""

print("============== Importing dependencies ============= ")

import warnings 
warnings.filterwarnings("ignore")

import tool_functions
import random

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

import iupred_funcs
import numpy as np

from pyHCA import HCA
from pyHCA.core.annotateHCA import _annotation_aminoacids
from pyHCA.core.classHCA import compute_disstat
from MMseqs import MMseqs2API

import multiprocessing

import argparse
import polars as pl  

import os



# ------------------------ #
#   AA index descriptors   #
# ------------------------ #

AA = { "A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V" }

with open("selected_aaindex1_reformated_58_Van_Westen.dms","r") as f:
    aa_indexes_names = []
    for line in f:
        if line.startswith("#"):
            aa = line.split()[1:]
            continue
        aa_indexes_names.append(line.split()[0].strip())

with open("AA_index_scaled.txt","r") as f:
    aa_indexes = {}
    for line in f:
        aa = line.split()[0].strip()
        for i,aaj in enumerate(line.split()[1:]):
            if aa_indexes_names[i] not in aa_indexes:
                aa_indexes[aa_indexes_names[i]] = {}
            aa_indexes[aa_indexes_names[i]][aa] = float(aaj)
            

def generate_random_protein_sequences(n, min_length, max_length) -> list[str]:
    
    sequences = []
    aa = 'ACDEFGHIKLMNPQRSTVWY'
    
    for _ in range(n):
        length = random.randint(min_length, max_length)
        sequence = ''.join(random.choices(aa, k=length))
        sequences.append(sequence)
    
    return sequences
    

def calculate_PDT(seq , aa_indexes , lamda):

    PDT = {}
    lamda = 2
    for aaj in aa_indexes:
        Ds = []
        for x in range(len(seq) - lamda):
            aa1 = seq[x]
            aa2 = seq[x + lamda]
            D = np.square(tool_functions.I(aa=aa1, aaj=aaj, aa_indexes=aa_indexes) - tool_functions.I(aa=aa2, aaj=aaj, aa_indexes=aa_indexes))
            Ds.append(D)

        PDT[aaj] = str(round(np.sum(Ds) / (len(seq) - lamda), 4))
        
    return(PDT)


def compute_sequence_metrics(record : SeqRecord):

    seq = str(record.seq)
    name = record.id

    try:

        # Calculate the mean of 58 aa index descriptors
        mean_aa_index = tool_functions.calculate_the_mean(seq, aa_indexes)

        # ==================================================================== #
        # Calculate the PDTn of 58 aa index descriptors
        #PDT2 = calculate_PDT(seq =seq , aa_indexes = aa_indexes, lamda =2)
        #PDT3 = calculate_PDT(seq =seq , aa_indexes = aa_indexes, lamda =3)
        #PDT4 = calculate_PDT(seq =seq , aa_indexes = aa_indexes, lamda =4)
        #PDT5 = calculate_PDT(seq =seq , aa_indexes = aa_indexes, lamda =5)
        #PDT8 = calculate_PDT(seq =seq , aa_indexes = aa_indexes, lamda =8)



        # ------ Protein Sequence Based Descriptors:   ----------------------------------- #

        hydrophobic_portion = tool_functions.hydrophobic_perc(seq)
        iupred_score   = iupred_funcs.iupred(seq=seq, mode="short")[0]
        iupred_portion = tool_functions.calculate_proportion_of_seq_disordered(iupred_score)
        iupred_mean    = round(sum(iupred_score) / len(iupred_score), 3)
        anchor_score   = iupred_funcs.anchor2(seq=seq,iupred_scores=iupred_score)
        anchor_portion = tool_functions.calculate_proportion_of_seq_disordered(anchor_score)
        score, pvalue = compute_disstat(0, len(seq), _annotation_aminoacids(seq=seq, method="domain", verbose=False)["cluster"])
        HCA_score = round(score, 2)
        
        aa_frequency = { aa : seq.count(aa) / len(seq) for aa in AA }

        #hca = HCA(seq=seq, querynames=name)
        #bc, bc_SS = tool_functions.get_hca_barcode(hca=hca, orf=name)

        return { 

                "id" : name,
                **mean_aa_index,
                "hydrophobic_portion" : hydrophobic_portion,
                "iupred_portion" : iupred_portion,
                "iupred_mean" : iupred_mean,
                "anchor_portion" : anchor_portion,
                "HCA_score" : HCA_score,
                **aa_frequency,
                "seq_len" : len(seq)
                #"bc" : bc,
                #"bc_SS" : bc_SS

            }
    
    except Exception as e:
        print(f"Exception while parsing seq {id} : {e}")
        return None


def pool_process(records, num_processes):
    
    with multiprocessing.Pool(processes=num_processes) as pool:
        print(f"Processing {len(records)} sequences using {num_processes} processes")
        results = pool.map(compute_sequence_metrics, records)

    return [ res for res in results if res is not None ]

def main(records, num_processes, category):

    # Filter out bad records 
    records = [ record for record in records if len(record.seq) > 0 and set(record.seq) <= AA ]
    
    results = pool_process(records, num_processes)

    columns = [
        "id",
        "seq_len",
        "hydrophobic_portion",
        "iupred_portion",
        "iupred_mean",
        "anchor_portion",
        "HCA_score",
        "A",
        "R",
        "N",
        "D",
        "C",
        "Q",
        "E",
        "G",
        "H",
        "I",
        "L",
        "K",
        "M",
        "F",
        "P",
        "S",
        "T",
        "W",
        "Y",
        "V",
        "ARGP820103",
        "BHAR880101",
        "CHAM810101",
        "CHAM820101",
        "CHAM830101",
        "CHAM830107",
        "CHAM830108",
        "CHOP780201",
        "CHOP780202",
        "CHOP780203",
        "CIDH920105",
        "FASG760101",
        "FAUJ880102",
        "FAUJ880103",
        "FAUJ880104",
        "FAUJ880105",
        "FAUJ880106",
        "FAUJ880109",
        "FAUJ880110",
        "FAUJ880111",
        "FAUJ880112",
        "FAUJ880113",
        "GRAR740102",
        "JANJ780102",
        "JANJ780103",
        "JOND920102",
        "JUNJ780101",
        "KLEP840101",
        "KRIW790101",
        "KYTJ820101",
        "LEVM760102",
        "LEVM760103",
        "LEVM760104",
        "LEVM760105",
        "LEVM760106",
        "LEVM760107",
        "NISK800101",
        "NISK860101",
        "PONP800101",
        "RACS770103",
        "RADA880108",
        "ROSG850101",
        "ROSG850102",
        "ROSM880102",
        "WARP780101",
        "WOLR810101",
        "VINM940101",
        "TAKK010101",
        "MONM990201",
        "KOEP990101",
        "KOEP990102",
        "MITS020101",
        "COSI940101",
        "PONP930101",
        "ZHOH040102",
        "ZHOH040103",
        "BAEK050101",
        "CASG920101",
        "category"
    ]


    special_columns = {

        "id": pl.Utf8,
        "category": pl.Utf8,
        "seq_len": pl.Int32,
    }

    schema = {col: special_columns.get(col, pl.Float32) for col in columns}
        
    df = pl.DataFrame(results, schema=schema)

    df = df.with_columns(category = pl.lit(category)).select(columns)

    print(f"Dataframe shape : {df.shape}")

    print(df.head())

    return df

def validate_fasta(fasta):

    if os.path.exists(fasta):
        return True
    else:
        return False
    

if __name__ == '__main__':

    categories = {

        "Small" : "fastas/Small_christos.fasta",
        #"Small" : "../christos/fastas/Small_christos.fasta",
        "Folded" : "fastas/Folded_christos.fasta",
        "DIBS" : "fastas/DIBS_christos.fasta",
        "Disordered" : "fastas/Disordered_christos.fasta",
        "Bound" : "/store/EQUIPES/BIM/MEMBERS/simon.herman/Peptides/save/fasta/test.fasta",
        "Bitopic" : "fastas/Transmembrane_christos.fasta", 
        #"Polytopic" : "../save/fasta/poly_shorts.fasta",

    }
    
    dfs = {}

    mmseqs = MMseqs2API(threads = 50, cleanup = True)

    assert all([validate_fasta(fasta) for fasta in categories.values()]), "One or more files are missing"
    
    random_sequences = generate_random_protein_sequences(500,20,70)
    
    random_records = [ SeqRecord(seq = Seq(random_seq), id = str(i)) for i, random_seq in enumerate(random_sequences) ]
    
    data = main(random_records, num_processes = 50, category  = "random")
    
    dfs["random"] = data
    
    for category, fasta in categories.items():

        print(f"Processing category : {category}")

        records = list(SeqIO.parse(fasta, "fasta"))

        records = [SeqRecord(seq = Seq(record.seq.upper()), id = record.id, description = "") for record in records if 20 <= len(record.seq) <= 70 and set(record.seq.upper()).issubset(AA)]

        if len(records) > 500: 

            records = np.random.choice(records, 500, replace = False)

        #SeqIO.write(records, f"sample.fasta", "fasta")

        #mmseqs.createdb(fasta = "sample.fasta")
        #mmseqs.cluster(coverage=0.7, identity=0.3, cov_mode=0)
        #mmseqs.createtsv(tsv_file = f"{category}_clu.tsv")

        #os.remove("sample.fasta")
        
        """
        buff = pl.read_csv(f"mmseqs_temp/{category}_clu.tsv", separator = "\t", has_header = False, new_columns = ["representative", "id"])

        clusters = (
            buff
            .join(buff.group_by("representative").agg(pl.count("id").alias("cluster_size")), on="representative")
            .with_columns(weight = 1 / pl.col("cluster_size"))
            .select(["id", "weight"])
            .cast({"id": pl.Utf8, "weight": pl.Float32})
        )
        """
        data = main(records, num_processes = 50, category  = category)
        
        dfs[category] = data

        #data = data.join(clusters, on = "id")
        
        
        
    for df in dfs.values(): 
        
        print(df.shape)

    pl.concat(dfs.values()).write_parquet("dry_run.parquet")




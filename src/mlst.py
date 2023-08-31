# -*- coding: utf-8 -*-

import os

from Bio import SeqIO
import matplotlib
import matplotlib.pyplot as plt

from src import dataframe as dfio

matplotlib.use("Agg")


def _ensure_fasta_present(sample_name):
    gene_file = "output/genes/H.pylori.J99.2001/" + sample_name + "/genes.fa"
    if not os.path.isfile(gene_file):
        raise FileNotFoundError(f"File {gene_file} was not found.")
    return gene_file


def _concatenate_pylori_mlst(sample_name):
    gene_file = _ensure_fasta_present(sample_name)
    mlst_genes = ("atpA", "efp", "mutY", "ppa", "trpC", "ureI", "yphC")
    reverse_complement = (
        "atpA", "mutY", "trpC", "ureI"
    )
    sequences = []
    concatenated_sequence = ""
    for record in SeqIO.parse(gene_file, "fasta"):
        if record.id in mlst_genes:
            if record.id in reverse_complement:
                record.seq = record.reverse_complement().seq
                print(record.id)
            sequences.append(record)
    
    concatenated_sequence = "".join([f">{seq.id}\n{seq.seq}\n" for seq in sequences])
    return concatenated_sequence


def _check_cagA(sample_name):
    gene_file = _ensure_fasta_present(sample_name)
    cagA: SeqIO.SeqRecord = None
    for record in SeqIO.parse(gene_file, "fasta"):
        if record.id == "cagA":
            cagA = record
    
    
    dataframe = dfio.get_or_create('pylori')
    if cagA.count('N') > .75*len(cagA):
        dataframe.at[sample_name, 'cagA'] = "cagA-"
    else:
        dataframe.at[sample_name, 'cagA'] = "cagA+"

    dfio.save(dataframe)
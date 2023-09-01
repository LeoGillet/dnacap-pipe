# -*- coding: utf-8 -*-

import operator
import os
import threading
import time

from Bio import SeqIO, pairwise2
import matplotlib
import matplotlib.pyplot as plt
from tqdm import tqdm

from src import dataframe as dfio
from src import common

matplotlib.use("Agg")


def _ensure_fasta_present(sample_name):
    gene_file = "output/genes/H.pylori.J99.2001/" + sample_name + "/genes.fa"
    if not os.path.isfile(gene_file):
        raise FileNotFoundError(f"File {gene_file} was not found.")
    return gene_file


def _concatenate_pylori_mlst(sample_name):
    gene_file = _ensure_fasta_present(sample_name)
    mlst_genes = ("atpA", "efp", "mutY", "ppa", "trpC", "ureI", "yphC")
    reverse_complement = ("atpA", "mutY", "trpC", "ureI")
    sequences = {}
    concatenated_sequence = ""
    for record in SeqIO.parse(gene_file, "fasta"):
        if record.id in mlst_genes:
            if record.id in reverse_complement:
                record.seq = record.reverse_complement().seq
            sequences[record.id] = record

    for gene in mlst_genes:
        concatenated_sequence += sequences[gene].seq
    return concatenated_sequence


def _load_falush_populations(path="typing/populations_falush.fasta"):
    seq_records = []
    for record in SeqIO.parse(path, "fasta"):
        seq_records.append(record)
    return seq_records


def _align_to_falush(concatenated_sequence, falush_sequences):
    alignments = {}
    for falush_seq in falush_sequences:
        alignment_score = pairwise2.align.globalxx(
            concatenated_sequence, falush_seq.seq, score_only=True
        )
        alignments[falush_seq.id] = alignment_score
    return dict(sorted(alignments.items(), key=operator.itemgetter(1), reverse=True)[:5])


def _thread_mlst_typing(sample_name):
    alignment_results[sample_name] = _align_to_falush(
        _concatenate_pylori_mlst(sample_name), _load_falush_populations()
    )


def start_mlst_typing(pipeline_type):
    dataframe = dfio.get_or_create(pipeline_type)
    global alignment_results
    alignment_results = {}

    falush_threads = []
    with tqdm(total=len(list(dataframe.index))) as pbar:
        for sample_name in list(dataframe.index):
            thread = threading.Thread(
                target=_thread_mlst_typing, args=(sample_name,), daemon=True
            )
            pbar.set_description(f"Typing {sample_name}")
            falush_threads.append(thread)
            if threading.active_count() >= common.get_max_threads("alignment"):
                time.sleep(1)
                pbar.refresh()
            thread.start()
            pbar.update()

        for thread in falush_threads:
            thread.join() 
        pbar.close()

    for sample_name, top_matches in alignment_results.items():
        pbar.set_description(f"Saving alignment results of {sample_name}")
        best_match, match2, match3, match4, match5 = top_matches.keys()
        dataframe.at[
            sample_name, "MLST best match"
        ] = f"{best_match} ({top_matches[best_match]})"
        dataframe.at[
            sample_name, "MLST match 2"
        ] = f"{match2} ({top_matches[match2]})"
        dataframe.at[
            sample_name, "MLST match 3"
        ] = f"{match3} ({top_matches[match3]})"
        dataframe.at[
            sample_name, "MLST match 4"
        ] = f"{match4} ({top_matches[match4]})"
        dataframe.at[
            sample_name, "MLST match 5"
        ] = f"{match5} ({top_matches[match5]})"
    dfio.save(dataframe)


def _check_cagA(sample_name):
    gene_file = _ensure_fasta_present(sample_name)
    cagA: SeqIO.SeqRecord = None
    for record in SeqIO.parse(gene_file, "fasta"):
        if record.id == "cagA":
            cagA = record

    dataframe = dfio.get_or_create("pylori")
    if cagA.count("N") > 0.75 * len(cagA):
        dataframe.at[sample_name, "cagA"] = "cagA-"
    else:
        dataframe.at[sample_name, "cagA"] = "cagA+"

    dfio.save(dataframe)

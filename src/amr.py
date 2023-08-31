# -*- coding: utf-8 -*-

from Bio import SeqIO
import pandas as pd

from src import common
from src import dataframe as dfio


def _parse_amr_tab(tab_file) -> list[tuple[str, str, str, str, str]]:
    with open(tab_file, "r", encoding="UTF-8") as tabf:
        lines = [line.strip() for line in tabf.readlines()]

    mutations = []
    for line in lines:
        if line.startswith("#") or line.startswith(">"):
            continue
        mutations.append(line.split("\t"))
    return mutations


def _variant_detection(nuc_seq: SeqIO.SeqRecord, known_mutations, type):
    detected_mutations = set()
    unknown_char = "N" if type == "nuc" else "X"
    for _, start_pos, end_pos, ref_char, alt_char in known_mutations:
        start_pos, end_pos = int(start_pos), int(end_pos)
        if start_pos == end_pos:
            try:
                seq_char = nuc_seq.seq[start_pos - 1].upper()
            except IndexError:
                print(start_pos, len(nuc_seq.seq))
            mutation_code = f"{ref_char}{start_pos}{alt_char}"
        else:
            seq_char = nuc_seq.seq[start_pos - 1 : end_pos].upper()
            mutation_code = f"{ref_char}{start_pos}-{end_pos}{alt_char}"
            print(seq_char, alt_char, mutation_code)
        if unknown_char in seq_char:
            detected_mutations.add("?")
            continue
        if seq_char == alt_char:
            detected_mutations.add(mutation_code)
    return detected_mutations


def _fetch_all_gene(file_path):
    sequences = []
    for record in SeqIO.parse(file_path, "fasta"):
        sequences.append(record)
    return sequences


def _rename_seqs_no_gene(sequences: list[SeqIO.SeqRecord], gene):
    for seq in sequences:
        seq.id = seq.id.replace(gene + "|", "")
        seq.name = seq.name.replace(gene + "|", "")
    return sequences


def _export_mutations_to_dataframe(mutations: dict, pipeline_type):
    dataframe = dfio.get_or_create(pipeline_type)
    for sample_name, genes in mutations.items():
        for gene, mutations in genes.items():
            dataframe.at[sample_name, gene + " Mutations"] = ";".join(mutations)
    dfio.save(dataframe)


def start_pylori_amr_analysis():
    per_gene_dir = "output/genes/H.pylori.J99.2001/per_gene/"
    detected_mutations = {}

    # 23S - Chlarithromycine
    seqs_23S_1 = _fetch_all_gene(per_gene_dir + "23S_rRNA_1.fa")
    seqs_23S_1 = _rename_seqs_no_gene(seqs_23S_1, "23S_rRNA_1")
    for seq in seqs_23S_1:
        if seq.id not in detected_mutations:
            detected_mutations[seq.id] = {}
        detected_mutations[seq.id]["23S_1"] = _variant_detection(
            seq.reverse_complement(),
            _parse_amr_tab(common.get_mutations("pylori", "23S")),
            "nuc",
        )
    seqs_23S_2 = _fetch_all_gene(per_gene_dir + "23S_rRNA_2.fa")
    seqs_23S_2 = _rename_seqs_no_gene(seqs_23S_2, "23S_rRNA_2")
    for seq in seqs_23S_2:
        if seq.id not in detected_mutations:
            detected_mutations[seq.id] = {}
        detected_mutations[seq.id]["23S_2"] = _variant_detection(
            seq.reverse_complement(),
            _parse_amr_tab(common.get_mutations("pylori", "23S")),
            "nuc",
        )

    # 16S - Tétracycline
    seqs_16S_1 = _fetch_all_gene(per_gene_dir + "16S_rRNA_1.fa")
    seqs_16S_1 = _rename_seqs_no_gene(seqs_16S_1, "16S_rRNA_1")
    for seq in seqs_16S_1:
        if seq.id not in detected_mutations:
            detected_mutations[seq.id] = {}
        detected_mutations[seq.id]["16S_1"] = _variant_detection(
            seq.reverse_complement(),
            _parse_amr_tab(common.get_mutations("pylori", "16S")),
            "nuc",
        )
    seqs_16S_2 = _fetch_all_gene(per_gene_dir + "16S_rRNA_2.fa")
    seqs_16S_2 = _rename_seqs_no_gene(seqs_16S_2, "16S_rRNA_2")
    for seq in seqs_16S_2:
        if seq.id not in detected_mutations:
            detected_mutations[seq.id] = {}
        detected_mutations[seq.id]["16S_2"] = _variant_detection(
            seq.reverse_complement(),
            _parse_amr_tab(common.get_mutations("pylori", "16S")),
            "nuc",
        )

    # gyrA - Lévofloxacine
    seqs_gyrA = _fetch_all_gene(per_gene_dir + "gyrA.fa")
    seqs_gyrA = _rename_seqs_no_gene(seqs_gyrA, "gyrA")
    for seq in seqs_gyrA:
        if seq.id not in detected_mutations:
            detected_mutations[seq.id] = {}
        detected_mutations[seq.id]["gyrA"] = _variant_detection(
            seq.translate(),
            _parse_amr_tab(common.get_mutations("pylori", "gyrA")),
            "prot",
        )

    # rpoB - Rifampicine
    # seqs_rpoB = _fetch_all_gene(per_gene_dir + "rpoB.fa")
    # seqs_rpoB = _rename_seqs_no_gene(seqs_16S, "rpoB")
    # for seq in seqs_rpoB:
    #     if seq.id not in detected_mutations:
    #         detected_mutations[seq.id] = {}
    #     detected_mutations[seq.id]["rpoB"] = _variant_detection(
    #         seq.reverse_complement().translate(), _parse_amr_tab(common.get_mutations("pylori", "rpoB")), 'prot'
    #     )

    _export_mutations_to_dataframe(detected_mutations, "pylori")
    return detected_mutations

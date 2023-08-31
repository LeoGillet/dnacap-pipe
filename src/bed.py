# -*- coding: utf-8 -*-
"""
This module contains all functions related to the analysis of mapping results.
"""
import os
import subprocess
import threading
import time

from Bio import SeqIO
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from tqdm import tqdm

from src import common
from src import dataframe as dfio


matplotlib.use("agg")
# matplotlib.use("qtagg")


def _bedtools_coverage_command(bam_file, target_regions_bed_file, output_txt_file):
    return (
        f"{common.get_tool_path('bedtools')} coverage -abam {bam_file} "
        + f"-b {target_regions_bed_file} -mean > {output_txt_file}"
    )


def _mosdepth_command(bam_file, target_regions_bed_file, output_txt_file):
    return (
        f"{common.get_tool_path('mosdepth')} -b {target_regions_bed_file} "
        + f"--no-per-base {output_txt_file} {bam_file}"
    )


def _genomecov_command():
    ...


def _run_coverage(bam_file, target_regions_bed_file, output_txt_file):
    subprocess.run(
        _bedtools_coverage_command(bam_file, target_regions_bed_file, output_txt_file),
        shell=True,
        check=True,
        capture_output=True,
    )


def start_coverage_analysis(bam_files, genome_on_probe):
    """
    Creates subprocess threads for the execution of bedtools coverage
    on every genome and for each input sequence
    """
    bedtools_threads = []
    with tqdm(total=len(bam_files) + 1) as pbar:
        for file_path, basename, genome in bam_files:
            if genome not in genome_on_probe:
                continue
            pbar.set_description(
                f"Running bedtools coverage for sequence {basename} on {genome}"
            )
            if not os.path.isdir(f"output/coverage/{genome}"):
                os.makedirs(f"output/coverage/{genome}")
            output_txt_file = f"output/coverage/{genome}/{basename}.txt"
            thread = threading.Thread(
                target=_run_coverage,
                args=(file_path, genome_on_probe[genome], output_txt_file),
                daemon=True,
            )
            bedtools_threads.append(thread)
            while threading.active_count() >= common.get_max_threads("coverage"):
                pbar.refresh()
                time.sleep(1)
            thread.start()
            pbar.update()
        pbar.close()


def _depth_command(bam_file, output_txt_file):
    return (
        f"{common.get_tool_path('samtools')} depth -a {bam_file} "
        + f"> {output_txt_file}"
    )


def _run_depth(bam_file, output_txt_file):
    subprocess.run(
        _depth_command(bam_file, output_txt_file),
        shell=True,
        check=True,
        capture_output=True,
    )


def start_depth_analysis(bam_files):
    """
    Creates subprocess threads for the execution of samtools depth
    on every genome mapping
    """
    bedtools_threads = []
    with tqdm(total=len(bam_files)) as pbar:
        for file_path, basename, genome in bam_files:
            pbar.set_description(
                f"Running samtools depth for sequence {basename} on {genome}"
            )
            if not os.path.isdir(f"output/depth/{genome}"):
                os.makedirs(f"output/depth/{genome}")
            output_txt_file = f"output/depth/{genome}/{basename}.txt"
            thread = threading.Thread(
                target=_run_depth,
                args=(file_path, output_txt_file),
                daemon=True,
            )
            bedtools_threads.append(thread)
            while threading.active_count() >= common.get_max_threads("depth"):
                pbar.refresh()
                time.sleep(1)
            thread.start()
            pbar.update()
        for thread in bedtools_threads:
            thread.join()
        pbar.close()


def _txt_to_array(depth_file):
    x_val = []
    y_val = []
    try:
        with open(depth_file, "r", encoding="UTF-8") as depths:
            for row in depths:
                _, x_str, y_str = row.strip().split("\t")
                x_val.append(int(x_str))
                y_val.append(int(y_str))
    except UnicodeDecodeError:
        print("UnicodeDecodeError occurred on file", depth_file)
    return np.array(x_val), np.array(y_val)


def _plot_depth(x_val, y_val, regions, export_):
    fig_name = export_.replace("/", "_")
    fig = plt.figure(fig_name, figsize=(20, 20), dpi=300)
    axis = fig.add_subplot(111)
    axis.scatter(x_val, y_val, s=0.5, color="#0092bc")
    max_y = max(y_val)
    for region in regions:
        x_dot = int((region["start_pos"] + region["end_pos"]) / 2)
        value_at_pos = max(y_val[region["start_pos"] - 1 : region["end_pos"] - 1])
        marker = "g+"
        if value_at_pos < 30.0:
            marker = "y+"
        if value_at_pos < 10.0:
            marker = "r+"
        axis.plot(x_dot, value_at_pos, marker, markersize=15)
        axis.annotate(
            region["name"],
            (x_dot, value_at_pos),
            textcoords="offset points",
            xytext=(0, 20),
            ha="center",
            fontsize="xx-large",
        )
    plt.xticks(fontsize="25")
    plt.yticks(fontsize="25")
    plt.ylim([-1, max_y])
    plt.xlabel("Position (pb)", fontsize=30, labelpad=40)
    plt.ylabel("Depth at position", fontsize=30, labelpad=40)
    fig.savefig(export_)
    plt.close(fig_name)


def depth_to_dataframe(depth_files):
    local_df = dfio.get_or_create(
        "pylori" if "pylori" in depth_files[0][0].lower() else "genitalium"
    )
    with tqdm(total=len(depth_files), desc="Saving to dataframe") as pbar:
        for file_, basename, _, regions in depth_files:
            for region in regions:
                _, y_val = _txt_to_array(file_)
                subarray = y_val[region["start_pos"] - 1 : region["end_pos"] - 1]
                local_df.at[basename, region["name"] + " mean depth"] = sum(
                    subarray
                ) / len(subarray)
            pbar.update()
        pbar.close()
    dfio.save(local_df)


def _plot_thread(file_, basename, genome, regions, pbar):
    pbar.set_description(f"Reading and parsing depth of {basename} on {genome}")
    x_val, y_val = _txt_to_array(file_)
    pbar.set_description(f"Plotting depth of {basename} on {genome}")
    _plot_depth(x_val, y_val, regions, f"output/depth/{genome}/{basename}.png")


def bed_to_list(bed_file: str) -> list[dict[str, str, str, str]]:
    regions = []
    if not os.path.isfile(bed_file):
        raise FileNotFoundError("BED file", bed_file, "does not exist. Stopping...")
    with open(bed_file, "r", encoding="UTF-8") as bedf:
        file_entries = bedf.readlines()
    for entry in file_entries:
        try:
            scaffold, start_pos, end_pos, name, *_ = entry.strip().split("\t")
        except ValueError:
            if not entry == "\n":
                raise ValueError(
                    f"Couldn't unpack entry {entry} of BED file {bed_file}"
                )
        try:
            start_pos = int(start_pos)
            end_pos = int(end_pos)
        except ValueError as exc:
            raise ValueError(
                f"{exc}\nInvalid positions in entry {name} of BED file {bed_file}"
            )
        regions.append(
            {
                "scaffold": scaffold,
                "start_pos": start_pos,
                "end_pos": end_pos,
                "name": name,
            }
        )
    return regions


def generate_depth_plots(depth_files):
    """
    Generates PNG plots of depth per base for every mapping
    """
    plot_threads = []
    with tqdm(total=len(depth_files)) as pbar:
        for file_, basename, genome, regions in depth_files:
            thread = threading.Thread(
                target=_plot_thread,
                args=(file_, basename, genome, regions, pbar),
                daemon=True,
            )
            plot_threads.append(thread)
            pbar.set_description(f"Plotting depth of {basename} on {genome}")
            while threading.active_count() >= common.get_max_threads("plotting"):
                pbar.refresh()
                time.sleep(1)
            thread.start()
            pbar.update()
        pbar.close()
    for thread in plot_threads:
        thread.join()


def _create_roi_dirs(basename, genome):
    output_dir = f"output/genes/{genome}/{basename}"
    try:
        os.makedirs(output_dir)
    except FileExistsError:
        pass
    return output_dir


def _extract_whole_genome(bam_file, reference_file, output_dir):
    command = (
        f"{common.get_tool_path('samtools')} view -h {bam_file} | "
        + f"{common.get_tool_path('bcftools')} mpileup -d 0 -f {reference_file} - > {bam_file}_pileup.txt && "
        + f"{common.get_tool_path('bcftools')} call -c {bam_file}_pileup.txt | "
        + f"{common.get_tool_path('vcfutils')} vcf2fq > {output_dir}/genome.fastq"
    )
    subprocess.run(command, shell=True, check=True, capture_output=True)
    fasta_command = (
        f"{common.get_tool_path('seqtk')} seq -a -q30 -n N {output_dir}/genome.fastq > "
        + f"{output_dir}/genome.fasta"
    )
    subprocess.run(fasta_command, shell=True, check=True, capture_output=True)
    return f"{output_dir}/genome.fasta"


def _load_names(names_file):
    with open(names_file, "r", encoding="UTF-8") as namesf:
        raw_names = namesf.readlines()
    return [entry.strip() for entry in raw_names]


def _extract_roi(assembly_fasta, bed_file, names_file, output_dir):
    final_file = (
        f"{output_dir}genes.fa" if names_file else f"{output_dir}genes_unnamed.fa"
    )
    command = (
        f"{common.get_tool_path('seqtk')} subseq -l80 {assembly_fasta} {bed_file} > "
        + f"{output_dir}genes_unnamed.fa"
    )
    subprocess.run(command, shell=True, check=True)
    if names_file:
        command = (
            f"{common.get_tool_path('seqtk')} rename {output_dir}genes_unnamed.fa ROI_ > "
            + final_file
        )
        subprocess.run(command, shell=True, check=True)
        os.remove(f"{output_dir}genes_unnamed.fa")
        names = _load_names(names_file)
        with open(final_file, "r", encoding="UTF-8") as asf:
            data = asf.read()
            for i, name in enumerate(names):
                data = data.replace(f"ROI_{i+1}\n", name + "\n")
        with open(final_file, "w", encoding="UTF-8") as asf:
            asf.write(data)
    return final_file


def start_genome_extraction(bam_files: list[tuple[str, str, str, str]]):
    """
    Extracts all genomes from mapping results using
    samtools view, bcftools mpileup & call, vcfutils vcf2fq and seqtk seq
    with a base score minimum of 20
    """
    threads = []
    with tqdm(total=len(bam_files)) as pbar:
        for bam_file, basename, genome, reference_path in bam_files:
            output_dir = _create_roi_dirs(basename, genome)
            thread = threading.Thread(
                target=_extract_whole_genome,
                args=(
                    bam_file,
                    reference_path,
                    output_dir,
                ),
                daemon=True,
            )
            pbar.set_description(f"Extracting consensus genome from {basename}")
            threads.append(thread)
            while threading.active_count() >= common.get_max_threads(
                "genome_extraction"
            ):
                pbar.refresh()
                time.sleep(1)
            thread.start()
            pbar.update()
        pbar.close()

    for thread in threads:
        thread.join()


def _load_gene(file_, name):
    with open(file_, "r", encoding="UTF-8") as genef:
        data = genef.readlines()
    for i, line in enumerate(data):
        if name in line.strip():
            return data[i + 1].strip()
    return


def _get_genes_per_sample(genome, pipeline_type):
    sample_genes_dir = f"output/genes/{genome}/"
    samples = [
        sample
        for sample in os.listdir(sample_genes_dir)
        if os.path.isdir(sample_genes_dir + sample) and sample != 'per_gene'
    ]
    gene_names = _load_names(common.get_genome_regions(pipeline_type) + ".names.list")
    genes = {}
    for sample in samples:
        dir_ = sample_genes_dir + sample + "/"
        genes_file = dir_ + "genes.fa"
        for gene_name in gene_names:
            gene = _load_gene(genes_file, gene_name)
            if not gene:
                continue
            if gene_name not in genes:
                genes[gene_name] = {}
            genes[gene_name][sample.replace("\n", "")] = gene
    return genes


def write_genes_fasta(genome, pipeline_type):
    genes = _get_genes_per_sample(genome, pipeline_type)
    output_dir = f"output/genes/{genome}/per_gene/"
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)
    for gene_name, sequence in genes.items():
        lines = []
        for sample, gene_seq in sequence.items():
            lines.extend([">" + gene_name + "|" + sample + "\n", gene_seq + "\n"])
        with open(output_dir + gene_name + ".fa", "w", encoding="UTF-8") as genef:
            genef.writelines(lines)


def start_genes_grouping(pipeline_type):
    for genome in os.listdir("output/genes/"):
        if genome in common.get_ignored_genomes(pipeline_type) or not os.path.isdir(
            f"output/genes/{genome}"
        ):
            continue
        write_genes_fasta(genome, pipeline_type)


def start_genes_extractions(genes_files):
    """
    Extracts regions from genomes previously extracted
    """
    threads = []
    with tqdm(total=len(genes_files)) as pbar:
        for (
            assembly_fasta,
            bed_file,
            names_file,
            output_dir,
            basename,
        ) in genes_files:
            thread = threading.Thread(
                target=_extract_roi,
                args=(
                    assembly_fasta,
                    bed_file,
                    names_file,
                    output_dir,
                ),
                daemon=True,
            )
            pbar.set_description(f"Extracting regions of interest from {basename}")
            threads.append(thread)
            while threading.active_count() >= common.get_max_threads(
                "genome_extraction"
            ):
                pbar.refresh()
                time.sleep(1)
            thread.start()
            pbar.update()
        pbar.close()

    for thread in threads:
        thread.join()


def start_per_gene_coverage(pipeline_type):
    if pipeline_type == "pylori":
        per_gene_dir = "output/genes/H.pylori.J99.2001/per_gene/"
    else:
        per_gene_dir = "output/genes/M.genitalium.G37/per_gene/"

    dataframe = dfio.get_or_create(pipeline_type)
    for file_ in os.listdir(per_gene_dir):
        if not file_.endswith(".fa"):
            continue
        sequences = []
        for record in SeqIO.parse(per_gene_dir + file_, "fasta"):
            sequences.append(record)
        print(len(sequences))
        for seq in sequences:
            if '|' not in seq.id:
                continue
            gene_name, sample_name = seq.id.split('|')
            total_len = len(seq.seq)
            if total_len == 0:
                dataframe.at[sample_name, gene_name+' coverage'] = -1.0
                continue
            n_count = seq.count('N')
            dataframe.at[sample_name, gene_name+' coverage'] = (total_len-n_count)/total_len
    dfio.save(dataframe)
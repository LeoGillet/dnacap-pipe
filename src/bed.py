# -*- coding: utf-8 -*-
"""
This module contains all functions related to the analysis of mapping results.
"""
import os
import subprocess
import threading
import time

import matplotlib
import matplotlib.pyplot as plt
from tqdm import tqdm

from src import common


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
        f"{common.get_tool_path('samtools')} depth {bam_file} " + f"> {output_txt_file}"
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
    with open(depth_file, "r", encoding="UTF-8") as depths:
        for row in depths:
            _, x_str, y_str = row.strip().split("\t")
            x_val.append(int(x_str))
            y_val.append(int(y_str))
    return x_val, y_val


def _plot_depth(x_val, y_val, export_):
    fig_name = export_.replace("/", "_")
    fig = plt.figure(figsize=(20, 20), dpi=300)
    axis = fig.add_subplot(111)
    axis.scatter(x_val, y_val, s=0.5, color="#0092bc")
    plt.xticks(fontsize="25")
    plt.yticks(fontsize="25")
    plt.ylim([-1, max(y_val)])
    plt.xlabel("Position (pb)", fontsize=30, labelpad=40)
    plt.ylabel("Depth at position", fontsize=30, labelpad=40)
    fig.savefig(export_)
    plt.close(fig_name)


def _plot_thread(file_, basename, genome, pbar):
    pbar.set_description(f"Reading and parsing depth of {basename} on {genome}")
    x_val, y_val = _txt_to_array(file_)
    pbar.set_description(f"Plotting depth of {basename} on {genome}")
    _plot_depth(x_val, y_val, f"output/depth/{genome}/{basename}.png")


def generate_depth_plots(depth_files):
    """
    Generates PNG plots of depth per base for every mapping
    """
    plot_threads = []
    with tqdm(total=len(depth_files)) as pbar:
        for file_, basename, genome in depth_files:
            thread = threading.Thread(
                target=_plot_thread,
                args=(file_, basename, genome, pbar),
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
        + f"{common.get_tool_path('bcftools')} mpileup -d 0 -f {reference_file} - | "
        + f"{common.get_tool_path('bcftools')} call -c | "
        + f"{common.get_tool_path('vcfutils')} vcf2fq > {output_dir}/genome.fastq && "
        + f"{common.get_tool_path('seqtk')} seq -a -q20 -n N {output_dir}/genome.fastq > "
        + f"{output_dir}/genome.fasta"
    )
    subprocess.run(command, shell=True, check=True, capture_output=True)
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
        f"{common.get_tool_path('seqtk')} subseq {assembly_fasta} {bed_file} > "
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


def start_genes_extractions(genes_files):
    """
    Extracts regions from genomes previously extracted
    """
    threads = []
    with tqdm(total=len(genes_files)) as pbar:
        for assembly_fasta, bed_file, names_file, output_dir, basename in genes_files:
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

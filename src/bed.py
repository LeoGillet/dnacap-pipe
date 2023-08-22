# -*- coding: utf-8 -*-
"""
This module contains all functions related to the analysis of mapping results.
"""
import os
import matplotlib

matplotlib.use("qtagg")
import matplotlib.pyplot as plt
import numpy as np
import subprocess
import threading
import time
from tqdm import tqdm

from src import common


def _coverage_command(bam_file, target_regions_bed_file, output_txt_file):
    print(target_regions_bed_file)
    return (
        f"{common.get_tool_path('bedtools')} coverage -hist  -abam {bam_file} "
        + f"-b {target_regions_bed_file} > {output_txt_file}"
    )


def _run_coverage(bam_file, target_regions_bed_file, output_txt_file):
    subprocess.run(
        _coverage_command(bam_file, target_regions_bed_file, output_txt_file),
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
            output_txt_file = f"output/coverage/{genome}/{basename}_hist.all.txt"
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
            _, x, y = row.strip().split("\t")
            x_val.append(x)
            y_val.append(y)
    return x_val, y_val


def _plot_depth(x_val, y_val, export_):
    fig_name = export_.replace("/", "_")
    fig = plt.figure(figsize=(20, 20), dpi=300)
    print("created figure")
    ax = fig.add_subplot(111)
    ax.scatter(x_val, y_val, s=0.5, color="#0092bc")
    print("added scatter")
    plt.ylim([-1, max(y_val)])
    plt.xlabel("Depth at position")
    plt.ylabel("Position (pb)")
    print("added legends")

    print(export_)
    fig.savefig(export_)
    print("saved figure")
    plt.close("all")
    print("closed figure")


def generate_depth_plots(depth_files):
    with tqdm(total=len(depth_files)) as pbar:
        for file_, basename, genome in depth_files:
            start = time.time()
            pbar.set_description(f"Reading and parsing depth of {basename} on {genome}")
            x_val, y_val = _txt_to_array(file_)
            checkpoint = time.time()
            pbar.write(
                f"Parsed depth of {basename} on {genome} in {checkpoint-start:.2f}s"
            )
            pbar.set_description(f"Plotting depth of {basename} on {genome}")
            _plot_depth(x_val, y_val, f"output/depth/{genome}/{basename}.png")
            end = time.time()
            pbar.write(
                f"Plotted depth of {basename} on {genome} in {end-checkpoint:.2f}s"
            )
            pbar.write(f"Total time of {basename} on {genome} was {end-start:.2f}s")
            pbar.update()
        pbar.close()

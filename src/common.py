# -*- coding: utf-8 -*-
"""
This module stores any miscellaneous function that is used throughout multiple modules
"""
import json
import os
import platform

from src.logger import log
from src import bed


def clear_stdout():
    """
    Clears terminal
    """
    if platform.system() == "Windows":
        os.system("cls")
    else:
        os.system("clear")


def get_tool_path(tool: str):
    """
    Reads config.json file to fetch tool paths
    :param tool: name of the tool - key in 'paths' object of config.json
    :return: string containing path of the tool
    """
    with open("config.json", "r", encoding="UTF-8") as config_file:
        config = json.load(config_file)
    if tool not in config["paths"]:
        log(f"Tool {tool} not found in config file. Stopping", level="ERROR")
        raise KeyError(f"Tool {tool} not found in config file.")
    return config["paths"][tool]


def get_genome_path(genome_key, tool) -> str:
    """
    Reads config.json file to fetch genome paths
    :param genome_key: key corresponding to specific genome
    :param tool: key corresponding to mapping tool used to index genome
    """
    with open("config.json", "r", encoding="UTF-8") as config_file:
        config = json.load(config_file)
    if tool not in ("bwa", "bowtie", "no_tool"):
        log(f"Tool {tool} not supported.", level="ERROR")
        raise KeyError(f"Tool {tool} not supported. Supported mappers: 'bwa', 'bowtie'")
    suffix = ""
    folder = "refs/"
    if tool == "bowtie":
        folder += "bowtie/"
        suffix = ".indexed"
    elif tool == "bwa":
        folder += "bwa/"

    if genome_key not in config["genomes"]:
        log(
            f"Genome {genome_key} path not found in config file",
            level="ERROR",
        )
        raise KeyError(f"Genome {genome_key} path not found in config file")
    return folder + config["genomes"][genome_key] + suffix


def get_max_threads(operation):
    """
    Reads config.json file to fetch max threads for supported operations
    :param operation: specified operation ('mapping', 'fastqc', 'cleaning', ...)
    """
    with open("config.json", "r", encoding="UTF-8") as config_file:
        config = json.load(config_file)
    if operation not in config["threads"]:
        log(
            f"Operation {operation} unknown when searching for max threads",
            level="ERROR",
        )
        raise KeyError(f"Operation {operation} unknown when searching for max threads")
    if isinstance(config["threads"][operation], int):
        return config["threads"][operation]
    log(
        f"Max threads for operation {operation} in config file has an invalid format",
        level="ERROR",
    )
    raise TypeError(
        f"Max threads for operation {operation} in config file has an invalid format"
    )


def choose_genomes(pipeline_type) -> set:
    """
    Imports genomes selected for the mapping of input sequences from config file
    """
    with open("config.json", "r", encoding="UTF-8") as config_file:
        config = json.load(config_file)
    return set(config[pipeline_type]["reference_genomes"])


def get_ignored_genomes(pipeline_type) -> set:
    """
    Returns set of genomes that should be ignored for depth and genes extraction
    """
    with open("config.json", "r", encoding="UTF-8") as config_file:
        config = json.load(config_file)
    return config[pipeline_type]["ignored_genomes"]


def get_bam_files() -> list[tuple[str, str, str]]:
    """
    Returns a list of all constructed assemblies detected in output folder
    """
    bam_files = []
    mapping_dir = "output/mapping/"
    ignore_dirs = [
        mapping_dir + subfolder
        for subfolder in (
            "human_unmapped",
            "human_mapped",
            "human_flagstat",
            "human_bedcov",
        )
    ]

    dirs_in_folder = [
        mapping_dir + dir_ + "/"
        for dir_ in os.listdir(mapping_dir)
        if mapping_dir + dir_ not in ignore_dirs
    ]

    for genome_dir in dirs_in_folder:
        files_in_dir = os.listdir(genome_dir)
        for file_ in files_in_dir:
            if not file_.endswith(".bam") or "_unsorted" in file_:
                continue
            bam_files.append(
                (
                    genome_dir + file_,
                    file_.replace(".bam", ""),
                    genome_dir.split("/")[2],
                )
            )
    return bam_files


def get_probes_bed_files() -> list:
    """
    Returns a list of all constructed assemblies detected in probes folder
    """
    probes = []
    probes_folder = "probes/"
    for file_ in os.listdir(probes_folder):
        if not file_.endswith(".bed"):
            continue
        probes.append(probes_folder + file_)
    return probes


def get_genome_regions(pipeline_type) -> dict:
    """
    Returns dictionary of probes to search for depending on the genome
    the input sequence has been mapped to
    """
    with open("config.json", "r", encoding="UTF-8") as config_file:
        config = json.load(config_file)
    return config[pipeline_type]["regions"]


def get_depth_files(pipeline_type) -> list[tuple[str, str, str]]:
    """
    Gets all depth files created by previous steps and returns ones that shouldn't be ignored
    """
    depth_directory = "output/depth/"
    depth_files = []
    ignore_genomes = get_ignored_genomes(pipeline_type)
    for genome in os.listdir(depth_directory):
        if genome in ignore_genomes:
            continue
        for file_ in os.listdir(depth_directory + genome):
            if not file_.endswith(".txt"):
                continue
            depth_files.append(
                (
                    depth_directory + genome + "/" + file_,
                    file_.replace(".txt", ""),
                    genome,
                    bed.bed_to_list(get_genome_regions(pipeline_type)),
                )
            )
    return depth_files


def get_bam_files_extraction(ignore_genomes) -> list[tuple[str, str, str, str]]:
    """
    Gets all file required for the extraction of the genome
    """
    mapping_directory = "output/mapping/"
    bam_files = []
    for genome in os.listdir(mapping_directory):
        if genome in ignore_genomes:
            continue
        for file_ in os.listdir(mapping_directory + genome):
            if not file_.endswith(".bam") or file_.endswith("_unsorted.bam"):
                continue
            bam_files.append(
                (
                    mapping_directory + genome + "/" + file_,
                    file_.replace(".bam", ""),
                    genome,
                    get_genome_path(genome, "no_tool"),
                )
            )
    return bam_files


def get_genes_files(pipeline_type):
    """
    Gets all files required for the extraction of regions of interest
    """
    genomes_dir = "output/genes/"
    genomes = []
    for genome in os.listdir(genomes_dir):
        if genome in get_ignored_genomes(pipeline_type) or not os.path.isdir(genomes_dir + genome):
            continue
        for basename in os.listdir(genomes_dir + genome):
            if not os.path.isdir(genomes_dir + genome + "/" + basename):
                continue
            if os.path.isfile(
                genomes_dir + genome + "/" + basename + "/" + "genome.fasta"
            ):
                bed_names = None
                if os.path.isfile(get_genome_regions(pipeline_type) + ".names.list"):
                    bed_names = get_genome_regions(pipeline_type) + ".names.list"
                genomes.append(
                    (
                        genomes_dir + genome + "/" + basename + "/" + "genome.fasta",
                        get_genome_regions(pipeline_type),
                        bed_names,
                        genomes_dir + genome + "/" + basename + "/",
                        basename,
                    )
                )
    return genomes


def get_adapters_path(adapters):
    with open("config.json", "r", encoding="UTF-8") as config_file:
        config = json.load(config_file)
    return config["adapters"][adapters]


def get_mutations(pipeline_type, gene):
    with open("config.json", "r", encoding="UTF-8") as config_file:
        config = json.load(config_file)
    return config[pipeline_type]["mutations"][gene]

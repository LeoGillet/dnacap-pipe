# -*- coding: utf-8 -*-
"""
This module stores any miscellaneous function that is used throughout multiple modules
"""
import json
import os
import platform

from src.logger import log


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
    if tool not in ("bwa", "bowtie"):
        log(f"Tool {tool} not supported.", level="ERROR")
        raise KeyError(f"Tool {tool} not supported. Supported mappers: 'bwa', 'bowtie'")
    tool = "refs/bowtie/" if tool == "bowtie" else "refs/bwa/"
    suffix = ".indexed" if tool == "bowtie" else ""
    if genome_key not in config["genomes"]:
        log(
            f"Genome {genome_key} path not found in config file",
            level="ERROR",
        )
        raise KeyError(
            f"Genome {genome_key} path not found in config file"
        )
    return tool + config["genomes"][genome_key] + suffix


def get_max_threads(operation):
    """
    Reads config.json file to fetch max threads for supported operations
    :param operation: specified operation ('mapping', 'fastqc', 'sickle', ...)
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


def choose_genomes() -> set:
    """
    Imports genomes selected for the mapping of input sequences from config file
    """
    with open("config.json", "r", encoding="UTF-8") as config_file:
        config = json.load(config_file)
    return set(config["globals"]["loaded_genomes"])


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


def get_genome_regions() -> dict:
    """
    Returns dictionary of probes to search for depending on the genome
    the input sequence has been mapped to
    """
    with open("config.json", "r", encoding="UTF-8") as config_file:
        config = json.load(config_file)
    return config["globals"]["regions"]


def get_depth_files() -> list[tuple[str, str, str]]:
    depth_directory = "output/depth/"
    depth_files = []
    for genome in os.listdir(depth_directory):
        for file_ in os.listdir(depth_directory + genome):
            depth_files.append(
                (
                    depth_directory + genome + "/" + file_,
                    file_.replace(".txt", ""),
                    genome,
                )
            )
    return depth_files
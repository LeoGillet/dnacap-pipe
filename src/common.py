# -*- coding: utf-8 -*-
"""
This module stores any miscellaneous function that is used throughout multiple modules
"""
import json

from src.logger import log


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
        if tool not in config["genomes"]:
            log(f"Tool {tool} unknown when finding genome path", level="ERROR")
            raise KeyError(f"Tool {tool} unknown when finding genome path")
        if genome_key not in config["genomes"][tool]:
            log(
                f"Error occurred when finding path of genome {genome_key} with tool {tool}.",
                level="ERROR",
            )
            raise KeyError(
                f"Error occurred when finding path of genome {genome_key} with tool {tool}."
            )
        return config["genomes"][tool][genome_key]


def choose_genomes() -> set:
    """
    Imports genomes selected for the mapping of input sequences from config file
    """
    with open("config.json", "r", encoding="UTF-8") as config_file:
        config = json.load(config_file)
    return set(config["globals"]["loaded_genomes"])

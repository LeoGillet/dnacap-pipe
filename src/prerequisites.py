# -*- coding: utf-8 -*-

import glob
import json
import os
import sys
import shutil
from datetime import datetime

from src.logger import *


# ---------------------------------------------

def ensure_paths():
    with open('config.json', 'r') as config_j:
        config = json.load(config_j)
        not_found_tools = []
        for tool in config['paths'].values():
            if not os.path.isfile(tool):
                not_found_tools.append(tool)
                log("Tool {} path invalid.".format(tool), level='ERROR')
    if not_found_tools:
        raise FileNotFoundError(
            "Tool(s) {} were not found using paths located in config file. Stopping.".format(not_found_tools))
    else:
        log("All tools were found using paths in config file.")

# ---------------------------------------------

def rename_old_log():
    """
    Ensures 'log' directory exists and renames previous 'latest.log' to dated file name
    """
    if not os.path.isdir('log'):
        os.mkdir('log')
    if os.path.isfile('log/latest.log'):
        created_seconds = os.path.getctime('log/latest.log')
        created_timestamp = datetime.fromtimestamp(created_seconds).strftime("%Y%m%d_%H%M%S")
        os.rename('log/latest.log', f'log/log_{created_timestamp}.log')
        log("Previous log file renamed to {}".format(f'log_{created_timestamp}'))


def delete_older_logs():
    """
    Limits number of log files present in 'log' directory (max 10 archived logs) by deleting older files
    """
    logs = list(filter(os.path.isfile, glob.glob('log/*')))
    if len(logs) < 11: return
    logs.sort(key=os.path.getmtime)
    logs_to_delete = logs[:len(logs) - 10]
    log("Deleting older log files ({})".format(logs_to_delete))
    for logfile in logs_to_delete:
        os.remove(logfile)


# ---------------------------------------------

def _check_current_wd():
    """
    Ensures the program is running from its root directory
    """
    if 'main.py' not in os.listdir('.'):
        log("Main script was not executed from its root folder. Stopping.", level='ERROR')
        raise Exception("This script needs to be executed from its root folder (main.py).")


def prepare_dirs():
    """
    Ensures input & output folders are created before continuing
    """
    _check_current_wd()
    if not os.path.isdir('input'):
        os.mkdir('input')
        log("Created input folder")
    if not os.path.isdir('output'):
        os.mkdir('output')
        log("Created output folder")


# ---------------------------------------------

def _check_dir_contents(dir, filter_extension=None) -> bool:
    """
    Ensures folder contains a file. If 'filter_extension' is specified, ensures files with extensions are contained in directory
    :param dir: directory to scan
    :param filter_extension: extension to filter files with
    :return: boolean
    """
    if not filter_extension: return len(os.listdir(dir)) > 0
    return any(filter_extension in filename for filename in os.listdir(dir))


def _is_input_present() -> bool:
    """
    Ensures supported sequence file formats are found in input directory
    """
    return any(_check_dir_contents('input', filter_extension=ext) for ext in ('fastq.gz', 'fq.gz'))


def _is_output_empty():
    """
    Checks whether the output directory is empty before populating
    """
    return not _check_dir_contents('output')


def _archive_old_output():
    """
    Archives non-empty output folder contents in ZIP file
    """
    archive_date = datetime.fromtimestamp(os.path.getctime('output/' + os.listdir('output')[0])).strftime(
        "%Y%m%d_%H%M%S")
    shutil.make_archive(f'output_{archive_date}', 'zip', 'output')
    log("Archived previous results in file {}".format(f'output_{archive_date}.zip'), level='WARN')


def _empty_output():
    """
    Deletes files in output directory
    """
    os.system("cd output; rm -rf *")
    log("Erased previous results", level='WARN')


def _ask_overwrite_behaviour() -> bool:
    """
    Asks user for decision regarding previous output folder archiving or deletion
    """
    choice = input(str_warn(
        "Output folder is not empty. (A)rchive previous results, (D)elete previous results, (C)ancel? [A/d/c] > "))
    while choice not in ('', 'A', 'a', 'D', 'd', 'C', 'c'):
        choice = input(str_warn("Invalid user input. [A/d/c] > "))
    match choice:
        case '' | 'A' | 'a':
            _archive_old_output()
            _empty_output()
            return False
        case 'D' | 'd':
            _empty_output()
            return False
        case 'C' | 'c':
            return True
        case _:
            raise NotImplementedError


def check_folders():
    """
    Pre-processing input and output folders checks routine
    """
    if not _is_input_present():
        log("No file with <fastq.gz or fq.gz> extensions found in <input> folder. Exiting.", level='ERROR')
        raise Exception(
            str_err("No file with <fastq.gz or fq.gz> extensions found in <input> folder. Exiting."))
    while not _is_output_empty():
        if _ask_overwrite_behaviour(): sys.exit(1)
    os.makedirs("output/intermediate_files/sickled")
    os.makedirs("output/denovo")
    os.makedirs("output/mapping")


# ---------------------------------------------

def ask_sequencer_type():
    """
    Asks user for type of sequences outputted by sequencer
    :return: format
    """
    choice = input(
        str_info("Which sequencer type was used?\n\t1. Short-read Illumina sequencing (default)\n? "))
    while choice not in ('', '1'):
        choice = input(str_warn("Invalid user input. [1] "))
    match choice:
        case '' | '1':
            log("Short-read Illumina sequencer type selected by User.")
            return 'illumina'
        case _:
            raise NotImplementedError("Unknown sequencer type. Exiting.")

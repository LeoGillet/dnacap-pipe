# -*- coding: utf-8 -*-
"""
This module contains multiple functions related to the Quality Control of the input sequences
    - FastQC functions generate QC reports on all input sequences
    - Sickle trims sequences based on base quality
    - MultiQC scans output folder for logs and reports in order to create a run summary
"""
import subprocess
import threading
import time

from src import common
from src.logger import log, str_success


# ---------------------------------------------
#   FastQC
# ---------------------------------------------

def _exec_fastqc(fastqfile):
    completed = subprocess.run(
        f"{common.get_tool_path('fastqc')} {fastqfile} -o output/intermediate_files/fastqc_reports",
        shell=True, capture_output=True, check=True
    )
    if completed.returncode != 0:
        log(f"Error occurred when executing FastQC on sequence {fastqfile}", level='ERROR')
        log(f"Stacktrace : {completed.stderr}", level='DEBUG')
        raise RuntimeError(completed.stderr)
    log(f"Created FastQC report for sequence {fastqfile}")


def run_all_fastqc(sequence_pairs):
    """
    Runs all FastQC subprocesses in individual threads
    """
    fastqc_threads = []
    for i, (basename, (fwd, rev, ext)) in enumerate(sequence_pairs.items()):
        t_fwd = threading.Thread(target=_exec_fastqc, args=('input/' + fwd + ext,))
        t_rev = threading.Thread(target=_exec_fastqc, args=('input/' + rev + ext,))
        while threading.active_count() >= 8:
            log(f"{threading.active_count()} active threads. Sleeping for a second...")
            time.sleep(1)
        log(
            f"[{i + 1}/{len(sequence_pairs.keys())}] " +
            f"Starting FastQC thread for paired sequences {basename}",
            level='WARN'
        )
        print(
            f"[{i + 1}/{len(sequence_pairs.keys())}] " +
            f"Starting FastQC thread for paired sequences {basename}"
        )
        t_fwd.start()
        t_rev.start()
        fastqc_threads.extend([t_fwd, t_rev])

    for thread in fastqc_threads:
        thread.join()
    log(f"Done creating {len(sequence_pairs.keys()) * 2} FastQC reports.", level='SUCCESS')
    print(str_success(
        f"Done creating {len(sequence_pairs.keys()) * 2} FastQC reports."
    ))


# ---------------------------------------------
#   MultiQC
# ---------------------------------------------

def run_multiqc():
    """
    Runs MultiQC subprocess
    """
    completed = subprocess.run(
        "multiqc -c multiqc_config.yaml -o output/multiqc -ip --profile-runtime output",
        shell=True, capture_output=True, check=True)
    if completed.returncode != 0:
        log("Error occurred when creating MultiQC report", level='ERROR')
        log(f"Stacktrace : {completed.stderr}", level='DEBUG')
        raise RuntimeError(completed.stderr)
    log("Done creating MultiQC report", level='SUCCESS')
    print(str_success("Created MultiQC report"))


# ---------------------------------------------
#   Sickle
# ---------------------------------------------

def _sickle_command(pe_file1, pe_file2, basename):
    command = f"{common.get_tool_path('sickle')} pe " \
              f"-f input/{pe_file1} -r input/{pe_file2} -t sanger " \
              f"-o output/intermediate_files/sickled_{pe_file1} " \
              f"-p output/intermediate_files/sickled_{pe_file2} " \
              f"-s output/intermediate_files/sickled_singles_{basename}.fastq " \
              f"2>output/logs/sickle_log_{basename}.log"
    log("subprocess run:" + command, level='DEBUG')
    return command


def _run_sickle(forward_sequence: str, reverse_sequence: str, basename: str):
    completed = subprocess.run(
        _sickle_command(forward_sequence, reverse_sequence, basename),
        shell=True, capture_output=True, check=True)
    if completed.returncode != 0:
        log(
            f"Error occurred when sickling sequences {forward_sequence}+{reverse_sequence}",
            level='ERROR'
        )
        raise RuntimeError(completed.stderr)


def batch_sickle(sequence_pairs):
    """
    Runs sickle on all input sequences in individual threads
    """
    sickle_threads = []
    for i, (basename, (fwd, rev, ext)) in enumerate(sequence_pairs.items()):
        thread = threading.Thread(target=_run_sickle, args=(fwd + ext, rev + ext, basename))
        while threading.active_count() >= 8:
            log(f"{threading.active_count()} active threads. Sleeping for a second...")
            time.sleep(1)
        log(f"[{i + 1}/{len(sequence_pairs.keys())}] Sickling {basename}", level='WARN')
        print(f"[{i + 1}/{len(sequence_pairs.keys())}] Sickling {basename}")
        thread.start()

    for thread in sickle_threads:
        thread.join()
    log(f"Done sickling {len(sequence_pairs.keys())} paired sequences", level='SUCCESS')
    print(str_success(f"Done sickling {len(sequence_pairs.keys())} paired sequences"))

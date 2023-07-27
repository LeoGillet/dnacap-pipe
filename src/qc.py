# -*- coding: utf-8 -*-
import subprocess
import threading
import time

from src import common
from src.logger import *


# ---------------------------------------------
#   FastQC
# ---------------------------------------------

def _exec_fastqc(fastq_file):
    completed = subprocess.run(
        f"{common.get_tool_path('fastqc')} {fastq_file} -o output/intermediate_files/fastqc_reports",
        shell=True, capture_output=True
    )
    if completed.returncode != 0:
        log("Error occurred when executing FastQC on sequence {}".format(fastq_file), level='ERROR')
        log("Stacktrace : {}".format(completed.stderr), level='DEBUG')
        raise Exception(completed.stderr)
    else:
        log("Created FastQC report for sequence {}".format(fastq_file))


def run_all_fastqc(sequence_pairs):
    fastqc_threads = []
    for i, (bn, (fwd, rev, ext)) in enumerate(sequence_pairs.items()):
        t_fwd = threading.Thread(target=_exec_fastqc, args=('input/' + fwd + ext,))
        t_rev = threading.Thread(target=_exec_fastqc, args=('input/' + rev + ext,))
        while threading.active_count() >= 8:
            log(f"{threading.active_count()} active threads. Sleeping for a second...")
            time.sleep(1)
        log(f"[{i + 1}/{len(sequence_pairs.keys())}] Starting FastQC thread for paired sequences {bn}", level='WARN')
        print(f"[{i + 1}/{len(sequence_pairs.keys())}] Starting FastQC thread for paired sequences {bn}")
        t_fwd.start()
        t_rev.start()
        fastqc_threads.extend([t_fwd, t_rev])

    for thread in fastqc_threads:
        thread.join()
    log("Done creating {} FastQC reports.".format(len(sequence_pairs.keys()) * 2), level='SUCCESS')
    print(str_success("Done creating {} FastQC reports.".format(len(sequence_pairs.keys()) * 2)))


# ---------------------------------------------
#   MultiQC
# ---------------------------------------------

def run_multiqc():
    completed = subprocess.run(f"multiqc -c multiqc_config.yaml -o output/multiqc -ip --profile-runtime output",
                               shell=True, capture_output=True)
    if completed.returncode != 0:
        log("Error occurred when creating MultiQC report", level='ERROR')
        log("Stacktrace : {}".format(completed.stderr), level='DEBUG')
        raise Exception(completed.stderr)
    log("Done creating MultiQC report", level='SUCCESS')
    print(str_success("Created MultiQC report"))


# ---------------------------------------------
#   Sickle
# ---------------------------------------------

def _sickle_command(pe_file1, pe_file2, basename):
    command = f"{common.get_tool_path('sickle')} pe -f input/{pe_file1} -r input/{pe_file2} -t sanger " + \
              f"-o output/intermediate_files/sickled_{pe_file1} -p output/intermediate_files/sickled_{pe_file2} " + \
              f"-s output/intermediate_files/sickled_singles_{basename}.fastq 2>output/logs/sickle_log_{basename}.log"
    log("subprocess run:" + command, level='DEBUG')
    return command


def _run_sickle(forward_sequence: str, reverse_sequence: str, basename: str):
    completed = subprocess.run(_sickle_command(forward_sequence, reverse_sequence, basename), shell=True,
                               capture_output=True)
    if completed.returncode != 0:
        log("Error occurred when sickling sequences {}+{}".format(forward_sequence, reverse_sequence), level='ERROR')
        raise Exception(completed.stderr)


def batch_sickle(sequence_pairs):
    sickle_threads = []
    for i, (bn, (fwd, rev, ext)) in enumerate(sequence_pairs.items()):
        t = threading.Thread(target=_run_sickle, args=(fwd + ext, rev + ext, bn))
        while threading.active_count() >= 8:
            log(f"{threading.active_count()} active threads. Sleeping for a second...")
            time.sleep(1)
        log(f"[{i + 1}/{len(sequence_pairs.keys())}] Sickling {bn}", level='WARN')
        print(f"[{i + 1}/{len(sequence_pairs.keys())}] Sickling {bn}")
        t.start()

    for t in sickle_threads:
        t.join()
    log(f"Done sickling {len(sequence_pairs.keys())} paired sequences", level='SUCCESS')
    print(str_success(f"Done sickling {len(sequence_pairs.keys())} paired sequences"))

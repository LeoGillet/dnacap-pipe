# -*- coding: utf-8 -*-
"""
This module contains multiple functions related to the Quality Control of the input sequences
    - FastQC functions generate QC reports on all input sequences
    - FastP trims sequences based on base quality
    - MultiQC scans output folder for logs and reports in order to create a run summary
"""
import subprocess
import threading
import time
from tqdm import tqdm

from src import common
from src.logger import log, str_success


# ---------------------------------------------
#   FastQC
# ---------------------------------------------


def _exec_fastqc(fastqfile, mode):
    completed = subprocess.run(
        f"{common.get_tool_path('fastqc')} {fastqfile} -o output/intermediate_files/fastqc_reports/{mode}",
        shell=True,
        capture_output=True,
        check=True,
    )
    if completed.returncode != 0:
        log(
            f"Error occurred when executing FastQC on sequence {fastqfile}",
            level="ERROR",
        )
        log(f"Stacktrace : {completed.stderr}", level="DEBUG")
        raise RuntimeError(completed.stderr)
    log(f"Created FastQC report for sequence {fastqfile}")


def run_all_fastqc(sequence_pairs):
    """
    Runs all FastQC subprocesses in individual threads
    """
    fastqc_threads = []
    max_threads = common.get_max_threads("fastqc")
    with tqdm(total=len(sequence_pairs) * 2) as pbar:
        for mode, input_dir in (
            ("pretrim", "input/"),
            ("trimmed", "output/intermediate_files/trimmed/"),
        ):
            for basename, (fwd, rev, ext) in sequence_pairs.items():
                if fwd:
                    t_fwd = threading.Thread(
                        target=_exec_fastqc,
                        args=(
                            input_dir + fwd + ext,
                            mode,
                        ),
                        daemon=True,
                    )
                if rev:
                    t_rev = threading.Thread(
                        target=_exec_fastqc,
                        args=(
                            input_dir + rev + ext,
                            mode,
                        ),
                        daemon=True,
                    )
                pbar.set_description(
                    f"Generating FastQC report for paired sequences {basename} ({mode})"
                )
                while threading.active_count() >= max_threads:
                    pbar.refresh()
                    time.sleep(1)
                pbar.set_description(
                    f"Generating FastQC report for paired sequences {basename}"
                )
                if fwd:
                    t_fwd.start()
                    fastqc_threads.append(t_fwd)
                if rev:
                    t_rev.start()
                    fastqc_threads.append(t_rev)
                pbar.update()

            for thread in fastqc_threads:
                thread.join()
            pbar.update()
            log(
                f"Done creating {mode} {len(sequence_pairs.keys()) * 2} FastQC reports.",
                level="SUCCESS",
            )
            pbar.write(
                str_success(
                    f"Done creating {mode} {len(sequence_pairs.keys()) * 2} FastQC reports."
                )
            )
        pbar.close()


# ---------------------------------------------
#   Read cleaning & trimming
# ---------------------------------------------


def _trimmomatic_command_paired_end(pe_file1, pe_file2, basename):
    command = (
        f"java -jar {common.get_tool_path('trimmomatic')} PE "
        f"input/{pe_file1} input/{pe_file2} "
        f"output/intermediate_files/trimmed/{pe_file1} "
        f"output/intermediate_files/trimmed/unpaired_{pe_file1} "
        f"output/intermediate_files/trimmed/{pe_file2} "
        f"output/intermediate_files/trimmed/unpaired_{pe_file2} "
        f"ILLUMINACLIP:{common.get_adapters_path('TruSeq3-PE')}:2:30:10:2:True LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    )
    return command


def _trimmomatic_command_single_end(pe_file1, basename):
    return (
        f"java -jar {common.get_tool_path('trimmomatic')} SE "
        f"input/{pe_file1} "
        f"output/intermediate_files/trimmed/{pe_file1} "
        f"ILLUMINACLIP:{common.get_adapters_path('TruSeq3-SE')}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
    )

def _fastp_single_end(in1):
    command = (
        f"{common.get_tool_path('fastp')} -l 36 -3 -5 -P 20 "
        f"--trim_poly_g --trim_poly_x --cut_mean_quality 3 --cut_window_size 4 "
        f"--adapter_fasta {common.get_adapters_path('TruSeq3-SE')} "
        f"--in1 input/{in1} --out1 output/intermediate_files/trimmed/{in1}"
    )
    subprocess.run(
        command,
        shell=True,
        capture_output=True,
        check=True,
    )


def _fastp_paired_end(in1, in2):
    command = (
        f"{common.get_tool_path('fastp')} -l 36 -3 -5 -P 20 "
        f"--trim_poly_g --trim_poly_x --cut_mean_quality 3 --cut_window_size 4 "
        f"--adapter_fasta {common.get_adapters_path('TruSeq3-PE')} "
        f"--in1 input/{in1} --in2 input/{in2} "
        f"--out1 output/intermediate_files/trimmed/{in1} "
        f"--out2 output/intermediate_files/trimmed/{in2}"
    )
    subprocess.run(
        command,
        shell=True,
        capture_output=True,
        check=True,
    )


def batch_cleaning(sequence_pairs):
    """
    Runs fastp on all input sequences in individual threads
    """
    cleaning_threads = []
    max_threads = common.get_max_threads("cleaning")
    with tqdm(total=len(sequence_pairs) + 1) as pbar:
        for i, (basename, (fwd, rev, ext)) in enumerate(sequence_pairs.items()):
            if fwd and rev:
                thread = threading.Thread(
                    target=_fastp_paired_end,
                    args=(fwd + ext, rev + ext,),
                    daemon=True,
                )
            else:
                seq = fwd if fwd else rev
                thread = threading.Thread(
                    target=_fastp_single_end,
                    args=(seq + ext,),
                    daemon=True,
                )
            pbar.set_description(f"Cleaning {basename}")
            while threading.active_count() >= max_threads:
                pbar.refresh()
                time.sleep(1)
            log(
                f"[{i + 1}/{len(sequence_pairs.keys())}] Cleaning {basename}",
                level="WARN",
            )
            cleaning_threads.append(thread)
            thread.start()
            pbar.update()

        for thread in cleaning_threads:
            thread.join()
        pbar.update()
        pbar.write(
            str_success(f"Done sickling {len(sequence_pairs.keys())} paired sequences")
        )
        pbar.close()
    log(f"Done sickling {len(sequence_pairs.keys())} paired sequences", level="SUCCESS")

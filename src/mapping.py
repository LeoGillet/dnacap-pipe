# -*- coding: utf-8 -*-
"""
Contains all functions and tasks for the mapping
of all input sequences on all genomes mentioned in config json
"""
import os
import subprocess
import threading
import time
from tqdm import tqdm

from src import common
from src.logger import log, str_success

# ---------------------------------------------
#   Human decontamination
# ---------------------------------------------


def _start_human_mapping(fwd, rev, ext, basename):
    fwd_output_path = f"output/mapping/human_unmapped/{fwd+ext}"
    rev_output_path = f"output/mapping/human_mapped/{rev+ext}"
    fwd_input_path = f"output/intermediate_files/trimmed/{fwd+ext}"
    rev_input_path = f"output/intermediate_files/trimmed/{rev+ext}"
    subprocess.run(
        f"{common.get_tool_path('bowtie2')} -p 4 "
        + f"-x {common.get_genome_path('H.sapiens.GRCh38', 'bowtie')} "
        + f"-1 {fwd_input_path} -2 {rev_input_path} "
        + f"--un-conc-gz {fwd_output_path} "
        + f"--al-conc-gz {rev_output_path} "
        + f"-S output/mapping/human_sam/{basename}.sam",
        check=True,
        shell=True,
        capture_output=True,
    )
    os.rename(
        f"output/mapping/human_unmapped/{basename}_R1.fastq.1.gz",
        f"output/mapping/human_unmapped/{fwd + ext}",
    )
    os.rename(
        f"output/mapping/human_unmapped/{basename}_R1.fastq.2.gz",
        f"output/mapping/human_unmapped/{rev + ext}",
    )
    subprocess.run(
        f"{common.get_tool_path('samtools')} view -h -b -@ 3 "
        + f"output/mapping/human_sam/{basename}.sam > "
        + f"output/mapping/human_sam/{basename}_unsorted.bam",
        check=True,
        shell=True,
    )
    subprocess.run(
        f"{common.get_tool_path('samtools')} sort --write-index "
        + f"-o output/mapping/human_sam/{basename}.bam "
        + f"output/mapping/human_sam/{basename}_unsorted.bam",
        check=True,
        shell=True,
    )
    subprocess.run(
        f"{common.get_tool_path('samtools')} flagstat output/mapping/human_sam/{basename}.bam "
        + f"> output/mapping/human_flagstat/{basename}.txt",
        check=True,
        shell=True,
    )


def map_human(paired_sequences):
    """
    Starts threads specific to the mapping on human genome for decontamination
    """
    mapping_threads = []
    max_threads = common.get_max_threads("mapping")
    with tqdm(total=len(paired_sequences) + 1) as pbar:
        for i, (basename, (fwd, rev, ext)) in enumerate(paired_sequences.items()):
            thread = threading.Thread(
                target=_start_human_mapping,
                args=(
                    fwd,
                    rev,
                    ext,
                    basename,
                ),
                daemon=True,
            )
            mapping_threads.append(thread)
            pbar.set_description(f"Mapping {basename} on human genome...")
            while threading.active_count() > max_threads:
                pbar.set_description(
                    f"Mapping {basename} on human genome... (waiting for threads)"
                )
                pbar.refresh()
                time.sleep(1)
            pbar.set_description(f"Mapping {basename} on human genome...")
            thread.start()
            log(
                f"[{i+1}/{len(paired_sequences.keys())}] Mapping {basename} on human genome..."
            )
            pbar.update()

        for thread in mapping_threads:
            thread.join()
        pbar.update()
        log(
            f"Done mapping {len(paired_sequences.keys())} paired sequences on human genome",
            level="SUCCESS",
        )
        pbar.write(
            str_success(
                f"Done mapping {len(paired_sequences.keys())} paired sequences on human genome"
            )
        )
        pbar.close()


# ---------------------------------------------
#   Other species assembly
# ---------------------------------------------


def _start_mapping(fr_files: (str, str, str), basename, genome, on_human=True):
    fwd, rev, ext = fr_files
    output_dir = f"output/mapping/{genome}"
    if on_human:
        input_path = "output/mapping/human_unmapped"
    else:
        input_path = "output/intermediate_files/trimmed/"
    try:
        os.makedirs(output_dir)
    except FileExistsError:
        pass
    if fwd and rev:
        command = (
            f"{common.get_tool_path('bwa')} mem -t 2 {common.get_genome_path(genome, 'bwa')} "
            + f"{input_path}/{fwd+ext} "
            + f"{input_path}/{rev+ext} > {output_dir}/{basename}.sam"
        )
    elif fwd:
        command = (
            f"{common.get_tool_path('bwa')} mem -t 2 {common.get_genome_path(genome, 'bwa')} "
            + f"{input_path}/{fwd+ext} > {output_dir}/{basename}.sam"
        )
    elif rev:
        command = (
            f"{common.get_tool_path('bwa')} mem -t 2 {common.get_genome_path(genome, 'bwa')} "
            + f"{input_path}/{rev+ext} > {output_dir}/{basename}.sam"
        )
    subprocess.run(
        command,
        check=True,
        shell=True,
        capture_output=True,
    )
    subprocess.run(
        f"{common.get_tool_path('samtools')} view -h -b -@ 3 "
        + f"{output_dir}/{basename}.sam > "
        + f"{output_dir}/{basename}_unsorted.bam",
        check=True,
        shell=True,
    )
    subprocess.run(
        f"{common.get_tool_path('samtools')} sort --write-index "
        + f"-o {output_dir}/{basename}.bam "
        + f"{output_dir}/{basename}_unsorted.bam",
        check=True,
        shell=True,
    )
    os.remove(f"{output_dir}/{basename}_unsorted.bam")
    subprocess.run(
        f"{common.get_tool_path('samtools')} flagstat {output_dir}/{basename}.bam "
        + f"> {output_dir}/{basename}.txt",
        check=True,
        shell=True,
    )


def _start_mapping_no_human(fr_files: (str, str, str), basename, genome):
    return _start_mapping(fr_files, basename, genome, on_human=False)


def complete_mapping(sequence_pairs, genomes, on_human=False):
    """
    Maps all input sequences (decontaminated or not) to all genomes
    specified in config json
    """
    mapping_threads = []
    max_threads = common.get_max_threads("mapping")
    target_func = _start_mapping if on_human else _start_mapping_no_human
    with tqdm(total=len(sequence_pairs) * len(genomes) + 1) as pbar:
        for basename, (fwd, rev, ext) in sequence_pairs.items():
            for genome in genomes:
                pbar.set_description(f"Mapping {basename} on {genome}...")
                thread = threading.Thread(
                    target=target_func,
                    args=((fwd, rev, ext), basename, genome),
                    daemon=True,
                )
                mapping_threads.append(thread)
                while threading.active_count() > max_threads:
                    pbar.set_description(
                        f"Mapping {basename} on {genome}"
                    )
                    pbar.refresh()
                    time.sleep(1)
                thread.start()
                pbar.set_description(f"Mapping {basename} on {genome}...")
                pbar.update()
            tqdm.write(
                str_success(f"Done mapping {basename} on {len(genomes)} references")
            )

    for thread in mapping_threads:
        thread.join()
    pbar.update()
    pbar.close()

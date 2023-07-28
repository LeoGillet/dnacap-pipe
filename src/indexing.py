# -*- coding: utf-8 -*-
"""
This module contains all functions related to reference genomes indexing
by both BWA & Bowtie2 aligners
"""
import os
import subprocess

from src.common import get_tool_path
from src.logger import log


def _bwa_already_indexed(genome):
    index_files = (
        genome + ".amb",
        genome + ".ann",
        genome + ".bwt",
        genome + ".pac",
        genome + ".sa"
    )
    if all(os.path.isfile(f'refs/bwa/{path}') for path in index_files):
        return True
    return False


def _bowtie_already_indexed(genome):
    # /!\  This is a pretty lazy check and doesn't guarantee the integrity of the indexed genome
    #      + is slow because of repetitive folder scans
    for file in os.listdir('refs/bowtie'):
        if genome + '.indexed' in file:
            return True
    return False


def _bwa_indexing(genome):
    print(f"Indexing genome {genome} with bwa index...")
    completed = subprocess.run(f'{get_tool_path("bwa")} index "refs/bwa/{genome}"', shell=True,
                               stdout=subprocess.DEVNULL,
                               stderr=subprocess.STDOUT,
                               check=True)
    if completed.returncode != 0:
        log(f"Error occurred while indexing genome {genome} with BWA. Stopping.", level='ERROR')
        log(f"Stacktrace: {completed.stderr}", level='ERROR')
    else:
        log(f"bwa indexing: Genome {genome} has been indexed")


def _bowtie_indexing(genome):
    print(f"Indexing genome {genome} with bowtie2-build...")
    completed = subprocess.run(
        f'{get_tool_path("bowtie2-build")} "refs/bowtie/{genome}" "refs/bowtie/{genome}.indexed"',
        shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT, check=True)
    if completed.returncode != 0:
        log(
            f"Error occurred while indexing genome {genome} with bowtie2-build. Stopping.",
            level='ERROR'
        )
        log(f"Stacktrace: {completed.stderr}", level='ERROR')
    else:
        log(f"bowtie2-build indexing: Genome {genome} has been indexed")


def index_everything():
    """
    Finds all compatible sequences (fasta) and indexes them firstly with BWA
    and then with bowtie2
    """
    genomes_fa = [fa for fa in os.listdir('refs/bwa') if fa.endswith('.fasta')]
    for genome in genomes_fa:
        if not _bwa_already_indexed(genome):
            _bwa_indexing(genome)
    genomes_fa = [fa for fa in os.listdir('refs/bowtie') if fa.endswith('.fasta')]
    for genome in genomes_fa:
        if not _bowtie_already_indexed(genome):
            _bowtie_indexing(genome)

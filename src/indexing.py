# -*- coding: utf-8 -*-
import os
import subprocess

from src.common import *


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
    # TODO: src.indexing._bowtie_already_indexed not optimal checking
    # /!\  This is a pretty lazy check and doesn't guarantee the integrity of the indexed genome
    #      + is slow because of repetitive folder scans
    for file in os.listdir('refs/bowtie'):
        if genome + '.indexed' in file:
            return True
    return False


def _bwa_indexing(genome):
    print("Indexing genome {} with bwa index...".format(genome))
    completed = subprocess.run(f'{get_tool_path("bwa")} index "refs/bwa/{genome}"', shell=True,
                               stdout=subprocess.DEVNULL,
                               stderr=subprocess.STDOUT)
    if completed.returncode != 0:
        log("Error occurred while indexing genome {} with BWA. Stopping.".format(genome), level='ERROR')
        log("Stacktrace: {}".format(completed.stderr), level='ERROR')
    else:
        log("bwa indexing: Genome {} has been indexed".format(genome))


def _bowtie_indexing(genome):
    print("Indexing genome {} with bowtie2-build...".format(genome))
    completed = subprocess.run(
        f'{get_tool_path("bowtie2-build")} "refs/bowtie/{genome}" "refs/bowtie/{genome}.indexed"',
        shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    if completed.returncode != 0:
        log("Error occurred while indexing genome {} with bowtie2-build. Stopping.".format(genome), level='ERROR')
        log("Stacktrace: {}".format(completed.stderr), level='ERROR')
    else:
        log("bowtie2-build indexing: Genome {} has been indexed".format(genome))


def index_everything():
    # TODO: Can probably be handled better, but works
    genomes_fa = [fa for fa in os.listdir('refs/bwa') if fa.endswith('.fasta')]
    for genome in genomes_fa:
        if not _bwa_already_indexed(genome):
            _bwa_indexing(genome)
    genomes_fa = [fa for fa in os.listdir('refs/bowtie') if fa.endswith('.fasta')]
    for genome in genomes_fa:
        if not _bowtie_already_indexed(genome):
            _bowtie_indexing(genome)

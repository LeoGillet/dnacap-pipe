#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main entry point to the program
"""
import platform

from src import common
from src import indexing as idx
from src import logger
from src import mapping as mp
from src import prerequisites as pre
from src import qc
from src import reports as rep
from src import sequences as seq

if __name__ == '__main__':
    pre.rename_old_log()
    pre.delete_older_logs()
    logger.log(
        "Started workflow on machine " +
        f"({platform.platform()} {platform.machine()}-{platform.version()})"
    )
    pre.ensure_paths()
    pre.prepare_dirs()
    pre.check_folders()
    # SEQUENCER_TYPE = pre.ask_sequencer_type()
    sequence_pairs = seq.pair_sequences()
    idx.index_everything()
    qc.run_all_fastqc(sequence_pairs)
    qc.batch_sickle(sequence_pairs)
    # mp.map_human(sequence_pairs)

    genomes = common.choose_genomes()
    mp.complete_mapping(sequence_pairs, genomes)
    rep.generate_fastqc()
    rep.generate_mapping_reports(genomes)

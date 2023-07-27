#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import platform
import sys

from src import prerequisites as pre
from src import logger
from src import sequences as seq
from src import indexing as idx
from src import qc

if __name__ == '__main__':
    pre.rename_old_log()
    pre.delete_older_logs()
    logger.log(f"Started workflow on machine ({platform.platform()} {platform.machine()}-{platform.version()})")
    pre.ensure_paths()
    pre.prepare_dirs()
    pre.check_folders()
    sequencer_type = pre.ask_sequencer_type()
    sequence_pairs = seq.pair_sequences()
    idx.index_everything()
    qc.run_all_fastqc(sequence_pairs)
    qc.batch_sickle(sequence_pairs)
    qc.run_multiqc()

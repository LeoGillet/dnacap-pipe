# -*- coding: utf-8 -*-
"""
This module contains all functions related to the scanning of the input directory.
It tries to pair all sequences by following the standard Illumina sanger file conventions
"""
import os

from src.logger import log


def _list_supported_sequences() -> list[str]:
    """
    Lists supported sequences in input directory
    :return: list containing file names
    """
    all_files = os.listdir("input")
    supported_files = [
        file
        for file in all_files
        if file.endswith(".fastq.gz") or file.endswith(".fq.gz")
    ]
    log(
        f"Input folder contains {len(all_files)} files ({len(supported_files)} supported, "
        + f"{len(all_files) - len(supported_files)} ignored)"
    )
    return supported_files


def _display_paired_sequences(sequences: dict[tuple[str, str, str]]):
    """
    Displays in stdout all paired sequences
    :param sequences: supported files
    """
    for i, (_, (fwd, rev, ext)) in enumerate(sequences.items()):
        print(f"{i+1}:\t{fwd}\t<->\t{rev}\t({ext})")


def pair_sequences(quiet=False) -> dict[tuple[str, str, str]]:
    """
    Pairs all supported files in input directory based on file names
    :return: dict containing tuples (forward, reverse, file extension) (keys are basenames)
    """
    all_sequences = {}
    _done_files = []
    for file in _list_supported_sequences():
        _tmp_f = file
        if file in _done_files:
            continue
        if ".fastq.gz" in _tmp_f:
            ext = ".fastq.gz"
            _tmp_f = _tmp_f.replace(".fastq.gz", "")
        if "_R1" in _tmp_f:
            _tmp_f = _tmp_f.replace("_R1", "")
            if os.path.isfile("input/" + _tmp_f + "_R2" + ext):
                all_sequences[_tmp_f] = (_tmp_f + "_R1", _tmp_f + "_R2", ext)
                _done_files.extend([file, _tmp_f + "_R2" + ext])
            else:
                all_sequences[_tmp_f] = (_tmp_f + "_R1", None, ext)
                _done_files.append(file)
        if "_R2" in _tmp_f:
            _tmp_f = _tmp_f.replace("_R2", "")
            if os.path.isfile("input/" + _tmp_f + "_R1" + ext):
                all_sequences[_tmp_f] = (_tmp_f + "_R1", _tmp_f + "_R2", ext)
                _done_files.extend([_tmp_f + "_R1" + ext, file])
            else:
                all_sequences[_tmp_f] = (_tmp_f + "_R2", None, ext)
    log(f"{len(all_sequences.keys())} paired sequences found")
    if not quiet:
        _display_paired_sequences(all_sequences)
    return all_sequences

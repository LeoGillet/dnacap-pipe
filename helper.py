#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys

from src import common
from src import indexing as idx
from src import logger
from src import mapping as mp
from src import prerequisites as pre
from src import qc
from src import reports as rep
from src import sequences as seq

if __name__ == '__main__':
    print("╔══════════════════════════════════╗")
    print("║   % Pipeline Actions Helper %    ║")
    print("║ 1. Summarize flagstats           ║")
    print("╚══════════════════════════════════╝")
    choice = input("Your choice ? > ")

    match choice:
        case '1':
            rep.group_flagstat_results(common.choose_genomes())
        case _:
            print("Invalid choice. Quitting...")
            sys.exit(1)
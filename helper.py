#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Entry point file for running supplementary tasks
"""
import sys

from src import common
from src import reports as rep

if __name__ == "__main__":
    print("╔════════════════════════════════════╗")
    print("║    % Pipeline Actions Helper %     ║")
    print("║                                    ║")
    print("║ 1. Regenerate end MultiQC reports  ║")
    print("║ 2. Summarize flagstats             ║")
    print("║ Q. Quit                            ║")
    print("╚════════════════════════════════════╝")
    choice = input("Your choice ? > ")

    match choice:
        case "1":
            rep.generate_fastqc()
            rep.generate_mapping_reports(common.choose_genomes())
        case "2":
            on_human = input("Did you map on GRCh38 ? [YES|no] > ")
            while on_human.lower() not in ("yes", "no", "y", "n", ""):
                on_human = input("Did you map on GRCh38 ? [YES|no] > ")
            if on_human.lower() in ("yes", "y", ""):
                rep.group_flagstat_results(common.choose_genomes(), on_human=True)
            else:
                rep.group_flagstat_results(common.choose_genomes(), on_human=False)
        case "Q":
            sys.exit(0)
        case _:
            print("Invalid choice. Quitting...")
            sys.exit(1)

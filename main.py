#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Entry point file for running supplementary tasks
"""
import platform
import sys

from src import bed
from src import common
from src import indexing as idx
from src import logger
from src import mapping as mp
from src import prerequisites as pre
from src import reports as rep
from src import qc
from src import sequences as seq
from src import status as st

if __name__ == "__main__":
    wf_status = st.Status()
    genomes = []
    sequence_pairs = []
    common.clear_stdout()
    while True:
        print()
        print("               % Pipeline Actions Helper %")
        print()
        print(
            f"\t0. [Required] Run prerequisites {'✓' if wf_status.get('prerequisites') else ''}"
        )
        print(f"\t1. Generate FastQC reports {'✓' if wf_status.get('fastqc_reports') else ''}")
        print(f"\t2. [Required] Clean reads {'✓' if wf_status.get('cleaning') else ''}")
        print(f"\t3. [Required] Start mapping {'✓' if wf_status.get('mapping') else ''}")
        print(f"\t4. Regenerate end MultiQC reports {'✓' if wf_status.get('mapping') else ''}")
        print(f"\t5. Summarize flagstats {'✓' if wf_status.get('flagstats') else ''}")
        print(f"\t6. Generate coverage files {'✓' if wf_status.get('coverage') else ''}")
        print(f"\t7. Generate depth files {'✓' if wf_status.get('depth_txt') else ''}")
        print(f"\t8. Create depth plots {'✓' if wf_status.get('depth_plots') else ''}")
        print(f"\t9. Extract genome {'✓' if wf_status.get('genome_extracted') else ''}")
        print(
            f"\t10. Extract regions of interests {'✓' if wf_status.get('genes_extracted') else ''}"
        )
        print("\tQ. Quit\n")
        if wf_status.get('prerequisites'):
            print("\tCommon run information :")
            print(f"\t\t- {wf_status.get('stat_seq_n'): <3} sequences")
            print(f"\t\t- {wf_status.get('stat_ref_gen_n'): <3} reference genomes")
            print("\n\n")
        choice = input("Your choice ? > ").upper()

        match choice:
            case "0":
                pre.rename_old_log()
                pre.delete_older_logs()
                logger.log(
                    "Started workflow on machine "
                    + f"({platform.platform()} {platform.machine()}-{platform.version()})"
                )
                pre.ensure_paths()
                pre.prepare_dirs()
                pre.check_folders()
                status = st.Status()
                common.clear_stdout()
                sequence_pairs = seq.pair_sequences()
                idx.index_everything()
                genomes = common.choose_genomes()
                wf_status.set("prerequisites", 1)
                wf_status.set("stat_seq_n", len(sequence_pairs))
                wf_status.set("stat_ref_gen_n", len(genomes))
            case "1":
                sequence_pairs = seq.pair_sequences(quiet=True)
                qc.run_all_fastqc(sequence_pairs)
                wf_status.set("fastqc_reports", 1)
            case "2":
                if wf_status.prerequisites:
                    sequence_pairs = seq.pair_sequences(quiet=True)
                    qc.batch_sickle(sequence_pairs)
                    wf_status.set("cleaning", 1)
                else:
                    common.clear_stdout()
                    print(
                        "You must run prerequisites before starting the mapping workflow!"
                    )

            case "3":
                if wf_status.get("prerequisites") and wf_status.get("cleaning"):
                    sequence_pairs = seq.pair_sequences(quiet=True)
                    genomes = common.choose_genomes()
                    mp.complete_mapping(sequence_pairs, genomes)
                    rep.generate_fastqc()
                    rep.generate_mapping_reports(genomes)
                    wf_status.set("mapping", 1)
                else:
                    common.clear_stdout()
                    print(
                        "You must run prerequisites and sequence cleaning",
                        "before starting the mapping workflow!",
                    )
            case "4":
                rep.generate_fastqc()
                rep.generate_mapping_reports(common.choose_genomes())
            case "5":
                on_human = input("Did you map on GRCh38 ? [YES|no] > ")
                while on_human.lower() not in ("yes", "no", "y", "n", ""):
                    on_human = input("Did you map on GRCh38 ? [YES|no] > ")
                if on_human.lower() in ("yes", "y", ""):
                    rep.group_flagstat_results(common.choose_genomes(), on_human=True)
                else:
                    rep.group_flagstat_results(common.choose_genomes(), on_human=False)
                wf_status.set("flagstats", 1)
            case "6":
                bed.start_coverage_analysis(
                    common.get_bam_files(), common.get_genome_regions()
                )
                wf_status.set("coverage", 1)
            case "7":
                bed.start_depth_analysis(common.get_bam_files())
                wf_status.set("depth_txt", 1)
            case "8":
                bed.generate_depth_plots(
                    common.get_depth_files(common.get_ignored_genomes_step3())
                )
                wf_status.set("depth_plots", 1)
            case "9":
                bed.start_genome_extraction(
                    common.get_bam_files_extraction(common.get_ignored_genomes_step3())
                )
                wf_status.set("genome_extracted", 1)
            case "10":
                bed.start_genes_extractions(
                    common.get_genes_files(common.get_ignored_genomes_step3())
                )
                wf_status.set("genes_extracted", 1)
            case "Q":
                sys.exit(0)
            case _:
                pass

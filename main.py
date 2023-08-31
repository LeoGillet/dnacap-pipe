#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Entry point file for running supplementary tasks
"""
import platform
import sys

from src import bed
from src import common
from src import dataframe as dfio
from src import indexing as idx
from src import logger
from src import mapping as mp
from src import menu_genitalium as mg
from src import menu_pylori as mp
from src import prerequisites as pre
from src import reports as rep
from src import qc
from src import sequences as seq
from src import status as st

if __name__ == "__main__":
    wf_status = st.Status()
    pipeline_type = wf_status.get("pipeline_type")
    genomes = []
    sequence_pairs = []
    common.clear_stdout()
    while True:
        print()
        if pipeline_type == "pylori":
            print("               % H. pylori pipeline %")
        if pipeline_type == "genitalium":
            print("               % M. genitalium pipeline %")
        print()
        print(
            f"\t0. [Required] Run prerequisites {'✓' if wf_status.get('prerequisites') else ''}"
        )
        print(f"\t1. [Required] Clean reads {'✓' if wf_status.get('cleaning') else ''}")
        print(
            f"\t2. [Required] Start mapping {'✓' if wf_status.get('mapping') else ''}"
        )
        print(
            f"\t3. Generate FastQC reports {'✓' if wf_status.get('fastqc_reports') else ''}"
        )
        print(
            f"\t4. Regenerate end MultiQC reports {'✓' if wf_status.get('mapping') else ''}"
        )
        print(f"\t5. Summarize flagstats {'✓' if wf_status.get('flagstats') else ''}")
        print(f"\t6. Generate depth files {'✓' if wf_status.get('depth_txt') else ''}")
        print(f"\t7. Create depth plots {'✓' if wf_status.get('depth_plots') else ''}")
        print(f"\t8. Extract genome {'✓' if wf_status.get('genome_extracted') else ''}")
        print(
            f"\t9. Extract regions of interests {'✓' if wf_status.get('genes_extracted') else ''}"
        )
        print(
            f"\t10. Generate coverage files {'✓' if wf_status.get('coverage') else ''}"
        )
        print(
            "\n",
            mg.extra_options()
            if wf_status.get("pipeline_type") == "genitalium"
            else mp.extra_options(),
        )
        print("\n\tA. Run everything")
        print(f"\tC. Export generated data to CSV (size: {dfio.get_size_of_pickled_dataframe()})")
        print(f"\tG. Generate HTML report")
        print("\n\tQ. Quit\n")
        if wf_status.get("prerequisites"):
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
                genomes = common.choose_genomes(wf_status.get("pipeline_type"))
                wf_status.set("prerequisites", 1)
                wf_status.set("stat_seq_n", len(sequence_pairs))
                wf_status.set("stat_ref_gen_n", len(genomes))
            case "1":
                if wf_status.get("prerequisites"):
                    sequence_pairs = seq.pair_sequences(quiet=True)
                    qc.batch_cleaning(sequence_pairs)
                    wf_status.set("cleaning", 1)
                else:
                    common.clear_stdout()
                    print(
                        "You must run prerequisites before starting the mapping workflow!"
                    )
            case "2":
                if wf_status.get("prerequisites") and wf_status.get("cleaning"):
                    sequence_pairs = seq.pair_sequences(quiet=True)
                    genomes = common.choose_genomes(wf_status.get("pipeline_type"))
                    mp.complete_mapping(sequence_pairs, genomes)
                    wf_status.set("mapping", 1)
                else:
                    common.clear_stdout()
                    print(
                        "You must run prerequisites and sequence cleaning",
                        "before starting the mapping workflow!",
                    )
            case "3":
                sequence_pairs = seq.pair_sequences(quiet=True)
                qc.run_all_fastqc(sequence_pairs)
                wf_status.set("fastqc_reports", 1)
            case "4":
                rep.generate_fastqc()
                rep.generate_mapping_reports(
                    common.choose_genomes(wf_status.get("pipeline_type"))
                )
            case "5":
                on_human = input("Did you map on GRCh38 ? [YES|no] > ")
                while on_human.lower() not in ("yes", "no", "y", "n", ""):
                    on_human = input("Did you map on GRCh38 ? [YES|no] > ")
                if on_human.lower() in ("yes", "y", ""):
                    rep.group_flagstat_results(
                        common.choose_genomes(wf_status.get("pipeline_type")),
                        on_human=True,
                    )
                else:
                    rep.group_flagstat_results(
                        common.choose_genomes(wf_status.get("pipeline_type")),
                        on_human=False,
                    )
                wf_status.set("flagstats", 1)
            case "6":
                bed.start_depth_analysis(common.get_bam_files())
                wf_status.set("depth_txt", 1)
                bed.depth_to_dataframe(
                    common.get_depth_files(wf_status.get("pipeline_type"))
                )
            case "7":
                bed.generate_depth_plots(
                    common.get_depth_files(wf_status.get("pipeline_type"))
                )
                wf_status.set("depth_plots", 1)
            case "8":
                bed.start_genome_extraction(
                    common.get_bam_files_extraction(
                        common.get_ignored_genomes(wf_status.get("pipeline_type"))
                    )
                )
                wf_status.set("genome_extracted", 1)
            case "9":
                bed.start_genes_extractions(
                    common.get_genes_files(wf_status.get("pipeline_type"))
                )
                bed.start_genes_grouping(wf_status.get("pipeline_type"))
                wf_status.set("genes_extracted", 1)
            case "10":
                bed.start_coverage_analysis(
                    common.get_bam_files(),
                    common.get_genome_regions(wf_status.get("pipeline_type")),
                )
                bed.start_per_gene_coverage(wf_status.get("pipeline_type"))
                wf_status.set("coverage", 1)
            case "11" | "12":
                if wf_status.get("pipeline_type") == "pylori":
                    mp.start(int(choice))
                else:
                    mg.start(int(choice))
            case "Q":
                sys.exit(0)
            case "A":
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
                genomes = common.choose_genomes(wf_status.get("pipeline_type"))
                wf_status.set("prerequisites", 1)
                wf_status.set("stat_seq_n", len(sequence_pairs))
                wf_status.set("stat_ref_gen_n", len(genomes))
                if wf_status.get("prerequisites"):
                    sequence_pairs = seq.pair_sequences(quiet=True)
                    qc.batch_cleaning(sequence_pairs)
                    wf_status.set("cleaning", 1)
                else:
                    common.clear_stdout()
                    print(
                        "You must run prerequisites before starting the mapping workflow!"
                    )
                if wf_status.get("prerequisites") and wf_status.get("cleaning"):
                    sequence_pairs = seq.pair_sequences(quiet=True)
                    genomes = common.choose_genomes(wf_status.get("pipeline_type"))
                    mp.complete_mapping(sequence_pairs, genomes)
                    wf_status.set("mapping", 1)
                else:
                    common.clear_stdout()
                    print(
                        "You must run prerequisites and sequence cleaning",
                        "before starting the mapping workflow!",
                    )
                sequence_pairs = seq.pair_sequences(quiet=True)
                qc.run_all_fastqc(sequence_pairs)
                wf_status.set("fastqc_reports", 1)
                rep.generate_fastqc()
                rep.generate_mapping_reports(
                    common.choose_genomes(wf_status.get("pipeline_type"))
                )
                on_human = input("Did you map on GRCh38 ? [YES|no] > ")
                while on_human.lower() not in ("yes", "no", "y", "n", ""):
                    on_human = input("Did you map on GRCh38 ? [YES|no] > ")
                if on_human.lower() in ("yes", "y", ""):
                    rep.group_flagstat_results(
                        common.choose_genomes(wf_status.get("pipeline_type")),
                        on_human=True,
                    )
                else:
                    rep.group_flagstat_results(
                        common.choose_genomes(wf_status.get("pipeline_type")),
                        on_human=False,
                    )
                wf_status.set("flagstats", 1)
                bed.start_coverage_analysis(
                    common.get_bam_files(),
                    common.get_genome_regions(wf_status.get("pipeline_type")),
                )
                wf_status.set("coverage", 1)
                bed.start_depth_analysis(common.get_bam_files())
                wf_status.set("depth_txt", 1)
                bed.generate_depth_plots(
                    common.get_depth_files(wf_status.get("pipeline_type"))
                )
                wf_status.set("depth_plots", 1)
                bed.start_genome_extraction(
                    common.get_bam_files_extraction(
                        common.get_ignored_genomes(wf_status.get("pipeline_type"))
                    )
                )
                wf_status.set("genome_extracted", 1)
                bed.start_genes_extractions(
                    common.get_genes_files(wf_status.get("pipeline_type"))
                )
                bed.start_genes_grouping(wf_status.get("pipeline_type"))
                wf_status.set("genes_extracted", 1)
            case "C":
                dfio.export_to_csv()
            case "G":
                rep._to_html('DC2_10_S19', 'pylori')
            case _:
                common.clear_stdout()

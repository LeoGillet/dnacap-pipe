# -*- coding: utf-8 -*-
"""
Compiles FastQC reports and other pipeline logs 
such as flagstat for post-run analysis
"""
import os
import subprocess

from src.logger import log, str_success

# ---------------------------------------------
#   MultiQC Reports
# ---------------------------------------------


def generate_fastqc():
    """
    Runs MultiQC subprocess in order to generate summary from FastQC reports
    """
    completed = subprocess.run(
        "multiqc -o output/multiqc -ip --profile-runtime output/intermediate_files/fastqc_reports",
        shell=True,
        capture_output=True,
        check=True,
    )
    if completed.returncode != 0:
        log("Error occurred when creating MultiQC report", level="ERROR")
        log(f"Stacktrace : {completed.stderr}", level="DEBUG")
        raise RuntimeError(completed.stderr)
    log("Done creating MultiQC report", level="SUCCESS")
    print(str_success("Created MultiQC report"))


def generate_mapping_reports(genomes):
    """
    Runs MultiQC on all flagstats per genome mapping
    """
    for genome in genomes:
        scan_dir = f"output/mapping/{genome}"
        os.makedirs(f"output/multiqc/{genome}")
        subprocess.run(
            f"multiqc -o output/multiqc/{genome} -ip --profile-runtime {scan_dir}",
            shell=True,
            capture_output=True,
            check=True,
        )


def _read_flagstat(file) -> {}:
    with open(file, "r", encoding="UTF-8") as flagstat:
        lines = flagstat.readlines()
    qc_reads = lines[0].strip().split(" ")
    qc_passed, qc_failed = qc_reads[0], qc_reads[2]
    mapped_reads = lines[6].strip().split(" ")
    mapped_qc_passed = mapped_reads[0]
    mapped_qc_passed_perc = mapped_reads[4].replace("(", "")
    paired_reads = lines[11].strip().split(" ")
    paired_qc_passed = paired_reads[0]
    paired_qc_passed_perc = paired_reads[5].replace("(", "")
    singletons = lines[13].strip().split(" ")
    singletons_qc_passed = singletons[0]
    singletons_qc_passed_perc = singletons[4].replace("(", "")
    return {
        "QC Passed": qc_passed,
        "QC Failed": qc_failed,
        "Mapped QC+": mapped_qc_passed,
        "Mapped QC+ (%)": mapped_qc_passed_perc,
        "Properly paired QC+": paired_qc_passed,
        "Properly paired QC+ (%)": paired_qc_passed_perc,
        "Singletons QC+": singletons_qc_passed,
        "Singletons QC+ (%)": singletons_qc_passed_perc,
    }


def _export_flagstat_results(flagstat_results, genomes):
    with open("summarized_flagstat.csv", "w", encoding="UTF-8") as summary:
        # Table title row
        summary.write(f",{',,,,,,,,'.join(genomes)}\n")
        row_titles = ["Sample"]
        row_titles.extend(
            [
                "QC Passed",
                "QC Failed",
                "Mapped QC+",
                "Mapped QC+ (%)",
                "Properly paired QC+",
                "Properly paired QC+ (%)",
                "Singletons QC+",
                "Singletons QC+ (%)",
            ]
            * len(genomes)
        )
        summary.write(",".join(row_titles) + "\n")
        for sample, _ in flagstat_results.items():
            row = [sample]
            for genome in genomes:
                row.extend(list(flagstat_results[sample][genome].values()))
            summary.write(",".join(row) + "\n")


def group_flagstat_results(genomes, on_human=True):
    """
    Summarizes all flagstat files generated from assemblies for comparison analysis
    """
    flagstat_results = {}
    scan_dirs = [(f"output/mapping/{genome}/", genome) for genome in genomes]
    if on_human:
        scan_dirs.append(("output/mapping/human_flagstat/", "Human.GRCh38"))
        genomes.add("Human.GRCh38")
    files_to_scan = []
    for dir_, genome in scan_dirs:
        files_to_scan.extend(
            [
                (dir_ + str(file), str(file).replace(".txt", ""), genome)
                for file in os.listdir(dir_)
                if file.endswith(".txt")
            ]
        )
    for file, basename, genome in files_to_scan:
        if basename not in flagstat_results:
            flagstat_results[basename] = {}
        flagstat_results[basename][genome] = _read_flagstat(file)

    _export_flagstat_results(flagstat_results, genomes)

# -*- coding: utf-8 -*-
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

# -*- coding: utf-8 -*-
"""
Compiles FastQC reports and other pipeline logs 
such as flagstat for post-run analysis
"""
import os
import shutil
import subprocess

import pandas as pd

from src.logger import log, str_success
from src import bed
from src import common
from src import dataframe as dfio

# ---------------------------------------------
#   MultiQC Reports
# ---------------------------------------------


def generate_fastqc():
    """
    Runs MultiQC subprocess in order to generate summary from FastQC reports
    """
    subprocess.run(
        "multiqc -o output/multiqc/pretrim -ip --profile-runtime output/intermediate_files/fastqc_reports/pretrim",
        shell=True,
        capture_output=True,
        check=True,
    )
    subprocess.run(
        "multiqc -o output/multiqc/trimmed -ip --profile-runtime output/intermediate_files/fastqc_reports/trimmed",
        shell=True,
        capture_output=True,
        check=True,
    )

    log("Done creating MultiQC report", level="SUCCESS")
    print(str_success("Created MultiQC report"))


def generate_mapping_reports(genomes):
    """
    Runs MultiQC on all flagstats per genome mapping
    """
    for genome in genomes:
        scan_dir = f"output/mapping/{genome}"
        try:
            os.makedirs(f"output/multiqc/{genome}")
        except FileExistsError:
            shutil.rmtree(f"output/multiqc/{genome}")
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
    qc_passed = qc_reads[0]
    mapped_reads = lines[6].strip().split(" ")
    mapped_qc_passed = mapped_reads[0]
    mapped_qc_passed_perc = mapped_reads[4].replace("(", "")
    return {
        "QC Passed": qc_passed,
        "Mapped QC+": mapped_qc_passed,
        "Mapped QC+ (%)": mapped_qc_passed_perc,
    }


def _export_flagstat_results(flagstat_results, genomes):
    with open("summarized_flagstat.csv", "w", encoding="UTF-8") as summary:
        # Table title row
        summary.write(f",{',,,'.join(genomes)}\n")
        row_titles = ["Sample"]
        row_titles.extend(
            [
                "QC Passed",
                "Mapped QC+",
                "Mapped QC+ (%)",
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

def _get_classes_from_value(value, hival, lowval):
    classes = ""
    if value >= hival:
        classes += "success-text "
    elif value >= lowval:
        classes += "warning-text "
    else:
        classes += "error-text "
    return classes

def _html_depth_coverage_table(dataframe: pd.DataFrame, sample_name, genes: list):
    rows = []
    for gene in genes:
        cover_value = dataframe.loc[sample_name][gene+' coverage']*100
        depth_value = dataframe.loc[sample_name][gene+' mean depth']
        cover_classes = _get_classes_from_value(cover_value, 95.0, 75.0)
        depth_classes = _get_classes_from_value(depth_value, 30.0, 10.0)
        
        row = "<tr>"
        row += f"<td>{gene}</td>"
        row += f"<td class={cover_classes}>{cover_value:.1f}%</td>"
        row += f"<td class={depth_classes}>{depth_value:.2f}</td>"
        row += "</tr>"
        rows.append(row)

    return f"""
    <table>
        <tr>
            <th>Gène</th>
            <th>Couverture</th>
            <th>Profondeur moyenne</th>
        </tr>
        {"".join(rows)}
    </table>
    """

def _html_pylori_amr(dataframe: pd.DataFrame, sample_name):
    # 16S - Tétracycline
    mutations_16s = ""
    if dataframe.loc[sample_name]["16S_1 Mutations"]:
        mutations_16s += ""
    if dataframe.loc[sample_name]["16S_2 Mutations"]:
        ...

    html_16s = f"""
    <span class="bold underline">ADNr 16S - Résistance à la tétracycline</span>
    <br>
    Mutations détectées: {}
    """

_HTML_STYLES = """
<style>
.success-text { color: #00473e; }
.warning-text { color: #faae2b; }
.error-text { color: #fa5246; }
</style>
"""

def _to_html(sample_name, pipeline_type):
    if pipeline_type == "pylori":
        genome = "H.pylori.J99.2001"
    else:
        raise NotImplementedError()

    dataframe = dfio.get_or_create(pipeline_type)
    gene_names = bed._load_names(common.get_genome_regions(pipeline_type) + ".names.list")
    depth_img_origin = f"output/depth/{genome}/{sample_name}.png"
    depth_img_dest = f"output/html_reports/{sample_name}/depth.png"

    table_html_str = _html_depth_coverage_table(dataframe, sample_name, gene_names)
    html_str = f"""
    <!DOCTYPE html>
    <html lang="fr">
        <head>
            <title>Rapport {sample_name}</title>
            {_HTML_STYLES}
        </head>
        <body>
            <h1>Rapport d'analyses DNA Capture</h1>
            <h3>Type de pipeline: {pipeline_type}</h3>
            <h3>Nom de souche: {sample_name}</h3>
            <hr style="padding-bottom=10px;"/>
            <a href="{sample_name}/depth.png" target="_blank">
                <img src="{sample_name}/depth.png" width="480px"/>
            </a>
            <hr style="padding-bottom=10px;"/>
            <h5>Analyse des couvertures et profondeurs moyennes par gènes</h5>
            {table_html_str}
            <hr style="padding-bottom=10px;"/>
            <h5>Analyse des résistances aux antibiotiques</h5>
            <hr style="padding-bottom=10px;"/>
            <h5>Typage MLST</h5>
            <hr style="padding-bottom=10px;"/>
            <h5>Analyse de la virulence</h5>
            <hr style="padding-bottom=10px;"/>
            <h5>Analyse phylogénétique</h5>
        </body>
    </html>\n
    """

    if not os.path.isdir(f"output/html_reports/{sample_name}"):
        os.makedirs(f"output/html_reports/{sample_name}")
    with open(f"output/html_reports/{sample_name}.html", "w", encoding="UTF-8") as html_file:
        html_file.write(html_str)

    shutil.copyfile(depth_img_origin, depth_img_dest)

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
        cover_value = dataframe.loc[sample_name][gene + " coverage"] * 100
        depth_value = dataframe.loc[sample_name][gene + " mean depth"]
        cover_classes = _get_classes_from_value(cover_value, 95.0, 75.0)
        depth_classes = _get_classes_from_value(depth_value, 30.0, 10.0)

        row = "<tr>"
        row += f"<td style='text-align: right;'>{gene}</td>"
        row += f"<td style='text-align: center;' class={cover_classes}>{cover_value:.1f}%</td>"
        row += f"<td style='text-align: center;' class={depth_classes}>{depth_value:.2f}</td>"
        row += "</tr>"
        rows.append(row)

    return f"""
    <table class="cover-depth-table">
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
    muts_16s1 = dataframe.loc[sample_name]["16S_1 Mutations"]
    muts_16s2 = dataframe.loc[sample_name]["16S_2 Mutations"]
    mutations_16s = ""
    if muts_16s1:
        mutations_16s += f"⚠️ Mutation(s) détectée(s) dans le 16S (1ère c.): <span class='mono bold'>{muts_16s1}</span><br>"
    else:
        mutations_16s += "✅ Aucune mutation détectée dans le 16S (1ère c.)<br>"
    if muts_16s2:
        mutations_16s += f"⚠️ Mutation(s) détectée(s) dans le 16S (2ème c.): <span class='mono bold'>{muts_16s2}</span>"
    else:
        mutations_16s += "✅ Aucune mutation détectée dans le 16S (2ème c.)<br>"

    mutations_16s += "<span class='bold underline'>Phénotype de résistance :</span> "
    if muts_16s1 or muts_16s2:
        mutations_16s += "⚠️ <span class='bold error-text'> probablement résistante</span> à la tétracycline<br>"
    else:
        mutations_16s += "✅ <span class='bold success-text'> probablement sensible</span> à la tétracycline<br>"

    # 23S - Chlarithromycine
    muts_23s1 = dataframe.loc[sample_name]["23S_1 Mutations"]
    muts_23s2 = dataframe.loc[sample_name]["23S_2 Mutations"]
    mutations_23s = ""
    if muts_23s1:
        mutations_23s += f"⚠️ Mutation(s) détectée(s) dans le 23S (1ère c.): <span class='mono bold'>{muts_23s1}</span><br>"
    else:
        mutations_23s += "✅ Aucune mutation détectée dans le 23S (1ère c.)<br>"
    if muts_23s2:
        mutations_23s += f"⚠️ Mutation(s) détectée(s) dans le 23S (2ème c.): <span class='mono bold'>{muts_23s2}</span>"
    else:
        mutations_23s += "✅ Aucune mutation détectée dans le 23S (2ème c.)<br>"

    mutations_23s += "<span class='bold underline'>Phénotype de résistance :</span> "
    if muts_23s1 or muts_23s2:
        mutations_23s += "⚠️ <span class='bold error-text'> probablement résistante</span> à la chlarithromycine<br>"
    else:
        mutations_23s += "✅ <span class='bold success-text'> probablement sensible</span> à la chlarithromycine<br>"

    # gyrA - Lévofloxacine
    muts_gyra = dataframe.loc[sample_name]["gyrA Mutations"]
    mutations_gyra = ""
    if muts_gyra:
        mutations_gyra += f"⚠️ Mutation(s) détectée(s) dans gyrA : <span class='mono bold'>{muts_gyra}</span><br>"
        mutations_gyra += (
            "<span class='bold underline'>Phénotype de résistance :</span> "
        )
        mutations_gyra += "⚠️ <span class='bold error-text'> probablement résistante</span> à la lévofloxacine<br>"
    else:
        mutations_gyra += "Aucune mutation détectée dans gyrA<br>"
        mutations_gyra += (
            "<span class='bold underline'>Phénotype de résistance :</span> "
        )
        mutations_gyra += "✅ <span class='bold success-text'> probablement sensible</span> à la lévofloxacine<br>"

    return f"""
    <span class="big bold underline">ADNr 16S - Résistance à la tétracycline</span><br>
    <span class="subscript">1ère copie : 1188028-1189528 | 2ème copie : 1463046-1464546</span>
    <br>
    {mutations_16s}
    <br>
    <span class="big bold underline">ADNr 23S - Résistance à la chlarithromycine</span><br>
    <span class="subscript">1ère copie : <span class="mono">1057558-1060528</span> | 2ème copie : <span class="mono">1427396-1430366</span></span>
    <br>
    {mutations_23s}
    <br>
    <span class="big bold underline">gyrA - Résistance à la lévofloxacine</span><br>
    {mutations_gyra}
    <br>
    <span class="big bold underline">rpoB - Résistance à la rifampicine</span><br>
    <span class="italic">Non implémenté</span>
    <br>
    """


def _parse_score(score_str):
    right = score_str.split("(")[1]
    return float(right.replace(")", ""))


def _html_pylori_mlst(dataframe, sample_name):
    best_match = dataframe.at[sample_name, "MLST best match"]
    match_2 = dataframe.at[sample_name, "MLST match 2"]
    print(_parse_score(best_match), type(_parse_score(best_match)))

    if _parse_score(best_match) <= 100.0 or best_match == match_2:
        return '<span class="italic">Typage MLST non conclusif. Gènes MLST non séquencés ou scores d\'alignement inférieurs à 100.</span>'

    return f"""
    <span class="big bold">Meilleurs matchs</span>
        <ul>
            <li class="bold">{dataframe.at[sample_name, "MLST best match"]}</li>
            <li>{dataframe.at[sample_name, "MLST match 2"]}</li>
            <li>{dataframe.at[sample_name, "MLST match 3"]}</li>
            <li>{dataframe.at[sample_name, "MLST match 4"]}</li>
            <li>{dataframe.at[sample_name, "MLST match 5"]}</li>
        </ul>
    <span class="subscript">Falush, D. (2003). Traces of Human Migrations in Helicobacter pylori Populations. Science, 299(5612), 1582–1585. doi:10.1126/science.1080857</span>
    """


def _html_pylori_cagA(dataframe, sample_name):
    mean_depth = dataframe.at[sample_name, "cagA mean depth"]
    coverage = dataframe.at[sample_name, "cagA coverage"] * 100
    if coverage <= 30.0:
        return """
            Couverture sur le gène <span class="bold">cagA</span> très faible ou nulle.<br>
            ✅ Souche probablement <span class="bold success-text">cagA négative</span>.
        """
    return f"""
            Gène <span class="bold">cagA</span> séquencé. (couverture&#x2248;{coverage:.0f}%, profondeur&#x2248;{mean_depth:.0f})<br>
            Souche probablement <span class="bold error-text">cagA positive</span>.<br>
            Motif EPIYA détecté: <span class="italic">non implémenté</span>.
        """

def _html_pylori_vacA(dataframe, sample_name):
    return """
        Sous-type <span class="bold">vacA</span>: <span class="italic">non implémenté</span>
    """

_HTML_STYLES = """
<style>
    body { font-family: "Space Grotesk", sans-serif; }
    .mono { font-family: "Space Mono", monospace; }
    .success-text { color: #00473e; }
    .warning-text { color: #faae2b; }
    .error-text { color: #fa5246; }
    .big { font-size: 1.25rem; }
    .subscript { font-size: 0.75rem; color: #444; }
    .bold { font-weight: bold; }
    .italic { font-style: italic; }
    .underline { text-decoration: underline; }
    
    .grid-1-2 {
        display: grid;
        grid-template-columns: repeat(2, 1fr);
        grid-template-rows: 1fr;
        grid-column-gap: 3px;
        grid-row-gap: 0px;
    } 
    .grid-1-2-1-2 { grid-area: 1 / 1 / 2 / 2; }
    .grid-1-2-1-1 { grid-area: 1 / 2 / 2 / 3; } 
    
    .cover-depth-table {
        width: 100%;
        background-color: #eee;
        border-collapse: collapse;
    }
    .cover-depth-table th {
        background-color: #000;
        color: white;
        width: 33%
    }
    .cover-depth-table td, .cover-depth-table th {
        padding: 5px; 
        border: 1px solid #000;
    }
    .cover-depth-table td {
        font-family: "Space Mono", monospace;
    }
</style>
"""


def to_html(sample_name, pipeline_type):
    if pipeline_type == "pylori":
        genome = "H.pylori.J99.2001"
    else:
        raise NotImplementedError()

    dataframe = dfio.get_or_create(pipeline_type)
    gene_names = bed._load_names(
        common.get_genome_regions(pipeline_type) + ".names.list"
    )
    depth_img_origin = f"output/depth/{genome}/{sample_name}.png"
    depth_img_dest = f"output/html_reports/{sample_name}/depth.png"

    table_html_str = _html_depth_coverage_table(dataframe, sample_name, gene_names)
    html_str = f"""
    <!DOCTYPE html>
    <html lang="fr">
        <head>
            <meta charset="UTF-8">
            <title>Rapport {sample_name}</title>
            <link rel="preconnect" href="https://fonts.googleapis.com">
            <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
            <link href="https://fonts.googleapis.com/css2?family=Space+Grotesk:wght@300;400;500;600;700&family=Space+Mono:ital,wght@0,400;0,700;1,400;1,700&display=swap" rel="stylesheet">
            {_HTML_STYLES}
        </head>
        <body>
            <h1>Rapport d'analyses génotypiques DNA Capture</h1>
            <h3>Type de pipeline: {pipeline_type}</h3>
            <h3>Nom de souche: {sample_name}</h3>
            <hr style="margin-bottom: 10px;"/>
            <h3>Analyse des couvertures et profondeurs moyennes par gènes</h3>
            <div class="grid-1-2">
                <div class="grid-1-2-1-2">
                    {table_html_str}
                </div>
                <div class="grid-1-2-1-1">
                    <a href="{sample_name}/depth.png" target="_blank">
                        <img src="{sample_name}/depth.png" width="100%"/>
                    </a>
                </div>
            </div>
            <hr style="margin-bottom: 10px;"/>
            <h3>Analyse des résistances aux antibiotiques</h3>
            {_html_pylori_amr(dataframe, sample_name)}
            <hr style="margin-bottom: 10px;"/>
            <h3>Typage MLST</h3>
            {_html_pylori_mlst(dataframe, sample_name)}
            <hr style="margin-bottom: 10px;"/>
            <h3>Analyse de la virulence</h3>
            {_html_pylori_cagA(dataframe, sample_name)}
            <hr style="margin-bottom: 10px;"/>
            <h3>Analyse phylogénétique</h3>
        </body>
    </html>\n
    """

    if not os.path.isdir(f"output/html_reports/{sample_name}"):
        os.makedirs(f"output/html_reports/{sample_name}")
    with open(
        f"output/html_reports/{sample_name}.html", "w", encoding="UTF-8"
    ) as html_file:
        html_file.write(html_str)

    shutil.copyfile(depth_img_origin, depth_img_dest)

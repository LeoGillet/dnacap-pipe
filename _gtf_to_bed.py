"""
"""

INPUT_FILE = "annotations/genitalium/G37.gff"
OUTPUT_FILE = "probes/genitalium_CDS.bed"

with open(INPUT_FILE, "r", encoding="UTF-8") as gff:
    lines = [line.strip() for line in gff.readlines()]

bed_lines = []

for line in lines:
    if line == "###":
        break
    if line.startswith("#") or "RefSeq\tgene" not in line:
        continue

    scaffold, _, _, start_pos, end_pos, _, _, _, _info = line.split("\t")
    for field in _info.split(";"):
        if "Name=" in field:
            gene_name = field.replace("Name=", "")
            break
    bed_lines.append(f"{scaffold}\t{start_pos}\t{end_pos}\t{gene_name}\n")

with open(OUTPUT_FILE, "w", encoding="UTF-8") as bedfile:
    bedfile.writelines(bed_lines)

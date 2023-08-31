---
author: "Léo Gillet <leo.gillet@chu-bordeaux.fr>"
---
⚠️ This project is currently in development and not all of its objectives have been implemented yet. 

# Python workflow for the analysis of DNA hybridization capture sequencing results

### Requirements :
- Python (>= 3.10)
- bcftools == 1.18
- samtools == 1.18
- bedtools == 2.31.0 
- MultiQC == 1.14
- Bowtie2 == 2.5.1
- fastp == 0.23.4
- FastQC == 0.12.1
- seqtk == 1.4-r130-dirty

### TODO :
- ✅ Logging
- ✅ Paths and folder checks
- ✅ Sequence parsing
- ✅ Automatic reference genome indexing
- ✅ FastQC report creation
- ✅ Sequence trimming (fastp)
- ❌ _De novo_ assembly
- ✅ Decontamination
- ✅ Mapping
- ✅ MultiQC report creation
- ✅ Mapping coverage and depth analysis
- ✅ Whole genome & region of intersect consensus extraction from mapping 

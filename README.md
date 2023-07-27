---
title: "dnacap-pipe Readme"
author: "Léo Gillet <leo.gillet@chu-bordeaux.fr>"
---
⚠️ This project is currently in development and not all of its objectives have been implemented yet. 

# Python workflow for the analysis of DNA hybridization capture sequencing results

### Requirements :
- Python (>= 3.10)
- MultiQC == 1.14
- Bowtie2 == 2.5.1
- BWA == 0.7.17-r1198-dirty
- FastQC == 0.12.1
- Sickle == 1.33

### TODO :
- ✅ Logging
- ✅ Paths and folder checks
- ✅ Sequence parsing
- ✅ Automatic reference genome indexing
- ✅ FastQC report creation
- ✅ Sequence trimming (sickle)
- ▶️ _De novo_ assembly
- ❌ Decontamination
- ❌ Mapping
- ✅ MultiQC report creation

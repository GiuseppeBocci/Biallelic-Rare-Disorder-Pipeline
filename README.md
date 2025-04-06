# Biallelic-Rare-Disorder-Pipeline
Genomics 2025/2026 Project a bit enhanced :)

Professor: [Chiara Matteo](https://www.unimi.it/it/ugov/person/matteo-chiara)

# MANUAL
-s => steps

-g => hg19

`./pipeline.sh -sg study` (use multiple threads)

for multiple studies (studies are elaborated one after another usign nultiple threads for each one)

`./pipeline.sh -sg $(ls | greap "study_")`

N.B. At the moement pipeline.sh = script2.sh

# THIS PROJECT IS NOT FINISHED
TODO:
- fix reports
- multiallelic sites
- add bed files
- pipeline in nexflow/snakemake

# Biallelic-Rare-Disorder-Pipeline
Genomics 2025/2026 Project a bit enhanced :)

Professor: [Chiara Matteo](https://www.unimi.it/it/ugov/person/matteo-chiara)

# PROS
- No need to install vep tool, rest-api used

# MANUAL
-s => steps

-g => hg19

`./pipeline.sh -sg study` (use multiple threads)

for multiple studies (studies are elaborated one after another usign multiple threads for each one)

`./pipeline.sh -sg $(ls | grep "study_")` => `./pipeline.sh -sg study_1 study_2 ...`

N.B. At the moment pipeline.sh = script2.sh

# THIS PROJECT IS NOT FINISHED
TODO:
- fix reports
- multiallelic sites
- add bed files
- pipeline in nexflow/snakemake

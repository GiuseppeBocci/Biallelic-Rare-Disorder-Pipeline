#!/bin/bash

# Default values
USE_STEPS=false
HG19=false
HTML=false # TODO: deprecated
CASE=""
TARGET="exons16_merged.bed"

# Parse options using getopts
while getopts ":sgH" opt; do
  case ${opt} in
    s )
      USE_STEPS=true
      ;;
    g )
      HG19=true
      ;;
    H )
      HTML=true
      ;;
    \? )
      echo "Usage: $0 [-s] [-g] case_id"
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))  # Shift to get positional arguments

# Ensure CASE is provided
if [[ -z "$1" ]]; then
    echo "Error: CASE_ID is required"
    echo "Usage: $0 [-s] [-g] [-H] case_id"
    exit 2
fi


# Trap Ctrl+C (SIGINT) and exit
trap 'echo "Interrupt received. Exiting script..."; exit 3' SIGINT
# Exit immediately if any command fails
set -e
trap 'echo "Error: Command failed on line $LINENO"; exit 4' ERR
# Stop if any command in a pipeline fails
set -o pipefail

for file_arg in $@; do

export CASE=$file_arg
STEPS_FILE="$CASE/steps"

echo -en "\n#####"
for i in $(seq 1 ${#CASE}); do echo -n '#'; done
echo "#####"
echo "#    $CASE    #"
echo -n "#####"
for i in $(seq 1 ${#CASE}); do echo -n '#'; done
echo -e "#####\n"


# Initialize steps file if needed
if [[ $USE_STEPS == true ]]; then
    touch "$STEPS_FILE"
else
    rm -f $STEPS_FILE
fi

# Function to check if a step is completed
is_step_done() {
    if [[ $USE_STEPS == true ]]; then
        grep -q "^$1$" "$STEPS_FILE"
    else
        return 1  # Always return false (run all steps)
    fi
}

# Function to mark a step as completed
mark_step_done() {
    echo "$1" >> "$STEPS_FILE"
}

order_priotirization() {
   local -n names=$1
   local parent=$2
   local child=$3
   local priotirization=""

   if [[ ${names[0]} == "child" ]]; then
       priotirization="\"$child.*$parent.*$parent\""
   elif [[ ${names[1]} == "child" ]]; then
       priotirization="\"$parent.*$child.*$parent\""
   else
       priotirization="\"$parent.*$parent.*$child\""
   fi

   echo "$priotirization"
}


###### Step 1: Alignment
STEP1_N="Alignment"
if ! is_step_done $STEP1_N; then
    echo "Starting Step 1: Alignment"

    ls $CASE/*.fq.gz | parallel 'base=$(basename {} .fq.gz); bowtie2 -x align/uni -U {} --rg-id "${base}" --rg "SM:${base}" | samtools view -Sb | samtools sort -o $CASE/${base}.bam && samtools index $CASE/${base}.bam'

    echo "Bam files created for $CASE"
    mark_step_done $STEP1_N
else
    echo "Step 1 already completed, skipping..."
fi


###### Step 2: Fastqc & Qualimap reports
STEP2_N="Reports"
if ! is_step_done $STEP2_N; then
	echo "Starting Step 2: Reports generation"

	ls $CASE/*.fq.gz | parallel "fastqc {}"
	ls $CASE/*.bam | parallel "qualimap bamqc -bam {} -outdir $CASE/{/.}"
	# Aggregate results
	multiqc -f $CASE -o $CASE

    mark_step_done $STEP2_N
else
    echo "Step 2 already completed, skipping..."
fi


###### Step 3: FreeBayes Variant Calling
STEP3_N="freebayes"
if ! is_step_done $STEP3_N; then
    echo "Starting Step 3: FreeBayes Variant Calling"

    freebayes -f align/universe.fasta -m 20 -C 5 -Q 10 --min-coverage 10 --targets $TARGET $CASE/child.bam $CASE/father.bam $CASE/mother.bam > $CASE/joint.vcf

    echo "FreeBayes ended for $CASE"
    mark_step_done $STEP3_N
else
    echo "Step 3 already completed, skipping..."
fi

###### Step 4: Filtering
STEP4_N="Filtering"
PRIORITIZATION=$(grep $CASE metadata.tsv | cut -f2)
if ! is_step_done $STEP4_N; then
    echo "Starting Step 4: Filtering by quality"

    if [[ $PRIORITIZATION == "AR" ]]; then
        bcftools filter -i 'QUAL > 30' $CASE/joint.vcf -o $CASE/joint.vcf.tmp && mv $CASE/joint.vcf.tmp $CASE/joint.vcf
    else
        echo "No filtering"
    fi
    mark_step_done $STEP4_N
else
    echo "Step 4 already completed, skipping..."
fi

##### Step 5: Priotirization
STEP5_N="Prioritization"
if ! is_step_done $STEP5_N; then
    echo "Starting Step 5: Prioritization"

    TRIO_NAMES=$(grep "^#CHR" $CASE/joint.vcf | rev | cut -f1-3 | rev)
    NO_PRIO=false
    COMMAND="grep "

    if [[ $PRIORITIZATION == "AR" ]]; then
        COMMAND+=$(order_priotirization TRIO_NAMES "0/1" "1/1") # 0/1.*0/1.*1/1
    elif [[ $PRIORITIZATION == "DN" ]]; then
        COMMAND+=$(order_priotirization TRIO_NAMES "0/0" "0/1")  # 0/0.*0/0.*0/1
    else
        echo "No prioritization strategy selected for $CASE"
        NO_PRIO=true
    fi

    COMMAND+=" $CASE/joint.vcf | cut -f 1-6 > $CASE/possible.vcf"

    if [[ $NO_PRIO == false ]]; then
        echo "Priotirization stategy is $PRIORITIZATION"
	echo $COMMAND
        eval "$COMMAND"
    else
        cp "$CASE/joint.vcf" "$CASE/possible.vcf"
    fi

    mark_step_done $STEP5_N
else
    echo "Step 5 already completed, skipping..."
fi

###### Step 6: Request
STEP6_N="Request"
if ! is_step_done $STEP6_N; then
	echo "Starting Step 6: Request to Ensembl"

    URL="/vep/homo_sapiens/region"
    if [[ $HG19 == true ]]; then
        URL="https://grch37.rest.ensembl.org$URL"
    else
        URL="https://rest.ensembl.org$URL"
    fi

    json_array=()

    # Read each line from the file
    while IFS= read -r line; do
        # Escape double quotes by replacing them with a backslash followed by a double quote
        # sanitized_line=$(echo "$line" | sed 's/"/\\\"/g')

        # Append the sanitized line to the array (with proper quotes)
        json_array+=("\"$line\"")
    done < "$CASE/possible.vcf"

    # Join the array elements into a single string (separated by commas)
    BODY="{\"variants\":["$(IFS=,; echo "${json_array[*]}")"]}"
    echo $BODY > $CASE/req.json.tmp

    URL+="?refseq=1&Phenotypes=1&hgvs=1"

    # echo "curl -X POST \"$URL\" -H \"Accept: application/json\" -H \"Content-Type: application/json\" -d '$BODY'"


  #  if [[ $HTML == true ]]; then
  #  	RESPONSE=$(curl -X POST "$URL" -H "Accept: text/html" -H "Content-Type: application/json" -d @$CASE/req.json.tmp)
  #  	rm $CASE/req.json.tmp
  #  	echo $RESPONSE > $CASE/res.html
  #  else
        RESPONSE=$(curl -X POST "$URL" -H "Accept: application/json" -H "Content-Type: application/json" -d @$CASE/req.json.tmp)
        rm $CASE/req.json.tmp
        echo $RESPONSE > $CASE/res.json
  # fi

    mark_step_done $STEP6_N
else
    echo "Step 6 already completed, skipping..."
fi


###### Step 7: Selection
STEP7_N="Selection"
if ! is_step_done $STEP7_N; then
    echo "Starting Step 7: Annotation Selection"

    python3 vep_json_parser.py $CASE/res.json > $CASE/selected &
    pid_py=$!

    wait $pid_py

	echo "Selected:"
	cat $CASE/selected

    #mark_step_done $STEP7_N
else
    echo "Step 7 already completed, skipping..."
fi

done

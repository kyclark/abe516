#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -p development
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J getfasta
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user kyclark@email.arizona.edu

module load launcher
module load sratoolkit

NUM=10000000
ACCS="stage4.txt"
OUT_DIR="$PWD/stage4"
PARAM=$(mktemp)

if [[ ! -f "$ACCS" ]]; then
    echo "Missing ACCS \"$ACCS\""
    exit 1
fi

[[ ! -d "$OUT_DIR" ]] && mkdir -p "$OUT_DIR"

while read -r ACC; do
    #if [[ ! -s "$OUT_DIR/$ACC.fasta" ]]; then
        echo "fastq-dump -X $NUM --fasta -O $OUT_DIR $ACC" >> "$PARAM"
    #fi
done < "$ACCS"

NUM_JOBS=$(wc -l "$PARAM" | awk '{print $1}')

echo "Launching NUM_JOBS \"$NUM_JOBS\""

PARAMRUN="$TACC_LAUNCHER_DIR/paramrun"

export LAUNCHER_PLUGIN_DIR="$TACC_LAUNCHER_DIR/plugins"
export LAUNCHER_WORKDIR="$PWD"
export LAUNCHER_RMI="SLURM"
export LAUNCHER_SCHED="interleaved"
export LAUNCHER_JOB_FILE="$PARAM"
export LAUNCHER_PPN=16
$PARAMRUN

rm "$PARAM"

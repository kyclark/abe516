#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -p normal
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J align
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user kyclark@email.arizona.edu

set -u

date
echo "Starting alignment"

module load launcher
module load star

OUT_DIR="$PWD/star-out-stage4"
FASTA=$(mktemp)

[[ ! -d "$OUT_DIR" ]] && mkdir -p "$OUT_DIR"

#find "$PWD/fasta" -type f > "$FASTA"
find "$PWD/stage4" -type f > "$FASTA"

NUM_FILES=$(wc -l "$FASTA" | awk '{print $1}')
echo "Found NUM_FILES \"$NUM_FILES\" in fasta dir"

PARAM=$(mktemp)

i=0
while read FILE; do
    let i++
    BASENAME=$(basename "$FILE" ".fasta")
    printf "%3d: %s\n" $i $BASENAME

    DIR="$OUT_DIR/$BASENAME"
    [[ ! -d "$DIR" ]] && mkdir -p "$DIR"

    SAM_FILE="$DIR/Aligned.out.sam"

    if [[ ! -f "$SAM_FILE" ]]; then
        echo STAR --runThreadN 16 \
            --genomeDir $PWD/sorghum \
            --sjdbGTFfile $PWD/genome/sorghum.gtf \
            --readFilesIn $FILE \
            --outFileNamePrefix "$DIR/" >> "$PARAM"
    fi
done < "$FASTA"

NUM_JOBS=$(wc -l "$PARAM" | awk '{print $1}')

echo "Launching NUM_JOBS \"$NUM_JOBS\""
if [[ $NUM_JOBS -lt 4 ]]; then
    LAUNCHER_PPN=$NUM_JOBS
else
    LAUNCHER_PPN=4
fi

export LAUNCHER_PPN
export LAUNCHER_PLUGIN_DIR="$TACC_LAUNCHER_DIR/plugins"
export LAUNCHER_WORKDIR="$PWD"
export LAUNCHER_RMI="SLURM"
export LAUNCHER_SCHED="interleaved"
export LAUNCHER_JOB_FILE="$PARAM"
$TACC_LAUNCHER_DIR/paramrun

rm "$FASTA"
rm "$PARAM"

date
echo "Done"

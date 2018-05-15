#!/bin/bash

module load samtools

set -u

TMP=$(mktemp)
find star-out-stage4 -name \*.sam -size +0c > "$TMP"

NUM=$(wc -l "$TMP" | awk '{print $1}')

echo "Found NUM \"$NUM\" SAM files"

BAM_DIR="bamfiles"

[[ ! -d "$BAM_DIR" ]] && mkdir -p "$BAM_DIR"

i=0
while read -r SAM_FILE; do
    ACC=$(basename $(dirname "$SAM_FILE"))
    let i++
    printf "%3d: %s\n" $i $ACC
    samtools view -bS -o "$BAM_DIR/$ACC.bam" "$SAM_FILE"
    #samtools sort "$SAM_FILE" "$SAM_FILE.sorted"
done < "$TMP"

rm "$TMP"

echo "Done."

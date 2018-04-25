#!/bin/bash

i=0
while read -r FILE; do
    BASENAME=$(basename $FILE)
    if [[ ! -s "$BASENAME" ]]; then
        let i++
        printf "%4d: %s\n" $i $BASENAME
        iget $FILE
    fi
done < mash-files

echo "Done."

#!/bin/bash

#SBATCH -A iPlant-Collabs
#SBATCH -p development
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -J mk-srghm
#SBATCH --mail-type BEGIN,END,FAIL
#SBATCH --mail-user kyclark@email.arizona.edu

set -u

OUT_DIR=$PWD/sorghum

if [[ -d "$OUT_DIR" ]]; then
    rm -rf $OUT_DIR/*
else
    mkdir -p "$OUT_DIR"
fi

date
echo "Making sorghum into $OUT_DIR"

module load star

STAR --runThreadN 16 \
    --runMode genomeGenerate \
    --genomeDir $OUT_DIR \
    --genomeFastaFiles $PWD/genome/Sorghum_bicolor.Sorghum_bicolor_NCBIv3.dna.toplevel.fa \
    --sjdbGTFfile $PWD/genome/sorghum.gtf

date
echo "Done"

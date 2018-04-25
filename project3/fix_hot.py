#!/usr/bin/env python3
"""docstring"""

import os
import re
import sys

args = sys.argv[1:]

if len(args) != 1:
    print('Usage: {} DIST > OUT'.format(os.path.basename(sys.argv[0])))
    sys.exit(1)

dist_file = args[0]

names = """
HOT194_500m,SRR1303812.fasta
HOT194_25m,SRR1303791.fasta
HOT198_25m,SRR1303795.fasta
HOT198_25m,SRR1303840.fasta
HOT195_25m,SRR1303792.fasta
HOT196_25m,SRR1303793.fasta
HOT199_25m,SRR1303796.fasta
HOT200_25m,SRR1303797.fasta
HOT201_25m,SRR1303798.fasta
HOT202_25m,SRR1303799.fasta
HOT203_25m,SRR1303800.fasta
HOT204_25m,SRR1303801.fasta
HOT205_25m,SRR1303802.fasta
HOT206_25m,SRR1303803.fasta
HOT208_25m,SRR1303804.fasta
HOT209_25m,SRR1303805.fasta
HOT210_25m,SRR1303806.fasta
HOT211_25m,SRR1304300.fasta
HOT212_25m,SRR1303808.fasta
HOT213_25m,SRR1303809.fasta
HOT214_25m,SRR1303810.fasta
HOT215_25m,SRR1303811.fasta
HOT197_25m,SRR1303794.fasta
HOT195_500m,SRR1303813.fasta
HOT196_500m,SRR1303814.fasta
HOT197_500m,SRR1303815.fasta
HOT197_500m,SRR1303839.fasta
HOT198_500m,SRR1303816.fasta
HOT199_500m,SRR1303817.fasta
HOT200_500m,SRR1303818.fasta
HOT201_500m,SRR1303819.fasta
HOT202_500m,SRR1303820.fasta
HOT203_500m,SRR1303821.fasta
HOT204_500m,SRR1303822.fasta
HOT205_500m,SRR1303823.fasta
HOT206_500m,SRR1303824.fasta
HOT208_500m,SRR1303825.fasta
HOT209_500m,SRR1303826.fasta
HOT210_500m,SRR1303827.fasta
HOT211_500m,SRR1303828.fasta
HOT212_500m,SRR1303829.fasta
HOT213_500m,SRR1303830.fasta
HOT214_500m,SRR1303831.fasta
HOT215_500m,SRR1303832.fasta
"""

hot_names = {}
for line in names.splitlines():
    if not line:
        continue

    depth, srr = line.rstrip().split(',')
    hot_names[srr] = depth

def goodname(file):
    return hot_names[os.path.basename(file)]

for i, line in enumerate(open(dist_file, 'rt')):
    flds = line.rstrip().split('\t')
    first = flds.pop(0)

    if i == 0:
        print('\t'.join([''] + list(map(goodname, flds))))
    else:
        print('\t'.join([goodname(first)] + flds))

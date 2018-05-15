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

species_file = '/Users/kyclark/work/fizkin-paper/icu/mash/srr-species.csv'
names = {}
for line in open(species_file):
    srr, species_name = line.rstrip().split(',')
    names[srr] = species_name.replace(' ', '_')

def goodname(srr):
    return names[srr] + '_' + srr

for i, line in enumerate(open(dist_file, 'rt')):
    flds = line.rstrip().split('\t')
    first = flds.pop(0)

    if i == 0:
        print('\t'.join([''] + list(map(goodname, flds))))
    else:
        print('\t'.join([goodname(first)] + flds))

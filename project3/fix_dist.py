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

biomes = {}
for line in open('biome.txt'):
    sample_id, biome = line.rstrip().split('\t')
    biomes[sample_id] = biome

regex = re.compile(r'\/samples\/(\d+)\/')
def idfromname(filename):
    match = regex.search(filename)
    if match:
        sid = match.groups()[0]
        biome = biomes[sid].replace(' ', '_')
        return biome + '_' + sid

for i, line in enumerate(open(dist_file, 'rt')):
    flds = line.rstrip().split('\t')
    first = flds.pop(0)
    if i == 0:
        sample_names = list(map(idfromname, flds))
        print('\t'.join([''] + sample_names))
    else:
        print('\t'.join([idfromname(first)] + flds))

#!/usr/bin/env python3

import re

ids = set(open('metagenomes.txt').read().splitlines())

regex = re.compile(r'\/samples\/(\d+)\/(.+)\n')
for line in open('mash/mash-files.txt'):
    match = regex.search(line)
    if match:
        sample_id, sketch = match.groups()
        if sample_id in ids:
            print('mv mash/{} mash/{}.msh'.format(sketch, sample_id))

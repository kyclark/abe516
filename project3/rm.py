#!/usr/bin/env python3

import re

ids = set(open('metagenomes.txt').read().splitlines())

regex = re.compile(r'\/samples\/(\d+)\/(.+)\n')
for line in open('mash/mash-files.txt'):
    match = regex.search(line)
    if match:
        sample_id, sketch = match.groups()
        if not sample_id in ids:
            print('rm mash/{}'.format(sketch))

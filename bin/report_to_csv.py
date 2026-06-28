#!/usr/bin/env python3
import os, sys

for fname in os.listdir('.'):
    if not fname.endswith('.report'):
        continue
    out = fname.rsplit('.report', 1)[0] + '.csv'
    with open(fname) as f:
        lines = f.readlines()

    with open(out, 'w') as fout:
        if 'bracken' in fname:
            fout.write(lines[0])
            for line in lines[1:]:
                fout.write(line)
        else:
            headers = ["percentage", "num_reads", "num_unique_reads", "rank_code", "taxonomy_id", "name"]
            fout.write('\t'.join(headers) + '\n')
            for line in lines:
                if not line.strip():
                    continue
                parts = line.rstrip('\n').split('\t')
                if len(parts) >= 5:
                    name = parts[-1]
                    if len(parts) > 5:
                        name = '\t'.join(parts[5:])
                    fout.write('\t'.join(parts[:5] + [name]) + '\n')

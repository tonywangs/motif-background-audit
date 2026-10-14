#!/usr/bin/env python3
"""Preregistered finite-space frequency diagnostic; not a uniformity proof."""
from collections import Counter
import itertools
import json
import math
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from motif_audit.core import canonical, digest
from motif_audit.dinucleotide import shuffle_dinucleotide


def signature(s):
    return s[0],s[-1],tuple(sorted(Counter(zip(s,s[1:])).items()))


def space(sequence):
    # Enumerate unique multiset permutations recursively, without shuffler code.
    counts=Counter(sequence); result=[]
    def visit(prefix):
        if len(prefix)==len(sequence):
            if signature(prefix)==signature(sequence): result.append(prefix)
            return
        for base in sorted(counts):
            if counts[base]:
                counts[base]-=1; visit(prefix+base); counts[base]+=1
    visit('')
    return result


def check():
    protocol=(ROOT/'experiments/dinucleotide-protocol.json').read_bytes()
    p=json.loads(protocol)['sampling_check']
    spaces={s:space(s) for s in p['distribution_sequences']}
    total=sum(map(len,spaces.values())); n=p['draws_per_sequence']
    threshold=math.sqrt(n/2*math.log(2*total/.001))
    rows=[]
    for s,valid in spaces.items():
        counts=Counter(shuffle_dinucleotide(s,p['seed_start']+i) for i in range(n))
        assert set(counts)<=set(valid)
        row={'sequence':s,'outcomes':{v:counts[v] for v in valid},'draws':n,
             'expected_per_outcome':n/len(valid),'absolute_deviation_limit':threshold,
             'maximum_absolute_deviation':max(abs(counts[v]-n/len(valid)) for v in valid)}
        row['passed']=row['maximum_absolute_deviation']<=threshold
        rows.append(row)
    return {'protocol_sha256':digest(protocol),'complete':True,'passed':all(r['passed'] for r in rows),
            'interpretation':'Empirical check only, not a proof of uniformity. Fixed seed streams are pseudorandom.',
            'checks':rows}


if __name__=='__main__':
    result=check()
    print(canonical(result).decode(),end='')
    if not result['passed']: sys.exit(1)

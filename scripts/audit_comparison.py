#!/usr/bin/env python3
"""Recompute every paired count and summary using independent window matching."""
from collections import Counter
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import re
import sys

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT),str(ROOT/'tests')]
from test_core import oracle
from motif_audit.core import shuffle_sequence
from motif_audit.dinucleotide import shuffle_dinucleotide


def sha(data): return hashlib.sha256(data).hexdigest()


def adjusted(ps):
    # Closed rank-prefix formula with rational arithmetic, independently coded.
    order=sorted(range(len(ps)),key=lambda i:(ps[i],i))
    return {i:float(min(Fraction(1),max(ps[j]*(len(ps)-rank)
                 for rank,j in enumerate(order[:position+1])))) for position,i in enumerate(order)}


def invariants(before,after):
    assert len(before)==len(after)
    assert [(i,c) for i,c in enumerate(before) if c not in 'ACGT']==[(i,c) for i,c in enumerate(after) if c not in 'ACGT']
    for m in re.finditer('[ACGT]+',before):
        a,b=m[0],after[m.start():m.end()]
        assert a[0]==b[0] and a[-1]==b[-1]
        assert Counter(a)==Counter(b)
        assert Counter(zip(a,a[1:]))==Counter(zip(b,b[1:]))


def close(a,b):
    if a is None or b is None: assert a==b
    elif isinstance(a,(float,int)) and isinstance(b,(float,int)):
        assert math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-12),(a,b)
    else: assert a==b,(a,b)


def audit(directory):
    directory=Path(directory)
    frozen=json.loads((ROOT/'experiments/dinucleotide-frozen.json').read_bytes())
    p=json.loads((ROOT/'experiments/dinucleotide-protocol.json').read_bytes())
    summary=json.loads((directory/'summary.json').read_bytes())
    assert summary['freeze_sha256']==sha((ROOT/'experiments/dinucleotide-frozen.json').read_bytes())
    complete=json.loads((directory/'completion.json').read_bytes())
    for name,meta in complete['files'].items():
        b=(directory/name).read_bytes(); assert meta=={'bytes':len(b),'sha256':sha(b)}
    actual_names={str(f.relative_to(directory)) for f in directory.rglob('*') if f.is_file()}
    assert actual_names==set(complete['files'])|{'completion.json'}
    expected_rows=[]; tests=controls=observations=0
    for d in frozen['datasets']:
        a=json.loads((directory/d['name']/'analysis.json').read_bytes())
        assert a['input_sha256']=={'fasta':d['sha256'],'manifest':frozen['manifest_sha256']}
        assert a['statistics']==a['backgrounds']['mononucleotide']['statistics']
        assert a['matches']==sorted(a['matches'],key=lambda m:(m['record'],m['start'],m['motif'],m['strand']))
        expected_matches=[]
        for ri,r in enumerate(a['records']):
            for mi,m in enumerate(a['motifs']):
                for start,end,strand in oracle(r['sequence'],m['pattern']):
                    expected_matches.append({'record':ri,'motif':mi,'start':start,'end':end,'strand':strand})
        assert sorted(expected_matches,key=lambda m:(m['record'],m['start'],m['motif'],m['strand']))==a['matches']
        pair=[]; pair_ps=[]
        for bg,b in a['backgrounds'].items():
            assert int(b['settings']['seed'])==d['control_seed']
            assert b['settings']['shuffles']==p['shuffles']
            observed=[sum(len(oracle(r['sequence'],m['pattern'])) for r in a['records']) for m in a['motifs']]
            counts=[[] for _ in a['motifs']]
            for rep,seeds in enumerate(b['control_seeds']):
                sequences=[]
                for ri,(record,seed) in enumerate(zip(a['records'],seeds)):
                    expected_seed=int.from_bytes(hashlib.sha256(f"motif-audit-v1:{d['control_seed']}:{rep}:{ri}".encode()).digest()[:8],'big')
                    assert str(expected_seed)==seed
                    s=record['sequence']
                    t=(shuffle_dinucleotide if bg=='dinucleotide' else shuffle_sequence)(s,int(seed))
                    assert Counter(s)==Counter(t)
                    if bg=='dinucleotide': invariants(s,t)
                    sequences.append(t)
                for mi,m in enumerate(a['motifs']):
                    counts[mi].append(sum(len(oracle(s,m['pattern'])) for s in sequences)); controls+=1
            ps=[]; group=[]
            for i,stat in enumerate(b['statistics']):
                obs=observed[i]; cc=counts[i]; assert len(cc)==p['shuffles']
                assert cc==stat['control_counts']
                mean=Fraction(sum(cc),len(cc)); ge=sum(c>=obs for c in cc)
                pv=Fraction(ge+1,len(cc)+1); ps.append(pv)
                expected={'observed':obs,'control_mean':float(mean),'exceedances':ge,
                          'p_numerator':ge+1,'p_denominator':len(cc)+1,'p_upper':float(pv),
                          'difference':float(obs-mean),'fold':float(obs/mean) if mean else None,
                          'fold_note':None if mean else 'undefined: control mean is zero'}
                for key,value in expected.items(): close(stat[key],value)
                row={'dataset':d['name'],'condition':d['condition'],'background':bg,'motif':a['motifs'][i]['id'],**expected}
                group.append(row); observations+=1
            adj=adjusted(ps)
            for i,row in enumerate(group):
                row['p_holm']=adj[i]; close(b['statistics'][i]['p_holm'],adj[i])
            pair+=group; pair_ps+=ps
        adj=adjusted(pair_ps)
        for i,row in enumerate(pair):
            row['p_holm_comparison']=adj[i]
            bg=a['backgrounds'][row['background']]
            mi=[m['id'] for m in a['motifs']].index(row['motif'])
            close(bg['statistics'][mi]['p_holm_comparison'],adj[i])
        expected_rows+=pair
    study=adjusted([Fraction(r['p_numerator'],r['p_denominator']) for r in expected_rows])
    # JSON canonical sorting changes background dictionary order; compare by key.
    for i,row in enumerate(expected_rows): row['p_holm_study']=study[i]
    key=lambda r:(r['dataset'],r['background'],r['motif'])
    assert len(expected_rows)==len(summary['rows'])==128
    for expected,actual in zip(sorted(expected_rows,key=key),sorted(summary['rows'],key=key)):
        assert set(expected)==set(actual)
        for field,value in expected.items(): close(actual[field],value)
    for condition,backgrounds in summary['rejections'].items():
        for bg,values in backgrounds.items():
            for field,count in values.items():
                assert count==sum(r[field]<=p['alpha'] for r in expected_rows if r['condition']==condition and r['background']==bg)
    return {'complete':True,'independent_observed_counts':observations,'independent_control_counts':controls,
            'summary_rows':len(expected_rows),'all_adjustments_recomputed':True,
            'limitation':'Replay uses production shufflers; invariants, enumeration, matching and arithmetic are independently checked.'}


if __name__=='__main__':
    print(json.dumps(audit(sys.argv[1] if len(sys.argv)>1 else ROOT/'results/dinucleotide'),sort_keys=True))

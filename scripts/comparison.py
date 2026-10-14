#!/usr/bin/env python3
"""Freeze or reproduce the paired-background experiment, entirely offline."""
import argparse
from collections import Counter
import json
from pathlib import Path
import resource
import shutil
import signal
import sys
import time

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from motif_audit.core import Budget, canonical, compare, digest, holm
from motif_audit.cli import write_output
from comparison_inputs import generate

PROTOCOL=ROOT/'experiments/dinucleotide-protocol.json'
FROZEN=ROOT/'experiments/dinucleotide-frozen.json'
INPUTS=ROOT/'experiments/dinucleotide-inputs'
SOURCES=['experiments/dinucleotide-protocol.json','scripts/comparison_inputs.py',
         'scripts/comparison.py','motif_audit/core.py','motif_audit/dinucleotide.py']


def hashes():
    return {name:digest((ROOT/name).read_bytes()) for name in SOURCES}


def freeze():
    if FROZEN.exists() or INPUTS.exists():
        raise ValueError('refusing to replace a freeze')
    p=json.loads(PROTOCOL.read_bytes())
    INPUTS.mkdir()
    manifest=canonical({'motifs':[{'id':s,'pattern':s} for s in p['motifs']]})
    (INPUTS/'motifs.json').write_bytes(manifest)
    datasets=[]
    for condition in p['conditions']:
        for repetition in range(p['datasets_per_condition']):
            name=f'{condition}-{repetition:02}'
            fasta,seed,plants=generate(p,condition,repetition)
            (INPUTS/(name+'.fasta')).write_bytes(fasta)
            datasets.append({'name':name,'condition':condition,'repetition':repetition,
                             'generation_seed':seed,'control_seed':p['control_seed']+len(datasets),
                             'planted_sites':plants,'sha256':digest(fasta)})
    frozen={'schema':'motif-background-audit/freeze/1','source_sha256':hashes(),
            'manifest_sha256':digest(manifest),'datasets':datasets}
    FROZEN.write_bytes(canonical(frozen))
    print('Frozen SHA256: '+digest(FROZEN.read_bytes()))


def evaluate(out):
    frozen=json.loads(FROZEN.read_bytes()); p=json.loads(PROTOCOL.read_bytes())
    assert frozen['source_sha256']==hashes(), 'source differs from preregistered freeze'
    manifest=(INPUTS/'motifs.json').read_bytes()
    assert digest(manifest)==frozen['manifest_sha256']
    for d in frozen['datasets']:
        assert digest((INPUTS/(d['name']+'.fasta')).read_bytes())==d['sha256']
        assert generate(p,d['condition'],d['repetition'])[0]==(INPUTS/(d['name']+'.fasta')).read_bytes()
    # CPU timer includes matching, shuffling, serialization, and file output.
    def timeout(signum,frame): raise TimeoutError('preregistered CPU budget exhausted')
    old=signal.signal(signal.SIGPROF,timeout)
    timer=signal.setitimer(signal.ITIMER_PROF,p['cpu_budget_seconds'])
    start=time.monotonic(); cpu=time.process_time()
    out.mkdir()
    try:
        rows=[]
        for d in frozen['datasets']:
            result=compare((INPUTS/(d['name']+'.fasta')).read_bytes(),manifest,
                           seed=d['control_seed'],shuffles=p['shuffles'])
            write_output(out/d['name'],result,Budget())
            for bg,analysis in result['backgrounds'].items():
                for motif,stat in zip(result['motifs'],analysis['statistics']):
                    rows.append({'dataset':d['name'],'condition':d['condition'],'background':bg,
                                 'motif':motif['id'], **{k:v for k,v in stat.items() if k!='control_counts'}})
        for row,q in zip(rows,holm([r['p_upper'] for r in rows])):
            row['p_holm_study']=q
        counts={c:{bg:{key:sum(r[key]<=p['alpha'] for r in rows if r['condition']==c and r['background']==bg)
                       for key in ['p_upper','p_holm','p_holm_comparison','p_holm_study']}
                   for bg in p['backgrounds']} for c in p['conditions']}
        summary={'freeze_sha256':digest(FROZEN.read_bytes()),'rows':rows,'rejections':counts,
                 'warning':'Synthetic diagnostic only; do not select background by significance.'}
        (out/'summary.json').write_bytes(canonical(summary))
        metrics={'elapsed_seconds':time.monotonic()-start,'cpu_seconds':time.process_time()-cpu,
                 'peak_process_rss_kib':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                 'python':sys.version,'runtime_dependencies':'Python standard library',
                 'datasets':len(frozen['datasets']),'tests':len(rows)}
        (out/'metrics.json').write_bytes(canonical(metrics))
        files={str(f.relative_to(out)):{'sha256':digest(f.read_bytes()),'bytes':f.stat().st_size}
               for f in sorted(out.rglob('*')) if f.is_file()}
        (out/'completion.json').write_bytes(canonical({'complete':True,'files':files}))
        print(json.dumps({'metrics':metrics,'rejections':counts},sort_keys=True))
    except BaseException:
        shutil.rmtree(out)
        raise
    finally:
        signal.setitimer(signal.ITIMER_PROF,*timer)
        signal.signal(signal.SIGPROF,old)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    group=parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--freeze',action='store_true')
    group.add_argument('--out',type=Path)
    args=parser.parse_args()
    if args.freeze: freeze()
    else: evaluate(args.out)

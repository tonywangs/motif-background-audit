#!/usr/bin/env python3
"""One-command reproduction and independent audit, optionally with offline Chromium."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time
import zipfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT))
from motif_audit.core import Budget, canonical, compare, digest
from motif_audit.cli import write_output
from audit_comparison import audit
from sampling_check import check as sampling_check
from verify import source_hashes, run


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--browser',action='store_true')
    parser.add_argument('--playwright-module')
    parser.add_argument('--browsers-path')
    parser.add_argument('--record',type=Path)
    args=parser.parse_args()
    start=time.monotonic(); sources=source_hashes()
    with tempfile.TemporaryDirectory(prefix='motif-compare-verify-') as tmp:
        scratch=Path(tmp); evidence=scratch/'evidence'; evidence.mkdir()
        # Preserve historic source/fixtures/results; fresh evidence goes elsewhere.
        baseline=[sys.executable,ROOT/'scripts/verify.py','--record',scratch/'legacy']
        if args.browser:
            baseline+=['--browser']
            if args.playwright_module: baseline+=['--playwright-module',args.playwright_module]
            if args.browsers_path: baseline+=['--browsers-path',args.browsers_path]
        verified=run(baseline,cwd=ROOT); print(verified.stdout.decode(),end='',flush=True)
        legacy=json.loads((scratch/'legacy/verification.json').read_bytes())
        (evidence/'legacy.json').write_bytes(canonical(legacy))
        # Store test output as structured evidence, not another .log file.
        (evidence/'tests.json').write_bytes(canonical({'complete':True,'unittest_output':(scratch/'legacy/tests.log').read_text()}))
        sampling=sampling_check(); assert sampling['passed']
        (evidence/'sampling.json').write_bytes(canonical(sampling))
        audited=audit(ROOT/'results/dinucleotide')
        (evidence/'audit.json').write_bytes(canonical(audited))
        reproduced=scratch/'experiment'
        run([sys.executable,ROOT/'scripts/comparison.py','--out',reproduced],cwd=ROOT)
        repeat={}
        for p in sorted((ROOT/'results/dinucleotide').rglob('*')):
            if not p.is_file() or p.parent==ROOT/'results/dinucleotide' and p.name in ('metrics.json','completion.json'): continue
            rel=p.relative_to(ROOT/'results/dinucleotide')
            assert p.read_bytes()==(reproduced/rel).read_bytes(),f'reproduction mismatch: {rel}'
            repeat[str(rel)]={'bytes':p.stat().st_size,'sha256':digest(p.read_bytes())}
        (evidence/'reproduction.json').write_bytes(canonical({'complete':True,'identical_files':repeat,
            'repeat_metrics':json.loads((reproduced/'metrics.json').read_bytes())}))
        run([sys.executable,ROOT/'scripts/install.py','--prefix',scratch/'installed'],cwd=ROOT)
        exe=scratch/'installed/bin/motif-audit'
        with zipfile.ZipFile(exe) as archive:
            assert archive.read('motif_audit/USHUFFLE_LICENSE.txt')==(ROOT/'motif_audit/USHUFFLE_LICENSE.txt').read_bytes()
        isolated=scratch/'isolated'; isolated.mkdir()
        guard=scratch/'guard'; guard.mkdir()
        (guard/'sitecustomize.py').write_text('import sys\ndef guard(event,args):\n if event.startswith("socket."): raise RuntimeError("offline verification")\nsys.addaudithook(guard)\n')
        env={**os.environ,'PYTHONPATH':str(guard),'PYTHONDONTWRITEBYTECODE':'1'}
        fasta=isolated/'input.fasta'; motifs=isolated/'motifs.json'
        fasta.write_bytes((ROOT/'examples/sequences.fasta').read_bytes()); motifs.write_bytes((ROOT/'examples/motifs.json').read_bytes())
        before={p.name:digest(p.read_bytes()) for p in (fasta,motifs)}
        base=[exe,'--fasta',fasta,'--motifs',motifs,'--background','both','--shuffles','199']
        timing=scratch/'timing.json'
        run(['/usr/bin/time','-f','{"elapsed_seconds":%e,"peak_process_rss_kib":%M}','-o',timing,*base,'--out',scratch/'one'],cwd=isolated,env=env)
        run([*base,'--out',scratch/'two'],cwd=isolated,env=env)
        outputs={}
        for p in sorted((scratch/'one').iterdir()):
            assert p.read_bytes()==(scratch/'two'/p.name).read_bytes()
            outputs[p.name]={'bytes':p.stat().st_size,'sha256':digest(p.read_bytes())}
        assert before=={p.name:digest(p.read_bytes()) for p in (fasta,motifs)}
        failure_checks=[]
        for extra in [[],['--background','bad'],['--shuffles','0'],['--seconds','0.000001']]:
            dest=scratch/'one' if not extra else scratch/'failure'
            result=subprocess.run([str(x) for x in [*base,'--out',dest,*extra]],capture_output=True,env=env,cwd=isolated,timeout=10)
            assert result.returncode==2
            if extra: assert not dest.exists()
            failure_checks.append({'arguments':extra,'exit_code':result.returncode})
        assert all(digest((scratch/'one'/n).read_bytes())==v['sha256'] for n,v in outputs.items())
        installed={'complete':True,'offline_socket_audit_hook':True,'isolated_cwd':True,
                   'unchanged_inputs':before,'repeated_outputs':outputs,'license_included':True,
                   'failures':failure_checks,**json.loads(timing.read_bytes())}
        (evidence/'installed.json').write_bytes(canonical(installed))
        shutil.copytree(scratch/'one',evidence/'sample')
        browser=None
        if args.browser:
            hostile='</script><img src="https://invalid.test/p" onerror="window.pwned=1"> @@JS@@ & <svg/onload=window.pwned=1>'
            raw=(f'>=hostile {hostile}\n'+'ACGTAC'*100+'\n>reverse\n'+'GTACGT'*100+'\n>palindrome\n'+'ATAT'*200+'\n>empty\nNNRYSW\n').encode()
            manifest=canonical({'motifs':[{'id':'target','pattern':'ACGTAC'},{'id':'palindrome','pattern':'ATAT'},{'id':'any_concrete','pattern':'N'}]})
            fixture=scratch/'fixture'; empty=scratch/'empty'
            write_output(fixture,compare(raw,manifest,shuffles=7),Budget())
            write_output(empty,compare(b'>ambiguous\nNNRYSW',manifest,shuffles=7),Budget())
            benv=dict(os.environ)
            if args.playwright_module: benv['MOTIF_PLAYWRIGHT_MODULE']=str(Path(args.playwright_module).resolve())
            if args.browsers_path: benv['PLAYWRIGHT_BROWSERS_PATH']=str(Path(args.browsers_path).resolve())
            result=run(['node',ROOT/'scripts/browser.cjs',fixture,empty,scratch/'one',evidence],cwd=ROOT,env=benv)
            print(result.stdout.decode(),end='',flush=True)
            browser=json.loads((evidence/'browser.json').read_bytes())
        assert sources==source_hashes(),'verification changed source files'
        overview={'complete':True,'elapsed_seconds':time.monotonic()-start,'source_sha256':sources,
                  'independent_audit':audited,'sampling_passed':sampling['passed'],
                  'reproduced_files':len(repeat),'installed':installed,'browser':browser,
                  'legacy_complete':legacy['complete']}
        (evidence/'verification.json').write_bytes(canonical(overview))
        if args.record:
            args.record.mkdir(parents=True,exist_ok=True)
            shutil.copytree(evidence,args.record,dirs_exist_ok=True)
        print(json.dumps({'complete':True,'elapsed_seconds':overview['elapsed_seconds'],
                          'reproduced_files':len(repeat),'audit':audited,'browser':bool(browser)},sort_keys=True))


if __name__=='__main__': main()

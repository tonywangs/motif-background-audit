#!/usr/bin/env python3
"""One-command verification, including an isolated installed CLI and optional Chromium."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from motif_audit.core import analyze, canonical, digest, Budget
from motif_audit.cli import write_output
from audit_experiment import audit


def run(args, **kwargs):
    result = subprocess.run([str(x) for x in args], capture_output=True, timeout=180, **kwargs)
    if result.returncode:
        raise RuntimeError(f"command failed ({result.returncode}): {args}\n{result.stdout.decode()}\n{result.stderr.decode()}")
    return result


def source_hashes():
    paths = [p for base in ['motif_audit','scripts','tests','docs','examples','experiments']
             for p in (ROOT/base).rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc']
    paths += [ROOT/p for p in ['README.md','package.json','package-lock.json','.gitignore'] if (ROOT/p).exists()]
    return {str(p.relative_to(ROOT)):digest(p.read_bytes()) for p in sorted(paths)}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--browser',action='store_true',help='require Chromium acceptance tests (no silent skip)')
    parser.add_argument('--playwright-module',help='optional absolute path to an existing Playwright installation')
    parser.add_argument('--browsers-path',help='optional Playwright browser cache path')
    parser.add_argument('--record',type=Path,help='copy honest test evidence here; existing files may be replaced')
    args=parser.parse_args()
    started=time.monotonic()
    hashes=source_hashes()
    tests=run([sys.executable,'-m','unittest','discover','-s','tests','-v'],cwd=ROOT)
    log=tests.stdout+tests.stderr
    print(log.decode(),end='')
    experiment=audit(ROOT/'results/experiment')
    with tempfile.TemporaryDirectory(prefix='motif-verify-') as temp:
        scratch=Path(temp)
        evidence=scratch/'evidence'; evidence.mkdir()
        reproduced=scratch/'experiment'
        run([sys.executable,ROOT/'scripts/experiment.py','--out',reproduced],cwd=ROOT)
        original=ROOT/'results/experiment'
        compared=0
        for p in original.rglob('*'):
            if p.is_file() and p.name not in ('metrics.json','completion.json'):
                assert p.read_bytes()==(reproduced/p.relative_to(original)).read_bytes(), 'experiment reproduction mismatch: '+str(p)
                compared+=1
        experiment['reproduced_deterministic_files']=compared
        experiment['reproduction_metrics']=json.loads((reproduced/'metrics.json').read_bytes())
        run([sys.executable,ROOT/'scripts/install.py','--prefix',scratch/'installed'],cwd=ROOT)
        executable=scratch/'installed/bin/motif-audit'
        isolated=scratch/'isolated'; isolated.mkdir()
        # Prevent networking in the installed application even if future code adds it.
        guard=scratch/'guard'; guard.mkdir()
        (guard/'sitecustomize.py').write_text('''import sys

def offline(event,args):
    if event.startswith("socket."):
        raise RuntimeError("network disabled by offline verification")
sys.addaudithook(offline)
''')
        env={**os.environ,'PYTHONPATH':str(guard),'PYTHONDONTWRITEBYTECODE':'1'}
        fasta=isolated/'sequences.fasta'; motifs=isolated/'motifs.json'
        fasta.write_bytes((ROOT/'examples/sequences.fasta').read_bytes())
        motifs.write_bytes((ROOT/'examples/motifs.json').read_bytes())
        before={p.name:digest(p.read_bytes()) for p in [fasta,motifs]}
        timing=scratch/'timing.json'
        base=[executable,'--fasta',fasta,'--motifs',motifs,'--seed','1729','--shuffles','199']
        run(['/usr/bin/time','-f','{"elapsed_seconds":%e,"peak_process_rss_kib":%M}','-o',timing,*base,'--out',scratch/'one'],cwd=isolated,env=env)
        run([*base,'--out',scratch/'two'],cwd=isolated,env=env)
        cli_metrics=json.loads(timing.read_bytes())
        one,two=scratch/'one',scratch/'two'
        outputs={}
        for p in sorted(one.iterdir()):
            assert p.read_bytes()==(two/p.name).read_bytes(), 'nondeterministic output: '+p.name
            outputs[p.name]={'sha256':digest(p.read_bytes()),'bytes':p.stat().st_size}
        assert before=={p.name:digest(p.read_bytes()) for p in [fasta,motifs]},'input changed'
        complete=json.loads((one/'completion.json').read_bytes())
        for name,meta in complete['artifacts'].items():
            assert meta==outputs[name]
        # Installed entry point must preserve nonzero exit codes and existing outputs.
        collision=subprocess.run([str(x) for x in [*base,'--out',one]],capture_output=True,cwd=isolated,env=env,timeout=10)
        assert collision.returncode==2
        assert all(digest((one/n).read_bytes())==m['sha256'] for n,m in outputs.items())
        invalid=subprocess.run([str(x) for x in [*base,'--out',scratch/'invalid','--shuffles','0']],capture_output=True,cwd=isolated,env=env,timeout=10)
        assert invalid.returncode==2 and not (scratch/'invalid').exists()
        installed={'complete':True,'offline_socket_audit_hook':True,'isolated_cwd':True,
                   'identical_repeated_outputs':True,'unchanged_input_hashes':before,
                   'artifact_metadata':outputs,'installed_zipapp_bytes':executable.stat().st_size,**cli_metrics}
        (evidence/'installed.json').write_bytes(canonical(installed))
        shutil.copytree(one,evidence/'sample')
        (evidence/'tests.log').write_bytes(log)
        (evidence/'experiment-audit.json').write_bytes(canonical(experiment))
        browser=None
        if args.browser:
            hostile='</script><img src="https://invalid.test/p" onerror="window.pwned=1"> @@JS@@ & <svg/onload=window.pwned=1>'
            raw=(f'>=hostile {hostile}\n'+ 'ACGTAC'*100 +'\n>reverse\n'+'GTACGT'*100+'\n>palindrome\n'+'ATAT'*200+'\n>empty\nNNRYSW\n').encode()
            manifest=canonical({'motifs':[{'id':'target','pattern':'ACGTAC'},{'id':'palindrome','pattern':'ATAT'},{'id':'any_concrete','pattern':'N'}]})
            fixture=scratch/'fixture'; empty=scratch/'empty'
            write_output(fixture,analyze(raw,manifest,shuffles=7),Budget())
            write_output(empty,analyze(b'>ambiguous\nNNRYSW',manifest,shuffles=7),Budget())
            browser_env=dict(os.environ)
            if args.playwright_module:
                browser_env['MOTIF_PLAYWRIGHT_MODULE']=str(Path(args.playwright_module).resolve())
            if args.browsers_path:
                browser_env['PLAYWRIGHT_BROWSERS_PATH']=str(Path(args.browsers_path).resolve())
            result=run(['node',ROOT/'scripts/browser.cjs',fixture,empty,one,evidence],cwd=ROOT,env=browser_env)
            print(result.stdout.decode(),end='')
            browser=json.loads((evidence/'browser.json').read_bytes())
        assert hashes==source_hashes(),'verification changed source/input files'
        overview={'complete':True,'elapsed_seconds':time.monotonic()-started,'python':sys.version,
                  'platform':platform.platform(),'runtime_dependencies':'Python standard library only',
                  'source_sha256':hashes,'unit_tests':'passed; see results/tests.log for repository evidence',
                  'installed_cli':installed,'experiment_audit':experiment,
                  'browser':browser if browser else {'status':'not requested; not verified by this invocation'}}
        (evidence/'verification.json').write_bytes(canonical(overview))
        if args.record:
            target=args.record.resolve()
            target.mkdir(parents=True,exist_ok=True)
            for p in evidence.iterdir():
                if p.is_dir():
                    shutil.copytree(p,target/p.name,dirs_exist_ok=True)
                else:
                    destination = ROOT/'results/tests.log' if p.name=='tests.log' and target==ROOT/'results/verification' else target/p.name
                    shutil.copy2(p,destination)
        print(json.dumps({'complete':True,'elapsed_seconds':overview['elapsed_seconds'],
                          'installed_cli':cli_metrics,'browser':bool(browser),'independent_control_counts':experiment['control_motif_checks']},sort_keys=True))


if __name__=='__main__':
    main()

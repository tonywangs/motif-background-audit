#!/usr/bin/env python3
"""Check the proposed tracked tree against repository artifact limits (read-only)."""
import json
from pathlib import Path
import subprocess

ROOT = Path(__file__).resolve().parents[1]


def check():
    names = subprocess.check_output(['git','ls-files','--cached','--others','--exclude-standard','-z'],cwd=ROOT).decode().split('\0')
    files = sorted({n for n in names if n and (ROOT/n).is_file()})
    assert len(files) <= 1000, 'file count exceeds 1,000'
    sizes = {n:(ROOT/n).stat().st_size for n in files}
    for name,size in sizes.items():
        assert size <= 10*1024*1024, f'file exceeds 10 MiB: {name}'
        assert not name.endswith('.log') or name=='results/tests.log', f'disallowed log: {name}'
        assert not (ROOT/name).is_symlink(), f'unexpected publication symlink: {name}'
    total = sum(sizes.values())
    assert total <= 32*1024*1024, 'tree exceeds 32 MiB'
    return {'complete':True,'files':len(files),'total_bytes':total,'largest_file_bytes':max(sizes.values(),default=0)}


if __name__=='__main__':
    print(json.dumps(check(),sort_keys=True))

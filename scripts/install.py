#!/usr/bin/env python3
"""Install a dependency-free executable zipapp without pip or network access."""
import argparse
from pathlib import Path
import shutil
import tempfile
import zipapp


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prefix", type=Path, required=True)
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    target = args.prefix / "bin" / "motif-audit"
    target.parent.mkdir(parents=True, exist_ok=True)
    if target.exists() or target.is_symlink():
        parser.error(f"refusing to overwrite {target}")
    with tempfile.TemporaryDirectory() as directory:
        stage = Path(directory)
        shutil.copytree(root / "motif_audit", stage / "motif_audit", ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        with tempfile.TemporaryDirectory() as output:
            archive = Path(output) / "motif-audit"
            zipapp.create_archive(stage, archive, interpreter="/usr/bin/env python3", main="motif_audit.cli:entrypoint", compressed=True)
            with target.open("xb") as stream:
                stream.write(archive.read_bytes())
    target.chmod(0o755)
    print(target)


if __name__ == "__main__":
    main()

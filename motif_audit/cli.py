"""Local-only command line with bounded work and fail-closed output handling."""
import argparse
from contextlib import contextmanager
from dataclasses import replace
import os
from pathlib import Path
import shutil
import signal
import stat
import sys
from . import __version__
from .core import AuditError, Budget, DEFAULT_LIMITS, analyze, compare, canonical, digest
from .report import match_tsv, render


def read_bounded(path, limit):
    # Reject special files (including named pipes) before potentially blocking reads.
    flags = os.O_RDONLY | getattr(os, "O_NONBLOCK", 0)
    fd = os.open(path, flags)
    with os.fdopen(fd, "rb") as stream:
        if not stat.S_ISREG(os.fstat(stream.fileno()).st_mode):
            raise AuditError("inputs must be regular local files")
        data = stream.read(limit + 1)
    if len(data) > limit:
        raise AuditError(f"input byte limit exceeded for {path}")
    return data


@contextmanager
def deadline(seconds):
    """POSIX hard deadline covers parsing, scanning, serialization and disk writes."""
    previous = {}
    old_timer = None

    def timeout(signum, frame):
        raise AuditError("runtime limit exceeded; analysis incomplete")

    def cancel(signum, frame):
        raise KeyboardInterrupt

    try:
        for sig in (signal.SIGINT, signal.SIGTERM):
            previous[sig] = signal.signal(sig, cancel)
        if hasattr(signal, "setitimer"):
            previous[signal.SIGALRM] = signal.signal(signal.SIGALRM, timeout)
            old_timer = signal.setitimer(signal.ITIMER_REAL, seconds)
        yield
    finally:
        if old_timer is not None:
            signal.setitimer(signal.ITIMER_REAL, *old_timer)
        for sig, handler in previous.items():
            signal.signal(sig, handler)


def write_output(output, result, budget):
    """Reserve an exclusive directory; a completion manifest is published last.

    On any Python exception or handled signal, remove only this run's new directory.
    SIGKILL/power loss may leave an incomplete directory, never a completion marker.
    """
    output = Path(output)
    output.mkdir()  # Deliberately refuses any existing directory, file or symlink.
    try:
        artifacts = {"analysis.json": canonical(result), "matches.tsv": match_tsv(result),
                     "report.html": render(result)}
        max_bytes = result["settings"]["limits"]["artifact_bytes"]
        for name, content in artifacts.items():
            budget.check()
            if len(content) > max_bytes:
                raise AuditError(f"artifact-size limit exceeded for {name}")
            with (output / (name + ".part")).open("xb") as stream:
                stream.write(content)
                stream.flush()
                os.fsync(stream.fileno())
        complete = canonical({"schema": "motif-background-audit/completion/1", "complete": True,
                              "artifacts": {n: {"bytes": len(v), "sha256": digest(v)}
                                            for n, v in artifacts.items()}})
        for name in artifacts:
            budget.check()
            (output / (name + ".part")).rename(output / name)
        with (output / "completion.json.part").open("xb") as stream:
            stream.write(complete)
            stream.flush()
            os.fsync(stream.fileno())
        budget.check()
        (output / "completion.json.part").rename(output / "completion.json")
    except BaseException:
        shutil.rmtree(output)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description="Audit exact IUPAC DNA motifs against within-record shuffled controls, entirely offline.")
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("--fasta", type=Path, required=True, help="local DNA FASTA, UTF-8")
    parser.add_argument("--motifs", type=Path, required=True, help="explicit JSON motif manifest")
    parser.add_argument("--out", type=Path, required=True, help="new output directory; parent must exist")
    parser.add_argument("--background", choices=["mononucleotide", "dinucleotide", "both"], default="mononucleotide")
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--shuffles", type=int, default=199)
    parser.add_argument("--seconds", type=float, default=DEFAULT_LIMITS.seconds,
                        help="wall-clock ceiling, >0 and at most 120 seconds")
    args = parser.parse_args(argv)
    if not 0 < args.seconds <= DEFAULT_LIMITS.seconds:
        parser.error("--seconds must be >0 and at most 120")
    limits = replace(DEFAULT_LIMITS, seconds=args.seconds)
    try:
        with deadline(limits.seconds):
            budget = Budget(limits.seconds)
            # Early collision check is advisory. Exclusive mkdir below handles races.
            if os.path.lexists(args.out):
                raise AuditError("output already exists; choose a new directory")
            fasta = read_bounded(args.fasta, limits.input_bytes)
            manifest = read_bounded(args.motifs, limits.manifest_bytes)
            function = compare if args.background == "both" else analyze
            options = {} if args.background == "both" else {"background": args.background}
            result = function(fasta, manifest, seed=args.seed, shuffles=args.shuffles,
                              limits=limits, budget=budget, **options)
            write_output(args.out, result, budget)
    except KeyboardInterrupt:
        print("cancelled; no completed output was published", file=sys.stderr)
        return 130
    except (AuditError, OSError, ValueError, TypeError, RecursionError) as e:
        print(f"motif-audit: {e}", file=sys.stderr)
        return 2
    print(f"Complete: {args.out / 'report.html'} ({len(result['matches'])} occurrences)")
    return 0


def entrypoint():
    raise SystemExit(main())

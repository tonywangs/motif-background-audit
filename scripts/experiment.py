#!/usr/bin/env python3
"""Execute the frozen synthetic protocol; refuse altered protocol/runner hashes."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import random
import resource
import signal
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from motif_audit.core import IUPAC, analyze, canonical, digest, holm, orientations


def generate(protocol, regime, repetition):
    seed = protocol["generator_seed_base"] + protocol["regimes"].index(regime)*100 + repetition
    rng = random.Random(seed)
    lines, plants = [], []
    for ri in range(protocol["records_per_dataset"]):
        n = protocol["length_per_record"]
        if regime == "locally_structured":
            sequence = [rng.choice("ACGT")]
            for _ in range(n-1):
                previous = sequence[-1]
                sequence.append(previous if rng.random() < .9 else rng.choice("".join(c for c in "ACGT" if c != previous)))
        else:
            weights = ([.7,.1,.1,.1] if ri == 0 else [.1,.7,.1,.1]) if regime == "composition_biased" else [1,1,1,1]
            sequence = rng.choices("ACGT", weights=weights, k=n)
        if regime == "planted":
            for j in range(12):
                start = 16 + 80*j
                sequence[start:start+6] = "ACGTAC"
                plants.append({"record": ri, "start": start, "end": start+6, "pattern": "ACGTAC"})
        lines.extend([f">{regime}_{repetition}_{ri} Synthetic protocol input; not biological data", "".join(sequence)])
    return ("\n".join(lines) + "\n").encode(), seed, plants


def uniform_baseline(result, motif):
    expectation = Fraction(0)
    for _, pattern in orientations(motif):
        probability = Fraction(1)
        for letter in pattern:
            probability *= Fraction(len(IUPAC[letter]), 4)
        windows = sum(sum(all(c in "ACGT" for c in r["sequence"][i:i+len(pattern)])
                          for i in range(max(0, r["length"]-len(pattern)+1))) for r in result["records"])
        expectation += probability * windows
    return float(expectation)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True, help="new directory only")
    args = parser.parse_args()
    protocol_bytes = (ROOT / "experiments/protocol.json").read_bytes()
    frozen = json.loads((ROOT / "experiments/frozen.json").read_bytes())
    assert digest(protocol_bytes) == frozen["protocol_sha256"], "protocol changed after freeze"
    assert digest(Path(__file__).read_bytes()) == frozen["runner_sha256"], "runner changed after freeze"
    protocol = json.loads(protocol_bytes)
    args.out.mkdir()
    (args.out / "inputs").mkdir()
    (args.out / "analyses").mkdir()
    manifest = canonical({"motifs": protocol["motifs"]})
    (args.out / "motifs.json").write_bytes(manifest)
    start = time.monotonic()
    # Independent global budget, including all datasets and serialization.
    def timeout(signum, frame):
        raise TimeoutError("frozen experiment budget exhausted; no completion manifest")
    signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, protocol["maximum_elapsed_seconds"])
    rows = []
    for regime in protocol["regimes"]:
        for repetition in range(protocol["repetitions"]):
            name = f"{regime}-{repetition:02}"
            fasta, seed, plants = generate(protocol, regime, repetition)
            (args.out / "inputs" / (name + ".fasta")).write_bytes(fasta)
            result = analyze(fasta, manifest, seed=seed+1000000, shuffles=protocol["shuffles"])
            (args.out / "analyses" / (name + ".json")).write_bytes(canonical(result))
            rows.append({"name": name, "regime": regime, "repetition": repetition,
                         "generator_seed": seed, "control_seed": seed+1000000,
                         "input_sha256": digest(fasta), "planted_sites": plants,
                         "motifs": [{"id": m["id"], "uniform_expected": uniform_baseline(result, m),
                                     **{k:v for k,v in stat.items() if k != "control_counts"}}
                                    for m, stat in zip(result["motifs"], result["statistics"])]})
            print(f"{name}: complete", flush=True)
    all_stats = [s for row in rows for s in row["motifs"]]
    for s,p in zip(all_stats, holm([s["p_upper"] for s in all_stats])):
        s["p_holm_study"] = p
    summary = {"protocol_sha256": digest(protocol_bytes), "runner_sha256": frozen["runner_sha256"],
               "datasets": rows, "rejection_counts": {
                   regime: {m["id"]: sum(s["p_holm"] <= .05 for row in rows if row["regime"] == regime
                                        for s in row["motifs"] if s["id"] == m["id"])
                            for m in protocol["motifs"]} for regime in protocol["regimes"]}}
    (args.out / "summary.json").write_bytes(canonical(summary))
    metrics = {"elapsed_seconds": time.monotonic()-start,
               "peak_process_rss_kib": resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
               "python": sys.version, "datasets": len(rows), "shuffles_per_dataset": protocol["shuffles"],
               "complete": True}
    (args.out / "metrics.json").write_bytes(canonical(metrics))
    hashes = {str(p.relative_to(args.out)): {"sha256": digest(p.read_bytes()), "bytes": p.stat().st_size}
              for p in sorted(args.out.rglob("*")) if p.is_file()}
    (args.out / "completion.json").write_bytes(canonical({"complete": True, "files": hashes}))
    signal.setitimer(signal.ITIMER_REAL, 0)
    print(json.dumps(metrics, sort_keys=True))


if __name__ == "__main__":
    main()

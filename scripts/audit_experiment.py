#!/usr/bin/env python3
"""Independent audit of all saved observed and shuffled counts, without core imports."""
from collections import Counter
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import random
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
SETS = {"A":"A", "C":"C", "G":"G", "T":"T", "R":"AG", "Y":"CT", "S":"CG", "W":"AT", "K":"GT", "M":"AC", "B":"CGT", "D":"AGT", "H":"ACT", "V":"ACG", "N":"ACGT"}
COMP = dict(zip("ATGCRYMKSWBDHVN", "TACGYRKMSWVHDBN"))


def sha(data):
    return hashlib.sha256(data).hexdigest()


def rev(s):
    return "".join(COMP[c] for c in reversed(s))


def count_and_matches(sequence, pattern):
    # Expand small frozen motifs into concrete strings, use literal searches on
    # both full reference orientations. No regex or production matcher code.
    words = ["".join(w) for w in itertools.product(*(SETS[c] for c in pattern))]
    palindrome = pattern == rev(pattern)
    found = set()
    for strand, text in [("+", sequence), ("-", rev(sequence))]:
        for word in words:
            i = text.find(word)
            while i >= 0:
                start = i if strand == "+" else len(sequence)-i-len(word)
                found.add((start,start+len(word),"both" if palindrome else strand))
                i = text.find(word,i+1)
    return found


def permutation(s, seed):
    rng = random.Random(seed)
    positions = [i for i,c in enumerate(s) if c in "ACGT"]
    result = list(s)
    for stop in range(len(positions),1,-1):
        pick = rng.getrandbits(stop.bit_length())
        while pick >= stop:
            pick = rng.getrandbits(stop.bit_length())
        a,b = positions[stop-1],positions[pick]
        result[a],result[b] = result[b],result[a]
    result = "".join(result)
    assert len(result) == len(s) and Counter(result) == Counter(s)
    assert all(result[i] == c for i,c in enumerate(s) if c not in "ACGT")
    return result


def adjusted(ps):
    ordered = sorted(set(ps))
    # Holm adjusted value at p: largest rank-weighted p' <= p.
    # Ties take their first rank; subsequent equal-p ranks cannot increase it.
    return [min(Fraction(1),max(q*sum(x>=q for x in ps) for q in ordered if q<=p)) for p in ps]


def audit(directory):
    started = time.monotonic()
    directory = Path(directory)
    freeze = json.loads((ROOT / "experiments/frozen.json").read_bytes())
    assert sha((ROOT / "experiments/protocol.json").read_bytes()) == freeze["protocol_sha256"]
    assert sha((ROOT / "scripts/experiment.py").read_bytes()) == freeze["runner_sha256"]
    protocol = json.loads((ROOT / "experiments/protocol.json").read_bytes())
    complete = json.loads((directory / "completion.json").read_bytes())
    assert complete["complete"] is True
    actual_files = {str(p.relative_to(directory)) for p in directory.rglob("*") if p.is_file()} - {"completion.json"}
    assert actual_files == set(complete["files"])
    for name, metadata in complete["files"].items():
        data = (directory / name).read_bytes()
        assert len(data) == metadata["bytes"] and sha(data) == metadata["sha256"], name
    report = json.loads((directory / "summary.json").read_bytes())
    assert report["protocol_sha256"] == freeze["protocol_sha256"]
    assert report["runner_sha256"] == freeze["runner_sha256"]
    assert len(report["datasets"]) == 32
    comparisons = observed_comparisons = 0
    all_ps, all_summaries = [], []
    for dataset in report["datasets"]:
        name = dataset["name"]
        result = json.loads((directory / "analyses" / (name + ".json")).read_bytes())
        raw = (directory / "inputs" / (name + ".fasta")).read_bytes()
        assert sha(raw) == result["input_sha256"]["fasta"] == dataset["input_sha256"]
        assert sha((directory / "motifs.json").read_bytes()) == result["input_sha256"]["manifest"]
        sequences = [part.split('\n',1)[1].replace('\n','') for part in raw.decode().split('>')[1:]]
        assert sequences == [r["sequence"] for r in result["records"]]
        assert len(sequences)==2 and all(len(s)==1024 for s in sequences)
        assert result["complete"] is True
        assert result["settings"]["shuffles"]==199
        assert int(result["settings"]["seed"]) == dataset["control_seed"] == dataset["generator_seed"]+1000000
        assert dataset["generator_seed"] == protocol["generator_seed_base"] + protocol["regimes"].index(dataset["regime"])*100 + dataset["repetition"]
        patterns = [m["pattern"] for m in result["motifs"]]
        assert patterns == [m["pattern"] for m in protocol["motifs"]]
        observed = []
        for mi,pattern in enumerate(patterns):
            expected = {(ri,*m) for ri,s in enumerate(sequences) for m in count_and_matches(s,pattern)}
            actual = {(m["record"],m["start"],m["end"],m["strand"]) for m in result["matches"] if m["motif"]==mi}
            assert actual == expected and len(actual) == sum(m["motif"]==mi for m in result["matches"])
            observed.append(len(expected)); observed_comparisons += 1
        controls = [[] for _ in patterns]
        for b,seeds in enumerate(result["control_seeds"]):
            shuffled = []
            for ri,(s,seed) in enumerate(zip(sequences,seeds)):
                derived = int(sha(f'motif-audit-v1:{dataset["control_seed"]}:{b}:{ri}'.encode())[:16],16)
                assert int(seed) == derived
                shuffled.append(permutation(s,int(seed)))
            for mi,p in enumerate(patterns):
                n = sum(len(count_and_matches(s,p)) for s in shuffled)
                controls[mi].append(n)
                assert result["statistics"][mi]["control_counts"][b] == n, (name,b,mi)
                comparisons += 1
        assert len(result["control_seeds"])==199
        ps=[]
        for mi,(n,counts,stat,short) in enumerate(zip(observed,controls,result["statistics"],dataset["motifs"])):
            assert len(counts)==199 and stat["observed"]==n
            ge = len([x for x in counts if x>=n])
            p,mean = Fraction(ge+1,200),Fraction(sum(counts),199)
            ps.append(p)
            assert stat["p_numerator"]==ge+1 and stat["p_denominator"]==200
            assert stat["exceedances"]==ge and stat["p_upper"]==float(p)
            assert stat["control_mean"]==float(mean)
            assert abs(stat["difference"]-float(n-mean)) < 1e-10
            if mean:
                assert abs(stat["fold"]-float(n/mean)) < 1e-10
            else:
                assert stat["fold"] is None
            for key,value in stat.items():
                if key != "control_counts":
                    assert short[key]==value
            degeneracy = 1
            for c in patterns[mi]:
                degeneracy *= len(SETS[c])
            expected_uniform = (1 if patterns[mi]==rev(patterns[mi]) else 2)*2*(1024-len(patterns[mi])+1)*degeneracy/(4**len(patterns[mi]))
            assert short["uniform_expected"]==expected_uniform
        for stat,p in zip(result["statistics"],adjusted(ps)):
            assert abs(stat["p_holm"]-float(p))<1e-12
        if dataset["regime"]=="planted":
            assert len(dataset["planted_sites"])==24
            for site in dataset["planted_sites"]:
                assert sequences[site["record"]][site["start"]:site["end"]]=="ACGTAC"
        else:
            assert dataset["planted_sites"]==[]
        all_ps.extend(ps); all_summaries.extend(dataset["motifs"])
    for stat,p in zip(all_summaries,adjusted(all_ps)):
        assert abs(stat["p_holm_study"]-float(p))<1e-12
    for regime,motifs in report["rejection_counts"].items():
        for motif,n in motifs.items():
            assert n==sum(s["p_holm"]<=.05 for row in report["datasets"] if row["regime"]==regime for s in row["motifs"] if s["id"]==motif)
    return {"complete":True,"datasets":32,"observed_motif_checks":observed_comparisons,
            "control_motif_checks":comparisons,"shuffle_invariant_checks":32*199*2,
            "elapsed_seconds":time.monotonic()-started,"independent_matcher":"concrete-word expansion and literal search on both reference orientations",
            "independent_statistics":"Fraction arithmetic and rank-based Holm; no production imports"}


if __name__ == "__main__":
    print(json.dumps(audit(sys.argv[1] if len(sys.argv)>1 else ROOT / 'results/experiment'),sort_keys=True,indent=2))

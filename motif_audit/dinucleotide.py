"""k=2 uShuffle/Wilson adaptation, v1. See docs/licenses/ushuffle.txt.

Copyright (c) 2007 Minghui Jiang, James Anderson, Joel Gillespie, Martin Mayne.
BSD-3-Clause; full notice distributed in motif_audit/USHUFFLE_LICENSE.txt.
"""
import hashlib
import random
import re
from .core import AuditError, Budget

VERSION = "wilson-dinucleotide-python/1"
MAX_STEPS = 2_000_000


def run_seed(seed, index):
    return int.from_bytes(hashlib.sha256(
        f"motif-audit-dinucleotide-v1:{seed}:{index}".encode("ascii")).digest()[:8], "big")


def _permute(values, rng):
    for i in range(len(values) - 1, 0, -1):
        j = rng.randrange(i + 1)
        values[i], values[j] = values[j], values[i]


def shuffle_run(sequence, seed, budget=None, *, max_steps=MAX_STEPS):
    budget = budget or Budget()
    budget.check()
    if len(sequence) <= 2:
        return sequence
    edges = {c: [] for c in dict.fromkeys(sequence)}
    for a, b in zip(sequence, sequence[1:]):
        edges[a].append(b)
    root = sequence[-1]
    tree = {root}
    last = {}
    rng = random.Random(seed)
    steps = 0
    for start in edges:
        current = start
        while current not in tree:
            steps += 1
            if steps > max_steps:
                raise AuditError("dinucleotide walk step limit exceeded; analysis incomplete")
            if steps % 1024 == 0:
                budget.check()
            last[current] = rng.randrange(len(edges[current]))
            current = edges[current][last[current]]
        current = start
        while current not in tree:
            tree.add(current)
            current = edges[current][last[current]]
    for base, outgoing in edges.items():
        if base != root:
            reserved = outgoing.pop(last[base])
            _permute(outgoing, rng)
            outgoing.append(reserved)
        else:
            _permute(outgoing, rng)
    positions = dict.fromkeys(edges, 0)
    result = [sequence[0]]
    for i in range(len(sequence) - 1):
        if i % 1024 == 0:
            budget.check()
        base = result[-1]
        result.append(edges[base][positions[base]])
        positions[base] += 1
    return "".join(result)


def shuffle_dinucleotide(sequence, seed, budget=None):
    budget = budget or Budget()
    index = 0
    def replace(match):
        nonlocal index
        value = shuffle_run(match[0], run_seed(seed, index), budget)
        index += 1
        return value
    budget.check()
    return re.sub("[ACGT]+", replace, sequence)

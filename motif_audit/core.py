"""Pure matching and statistics; no network, filesystem writes, or HTML."""
from collections import Counter
from dataclasses import asdict, dataclass
import hashlib
import json
import random
import re
import time

IUPAC = dict(zip("ACGTRYSWKMBDHVN", (
    "A", "C", "G", "T", "AG", "CT", "CG", "AT", "GT", "AC",
    "CGT", "AGT", "ACT", "ACG", "ACGT")))
COMPLEMENT = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")


class AuditError(ValueError):
    """An input, output, or resource contract could not be satisfied."""


@dataclass(frozen=True)
class Limits:
    input_bytes: int = 1_048_576
    manifest_bytes: int = 65_536
    records: int = 1_000
    sequence_length: int = 100_000
    total_bases: int = 200_000
    motifs: int = 32
    motif_length: int = 64
    matches: int = 30_000
    shuffles: int = 999
    work: int = 120_000_000
    seconds: float = 120.0
    artifact_bytes: int = 8_388_608


DEFAULT_LIMITS = Limits()


class Budget:
    def __init__(self, seconds=DEFAULT_LIMITS.seconds):
        self.deadline = time.monotonic() + seconds

    def check(self):
        if time.monotonic() >= self.deadline:
            raise AuditError("runtime limit exceeded; analysis incomplete")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def canonical(data):
    return (json.dumps(data, sort_keys=True, indent=2, ensure_ascii=True,
                       allow_nan=False) + "\n").encode("ascii")


def rc(pattern):
    return pattern.translate(COMPLEMENT)[::-1]


def _text(data, limit, kind):
    if len(data) > limit:
        raise AuditError(f"{kind} byte limit exceeded")
    try:
        text = data.decode("utf-8")
    except UnicodeDecodeError as e:
        raise AuditError(f"{kind} must be UTF-8") from e
    if any((ord(c) < 32 and c not in "\t\n\r") or ord(c) == 127 for c in text):
        raise AuditError(f"{kind} contains a forbidden control character")
    return text


def parse_fasta(data, limits=DEFAULT_LIMITS):
    text = _text(data, limits.input_bytes, "FASTA")
    records, seen, chunks = [], set(), []
    header = None
    length = total = 0

    def finish():
        if header is not None:
            if not chunks:
                raise AuditError(f"empty sequence for {header[0]!r}")
            sequence = "".join(chunks)
            records.append({"id": header[0], "description": header[1],
                            "sequence": sequence})

    for number, raw in enumerate(text.splitlines(), 1):
        if not raw.strip():
            continue
        if raw.startswith(">"):
            finish()
            if len(records) >= limits.records:
                raise AuditError("record count limit exceeded")
            parts = raw[1:].strip().split(maxsplit=1)
            if not parts or len(parts[0]) > 128 or len(raw) > 1024:
                raise AuditError(f"invalid or oversized FASTA header on line {number}")
            identifier = parts[0]
            if identifier in seen:
                raise AuditError(f"duplicate FASTA identifier {identifier!r}")
            seen.add(identifier)
            header = (identifier, parts[1] if len(parts) > 1 else "")
            chunks, length = [], 0
        else:
            if header is None:
                raise AuditError(f"sequence before FASTA header on line {number}")
            # Spaces and tabs on sequence lines are formatting, never coordinates.
            raw_sequence = raw.replace(" ", "").replace("\t", "")
            seq = raw_sequence.upper()
            if not raw_sequence.isascii() or not seq or any(c not in IUPAC for c in seq):
                raise AuditError(f"non-DNA IUPAC symbol on line {number}")
            length += len(seq)
            total += len(seq)
            if length > limits.sequence_length or total > limits.total_bases:
                raise AuditError("sequence length or total bases limit exceeded")
            chunks.append(seq)
    finish()
    if not records:
        raise AuditError("FASTA has no records")
    return records


def parse_manifest(data, limits=DEFAULT_LIMITS):
    text = _text(data, limits.manifest_bytes, "manifest")

    def unique(pairs):
        result = {}
        for key, value in pairs:
            if key in result:
                raise AuditError(f"duplicate JSON key {key!r}")
            result[key] = value
        return result

    try:
        value = json.loads(text, object_pairs_hook=unique)
    except (ValueError, RecursionError) as e:
        raise AuditError(f"invalid motif JSON: {e}") from e
    if not isinstance(value, dict) or set(value) != {"motifs"}:
        raise AuditError('manifest must contain only a "motifs" array')
    motifs = value["motifs"]
    if not isinstance(motifs, list) or not 1 <= len(motifs) <= limits.motifs:
        raise AuditError("motif count outside limits")
    seen, result = set(), []
    for m in motifs:
        if not isinstance(m, dict) or set(m) != {"id", "pattern"}:
            raise AuditError("each motif must have exactly id and pattern")
        identifier, pattern = m["id"], m["pattern"]
        if (not isinstance(identifier, str) or not re.fullmatch(r"[A-Za-z0-9_.:-]{1,64}", identifier)
                or identifier in seen):
            raise AuditError("motif IDs must be unique safe identifiers (1–64 characters)")
        if not isinstance(pattern, str) or not 1 <= len(pattern) <= limits.motif_length:
            raise AuditError("motif length outside limits")
        if not pattern.isascii():
            raise AuditError("motifs accept ASCII DNA IUPAC letters only")
        pattern = pattern.upper()
        if any(c not in IUPAC for c in pattern):
            raise AuditError("motifs accept DNA IUPAC letters only")
        seen.add(identifier)
        result.append({"id": identifier, "pattern": pattern,
                       "reverse_complement": rc(pattern), "palindrome": rc(pattern) == pattern})
    return result


def orientations(motif):
    pattern = motif["pattern"]
    yield ("both" if motif["palindrome"] else "+"), pattern
    if not motif["palindrome"]:
        yield "-", motif["reverse_complement"]


def scan(records, motifs, *, keep=True, limits=DEFAULT_LIMITS, budget=None):
    """Lookahead regex includes overlaps. Ambiguous sequence letters never match."""
    counts = [0] * len(motifs)
    matches = []
    total = 0
    for ri, record in enumerate(records):
        seq = record["sequence"]
        for mi, motif in enumerate(motifs):
            if budget:
                budget.check()
            for strand, pattern in orientations(motif):
                regex = re.compile("(?=" + "".join("[" + IUPAC[c] + "]" for c in pattern) + ")")
                for found in regex.finditer(seq):
                    total += 1
                    if total > limits.matches:
                        raise AuditError("match count limit exceeded (observed or control); no truncation")
                    counts[mi] += 1
                    if total % 1024 == 0 and budget:
                        budget.check()
                    if keep:
                        start = found.start()
                        matches.append({"record": ri, "motif": mi, "start": start,
                                        "end": start + len(pattern), "strand": strand})
    if keep:
        matches.sort(key=lambda m: (m["record"], m["start"], m["motif"], m["strand"]))
    return counts, matches


def control_seed(seed, repetition, record):
    # Independent per-record/per-repetition stream; indices are zero based.
    key = f"motif-audit-v1:{seed}:{repetition}:{record}".encode("ascii")
    return int.from_bytes(hashlib.sha256(key).digest()[:8], "big")


def shuffle_sequence(sequence, seed):
    concrete = [c for c in sequence if c in "ACGT"]
    rng = random.Random(seed)
    # Explicit Fisher–Yates; getrandbits rejection avoids modulo bias.
    for i in range(len(concrete) - 1, 0, -1):
        j = rng.getrandbits((i + 1).bit_length())
        while j > i:
            j = rng.getrandbits((i + 1).bit_length())
        concrete[i], concrete[j] = concrete[j], concrete[i]
    bases = iter(concrete)
    return "".join(next(bases) if c in "ACGT" else c for c in sequence)


def summary(observed, controls):
    mean = sum(controls) / len(controls)
    ge = sum(c >= observed for c in controls)
    return {"observed": observed, "control_counts": controls, "control_mean": mean,
            "exceedances": ge, "p_numerator": ge + 1, "p_denominator": len(controls) + 1,
            "p_upper": (ge + 1) / (len(controls) + 1),
            "difference": observed - mean,
            "fold": observed / mean if mean else None,
            "fold_note": None if mean else "undefined: control mean is zero"}


def holm(pvalues):
    """Family-wise error adjustment, valid under arbitrary motif dependence."""
    result = [0.0] * len(pvalues)
    running = 0.0
    for rank, index in enumerate(sorted(range(len(pvalues)), key=lambda i: (pvalues[i], i))):
        running = max(running, min(1.0, (len(pvalues) - rank) * pvalues[index]))
        result[index] = running
    return result


def estimated_work(records, motifs, shuffles):
    scan_work = sum(max(0, len(r["sequence"]) - len(m["pattern"]) + 1)
                    * len(m["pattern"]) * (1 if m["palindrome"] else 2)
                    for r in records for m in motifs)
    return (shuffles + 1) * scan_work + shuffles * sum(len(r["sequence"]) for r in records)


def analyze(fasta, manifest, *, seed=1729, shuffles=199, limits=DEFAULT_LIMITS, budget=None):
    if type(seed) is not int or not 0 <= seed < 2**64:
        raise AuditError("seed must be an unsigned 64-bit integer")
    if type(shuffles) is not int or not 1 <= shuffles <= limits.shuffles:
        raise AuditError("shuffle count outside limits")
    budget = budget or Budget(limits.seconds)
    records = parse_fasta(fasta, limits)
    motifs = parse_manifest(manifest, limits)
    work = estimated_work(records, motifs, shuffles)
    if work > limits.work:
        raise AuditError(f"shuffle-work limit exceeded: {work} > {limits.work}")
    observed, matches = scan(records, motifs, limits=limits, budget=budget)
    controls = [[] for _ in motifs]
    seeds = []
    for b in range(shuffles):
        budget.check()
        row = [control_seed(seed, b, ri) for ri in range(len(records))]
        seeds.append([str(s) for s in row])  # JS cannot exactly represent uint64 numbers.
        shuffled = [{**r, "sequence": shuffle_sequence(r["sequence"], row[ri])}
                    for ri, r in enumerate(records)]
        counts, _ = scan(shuffled, motifs, keep=False, limits=limits, budget=budget)
        for bucket, count in zip(controls, counts):
            bucket.append(count)
    stats = [summary(n, c) for n, c in zip(observed, controls)]
    for stat, adjusted in zip(stats, holm([s["p_upper"] for s in stats])):
        stat["p_holm"] = adjusted
    budget.check()
    return {"schema": "motif-background-audit/1", "complete": True,
            "software_version": "1.0.0", "input_sha256": {"fasta": digest(fasta), "manifest": digest(manifest)},
            "settings": {"seed": str(seed), "shuffles": shuffles, "estimated_work": work,
                         "limits": asdict(limits), "coordinates": "0-based half-open, forward reference",
                         "null": "within-record concrete-base permutation; ambiguous positions fixed",
                         "rng": "SHA256-derived uint64 streams; Python MT19937 getrandbits Fisher-Yates v1",
                         "multiplicity": "Holm across all motifs in this analysis"},
            "records": [{**r, "length": len(r["sequence"]), "counts": dict(sorted(Counter(r["sequence"]).items()))}
                        for r in records],
            "motifs": motifs, "matches": matches, "statistics": stats, "control_seeds": seeds}

import collections
from dataclasses import replace
from fractions import Fraction
import itertools
import json
import random
import unittest
from motif_audit.core import (AuditError, Budget, DEFAULT_LIMITS, analyze, canonical,
                              control_seed, holm, parse_fasta, parse_manifest, rc,
                              scan, shuffle_sequence, summary)

# Deliberately independent definitions and algorithm: inspect every window on
# the entire forward AND reverse reference, remapping reverse coordinates.
SETS = {"A": {"A"}, "T": {"T"}, "G": {"G"}, "C": {"C"}, "R": {"A", "G"},
        "Y": {"C", "T"}, "M": {"A", "C"}, "K": {"G", "T"}, "S": {"G", "C"},
        "W": {"A", "T"}, "B": {"C", "G", "T"}, "D": {"A", "G", "T"},
        "H": {"A", "C", "T"}, "V": {"A", "C", "G"}, "N": {"A", "C", "G", "T"}}
COMP = dict(A="T", T="A", C="G", G="C", R="Y", Y="R", M="K", K="M", S="S", W="W", B="V", V="B", D="H", H="D", N="N")


def reverse(s):
    return "".join(COMP[c] for c in reversed(s))


def oracle(sequence, pattern):
    found = set()
    palindrome = pattern == reverse(pattern)
    for strand, text in [("+", sequence), ("-", reverse(sequence))]:
        for i in range(len(text) - len(pattern) + 1):
            if all(base in SETS[letter] for base, letter in zip(text[i:i+len(pattern)], pattern)):
                start = i if strand == "+" else len(sequence) - i - len(pattern)
                found.add((start, start + len(pattern), "both" if palindrome else strand))
    return found


def manifest(*patterns):
    return json.dumps({"motifs": [{"id": f"m{i}", "pattern": p} for i, p in enumerate(patterns)]}).encode()


class MatchingTests(unittest.TestCase):
    def test_independent_oracle_500_seeded_cases(self):
        rng = random.Random(20261012)
        for case in range(500):
            sequence = "".join(rng.choices("ACGTNRYSWKMBDHV", weights=[25]*4+[1]*11, k=rng.randrange(1, 150)))
            pattern = "".join(rng.choices(list(SETS), k=rng.randrange(1, 13)))
            if case % 7 == 0:
                pattern = rng.choice(["ATAT", "NN", "CGCG", "W", "AAA", "RY", "ACGTAC"])
            with self.subTest(case=case):
                motifs = parse_manifest(manifest(pattern))
                counts, matches = scan([{"sequence": sequence}], motifs)
                actual = {(m["start"], m["end"], m["strand"]) for m in matches}
                expected = oracle(sequence, pattern)
                self.assertEqual(actual, expected)
                self.assertEqual(counts, [len(expected)])

    def test_boundaries_overlaps_palindromes_reverse_and_empty(self):
        for sequence, pattern, expected in [
            ("AAAA", "AAA", {(0, 3, "+"), (1, 4, "+")}),
            ("TTT", "AAA", {(0, 3, "-")}),
            ("ATATAT", "ATAT", {(0, 4, "both"), (2, 6, "both")}),
            ("ACGT", "NN", {(0, 2, "both"), (1, 3, "both"), (2, 4, "both")}),
            ("NRY", "N", set()), ("A", "AA", set()), ("TTTT", "AAAA", {(0, 4, "-")}),
            ("AG", "RN", {(0, 2, "+")}),
            ("AT", "AN", {(0, 2, "+"), (0, 2, "-")})]:
            with self.subTest(sequence=sequence, pattern=pattern):
                _, matches = scan([{"sequence": sequence}], parse_manifest(manifest(pattern)))
                self.assertEqual({(m["start"], m["end"], m["strand"]) for m in matches}, expected)

    def test_control_overflow_is_an_error_not_a_truncated_null(self):
        # The observed AATAA has zero AAA hits, but some permutations have two.
        with self.assertRaisesRegex(AuditError, "match count limit"):
            analyze(b">x\nAATAA", manifest("AAA"), shuffles=20,
                    limits=replace(DEFAULT_LIMITS, matches=1))

    def test_all_iupac_complements(self):
        for letter in SETS:
            self.assertEqual(rc(letter), COMP[letter])
            self.assertEqual(rc(rc(letter)), letter)
            self.assertEqual({COMP[x] for x in SETS[letter]}, SETS[rc(letter)])


class ParsingTests(unittest.TestCase):
    def test_normalization(self):
        r = parse_fasta(b"\r\n>x description\r\nac gt\tNR\r\nyS\n\n>y\nT\n")
        self.assertEqual(r[0], {"id": "x", "description": "description", "sequence": "ACGTNRYS"})
        self.assertEqual(parse_manifest(manifest("aNr"))[0]["pattern"], "ANR")

    def test_invalid_fasta(self):
        for text in [b"", b"ACGT", b">x", b">x\n>x\nA", b">\nACG", b">x\nU", b">x\nA-", b">x\nA?", b">x\nA\n>x desc\nA", b">x\nA\x00", b">x\n\xff", b">x\n>y\nA", b">x\nA\vT"]:
            with self.subTest(text=text), self.assertRaises(AuditError):
                parse_fasta(text)

    def test_unicode_case_expansion_is_not_DNA(self):
        with self.assertRaises(AuditError):
            parse_fasta(">x\nß".encode())
        with self.assertRaises(AuditError):
            parse_manifest(manifest("ß"))

    def test_invalid_manifest(self):
        for text in [b"", b"[]", b"null", b'{"motifs":[]}', b'{"motifs":[],"motifs":[]}',
                     b'{"motifs":[{"id":"x","id":"y","pattern":"A"}]}',
                     b'{"motifs":[{"id":"x","pattern":"A","extra":1}]}',
                     b'{"motifs":[{"id":"x","pattern":"A"},{"id":"x","pattern":"T"}]}',
                     b'{"motifs":[{"id":"=x","pattern":"A"}]}',
                     b'{"motifs":[{"id":"x","pattern":null}]}',
                     b'{"motifs":[{"id":"x","pattern":".*"}]}', b'{"motifs":[true]}']:
            with self.subTest(text=text), self.assertRaises(AuditError):
                parse_manifest(text)

    def test_limits(self):
        cases = [("input_bytes", 2, b">x\nA", manifest("A")),
                 ("manifest_bytes", 2, b">x\nA", manifest("A")),
                 ("sequence_length", 2, b">x\nAAA", manifest("A")),
                 ("total_bases", 2, b">x\nAA\n>y\nA", manifest("A")),
                 ("records", 1, b">x\nA\n>y\nA", manifest("A")),
                 ("motifs", 1, b">x\nA", manifest("A", "T")),
                 ("motif_length", 1, b">x\nA", manifest("AT")),
                 ("matches", 2, b">x\nAAA", manifest("A")),
                 ("work", 1, b">x\nAAA", manifest("A")),
                 ("shuffles", 1, b">x\nA", manifest("A"))]
        for field, value, fasta, motifs in cases:
            with self.subTest(limit=field), self.assertRaises(AuditError):
                analyze(fasta, motifs, shuffles=2, limits=replace(DEFAULT_LIMITS, **{field:value}))
        with self.assertRaises(AuditError):
            analyze(b">x\nA", manifest("A"), budget=Budget(-1))
        for seed in [-1, 2**64, 1.5, True]:
            with self.assertRaises(AuditError):
                analyze(b">x\nA", manifest("A"), seed=seed)
        for count in [0, -1, 1000, 1.5, True]:
            with self.assertRaises(AuditError):
                analyze(b">x\nA", manifest("A"), shuffles=count)


class ShuffleStatisticsTests(unittest.TestCase):
    def test_independent_shuffle_invariants_500(self):
        rng = random.Random(1717)
        for case in range(500):
            s = "".join(rng.choices(list(SETS), k=rng.randrange(0, 1000)))
            shuffled = shuffle_sequence(s, case)
            self.assertEqual(len(s), len(shuffled))
            self.assertEqual(collections.Counter(s), collections.Counter(shuffled))
            self.assertEqual([(i,c) for i,c in enumerate(s) if c not in "ACGT"],
                             [(i,c) for i,c in enumerate(shuffled) if c not in "ACGT"])
            self.assertEqual(shuffled, shuffle_sequence(s, case))
        self.assertGreater(len({shuffle_sequence("AACCGGTT", i) for i in range(30)}), 20)
        self.assertEqual(shuffle_sequence("NNRY", 3), "NNRY")

    def test_statistical_arithmetic_independent(self):
        rng = random.Random(22)
        for _ in range(250):
            observed = rng.randrange(20)
            controls = [rng.randrange(20) for _ in range(rng.randrange(1, 100))]
            stat = summary(observed, controls)
            tail = Fraction(len([x for x in controls if x >= observed]) + 1, len(controls) + 1)
            mean = Fraction(sum(controls), len(controls))
            self.assertEqual(stat["p_upper"], float(tail))
            self.assertEqual(stat["p_numerator"], len([x for x in controls if x >= observed]) + 1)
            self.assertEqual(stat["p_denominator"], len(controls) + 1)
            self.assertAlmostEqual(stat["difference"], float(observed - mean))
            self.assertAlmostEqual(stat["fold"], float(observed / mean))
        self.assertIsNone(summary(1, [0, 0])["fold"])
        self.assertEqual(summary(0, [0, 0])["p_upper"], 1)
        self.assertEqual(summary(10, [0, 0])["p_upper"], 1/3)
        for ps in [[.01,.04,.03], [1,1,1], [.01,.01,.01], [.001,.2,.9,1]]:
            sorted_ps = sorted(ps)
            expected = [float(min(Fraction(1), max(Fraction(str(sorted_ps[j])) * (len(ps)-j)
                         for j in range(sorted_ps.index(p)+1)))) for p in ps]
            for a,b in zip(holm(ps), expected):
                self.assertAlmostEqual(a,b)

    def test_end_to_end_control_counts_against_oracle(self):
        fasta = b">first\nACGNTTATATACGGACGTAC\n>second\nAAAAANAAAATTTTT"
        result = analyze(fasta, manifest("ATAT", "AN", "AAA", "ACGTAC"), seed=17, shuffles=31)
        for b, seeds in enumerate(result["control_seeds"]):
            shuffled = [shuffle_sequence(r["sequence"], int(seed)) for r, seed in zip(result["records"], seeds)]
            for mi, m in enumerate(result["motifs"]):
                n = sum(len(oracle(s, m["pattern"])) for s in shuffled)
                self.assertEqual(result["statistics"][mi]["control_counts"][b], n)
        self.assertEqual(canonical(result), canonical(analyze(fasta, manifest("ATAT", "AN", "AAA", "ACGTAC"), seed=17, shuffles=31)))
        self.assertEqual(len(set(itertools.chain.from_iterable(result["control_seeds"]))), 62)
        self.assertEqual(control_seed(17, 0, 0), int(result["control_seeds"][0][0]))


if __name__ == "__main__":
    unittest.main()

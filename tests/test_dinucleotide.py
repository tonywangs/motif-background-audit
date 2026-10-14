"""Independent invariants, finite outcome enumeration, and failure regressions."""
from collections import Counter
from dataclasses import replace
import itertools
import random
import re
import unittest
from unittest.mock import patch
from motif_audit.core import (AuditError, Budget, DEFAULT_LIMITS, analyze, compare,
                              canonical, scan, parse_manifest)
from motif_audit.dinucleotide import shuffle_dinucleotide, shuffle_run
from test_core import oracle, manifest
import test_io


def signature(s):
    return s[:1], s[-1:], tuple(sorted(Counter(zip(s, s[1:])).items()))


def outcomes(s):
    # Independent multiset permutation enumeration; no graph traversal.
    return sorted({''.join(p) for p in itertools.permutations(s)
                   if signature(''.join(p)) == signature(s)})


class DinucleotideTests(unittest.TestCase):
    def test_invariants_matching_replay_500(self):
        rng = random.Random(48151)
        for case in range(500):
            s = ''.join(rng.choices('ACGTNRYS', weights=[12]*4+[1]*4, k=rng.randrange(1,500)))
            shuffled = shuffle_dinucleotide(s,case)
            self.assertEqual(shuffled,shuffle_dinucleotide(s,case))
            self.assertEqual(Counter(s),Counter(shuffled))
            self.assertEqual(signature(s),signature(shuffled))
            self.assertEqual([(i,c) for i,c in enumerate(s) if c not in 'ACGT'],
                             [(i,c) for i,c in enumerate(shuffled) if c not in 'ACGT'])
            for m in re.finditer('[ACGT]+',s):
                t=shuffled[m.start():m.end()]
                self.assertEqual(Counter(m[0]),Counter(t))
                self.assertEqual(signature(m[0]),signature(t))
            pattern = rng.choice(['AN','AAA','ATAT','NN','ACGTAC','RY'])
            counts,matches=scan([{'sequence':shuffled}],parse_manifest(manifest(pattern)))
            expected=oracle(shuffled,pattern)
            self.assertEqual(counts,[len(expected)])
            self.assertEqual({(m['start'],m['end'],m['strand']) for m in matches},expected)

    def test_all_binary_strings_through_length_eight(self):
        for n in range(1,9):
            groups={}
            for p in itertools.product('AC',repeat=n):
                s=''.join(p); groups.setdefault(signature(s),set()).add(s)
            for space in groups.values():
                for s in space:
                    for seed in range(16):
                        self.assertIn(shuffle_dinucleotide(s,seed),space)
        for s in ['', 'A', 'AC', 'NNRY', 'AAAAAA', 'ACACAC']:
            self.assertEqual(shuffle_dinucleotide(s,12),s)

    def test_background_validation_limits_and_replay(self):
        f=b'>x\nAACGTACNACACAAA'; m=manifest('AAA','ACG')
        with self.assertRaises(AuditError): analyze(f,m,background='typo')
        with self.assertRaises(AuditError): shuffle_run('ACACA',1,Budget(-1))
        with self.assertRaises(AuditError): shuffle_run('ACACA',1,max_steps=0)
        with self.assertRaises(AuditError): compare(f,m,limits=replace(DEFAULT_LIMITS,work=1))
        with self.assertRaises(AuditError): compare(f,m,budget=Budget(-1))
        with patch('motif_audit.dinucleotide.shuffle_run',side_effect=KeyboardInterrupt):
            with self.assertRaises(KeyboardInterrupt): compare(f,m,shuffles=2)
        r=compare(f,m,shuffles=7)
        self.assertEqual(canonical(r),canonical(compare(f,m,shuffles=7)))
        self.assertEqual(r['backgrounds']['mononucleotide']['settings'],analyze(f,m,shuffles=7)['settings'])
        self.assertEqual(r['backgrounds']['mononucleotide']['statistics'][0]['control_counts'],
                         analyze(f,m,shuffles=7)['statistics'][0]['control_counts'])


# Run the same malformed-input, collisions, serialization, signal, cleanup and
# artifact limit tests on an actual comparison result, not only the old format.
class ComparisonOutputTests(test_io.OutputTests):
    def setUp(self):
        super().setUp()
        self.result=compare(b'>=hostile </script><img src=https://invalid.test onerror=alert(1)> @@JS@@\nACGTACATAT',
                            manifest('ACGTAC','ATAT'),shuffles=2)

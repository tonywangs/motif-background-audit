# Frozen paired-background experiment

This is a modest synthetic diagnostic, not a biological validation or a claim
of universally better calibration. No paid inference, Evo2, GPUs, personal data,
or biological experiment was used. Background choice was fixed before evaluation.

The [protocol](../experiments/dinucleotide-protocol.json),
[freeze](../experiments/dinucleotide-frozen.json) and 16 frozen FASTA inputs define
four datasets per condition, three 240-base records each, four motifs, and 199
controls for each background. Each paired run uses identical inputs, motifs,
master seeds and per-record seeds; dinucleotide mode derives independent run
streams. Conditions are independent uniform bases; a first-order chain with
probability .8 of copying the previous base (otherwise a uniform draw); random
six-base blocks ACGTAC/GTACGT with longer-range structure; and uniform sequences
with five ACGTAC insertions per record. Plant coordinates are retained in the
freeze. Analysis includes reverse complements and overlapping occurrences.

## Results

Counts below are rejections among 16 motif/dataset tests per condition, after
Holm correction across eight motif/background tests within each dataset.

| Condition | Mononucleotide | Dinucleotide |
|---|---:|---:|
| Independent bases | 0 | 0 |
| First-order dependent | 4 | 0 |
| Longer-range structured | 8 | 0 |
| Planted motif | 8 | 8 |

The first-order example `first_order-00` contains 290 AAA/TTT occurrences.
The mononucleotide control mean is 32.29 (rounded), upper-tail p = .005 and
paired Holm p = .04. With dinucleotide controls the mean is 287.20 and the
upper-tail p = .18 (paired Holm p = 1). This demonstrates how the conditioning
assumption can change an apparent enrichment conclusion on a known synthetic
input. It does not establish a biological interpretation.

Both backgrounds reject the planted six-base motif and its embedded ACG motif
in every planted dataset under the paired correction. The 16 planted six-base
occurrences in `planted-00` include the 15 designed insertions plus a chance
occurrence; occurrence counts are not plant counts.

No comparison passes study-wide Holm correction across all 128 tests. With
199 controls the minimum empirical p is .005; Holm's first threshold is
.05/128, so study-wide rejection is impossible at this resolution. This is a
limitation of the preregistered budget/design, not evidence of equivalence.
The study was not retuned after seeing results.

Dinucleotide controls also yield no upper-tail rejection on the longer-range
condition. That does not validate a first-order model for six-base blocks:
those blocks impose constraints that the shuffle does not preserve. In
`longer_range-00`, ACGTAC/GTACGT observed count is 224 while the dinucleotide
control mean is 236.82. This analysis only tests enrichment, not depletion or
all forms of model misspecification. Conditioning can also remove signal
that is represented in the preserved counts. Four datasets per condition are
far too few to estimate general calibration rates.

## Reproducibility and verification

`results/dinucleotide/summary.json` retains all 128 rows, effect sizes, probability
fractions and three adjustment families. Each dataset directory contains every
control count and seed in `analysis.json`, a self-contained comparison report,
TSV matches and an artifact hash manifest. The input and implementation hashes
are checked before evaluation; no protocol, old result, old example or old
experiment input was replaced.

Initial evaluation took **17.33 elapsed seconds, 8.85 CPU seconds**, with
**23,292 KiB peak process RSS**, on Python 3.12.3/Linux. The preregistered CPU
ceiling was 1,800 seconds. Timing is process/environment dependent; repeat timing
and deterministic artifact hashes are recorded separately by the verifier.
See `results/dinucleotide/metrics.json` for original measurements and
`results/dinucleotide-verification/` for final verification evidence, runtime
versions, artifact sizes, repeated hashes and installed-CLI measurements.

The independent audit rescans all observed and shuffled inputs with an explicit
window matcher and reverse-reference coordinate mapping, checks invariants, and
recomputes means, effects, probability fractions and all Holm adjustments using
rational arithmetic. It reuses the production shufflers for seeded replay;
separate exhaustive enumeration checks reachable outcomes independently.

All binary sequences of lengths one through eight are checked against enumerated
endpoint/dinucleotide equivalence classes, using 16 seeds each. Another 500 seeded
mixed-IUPAC cases check per-run endpoints, composition, adjacent-pair counts,
ambiguity positions, replay and independent matching. Five fixed finite spaces
are each sampled 20,000 times. Their maximum absolute frequency deviations are
48, 125, 37.33, 81 and 0; all satisfy the preregistered 323.90 threshold. The final
space is a singleton. These checks do not prove uniformity; they can miss small
or rare biases and do not validate larger spaces exhaustively.

There is no biological, cross-platform RNG compatibility, screen-reader or
non-Chromium validation. The algorithm conditions on fixed endpoints, and its
bounded random walk may fail on difficult inputs. Failure aborts the analysis;
no failed replicate is silently replaced or removed. CPU and wall-clock guards
are implementation limits, not biological sampling assumptions.

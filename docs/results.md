# Synthetic experiment results

This experiment tests the software and diagnoses the consequences of its background assumptions. It does not validate biological activity, binding, regulation, or function. No Evo2 inference, paid inference, or laboratory experiment was performed.

## Protocol and evidence

The [protocol](../experiments/protocol.json) and [runner](../scripts/experiment.py) were hashed locally at **2026-10-12 05:09:31 UTC**, before the first evaluation. The [freeze record](../experiments/frozen.json) is a local integrity checkpoint, not an external preregistration or trusted timestamp.

Protocol SHA-256: `6fd2a6efd1cb9e4d468fd3ad88bacc14a2caf5b100f6343b0cc66910f2096ad6`.

There were four prespecified generators, eight repetitions each, two records of 1,024 bases per dataset, four prespecified motifs, and 199 shuffled controls per dataset. Every dataset used its specified seed; none were excluded or rerun based on outcomes. Locally structured sequences used a 0.9 probability of repeating the previous base. Planted datasets overwrote `ACGTAC` at 12 fixed positions in each record. The biased generator favored A in one record and C in the other. Full generation details and hypotheses are in the protocol.

The first evaluation completed all **32 datasets** in **14.436 seconds**, with peak process RSS **20,992 KiB** on Python 3.12.3/Linux ([actual metrics](../results/experiment/metrics.json)). The experiment-wide 1,800-second ceiling was not approached. These are small local workloads, not a genome-scale performance benchmark.

All [input FASTA files](../results/experiment/inputs), [full analyses](../results/experiment/analyses), [summary data](../results/experiment/summary.json), and [artifact hashes](../results/experiment/completion.json) are retained. Analyses contain all 199 counts for every motif, per-record control seeds, normalized sequences, observed coordinates, exact tail numerators/denominators, and effect sizes. The independent auditor uses concrete-word expansion and literal searches, reconstructs every shuffle without importing production code, checks its invariants, and recomputes statistics using rational arithmetic.

## Observed findings

Each cell below is the number of repetitions with **per-dataset Holm-adjusted upper-tail p ≤ 0.05**, out of eight. Every dataset's statistical family contains all four motifs.

| Generator | AAAAAA | ACGTAC | CGCGCG (palindrome) | TATAWA |
|---|---:|---:|---:|---:|
| Independent uniform bases | **1/8** | 0/8 | 0/8 | 0/8 |
| Independent composition-biased bases | 0/8 | 0/8 | 0/8 | 0/8 |
| Locally persistent bases, no planting | **8/8** | 0/8 | 0/8 | 0/8 |
| Independent bases with ACGTAC planted | 0/8 | **8/8** | 0/8 | 0/8 |

Means below average the eight datasets within a generator. “Control mean” averages each dataset's 199-control mean. They are descriptive summaries, not estimates with a validated uncertainty interval.

| Generator and motif | Mean observed count | Mean control count | Independent uniform expected count |
|---|---:|---:|---:|
| Independent / AAAAAA | 1.625 | 0.989 | 0.995 |
| Composition-biased / AAAAAA | 126.500 | 119.971 | 0.995 |
| Locally persistent / AAAAAA | 613.125 | 1.940 | 0.995 |
| Planted / ACGTAC | 28.500 | 1.003 | 0.995 |

**Composition matters.** Under the biased generator, the homopolymer count is much larger than the uniform expectation, but the within-record shuffled mean is similarly large. No motif rejects after per-dataset Holm correction in those eight datasets. This supports using a composition-aware comparison for this deliberately independent, biased generator; it does not establish broad calibration or prove absence of enrichment.

**Local structure still matters after composition is preserved.** All eight persistent datasets have an AAAAAA upper-tail probability of 1/200 = 0.005 (Holm 0.02), although no sites were planted. Their generator creates long runs; controls destroy those runs. The test correctly detects disagreement with its null. Interpreting the result as motif-specific biological function would be misleading.

**The planted example is detected, with a counting nuance.** ACGTAC also has p = 0.005 (Holm 0.02) in all eight planted datasets. There are 24 deliberate overwrites per dataset but the average observed strand-occurrence count is 28.5. Incidental occurrences and overlapping/reverse-strand windows also count. Observed count is therefore not the number of deliberately inserted sites.

**The null result is not perfect zero rejection.** One of the eight independent uniform datasets rejects for AAAAAA after per-dataset Holm adjustment (raw p = 0.01). Chance rejections remain possible under a correct null, especially across multiple dataset families. Eight repetitions are far too few for a reliable type-I error calibration study. This outcome is retained without seed replacement.

**Study-wide adjustment changes the conclusion.** Holm across all **128** prespecified comparisons rejects **none** at 0.05. With 199 controls the smallest possible raw p is 0.005, while the first study-wide Holm threshold is 0.05/128. Thus the study-wide non-rejection was forced by Monte Carlo resolution, even for strong planted effects. The protocol declared both the per-dataset descriptive summaries and the study-wide adjustment in advance. Per-dataset detections must not be presented as study-wide discoveries.

## Verification and limits

The [one-command verifier](../scripts/verify.py) includes:

- 500 seeded comparisons against a separate brute-force, full-reverse-reference matcher, plus explicit overlaps, ambiguous bases, exact boundaries, short records, dual-orientation degenerate patterns, and palindromes.
- 500 independent shuffle-invariant cases and 250 exact statistical arithmetic cases.
- A separate experimental audit of **128 observed motif counts**, **25,472 control motif counts**, and **12,736 record-level permutations**, with independent Holm arithmetic.
- A fresh execution of the frozen protocol, comparing every deterministic input/result artifact with the preserved experiment; timings are allowed to differ.
- An installed, dependency-free executable run from an isolated directory with socket access blocked, repeated for byte-identical artifacts and unchanged input hashes.
- Chromium 141.0.7390.37 / Playwright 1.56.1 / Node v24.20.0: 11 acceptance checks, 2,599 reachable fixture occurrences, at most 75 mounted match rows, desktop and 375px widths, keyboard/context navigation, filters, downloads, empty results, and hostile imported text. External requests, JavaScript errors, and CSP violations were all zero in the recorded run.

See [verification metrics and source hashes](../results/verification/verification.json), [plain test output](../results/tests.log), [browser evidence](../results/verification/browser.json), and [sample report](../results/verification/sample/report.html). CLI runtime, peak process RSS, exact artifact sizes, and repeated-run hashes are recorded in [installed workflow evidence](../results/verification/installed.json); they are measurements of the small bundled example, not throughput guarantees.

Remaining limits: synthetic data only, few generator repetitions, no dinucleotide-preserving alternative null, no probabilistic motif scores, no Monte Carlo confidence intervals, no biological confirmation, and no validated screen-reader or non-Chromium behavior. Reports retain full input sequence content. The hard resource ceilings deliberately refuse large analyses instead of sampling or truncating them.

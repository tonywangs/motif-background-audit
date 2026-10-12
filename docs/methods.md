# Methods and prior implementations

## Existing work reviewed

Primary sources were inspected before implementation on 2026-10-12. This project makes no novelty claim.

- [NCBI nucleotide FASTA conventions](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) describe a definition line beginning with `>`, a unique identifier, wrapped sequence data, and IUPAC nucleotide symbols. The local parser follows that structure while documenting its more permissive identifier policy and stricter rejection of gaps; it is not a GenBank submission validator.
- [MEME Suite FIMO](https://meme-suite.org/meme/doc/fimo.html) scores motif occurrences with log-odds and a background model, derives position-level probabilities, and supports multiple-testing adjustments. A position weight matrix can represent differing preferences at each base; an exact IUPAC set cannot represent those preferences. This project's probability refers to an aggregate count under shuffled sequences, **not** the FIMO probability of an individual scored occurrence. Their thresholds, coordinates, and outputs should not be equated.
- [Biopython `nt_search` documentation](https://biopython.org/docs/latest/api/Bio.SeqUtils.html#Bio.SeqUtils.nt_search) and the [pinned Biopython 1.86 implementation](https://github.com/biopython/biopython/blob/biopython-186/Bio/SeqUtils/__init__.py#L254) provide an existing exact IUPAC matcher. The implementation builds character classes and advances search positions to retain overlapping forward-strand hits. This workbench uses regular-expression lookahead, adds explicit strand/coordinate rules and control analyses, and validates against an independently written window-by-window oracle. Biopython is not a runtime dependency and its implementation is not copied into the repository.

## Conditional shuffled background

Let `T_m(S)` be the pooled strand-occurrence count of motif `m` across input records. For each control repetition, each record's concrete A/C/G/T bases are uniformly permuted over that record's concrete positions. Every ambiguous letter stays at its original position. Length and concrete-base counts are invariant, separately for every record. No bases move between records.

The null assumes concrete positions are exchangeable **conditional on per-record composition and the fixed ambiguity mask**. Under that null, observed order and independent shuffled orders are comparable. Overlapping matches are allowed and their dependence is reproduced by scanning each full shuffled sequence; no independence-of-windows assumption is needed for the permutation comparison.

The assumption is often biologically implausible: local GC changes, repeats, dinucleotide preferences, motif clustering, and longer-range sequence dependence violate it. Preserving total composition does not address those structures. A small tail probability then diagnoses disagreement with this null, not biological regulatory activity. Record boundaries are part of the model; splitting or joining records can change the result.

## RNG and reproducibility

Default master seed: 1729. For zero-based repetition `b` and record index `r`, compute SHA-256 of ASCII `motif-audit-v1:SEED:b:r`; interpret the first eight digest bytes as an unsigned big-endian integer. Save every resulting integer as a decimal string. Initialize Python's MT19937 with that integer. The implementation explicitly performs descending Fisher–Yates swaps, drawing a candidate with `getrandbits((i+1).bit_length())` and rejecting candidates greater than `i`. Equal bases are interchangeable. Identity permutations and repeated controls are valid draws, not discarded.

The precise seed derivation and shuffle algorithm are part of schema v1. Reordering input records changes their stream assignments. Use the recorded Python version for strict long-term reproduction; compatibility with future changes to Python's random implementation is not promised. Original input hashes detect changes to source bytes even if normalization produces the same sequence.

## Statistics

For `B` controls `C_1 ... C_B` and observed count `O`:

- Upper-tail Monte Carlo probability: `(1 + sum(C_b >= O)) / (B+1)`; ties count in the tail. It is never zero and its minimum is `1/(B+1)`.
- Control mean: `sum(C_b)/B`.
- Difference: `O - mean(C)`.
- Fold: `O/mean(C)`, or JSON `null` with an explanation if the mean is zero. No pseudocount or infinity is silently substituted.

Holm's step-down adjustment covers all `M` motifs in a single analysis: sort their probabilities ascending and take the running maximum of `(M-rank+1) * p`, capped at one. It does not require independence among motifs. Duplicate-pattern IDs remain in this family. Selecting inputs or patterns after seeing results would require accounting for that selection; this CLI cannot infer a user's broader search history.

A low count's upper-tail probability may be one; this tool does not perform a lower-tail depletion test. Monte Carlo resolution can be coarse and low shuffle counts are for exploration. No precision interval, power guarantee, or biological discovery threshold is implied.

The synthetic experiment additionally includes an **independent uniform-base expected-count baseline**, calculated from IUPAC set sizes and eligible windows. It is an expected count, not a baseline p-value or a fitted model. Composition bias alone can make that expectation inappropriate. Per-dataset Holm families and a separate study-wide adjustment were fixed in advance; all results, including non-rejections, are retained.

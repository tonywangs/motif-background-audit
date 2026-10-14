# Dinucleotide controls: specification before implementation

The selected implementation is a Python, k=2-only adaptation of the Wilson
last-exit-tree algorithm in MEME's uShuffle source, pinned to
[commit 8ca7292109d0943d4a3d1d8bca188fe6234b1652](https://github.com/cinquin/MEME/blob/8ca7292109d0943d4a3d1d8bca188fe6234b1652/src/ushuffle.c).
The retained upstream notice is BSD-3-Clause: [license](licenses/ushuffle.txt).
The reviewed wrapper is [fasta-shuffle-letters.c](https://github.com/cinquin/MEME/blob/8ca7292109d0943d4a3d1d8bca188fe6234b1652/src/fasta-shuffle-letters.c).
Also reviewed the [uShuffle Cython wrapper](https://github.com/guma44/ushuffle/tree/b54df34f60246ff236ec58bbccb40af3cea608ed), whose checkout lacks src/ushuffle.c;
we do not depend on that wrapper. This is an application of existing work,
not a new shuffling algorithm.

[Jiang et al. (2008)](https://pmc.ncbi.nlm.nih.gov/articles/PMC2375906/)
describe sampling Eulerian walks using random arborescences. Our specialization
uses four vertices, labeled parallel edges, a tree rooted at the last base,
and randomizes outgoing edges except each nonroot's reserved last exit.
Python MT19937 with getrandbits rejection replaces upstream modulo selection.
Uniformity is the algorithmic target over distinct endpoint-conditioned strings;
finite seeded tests are empirical checks, not a proof of this implementation.
No external library, compiler or network is required at runtime.

## Contract

Mononucleotide remains the unchanged default: concrete bases can cross ambiguity
gaps within a record. Dinucleotide mode splits each record into maximal A/C/G/T
runs, holds every ambiguous symbol at its original position, and independently
shuffles runs using SHA256-derived run streams. Both endpoints of every run are
fixed; each run preserves exact mono- and dinucleotide counts. This also preserves
all adjacent pairs across ambiguity boundaries. Empty runs, lengths one and two,
homopolymers and singleton shuffle spaces legitimately return the input. Identity
outputs are retained; duplicates are never rejected. Balanced graphs also retain
both endpoints (no random rotation). Different runs need not share composition.
This is a more strongly conditioned null than whole-record mononucleotide.
The [MEME utility](https://meme-suite.org/meme/doc/fasta-shuffle-letters.html)
supports general k and alphabet handling; our ambiguity/run contract is explicit
and is not a claim of command-line equivalence with MEME.

A walk exceeding the fixed step ceiling fails the entire analysis; it never
falls back to a biased shuffle or discards a replicate. Runtime checks also abort.
Completed runs must report all control counts, seeds and implementation version.
Do not select a background by significance. A comparison corrects across both
backgrounds and all motifs, alongside within-background corrections.
Dinucleotide conditioning does not model higher-order repeats, regional
heterogeneity, biology or function. Validation here is synthetic only.

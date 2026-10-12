# Input, matching, and output contracts

## FASTA

Local regular files only; no stdin, URL, compressed input, named pipe, or device input. Read at most the configured byte ceiling plus one byte before rejection. UTF-8 is required; ASCII control characters other than line breaks and tabs are rejected. A byte-order mark is not supported.

A header starts with `>` in column zero. Its first whitespace-separated token is a unique, case-sensitive record ID (1–128 characters). Duplicate IDs fail, even if descriptions differ. Descriptions are retained as data. Headers may be up to 1,024 characters. Unlike NCBI's GenBank submission rules, local IDs are not restricted to its 25-character alphabet; this tool does not prepare GenBank submissions.

Sequence data before a header, empty records, gaps, punctuation, RNA `U`, and non-ASCII sequence letters fail. DNA IUPAC letters `ACGTRYSWKMBDHVN` are accepted in either case. Spaces/tabs within sequence lines and blank lines are ignored; remaining letters are uppercased. Wrapping and CRLF/LF line endings do not affect coordinates. These normalization decisions do not alter input files; hashes describe original bytes.

## Motif manifest

UTF-8 JSON with exactly one top-level key, `motifs`, containing 1–32 objects with exactly `id` and `pattern`. Duplicate JSON keys and duplicate motif IDs fail. IDs match `[A-Za-z0-9_.:-]{1,64}`. Patterns contain 1–64 ASCII DNA IUPAC letters and normalize to uppercase. Different IDs may use identical patterns; these remain distinct hypotheses and are included separately in Holm adjustment. No automatic motif discovery, regular-expression syntax, PWM, or MEME-format input is accepted.

## Matching

For each possible start, a window matches when every sequence base is a concrete A/C/G/T member of the corresponding motif letter's set. A sequence ambiguity is never evidence for a match, including when the motif letter is `N`. Overlapping occurrences are retained. A motif longer than a record produces zero hits.

Both the pattern and its reverse complement are searched on the forward reference. `[start,end)` always slices the forward normalized sequence; `end-start` equals motif length. Reverse-strand oriented sequence is the reverse complement of that slice. Palindromic patterns are searched once and labeled `both`. Nonpalindromic degenerate patterns may match both orientations at one start and contribute two strand occurrences. This is intentionally an occurrence count, not a unique-locus count.

Records and motifs retain input order. Matches sort by record index, start, motif index, and strand. JSON uses zero-based array indices to refer to records/motifs. TSV uses their IDs. The browser highlights the forward-reference slice and separately shows the motif-oriented sequence, with at most 24 flanking bases on either side.

## Resource ceilings

| Resource | Hard default ceiling |
|---|---:|
| FASTA bytes | 1,048,576 |
| Manifest bytes | 65,536 |
| Records | 1,000 |
| Bases per record / across records | 100,000 / 200,000 |
| Motifs / bases per motif | 32 / 64 |
| Observed or any individual control scan occurrences | 30,000 |
| Shuffles | 999 (minimum 1; default 199) |
| Conservative work estimate | 120,000,000 |
| CLI wall-clock seconds | 120 |
| Each JSON, TSV, or HTML artifact | 8,388,608 bytes |
| Mounted browser match rows | 75 |

For motif length `k` and record length `n`, the scan-work estimate is `max(0,n-k+1) * k * orientations`. Sum across records/motifs, multiply by `B+1`, then add `B * total_bases` for shuffling. This is a deterministic admission rule, not a measured CPU instruction count. Runtime and output-size checks remain necessary after admission. Controls use the same match ceiling as observations and fail the entire analysis if exceeded.

`--seconds` may lower the 120-second ceiling. Other CLI ceilings cannot be raised. Pure Python functions accept a `Limits` object for tests and embedding; callers that relax those defaults assume responsibility for their own budgets. Library calls cooperatively check time; the CLI additionally uses a POSIX timer covering reading, parsing, analysis, serialization, and writes. Linux/Python 3.12 is validated. Pure library use and non-POSIX platforms do not have the CLI's hard timer guarantee.

## Output and failure handling

Successful output contains exactly:

- `analysis.json`: normalized inputs, original input SHA-256 hashes, settings/limits, all observed matches, every control count, seeds, and statistics.
- `matches.tsv`: full observed occurrence export with forward and motif-oriented matched bases.
- `report.html`: offline self-contained interactive report; full data embedded, no remote resources.
- `completion.json`: written last, with SHA-256 and byte lengths for the three artifacts.

A valid completed output has its completion manifest **and matching artifact hashes**. The `complete` field in a standalone intermediate JSON file is not a substitute for this directory-level commit marker. The output directory is exclusively reserved after analysis. Files stage under `.part` names and the completion marker is renamed last. Existing directories, files, and symlinks are refused. Signals and Python exceptions clean up only the newly reserved output. No overwrite, partial-match success, or silent truncation occurs. Disk writes are flushed and fsynced, but this is not a database or a power-loss-durable directory transaction. SIGKILL or power failure may leave incomplete files without a marker; inspect and remove that directory before retrying. The filesystem parent is assumed to be trusted; hostile concurrent replacement of the output directory is outside the threat model.

CLI success exits 0, input/resource/output errors exit 2, and handled SIGINT/SIGTERM cancellation exits 130. The installer preserves those exit codes. A 64-bit unsigned master seed is required. Seed values are strings in JSON to avoid JavaScript integer precision loss. Repeated runs with the same input bytes, software, settings, and supported Python generator produce byte-identical artifacts. Wall-clock measurements live separately in verification evidence, not deterministic analysis output.

## Browser and export safety

Imported strings enter the page as escaped JSON or escaped static text, and dynamic UI uses `textContent`/text nodes. No imported HTML, JavaScript, CSS, URL, or event handler is executed. A hash-based Content Security Policy disallows external resources, connections, objects, forms, and base URLs. No CDN, analytics, fonts, or service worker is used. Downloads are local Blob URLs.

Sequence/motif/strand selections compose. Filters and context navigation preserve access to all results through page entry, first/last controls, and previous/next controls. Alt+Left/Right changes pages; native Tab/Enter works on buttons. Numeric statistics have captions and column/row headers, with horizontal overflow contained in focusable table regions. Browser statistics always describe the full input family, regardless of view filters.

TSV tabs/newlines in text fields are backslash-escaped; a leading `=`, `+`, `-`, or `@` receives an apostrophe to reduce spreadsheet formula interpretation. This also affects exported `+`/`-` strand cells; machine consumers wanting exact typed values should use JSON. The browser's full TSV export is checked byte-for-byte against the CLI export. The browser's JSON export preserves data but may use different key order/escaping from the deterministic CLI file.

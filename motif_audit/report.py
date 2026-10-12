"""Deterministic, self-contained reports. Imported strings are always inert text."""
import base64
import hashlib
import html
from importlib.resources import files
from .core import AuditError, canonical, rc

TSV_COLUMNS = ("sequence_id", "motif_id", "start", "end", "strand", "forward_match", "oriented_match")


def tsv_cell(value):
    # Spreadsheet-safe text; normal IDs unchanged. JSON retains original IDs.
    value = str(value).replace("\t", "\\t").replace("\r", "\\r").replace("\n", "\\n")
    return "'" + value if value.startswith(("=", "+", "-", "@")) else value


def match_tsv(result):
    rows = ["\t".join(TSV_COLUMNS)]
    for m in result["matches"]:
        record, motif = result["records"][m["record"]], result["motifs"][m["motif"]]
        fragment = record["sequence"][m["start"]:m["end"]]
        values = (record["id"], motif["id"], m["start"], m["end"], m["strand"],
                  fragment, rc(fragment) if m["strand"] == "-" else fragment)
        rows.append("\t".join(tsv_cell(v) for v in values))
    return ("\n".join(rows) + "\n").encode("utf-8")


def _csp_hash(text):
    return "'sha256-" + base64.b64encode(hashlib.sha256(text.encode()).digest()).decode() + "'"


def render(result):
    assets = files("motif_audit").joinpath("assets")
    template = assets.joinpath("report.html").read_text(encoding="utf-8")
    css = assets.joinpath("report.css").read_text(encoding="utf-8")
    js = assets.joinpath("report.js").read_text(encoding="utf-8")
    data = canonical(result).decode().replace("&", "\\u0026").replace("<", "\\u003c").replace(">", "\\u003e")
    rows = []
    for m, s in zip(result["motifs"], result["statistics"]):
        vals = [m["id"], m["pattern"], s["observed"], f'{s["control_mean"]:.3f}',
                f'{s["difference"]:+.3f}', "undefined" if s["fold"] is None else f'{s["fold"]:.3f}',
                f'{s["p_numerator"]}/{s["p_denominator"]}', f'{s["p_upper"]:.4f}', f'{s["p_holm"]:.4f}']
        rows.append("<tr>" + "".join(f'<{("th scope=\"row\"" if i == 0 else "td")}>{html.escape(str(v))}</{("th" if i == 0 else "td")}>' for i, v in enumerate(vals)) + "</tr>")
    replacements = {"CSS": css, "JS": js, "DATA": data,
                    "CSP": f"default-src 'none'; script-src {_csp_hash(js)}; style-src {_csp_hash(css)}; connect-src 'none'; img-src 'none'; base-uri 'none'; form-action 'none'; object-src 'none'",
                    "SUMMARY": "".join(rows), "N_RECORDS": str(len(result["records"])),
                    "N_BASES": f'{sum(r["length"] for r in result["records"]):,}',
                    "N_MATCHES": f'{len(result["matches"]):,}', "N_SHUFFLES": str(result["settings"]["shuffles"]),
                    "SEED": result["settings"]["seed"], "FASTA_HASH": result["input_sha256"]["fasta"],
                    "MANIFEST_HASH": result["input_sha256"]["manifest"]}
    # Single-pass substitution: imported text cannot inject another template token.
    import re
    document = re.sub(r"@@([A-Z_]+)@@", lambda m: replacements[m[1]], template).encode("utf-8")
    if len(document) > result["settings"]["limits"]["artifact_bytes"]:
        raise AuditError("report-size limit exceeded")
    return document

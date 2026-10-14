"use strict";
(() => {
  const data = JSON.parse(document.getElementById("analysis-data").textContent);
  const byId = id => document.getElementById(id);
  const pageSize = 75;
  let filtered = [], page = 0, selected = -1;
  for (const [key, items] of [["sequence", data.records], ["motif", data.motifs]]) {
    items.forEach((item, i) => {
      const option = document.createElement("option");
      option.value = String(i); option.textContent = item.id;
      byId(key + "-filter").append(option);
    });
  }
  function cell(row, value) {
    const td = document.createElement("td"); td.textContent = String(value); row.append(td); return td;
  }
  function clearContext() {
    selected = -1;
    byId("context-info").textContent = "Choose a match to inspect 24 flanking bases on each side.";
    for (const id of ["context-bases", "oriented", "description"]) byId(id).textContent = "";
    byId("context-previous").disabled = byId("context-next").disabled = true;
  }
  function render() {
    const pages = Math.max(1, Math.ceil(filtered.length / pageSize));
    page = Math.max(0, Math.min(pages - 1, page));
    const start = page * pageSize;
    byId("match-rows").replaceChildren();
    for (let i = start; i < Math.min(start + pageSize, filtered.length); i++) {
      const m = filtered[i], row = document.createElement("tr");
      cell(row, data.records[m.record].id); cell(row, data.motifs[m.motif].id);
      cell(row, m.start); cell(row, m.end); cell(row, m.strand === "both" ? "both · palindrome" : m.strand);
      const button = document.createElement("button"); button.type = "button"; button.textContent = "Context";
      button.setAttribute("aria-label", `Context: ${data.records[m.record].id}, ${data.motifs[m.motif].id}, ${m.start}, ${m.strand}`);
      button.addEventListener("click", () => context(i)); cell(row, "").append(button);
      byId("match-rows").append(row);
    }
    if (!filtered.length) {
      const row = document.createElement("tr"), td = cell(row, "No matches for these filters.");
      td.colSpan = 6; td.className = "empty"; byId("match-rows").append(row);
    }
    byId("match-status").textContent = `${filtered.length.toLocaleString()} of ${data.matches.length.toLocaleString()} occurrences · showing ${filtered.length ? start + 1 : 0}–${Math.min(start + pageSize, filtered.length)}. Statistics above remain unchanged.`;
    byId("page-number").value = page + 1; byId("page-number").max = pages;
    byId("pages").textContent = `of ${pages}`;
    byId("first").disabled = byId("previous").disabled = page === 0;
    byId("last").disabled = byId("next").disabled = page === pages - 1;
  }
  function filter() {
    const r = byId("sequence-filter").value, m = byId("motif-filter").value, s = byId("strand-filter").value;
    filtered = data.matches.filter(x => (r === "" || x.record === Number(r)) && (m === "" || x.motif === Number(m)) && (s === "" || x.strand === s));
    page = 0; clearContext(); render();
  }
  function reverseComplement(s) {return s.split("").reverse().map(c => ({A:"T",C:"G",G:"C",T:"A"})[c]).join("");}
  function context(i) {
    if (i < 0 || i >= filtered.length) return;
    selected = i;
    const m = filtered[i], r = data.records[m.record], motif = data.motifs[m.motif];
    const left = Math.max(0, m.start - 24), right = Math.min(r.length, m.end + 24);
    const fragment = r.sequence.slice(m.start, m.end);
    byId("context-info").textContent = `${r.id} · ${motif.id} (${motif.pattern}) · [${m.start}, ${m.end}) · strand ${m.strand} · context [${left}, ${right}) · filtered match ${i+1}/${filtered.length}`;
    const mark = document.createElement("mark"); mark.textContent = fragment;
    byId("context-bases").replaceChildren(document.createTextNode(r.sequence.slice(left, m.start)), mark, document.createTextNode(r.sequence.slice(m.end, right)));
    byId("oriented").textContent = `Forward-reference context above. Motif-oriented match: ${m.strand === "-" ? reverseComplement(fragment) : fragment}`;
    byId("description").textContent = r.description;
    byId("context-previous").disabled = i === 0; byId("context-next").disabled = i === filtered.length - 1;
    const targetPage = Math.floor(i / pageSize);
    if (page !== targetPage) {page = targetPage; render();}
  }
  for (const id of ["sequence-filter", "motif-filter", "strand-filter"]) byId(id).addEventListener("change", filter);
  byId("reset").addEventListener("click", () => {for (const id of ["sequence-filter", "motif-filter", "strand-filter"]) byId(id).value = ""; filter();});
  for (const [id, action] of [["first", () => 0], ["previous", () => page-1], ["next", () => page+1], ["last", () => Math.ceil(filtered.length/pageSize)-1]]) byId(id).addEventListener("click", () => {page = action(); render();});
  byId("page-number").addEventListener("change", () => {const n = Number(byId("page-number").value); page = Number.isFinite(n) ? Math.floor(n)-1 : 0; render();});
  document.addEventListener("keydown", e => {if (e.altKey && ["ArrowLeft", "ArrowRight"].includes(e.key)) {e.preventDefault(); page += e.key === "ArrowLeft" ? -1 : 1; render();}});
  byId("context-previous").addEventListener("click", () => context(selected - 1));
  byId("context-next").addEventListener("click", () => context(selected + 1));
  function download(name, content, mime) {
    const url = URL.createObjectURL(new Blob([content], {type:mime}));
    const link = document.createElement("a"); link.href = url; link.download = name;
    document.body.append(link); link.click(); link.remove(); setTimeout(() => URL.revokeObjectURL(url), 1000);
    byId("export-status").textContent = `Exported ${name}`;
  }
  function safeCell(x) {let s = String(x).replaceAll("\t", "\\t").replaceAll("\r", "\\r").replaceAll("\n", "\\n"); return /^[=+\-@]/.test(s) ? "'"+s : s;}
  byId("export-tsv").addEventListener("click", () => {
    const rows = [["sequence_id", "motif_id", "start", "end", "strand", "forward_match", "oriented_match"].join("\t")];
    for (const m of filtered) {
      const r = data.records[m.record], fragment = r.sequence.slice(m.start, m.end);
      rows.push([r.id, data.motifs[m.motif].id, m.start, m.end, m.strand, fragment, m.strand === "-" ? reverseComplement(fragment) : fragment].map(safeCell).join("\t"));
    }
    download("filtered-matches.tsv", rows.join("\n")+"\n", "text/tab-separated-values;charset=utf-8");
  });
  byId("export-json").addEventListener("click", () => download("analysis.json", JSON.stringify(data, null, 2)+"\n", "application/json"));
  const backgrounds = data.backgrounds || {[data.settings.background || "mononucleotide"]: data};
  data.motifs.forEach((m, i) => {
    const option = document.createElement("option"); option.value = String(i); option.textContent = m.id;
    byId("distribution-motif").append(option);
  });
  function distribution() {
    const stats = backgrounds[byId("background-select").value].statistics;
    const s = stats[Number(byId("distribution-motif").value)];
    const lo = Math.min(...s.control_counts), hi = Math.max(...s.control_counts);
    const width = Math.max(1, Math.ceil((hi-lo+1)/40));
    const bins = Array(Math.floor((hi-lo)/width)+1).fill(0);
    s.control_counts.forEach(n => bins[Math.floor((n-lo)/width)]++);
    byId("distribution-rows").replaceChildren();
    bins.forEach((n, i) => {
      const row = document.createElement("tr"), a = lo+i*width, b = Math.min(hi,a+width-1);
      cell(row, a===b ? a : `${a}–${b}`); cell(row,n); byId("distribution-rows").append(row);
    });
    byId("distribution-description").textContent = `Observed: ${s.observed}; ${s.control_counts.length} controls; range ${lo}–${hi}. Identity and repeated controls are retained.`;
  }
  function background() {
    const name = byId("background-select").value, analysis = backgrounds[name];
    byId("background-description").textContent = `${name}: ${analysis.settings.null}. Background choice is an assumption, not a significance filter.`;
    const body = byId("summary").querySelector("tbody"); body.replaceChildren();
    analysis.statistics.forEach((s,i) => {
      const row = document.createElement("tr"), motif = data.motifs[i];
      const th = document.createElement("th"); th.scope = "row"; th.textContent = motif.id; row.append(th);
      [motif.pattern,s.observed,s.control_mean.toFixed(3),s.difference.toFixed(3),s.fold===null ? "undefined" : s.fold.toFixed(3),
       `${s.p_numerator}/${s.p_denominator}`,s.p_upper.toFixed(4),s.p_holm.toFixed(4),
       s.p_holm_comparison===undefined ? "—" : s.p_holm_comparison.toFixed(4)].forEach(v => cell(row,v));
      body.append(row);
    });
    distribution();
  }
  byId("background-select").addEventListener("change", background);
  byId("distribution-motif").addEventListener("change", distribution);
  background();
  filter();
})();

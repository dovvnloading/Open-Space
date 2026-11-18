import { icon } from "./icons.js";

import "./observation-viewer.js";
import { esc } from "./html.js";
export { esc };

export function badge(status) {
  return `<span class="status" data-status="${esc(status)}">${esc(status)}</span>`;
}

export function sourceLink(link) {
  return `<a href="${esc(link.url)}" target="_blank" rel="noopener">${esc(link.label)}${icon("external")}</a>`;
}

const count = (value) => new Intl.NumberFormat("en-US").format(value);
const revisionLink = (catalog) => `${catalog.repository}/commit/${catalog.revision}`;
const definitions = {
  Measured: "The catalog reports a measurement from public observations.",
  Candidate: "The catalog identifies a promising signal that still needs independent confirmation.",
  Prediction: "The catalog describes a forecast that an observation can test.",
  Correction: "The catalog reports a correction to a data product, catalog or published result.",
  "Null result": "The search did not find the proposed signal within its tested sample and sensitivity.",
  "Research record":
    "This report or index appears in a research folder outside the classified catalog entries.",
};

function sourceList(links) {
  return `<div class="source-links">${links.map(sourceLink).join("")}</div>`;
}

function stamp(catalog) {
  return sourceLink({ url: revisionLink(catalog), label: "Revision " + catalog.revision.slice(0, 8) });
}

function statusLine(record) {
  return `<div class="disposition-line">${badge(record.status)}<span>${esc(definitions[record.status])}</span></div>`;
}

function archiveOrigin(record, catalog) {
  const origin = record.provenance;
  return `<div class="origin-strip">${origin ? sourceLink({ url: origin.url, label: origin.kind === "catalog" ? "Catalog entry" : "Research folder source" }) : ""}
    ${stamp(catalog)}</div>`;
}

function folderCoverage(record) {
  if (!record.evidenceChecks?.length) return "";
  const rows = record.evidenceChecks.map((row) => {
    const link = sourceLink({ url: row.url, label: row.folder });
    const checks = [
      `${count(row.matchedDigests)} matching digest assertions`,
      `${count(row.tracedNumbers)} of ${count(row.preciseNumbers)} precise numbers text-matched`,
      `${count(row.reproducedClaims)} of ${count(row.declaredClaims)} declared calculations reproduced`,
      `${count(row.openWarnings)} open verifier warnings`,
    ];
    return `<li>${link}<span>${esc(checks.join(" · "))}</span></li>`;
  });
  return `<details class="methods"><summary>Checks on the linked research folders</summary>
    <p>These counts cover entire folders, which may contain several findings. A number match shows that the value occurs in a machine-readable file.</p>
    <ul class="coverage-list">${rows.join("")}</ul></details>`;
}

export function renderRecord(record, catalog) {
  return `<article class="dialog-padding record-reader">
    <div class="eyebrow"><span class="small-line"></span>${esc(record.topic)}</div>
    <h2 id="dossier-title">${esc(record.title)}</h2>
    ${statusLine(record)}
    <p class="reader-lede">${esc(record.summary)}</p>
    <a class="image-study-link" href="#universe?target=${encodeURIComponent(catalog.destinations?.find(d => d.recordId === record.id)?.id || record.id)}">${icon("orbit")}Fly to investigation${icon("arrow")}</a>
    <h3>Source evidence</h3>${sourceList(record.links)}
    ${folderCoverage(record)}${archiveOrigin(record, catalog)}</article>`;
}

const tablePageSize = 8;

export function observationTable(table, {open = false, sourceRows = false} = {}) {
  if (!table?.rows.length) return "";
  const headings = table.columns.map((column) => `<th scope="col">${esc(column)}</th>`).join("");
  const body = table.rows
    .map(
      (row, i) =>
        `<tr ${i >= tablePageSize ? "hidden" : ""}>${row.map((value, column) => {
          const wrap = sourceRows && String(value).length > (column === 0 ? 18 : 32) && /\s/.test(value);
          return wrap ? `<td><span class="source-cell-text">${esc(value)}</span></td>` : `<td>${esc(value)}</td>`;
        }).join("")}</tr>`,
    )
    .join("");
  const rowKind = sourceRows ? "source" : "observation";
  return `<details class="methods table-reader" data-first-row="0" ${sourceRows ? "data-source-rows" : ""} ${open ? "open" : ""}><summary>${sourceRows ? "Source rows" : "Inspect observation rows"} (${count(table.rows.length)})</summary>
    <div class="table-scroll" tabindex="0" role="region" aria-label="${sourceRows ? "Source data" : "Observation data"}">
      <table><caption>${esc(table.caption)}</caption><thead><tr>${headings}</tr></thead><tbody>${body}</tbody></table>
    </div>
    <div class="table-nav">${sourceLink(table.source)}
      <span role="status">Rows 1–${Math.min(tablePageSize, table.rows.length)} of ${table.rows.length}</span>
      ${table.rows.length > tablePageSize ? `<button type="button" data-step="-1" aria-disabled="true" aria-label="Previous ${rowKind} rows">${icon("left")}</button><button type="button" data-step="1" aria-label="Next ${rowKind} rows">${icon("right")}</button>` : ""}
    </div></details>`;
}

export function moveObservationPage(button) {
  if (button.getAttribute("aria-disabled") === "true") return;
  const reader = button.closest(".table-reader");
  const rows = Array.from(reader.querySelectorAll("tbody tr"));
  const start = Number(reader.dataset.firstRow) + Number(button.dataset.step) * tablePageSize;
  const end = Math.min(start + tablePageSize, rows.length);
  reader.dataset.firstRow = String(start);
  rows.forEach((row, i) => {
    row.hidden = i < start || i >= end;
  });
  reader.querySelector('[role="status"]').textContent = `Rows ${start + 1}–${end} of ${rows.length}`;
  reader.querySelector('[data-step="-1"]').setAttribute("aria-disabled", String(start === 0));
  reader.querySelector('[data-step="1"]').setAttribute("aria-disabled", String(end === rows.length));
}

export function scienceFigure(figure, { inspect = false, label = "Observation figure" } = {}) {
  return `<figure class="science-plate">
      <div class="figure-head"><span>${icon("pulse")} ${esc(label)}</span>
      ${inspect ? `<button type="button" data-rfi-enlarge>Enlarge ${icon("focus")}</button>`
        : `<a href="./${esc(figure.src)}" target="_blank" rel="noopener">Full size ${icon("external")}</a>`}</div>
      <img src="./${esc(figure.src)}" alt="${esc(figure.alt)}" loading="lazy">
      <figcaption>${esc(figure.caption)} ${figure.credit ? esc(figure.credit) + ". " : ""}${sourceLink({ url: figure.url, label: "Figure source" })}</figcaption>
    </figure>`;
}

export function renderTarget(target, record, catalog) {
  const linked = new Set(target.links.map((link) => link.url));
  const sources = [...target.links, ...record.links.filter((link) => !linked.has(link.url))];
  const summary = target.measurements
    .map((m) => `<div class="measure"><strong>${esc(m.value)}</strong><span>${esc(m.label)}</span></div>`)
    .join("");
  const methods = target.measurements
    .filter((m) => m.method)
    .map(
      (m) => `<li><h4>${esc(m.label)}</h4>
    <p>${esc(m.method)}</p>${sourceList(m.sources || [])}</li>`,
    )
    .join("");
  const figure = target.figure;
  const interpretation = [
    ["What was checked", target.checked],
    ["What remains open", target.limitation],
    ["Next test", target.nextTest],
  ]
    .filter(([, value]) => value)
    .map(([label, value]) => `<div><h3>${esc(label)}</h3><p>${esc(value)}</p></div>`)
    .join("");
  return `<article class="dialog-padding target-reader">
    <div class="reader-intro"><div class="eyebrow"><span class="small-line"></span>${esc(target.instrument)}<span class="reader-epoch">${esc(target.epoch)}</span></div>
    <h2 id="dossier-title">${esc(target.title)}</h2>${statusLine(record)}
    <p class="reader-lede">${esc(target.body)}</p></div>
    <a class="image-study-link" href="#universe?target=${encodeURIComponent(target.id)}">${icon("orbit")}Fly to investigation${icon("arrow")}</a>
    <observation-viewer data-target="${esc(target.id)}">
    ${scienceFigure(figure)}
    </observation-viewer>
    ${observationTable(target.table)}
    <a class="image-study-link" href="#plate/${esc(target.id)}">${icon("focus")}${esc({
      mars: "Inspect the Viking mosaic", sun: "Watch the coronagraph image sequence",
      quasar: "Locate images A and B on the survey plate", dwarf: "View the credited illustration"
    }[target.id])}${icon("arrow")}</a>
    ${target.observingUnit ? `<dl class="observation-unit"><dt>Observing unit</dt><dd>${esc(target.observingUnit)}</dd></dl>` : ""}
    <div class="measures">${summary}</div>
    ${methods ? `<details class="methods"><summary>How these values were checked</summary><ol>${methods}</ol></details>` : ""}
    <div class="interpretation">${interpretation}</div>
    <h3>Source evidence</h3>${sourceList(sources)}${archiveOrigin(record, catalog)}
    </article>`;
}

export function renderAudit(catalog) {
  const audit = catalog.audit;
  if (!audit)
    return `<div class="dialog-padding"><h2 id="dossier-title">Archive checks</h2><p>This bundled snapshot predates the detailed audit. The source revision and record links are available in the catalog.</p>${stamp(catalog)}</div>`;
  const checks = audit.repositoryChecks;
  const stats = [
    [audit.catalogEntries, "Classified catalog entries"],
    [audit.supplementalRecords, "Additional research records"],
    [audit.recordLinks, "Links from records to evidence"],
  ]
    .map(
      ([value, label]) =>
        `<div class="measure"><strong>${count(value)}</strong><span>${esc(label)}</span></div>`,
    )
    .join("");
  const checksBody = checks
    ? `<dl class="audit-ledger">
    <div><dt>File integrity</dt><dd><b>${count(checks.matchedDigests)}</b> digest assertions matched files. The ledger separately records ${count(checks.externalReceipts)} upstream receipts.</dd></div>
    <div><dt>Catalog number tracing</dt><dd><b>${count(checks.tracedCatalogNumbers)} of ${count(checks.preciseCatalogNumbers)}</b> precise catalog numbers occurred in a linked machine-readable result. The trace covered ${count(checks.catalogTraceEntries)} catalog entries.</dd></div>
    <div><dt>Declared calculations</dt><dd><b>${count(checks.reproducedClaims)} of ${count(checks.declaredClaims)}</b> arithmetic checks reproduced their expected values across ${count(checks.claimFiles)} claim files.</dd></div>
    <div><dt>Report text</dt><dd><b>${count(checks.untracedReportNumbers)} of ${count(checks.reportNumbers)}</b> precise report numbers lacked a machine-readable match under the verifier’s rules. Some quantities require a manual reading of their source.</dd></div>
    </dl><p>A checksum establishes file identity. Text tracing locates a value. Confirming a physical interpretation requires appropriate controls or another observation.</p>
    ${sourceList([{ url: checks.methodUrl, label: "Verifier method" }])}`
    : "<p>This snapshot includes source-link checks. It does not include a repository verification ledger.</p>";
  return `<article class="dialog-padding audit-reader"><div class="eyebrow"><span class="small-line"></span>ARCHIVE PROVENANCE</div>
    <h2 id="dossier-title">How the archive was checked</h2>
    <p class="reader-lede">The archive preserves the catalog’s scientific dispositions. New research folders can appear before anyone has assigned a finding a disposition.</p>
    <div class="measures">${stats}</div>
    <p>The build resolves ${count(audit.uniqueEvidencePaths)} distinct evidence paths at one Git revision and checks linked Markdown sections. Each classified record links to the catalog line that supplied its disposition.</p>
    <h3>Replaying the featured observations</h3>
    <p>Each featured readout is recalculated from its source table. The replay notes document which original instrument files were recovered, which inputs came from the repository and which interpretations remain open.</p>
    ${sourceList([{ url: `${catalog.repository}/blob/${catalog.revision}/observatory/data/spotlight-replay.md`, label: "Source replay and limitations" }])}
    <h3>Repository checks</h3>${checksBody}
    <div class="origin-strip">${stamp(catalog)}<a href="./data/catalog.json" download="open-space-catalog.json">Download catalog JSON ${icon("arrow")}</a></div></article>`;
}

export function renderCredits(plates, catalog) {
  const images = Object.values(plates)
    .map(
      (p) => `<li><img src="./assets/${esc(p.file)}" alt="">
    <div><h3>${esc(p.name)} <small>${esc(p.kind)}</small></h3><p>${esc(p.description)}</p>
    ${sourceLink({ label: p.credit, url: p.source })}</div></li>`,
    )
    .join("");
  return `<article class="dialog-padding credits-reader"><div class="eyebrow"><span class="small-line"></span>IMAGE PROVENANCE</div>
    <h2 id="dossier-title">Image &amp; data credits</h2><p>The deck’s context images introduce the targets. Choose <b>Measured data</b> on a target to inspect its observation figure.</p>
    <ul class="credit-list">${images}</ul><h3>Simulation and type</h3>
    <p>The simulator generates its sky and planetary systems from integer addresses. Instrument Serif, DM Sans and IBM Plex Mono are distributed under the SIL Open Font License.</p>
    ${
      catalog
        ? sourceList([
            {
              label: "Image sources and reuse terms",
              url: `${catalog.repository}/blob/${catalog.revision}/observatory/assets/CREDITS.md`,
            },
            { label: "Repository", url: catalog.repository },
          ])
        : ""
    }
    </article>`;
}

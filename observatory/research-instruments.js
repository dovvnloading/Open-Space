import { badge, esc, moveObservationPage, observationTable, scienceFigure, sourceLink } from "./evidence-view.js";
import { openImageDesk } from "./image-desk.js";
import { icon } from "./icons.js";

let instrumentId = 0;
const plateIds = new Set(["mars", "sun", "quasar"]);

function sources(links) {
  const unique = new Map();
  for (const link of links) {
    if (link?.url && !unique.has(link.url)) unique.set(link.url, link);
  }
  return [...unique.values()].map(link => sourceLink({ url: link.url, label: link.label || link.path || "Source file" })).join("");
}

function reading(label, content) {
  return content ? `<section class="rfi-reading"><h3>${esc(label)}</h3><p>${esc(content)}</p></section>` : "";
}

function renderReport({ target, record, destination, catalog, figure, table, measurements }) {
  const description = target?.body || record?.summary || destination.summary;
  const excerpt = destination.excerpt;
  const quoted = typeof excerpt === "string" ? excerpt : excerpt?.text;
  const report = record?.provenance;
  const revision = catalog?.revision && catalog.repository ?
    { url: catalog.repository + "/commit/" + catalog.revision, label: "Snapshot " + catalog.revision.slice(0, 8) } : null;
  const position = destination.coordinates;
  const reportSources = sources([
    ...(target?.links || []), ...(record?.links || []), excerpt?.source,
    table?.source, figure?.url && {url: figure.url, label: "Original figure"},
    destination.image?.url && {url: destination.image.url, label: "Source image"},
    ...(measurements.flatMap(m => m.sources || [])), position?.source,
    report && {url: report.url, label: report.kind === "catalog" ? "Catalog entry" : "Research folder source"}, revision,
  ]);
  return `<div class="rfi-report">
    ${reading("Investigation", description)}
    ${quoted && quoted !== description ? `<blockquote>${esc(quoted)}
      ${destination.excerptSource ? `<cite>${sourceLink(destination.excerptSource)}</cite>` : ""}</blockquote>` : ""}
    ${reading("Observing unit", target?.observingUnit || destination.observingUnit)}
    ${position ? `<section class="rfi-reading"><h3>${esc(position.label || "Sky position")}</h3><dl class="rfi-coordinates">
      <div><dt>RA</dt><dd>${esc(position.raDeg)}°</dd></div><div><dt>Dec</dt><dd>${esc(position.decDeg)}°</dd></div></dl>
      ${position.frame ? `<small>${esc(position.frame)}</small>` : ""} ${position.source?.url ? sourceLink(position.source) : ""}</section>` : ""}
    ${reading("What was checked", target?.checked || destination.checked)}
    ${reading("What remains open", target?.limitation || destination.limitation)}
    ${reading("Next test", target?.nextTest || destination.nextTest)}
    ${reportSources ? `<section class="rfi-reading"><h3>Sources</h3><div class="rfi-sources">${reportSources}</div></section>` : ""}
  </div>`;
}

function renderData({target, record, figure, table, measurements}) {
  const viewer = Boolean(target && figure && table?.rows?.length);
  const measures = measurements.map(m => `<div><dt>${esc(m.label)}</dt><dd>${esc(m.value)}</dd></div>`).join("");
  const methods = measurements.filter(m => m.method || m.sources?.length).map(m => `<li>
    <h4>${esc(m.label)}</h4>${m.method ? `<p>${esc(m.method)}</p>` : ""}
    <div class="rfi-sources">${sources(m.sources || [])}</div></li>`).join("");
  const preview = table?.totalRows > table?.rows.length ?
    `<p class="rfi-preview">Preview of ${table.rows.length.toLocaleString("en-US")} of ${table.totalRows.toLocaleString("en-US")} source rows.</p>` : "";
  return `<article class="rfi-data">
    ${!target && record?.summary ? `<p class="rfi-context">${esc(record.summary)}</p>` : ""}
    ${viewer ? `<observation-viewer data-target="${esc(target.id)}">${scienceFigure(figure, {inspect: true})}</observation-viewer>` : ""}
    ${!viewer && measures ? `<dl class="rfi-values">${measures}</dl>` : ""}
    ${table?.rows?.length ? observationTable(table, {open: !viewer, sourceRows: !viewer}) : ""}
    ${preview}${viewer && measures ? `<dl class="rfi-values">${measures}</dl>` : ""}
    ${methods ? `<details class="rfi-methods"><summary>How the values were obtained</summary><ul>${methods}</ul></details>` : ""}
  </article>`;
}

function addColumnControls(panel) {
  const scroll = panel.querySelector("[data-source-rows] .table-scroll");
  if (!scroll) return;
  const controls = document.createElement("div");
  controls.className = "rfi-columns";
  controls.hidden = true;
  controls.innerHTML = `<span>Scroll columns</span><button type="button" aria-label="Scroll columns left">${icon("left")}</button>
    <button type="button" aria-label="Scroll columns right">${icon("right")}</button>`;
  scroll.before(controls);
  const [left, right] = controls.querySelectorAll("button");
  const update = () => {
    controls.hidden = !scroll.clientWidth || scroll.scrollWidth <= scroll.clientWidth + 1;
    left.setAttribute("aria-disabled", String(scroll.scrollLeft < 1));
    right.setAttribute("aria-disabled", String(scroll.scrollLeft + scroll.clientWidth >= scroll.scrollWidth - 1));
  };
  const move = direction => scroll.scrollBy({left: direction * Math.round(scroll.clientWidth * .75)});
  left.onclick = () => move(-1);
  right.onclick = () => move(1);
  const observer = new ResizeObserver(update);
  observer.observe(scroll);
  observer.observe(scroll.querySelector("table"));
  scroll.addEventListener("scroll", update, {passive: true});
  update();
  return () => observer.disconnect();
}

function inspectFigure(root, figure) {
  const dialog = document.createElement("dialog");
  dialog.className = "rfi-lightbox";
  dialog.setAttribute("aria-label", "Full-size source image");
  dialog.innerHTML = `<div class="rfi-light-tools"><button type="button" data-rfi-size aria-pressed="false">Actual pixels</button>
    <button type="button" data-rfi-close aria-label="Close full-size image">${icon("x")}</button></div>
    <div class="rfi-light-scroll" tabindex="0" role="region" aria-label="Source image"><img alt="${esc(figure.alt)}"></div>
    <footer>${esc(figure.caption)} ${sourceLink({url:figure.url, label:"Source"})}</footer>`;
  root.append(dialog);
  dialog.querySelector("img").src = "./" + figure.src;
  dialog.addEventListener("close", () => dialog.remove(), {once:true});
  dialog.addEventListener("cancel", event => event.stopPropagation());
  dialog.querySelector("[data-rfi-close]").onclick = () => dialog.close();
  dialog.querySelector("[data-rfi-size]").onclick = event => {
    const native = dialog.classList.toggle("rfi-actual-size");
    event.currentTarget.setAttribute("aria-pressed", String(native));
  };
  dialog.showModal();
}

export function mountResearchInstruments(host, { destination, record, target, plate, catalog }) {
  const figure = target?.figure || destination.figure;
  const table = target?.table || destination.table;
  const measurements = target?.measurements || destination.measurements || [];
  const id = target?.id || destination.id;
  const desk = plateIds.has(id) && plate;
  const image = destination.image;
  const hasData = measurements.length || table?.rows?.length;
  const sections = [
    ...(hasData ? [["data", "Measurements", "pulse"]] : []),
    ...(desk || figure || image ? [["images", id === "dwarf" ? "Radio images" : desk || image ? "Imagery" : "Figure", "focus"]] : []),
    ["report", "Report", "book"],
  ];
  const prefix = "rfi-" + ++instrumentId;
  const root = document.createElement("section");
  root.className = "rfi-console";
  root.setAttribute("aria-label", "Research instruments");
  const context = {target, record, destination, catalog, figure, table, measurements};
  const tabs = sections.map(([key,label,glyph]) => `<button type="button" id="${prefix}-${key}-tab" role="tab"
    aria-controls="${prefix}-${key}" data-rfi-tab="${key}">${icon(glyph)}<span>${label}</span></button>`).join("");
  const panels = sections.map(([key]) => `<div id="${prefix}-${key}" role="tabpanel" tabindex="0"
    aria-labelledby="${prefix}-${key}-tab" data-rfi-panel="${key}" hidden></div>`).join("");
  const metadata = [target?.instrument, target?.epoch].filter(Boolean).join(" · ") || record?.topic;
  root.innerHTML = `<h2 class="rfi-mobile-target">${esc(destination.name || target?.target || record?.title)}</h2>
    <header class="rfi-meta">${record?.status ? badge(record.status) : ""}
    ${metadata ? `<span>${esc(metadata)}</span>` : ""}</header>
    <div class="rfi-tools" role="tablist" aria-label="Evidence views">${tabs}</div>${panels}`;
  host.replaceChildren(root);
  let clearColumns;

  function activate(key) {
    root.querySelectorAll("[data-rfi-tab]").forEach(tab => {
      const active = tab.dataset.rfiTab === key;
      tab.setAttribute("aria-selected", String(active));
      tab.tabIndex = active ? 0 : -1;
    });
    for (const panel of root.querySelectorAll("[data-rfi-panel]")) {
      panel.hidden = panel.dataset.rfiPanel !== key;
      if (panel.hidden) {
        panel.querySelector("image-desk")?.stop();
        continue;
      }
      if (panel.childElementCount) continue;
      if (key === "data") {
        panel.innerHTML = renderData(context);
        clearColumns = addColumnControls(panel);
      }
      if (key === "report") panel.innerHTML = renderReport(context);
      if (key === "images") {
        if (desk) panel.append(openImageDesk(id, {[id]:plate}, {embedded:true}));
        if (!desk && figure) panel.innerHTML = scienceFigure(figure, {inspect:true, label:target ? "Observation figure" : "Repository figure"});
        if (!desk && image && image.src !== figure?.src) panel.insertAdjacentHTML("afterbegin", scienceFigure(image, {inspect:true, label:"Source image"}));
      }
    }
  }
  // Catch load failures before parsing the first panel, including a cached broken image.
  root.addEventListener("error", event => {
    const img = event.target;
    if (!img.matches?.(".science-plate img, .rfi-light-scroll img")) return;
    const frame = img.parentElement;
    if (frame.querySelector(".rfi-image-error")) return;
    img.hidden = true;
    const notice = document.createElement("div");
    notice.className = "rfi-image-error";
    notice.innerHTML = '<span role="status">Image unavailable</span><button type="button">Retry image</button>';
    notice.querySelector("button").onclick = () => {
      img.hidden = false;
      notice.remove();
      img.src = img.src;
    };
    img.after(notice);
  }, true);
  root.addEventListener("click", event => {
    const button = event.target.closest("button");
    if (!button) return;
    if (button.dataset.rfiTab) activate(button.dataset.rfiTab);
    if (button.matches(".table-reader [data-step]")) moveObservationPage(button);
    if (button.hasAttribute("data-rfi-enlarge")) {
      const src = button.closest("figure").querySelector("img").getAttribute("src");
      const source = [figure, image].find(item => item && "./" + item.src === src);
      if (source) inspectFigure(root, source);
    }
  });
  root.addEventListener("keydown", event => {
    if (!event.target.matches("[data-rfi-tab]")) return;
    const buttons = [...root.querySelectorAll("[data-rfi-tab]")];
    const index = buttons.indexOf(event.target);
    const nextIndex = {ArrowLeft: (index-1+buttons.length)%buttons.length, ArrowRight:(index+1)%buttons.length,
      Home:0, End:buttons.length-1}[event.key];
    if (nextIndex === undefined) return;
    event.preventDefault();
    event.stopPropagation();
    activate(buttons[nextIndex].dataset.rfiTab);
    buttons[nextIndex].focus();
  });
  activate(hasData ? "data" : "report");
  return () => {
    clearColumns?.();
    root.querySelector(".rfi-lightbox[open]")?.close();
    root.remove();
  };
}

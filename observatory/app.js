import { hydrateIcons, icon } from "./icons.js";
import {
  esc,
  badge,
  renderAudit,
  renderCredits,
  renderTarget,
  renderRecord,
  moveObservationPage,
} from "./evidence-view.js";

const $ = (selector) => document.querySelector(selector);
const repoURL = "https://github.com/dovvnloading/Open-Space";
let plates = {};
let plateView = "context";
const pageSize = 15;
const searchIndex = new Map();
let renderedArchiveHash;
let catalogFailed = false;
let catalog,
  records = [],
  targets = [],
  slide = 0,
  activeView = "deck";
let disposeUniverse,
  routeCounter = 0,
  dialogTrigger,
  lastEntryView,
  modalKey;
const filters = { q: "", status: "", topic: "", sort: "recent", page: 1 };

function text(id, value) {
  document.getElementById(id).textContent = value;
}
function dateLabel(iso) {
  if (!iso) return "";
  return new Intl.DateTimeFormat("en", {
    day: "2-digit",
    month: "short",
    year: "numeric",
    timeZone: "UTC",
  }).format(new Date(iso));
}
function starPlate() {
  const canvas = $(".star-plate");
  const context = canvas.getContext("2d");
  if (!context) return;
  function draw() {
    const { width, height } = canvas.getBoundingClientRect();
    if (!width || !height) return;
    const ratio = Math.min(window.devicePixelRatio || 1, 2);
    canvas.width = Math.round(width * ratio);
    canvas.height = Math.round(height * ratio);
    context.setTransform(ratio, 0, 0, ratio, 0, 0);
    let seed = 1807;
    function rand() {
      seed = (Math.imul(seed, 1664525) + 1013904223) >>> 0;
      return seed / 4294967296;
    }
    const count = Math.max(100, (width * height) / 1390);
    for (let i = 0; i < count; i++) {
      const x = rand() * width,
        y = rand() * height,
        bright = rand();
      const alpha = (0.23 + bright * 0.62) * (x < width * 0.4 ? 0.5 : 1);
      context.fillStyle = `rgba(218,213,206,${alpha})`;
      context.beginPath();
      context.arc(x, y, 0.28 + bright * 0.82, 0, 2 * Math.PI);
      context.fill();
      if (bright > 0.979) {
        context.fillStyle = `rgba(220,216,200,${alpha * 0.6})`;
        context.fillRect(x - 2.7, y - 0.36, 5.4, 0.72);
        context.fillRect(x - 0.36, y - 2.7, 0.72, 5.4);
      }
    }
    context.strokeStyle = "#babcbd3c";
    context.lineWidth = 0.65;
    for (const [px, py] of [
      [0.45, 0.17],
      [0.5, 0.82],
      [0.865, 0.81],
      [0.94, 0.2],
    ]) {
      const x = width * px,
        y = height * py;
      context.beginPath();
      context.moveTo(x - 4, y);
      context.lineTo(x + 4, y);
      context.moveTo(x, y - 4);
      context.lineTo(x, y + 4);
      context.stroke();
    }
  }
  new ResizeObserver(draw).observe(canvas);
}
function showPlate(index) {
  if (!targets.length) return;
  slide = (index + targets.length) % targets.length;
  const target = targets[slide],
    plate = plates[target.id];
  $(".hero").dataset.target = target.id;
  $("#hero-art").src = "./assets/" + plate.file;
  text("target-tag", plate.tag);
  updatePlateView(target, plate);
  text("hero-name", plate.name);
  text("hero-sub", target.subtitle);
  text("hero-summary", target.summary);
  text("instrument", target.instrument);
  text("epoch", target.epoch);
  text("disposition", target.status);
  $("#hero-detail").href = "#target/" + target.id;
  $("#slide-count").innerHTML =
    String(slide + 1).padStart(2, "0") + " <em>/ " + String(targets.length).padStart(2, "0") + "</em>";
}
function updatePlateView(target, plate) {
  const figureMode = plateView === "figure";
  $(".hero").dataset.view = plateView;
  $(".art-stage").hidden = figureMode;
  $("#plate-figure").hidden = !figureMode;
  for (const button of document.querySelectorAll("[data-plate-view]")) {
    button.setAttribute("aria-pressed", String(button.dataset.plateView === plateView));
  }
  const credit = $("#plate-credit");
  if (figureMode) {
    const fig = target.figure;
    const figureImage = $("#plate-figure-image");
    figureImage.src = "./" + fig.src;
    figureImage.alt = fig.alt;
    $("#plate-figure-link").href = "#target/" + target.id;
    credit.href = "#target/" + target.id;
    credit.textContent = "Inspect figure and method";
  } else {
    credit.href = "#credits";
    credit.textContent = plate.kind + " · " + plate.credit;
  }
}

function renderCard(target) {
  const p = plates[target.id];
  return `<a class="field-card" href="#target/${esc(target.id)}">
    <div class="card-art" data-plate="${esc(target.id)}"><img src="./assets/${esc(p.file)}" alt="" width="1280" height="533"><span class="card-origin">${esc(p.kind)}</span><span class="card-type" style="--ch:${p.color}">${esc(p.field.toUpperCase())}</span></div>
    <div class="card-copy"><div class="card-meta">${icon("antenna")}${esc(target.instrument)}</div><h3>${esc(target.title)}</h3><p>${esc(target.summary)}</p><div class="card-tail">${badge(target.status)}${icon("arrow")}</div></div></a>`;
}
function newest(a, b) {
  return (Date.parse(b.updatedAt) || 0) - (Date.parse(a.updatedAt) || 0) || a.title.localeCompare(b.title);
}
function resultRow(record, detail) {
  const spotlight = targets.find((target) => target.recordId === record.id);
  const url = spotlight ? "#target/" + spotlight.id : "#record/" + encodeURIComponent(record.id);
  return `<a class="result-row" href="${esc(url)}"><span class="row-icon">${icon("book")}</span><div class="row-content"><small>${esc(record.topic)}</small><h3>${esc(record.title)}</h3>${detail ? "<p>" + esc(record.summary) + "</p>" : ""}</div>${badge(record.status)}<span class="row-date">${esc(dateLabel(record.updatedAt))}</span>${icon("right")}</a>`;
}
function fillFilters(element, values) {
  for (const value of Array.from(new Set(values)).sort()) element.add(new Option(value, value));
}
function prepareCatalog() {
  records = catalog.records;
  targets = catalog.spotlights.filter((target) => Object.hasOwn(plates, target.id));
  for (const record of records) {
    const target = targets.find((item) => item.recordId === record.id);
    searchIndex.set(
      record.id,
      [
        record.title,
        record.summary,
        record.status,
        record.topic,
        record.path,
        ...record.links.map((link) => link.label),
        target?.target,
        target?.title,
        target?.instrument,
        target?.summary,
      ]
        .join(" ")
        .toLowerCase(),
    );
  }
  showPlate(0);
  $("#cards").innerHTML = targets
    .filter((t) => t.id !== "mars")
    .map(renderCard)
    .join("");
  $("#recent-list").innerHTML = [...records]
    .sort(newest)
    .slice(0, 4)
    .map((r) => resultRow(r, false))
    .join("");
  fillFilters(
    $("#status-filter"),
    records.map((r) => r.status),
  );
  fillFilters(
    $("#topic-filter"),
    records.map((r) => r.topic),
  );
  text(
    "archive-intro",
    `Explore ${records.length} public research records, from measured signals and candidate discoveries to forecasts, corrections, and null results.`,
  );
  text("snapshot-date", dateLabel(catalog.revisionDate) + " UTC");
  const revision = catalog.revision;
  $("#revision").innerHTML =
    `<a href="${esc(catalog.repository)}/commit/${esc(revision)}" rel="noopener" target="_blank">Revision ${esc(revision.slice(0, 7))}</a>`;
}
function archiveHash() {
  const params = new URLSearchParams();
  for (const [key, value] of Object.entries(filters)) {
    if (value && !(key === "page" && value === 1) && !(key === "sort" && value === "recent"))
      params.set(key, String(value));
  }
  return "#archive" + (params.size ? "?" + params : "");
}
function readFilters(search) {
  const params = new URLSearchParams(search);
  filters.q = params.get("q") || "";
  for (const [key, id] of [
    ["status", "status-filter"],
    ["topic", "topic-filter"],
  ]) {
    const value = params.get(key);
    filters[key] = [...document.getElementById(id).options].some((o) => o.value === value) ? value : "";
  }
  filters.sort = params.get("sort") === "title" ? "title" : "recent";
  filters.page = Math.max(1, parseInt(params.get("page"), 10) || 1);
  $("#query").value = filters.q;
  $("#status-filter").value = filters.status;
  $("#topic-filter").value = filters.topic;
  $("#sort").value = filters.sort;
}
function renderArchive() {
  const words = filters.q.toLowerCase().trim().split(/\s+/).filter(Boolean);
  const found = records
    .filter((r) => {
      if ((filters.status && r.status !== filters.status) || (filters.topic && r.topic !== filters.topic))
        return false;
      const content = searchIndex.get(r.id);
      return words.every((word) => content.includes(word));
    })
    .sort(filters.sort === "title" ? (a, b) => a.title.localeCompare(b.title) : newest);
  const pages = Math.max(1, Math.ceil(found.length / pageSize));
  filters.page = Math.min(filters.page, pages);
  renderedArchiveHash = archiveHash();
  const start = (filters.page - 1) * pageSize;
  text(
    "results-count",
    found.length +
      (found.length === 1 ? " record" : " records") +
      (found.length ? " · " + (start + 1) + "–" + Math.min(start + pageSize, found.length) : ""),
  );
  $("#result-list").innerHTML = found.length
    ? found
        .slice(start, start + pageSize)
        .map((r) => resultRow(r, true))
        .join("")
    : '<div class="empty">No investigations match these filters.<br><button type="button" id="reset-filters">Clear filters</button></div>';
  $("#pagination").innerHTML =
    pages > 1
      ? `<button class="pager" data-page="-1" ${filters.page === 1 ? "disabled" : ""}>Previous</button><span>Page ${filters.page} of ${pages}</span><button class="pager" data-page="1" ${filters.page === pages ? "disabled" : ""}>Next</button>`
      : "";
}
function updateFilter() {
  filters.q = $("#query").value;
  filters.status = $("#status-filter").value;
  filters.topic = $("#topic-filter").value;
  filters.sort = $("#sort").value;
  filters.page = 1;
  history.replaceState(history.state, "", archiveHash());
  renderArchive();
}
function showView(view) {
  activeView = view;
  $("#deck").hidden = view !== "deck";
  $("#archive").hidden = view !== "archive";
  for (const a of document.querySelectorAll("[data-nav]")) {
    if (a.dataset.nav === view) a.setAttribute("aria-current", "page");
    else a.removeAttribute("aria-current");
  }
  document.title =
    (view === "archive" ? "Research archive" : "Observation deck") + " · Open Space Observatory";
}
function closeOverlay() {
  if (history.state?.fromDeckSite) history.back();
  else location.hash = activeView === "archive" ? archiveHash() : "#deck";
}
function renderDossier(key, target, record) {
  const panel = $("#dossier");
  const body = $("#dossier-body");
  if (panel.open && modalKey === key) return;
  modalKey = key;
  text(
    "reader-label",
    key === "credits"
      ? "Observatory credits"
      : key === "audit"
        ? "Archive checks"
        : target
          ? "Featured investigation"
          : "Archive record",
  );
  if (key === "credits") body.innerHTML = renderCredits(plates, catalog);
  else if (key === "audit") body.innerHTML = renderAudit(catalog);
  else if (target)
    body.innerHTML = renderTarget(
      target,
      records.find((r) => r.id === target.recordId),
      catalog,
    );
  else if (record) body.innerHTML = renderRecord(record, catalog);
  else
    body.innerHTML =
      '<div class="dialog-padding"><h2 id="dossier-title">Record unavailable</h2><p>This identifier does not appear in this catalog revision.</p><a class="plain-link" href="#archive">Browse the research archive</a></div>';
  if (!panel.open) {
    dialogTrigger = document.activeElement;
    panel.showModal();
  }
  panel.scrollTop = 0;
  document.title = $("#dossier-title").textContent + " · Open Space Observatory";
}
async function route() {
  const generation = ++routeCounter;
  let restoreFocus;
  let [path, search] = location.hash.slice(1).split("?");
  if (!path || path === "main") path = "deck";
  if (disposeUniverse && path !== "universe") {
    disposeUniverse();
    disposeUniverse = null;
  }
  const [kind, key] = path.split("/");
  if (
    !catalog &&
    !catalogFailed &&
    ["universe", "archive", "target", "record", "audit", "credits"].includes(kind)
  )
    return;
  if (!["target", "record", "credits", "audit"].includes(kind) && $("#dossier").open) {
    $("#dossier").close();
    modalKey = null;
    restoreFocus = dialogTrigger;
    dialogTrigger = null;
  }
  if (kind === "universe") {
    if (disposeUniverse) return;
    try {
      const { mountUniverse } = await import("./universe.js");
      if (generation !== routeCounter) return;
      document.title = "Universe simulator · Open Space Observatory";
      disposeUniverse = mountUniverse($("#universe-host"), {
        targets,
        onClose: closeOverlay,
        onEvidence: (t) => {
          history.replaceState(history.state, "", "#target/" + encodeURIComponent(t.id));
          route();
        },
      });
    } catch (error) {
      showNotice("The universe simulator could not open.", error);
    }
    return;
  }
  if (["target", "record", "credits", "audit"].includes(kind)) {
    if (kind !== "credits" && !catalog) return;
    const target = kind === "target" ? targets.find((t) => t.id === key) : undefined;
    const record =
      kind === "record" ? records.find((r) => r.id === key || r.aliases?.includes(key)) : undefined;
    renderDossier(path, target, record);
    return;
  }
  const view = kind === "archive" ? "archive" : "deck";
  showView(view);
  if (view === "archive") {
    readFilters(search);
    if (renderedArchiveHash !== archiveHash()) renderArchive();
  }
  if (restoreFocus?.isConnected) restoreFocus.focus({ preventScroll: true });
  if (lastEntryView !== view) window.scrollTo(0, 0);
  lastEntryView = view;
}
function showNotice(message, error) {
  const notice = $("#notice");
  notice.hidden = false;
  notice.innerHTML =
    esc(message) +
    ' <button type="button">Reload</button> · <a class="plain-link" href="' +
    repoURL +
    '">Visit the repository</a>';
  notice.querySelector("button").onclick = () => location.reload();
  console.error(message, error);
}
hydrateIcons();
showView("deck");
starPlate();
$(".plate-controls").addEventListener("click", (event) => {
  const button = event.target.closest("[data-plate-view]");
  if (!button || !targets.length) return;
  plateView = button.dataset.plateView;
  updatePlateView(targets[slide], plates[targets[slide].id]);
});
$("#previous").onclick = () => showPlate(slide - 1);
$("#next").onclick = () => showPlate(slide + 1);
$("#filters").onsubmit = (e) => {
  e.preventDefault();
  updateFilter();
};
$("#query").addEventListener("input", updateFilter);
for (const id of ["#status-filter", "#topic-filter", "#sort"]) $(id).addEventListener("change", updateFilter);
$("#pagination").addEventListener("click", (e) => {
  const button = e.target.closest("[data-page]");
  if (!button || button.disabled) return;
  const nextPage = filters.page + Number(button.dataset.page);
  const params = new URLSearchParams(archiveHash().split("?")[1]);
  params.set("page", nextPage);
  history.pushState({ fromDeckSite: true }, "", "#archive?" + params);
  route();
  $("#results-count").setAttribute("tabindex", "-1");
  $("#results-count").focus({ preventScroll: true });
  $(".archive-toolbar").scrollIntoView({ block: "start" });
});
$("#result-list").addEventListener("click", (e) => {
  if (e.target.id === "reset-filters") {
    history.replaceState(history.state, "", "#archive");
    route();
    $("#query").focus();
  }
});
$(".close-dialog").onclick = closeOverlay;
$("#dossier").addEventListener("cancel", (e) => {
  e.preventDefault();
  closeOverlay();
});
$("#dossier").addEventListener("click", (e) => {
  const pageButton = e.target.closest(".table-reader [data-step]");
  if (pageButton) moveObservationPage(pageButton);
  if (e.target === $("#dossier")) closeOverlay();
});
document.addEventListener("click", (e) => {
  const a = e.target.closest("a[href^='#']");
  if (!a || e.defaultPrevented || e.button || e.ctrlKey || e.metaKey || e.shiftKey || e.altKey) return;
  const hash = a.getAttribute("href");
  if (hash === "#main") {
    e.preventDefault();
    $("#main").focus();
    $("#main").scrollIntoView();
    return;
  }
  e.preventDefault();
  if (location.hash === hash) return;
  history.pushState({ fromDeckSite: true }, "", hash);
  route();
});
window.addEventListener("hashchange", route);
async function start() {
  try {
    const loadJSON = async (path) => {
      const response = await fetch(path);
      if (!response.ok) throw Error(path + " HTTP " + response.status);
      return response.json();
    };
    [catalog, plates] = await Promise.all([
      loadJSON("./data/catalog.json"),
      loadJSON("./assets/plates.json"),
    ]);
    prepareCatalog();
    await route();
  } catch (error) {
    catalogFailed = true;
    showNotice("The research catalog could not load.", error);
    if (location.hash === "#universe") await route();
  }
}
start();

import("./universe.js")
  .then(({ drawUniversePreview, universeSize }) => {
    drawUniversePreview($("#universe-preview"));
    const total = (universeSize.galaxies * universeSize.systemsPerGalaxy).toLocaleString("en-US");
    text(
      "sim-scale",
      `Explore ${universeSize.galaxies} generated galaxies and ${total} star systems. Choose a star, enter orbit, and inspect its worlds.`,
    );
  })
  .catch((error) => console.warn("Universe preview unavailable", error));

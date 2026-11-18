import { esc } from "./html.js";
import { icon } from "./icons.js";

import { pointingArc } from "./flight-path.js";

const fixed = (d) => Number.isFinite(d.coordinates?.raDeg) && Number.isFinite(d.coordinates?.decDeg);
export const destinationHash = (id) => "#universe?target=" + encodeURIComponent(id);
export function flightImage(d, plates) {
  if (d.id === "dwarf" && (d.image || d.figure)) return (d.image || d.figure).src;
  return plates[d.id] ? "./assets/" + plates[d.id].file : d.image?.src || d.figure?.src || "";
}

function gridLines() {
  return (
    [0, 1, 2, 3, 4, 5, 6].map((i) => `<line x1="${i * 60}" x2="${i * 60}" y1="0" y2="180"/>`).join("") +
    [0, 45, 90, 135, 180].map((y) => `<line x1="0" x2="360" y1="${y}" y2="${y}"/>`).join("")
  );
}

function skyChart(items) {
  const pointings = items.filter(fixed);
  if (!pointings.length) return "";
  return `<details class="rf-sky-chart" open><summary>Sky pointings <span>${pointings.length}</span></summary>
    <div class="rf-chart-body">
    <div class="rf-chart-field"><svg aria-hidden="true" viewBox="0 0 360 180" preserveAspectRatio="none">${gridLines()}<path data-course fill="none" stroke="#f7bf8d" stroke-width="1.7" vector-effect="non-scaling-stroke"/></svg>
    ${pointings
      .map((d) => {
        const c = d.coordinates;
        return `<button type="button" data-sky-id="${esc(d.id)}" style="left:${100 * (1 - c.raDeg / 360)}%;top:${(100 * (90 - c.decDeg)) / 180}%" aria-label="Fly to ${esc(d.name)} at the recorded pointing" title="${esc(d.name)}"><i></i></button>`;
      })
      .join("")}
    <span class="rf-dec-pos">+90°</span><span class="rf-dec-neg">−90°</span></div>
    <div class="rf-axis"><span>24h</span><span>12h</span><span>0h</span></div>
    <p data-map-label>Select a recorded pointing to fly there.</p>
    <div class="rf-point-register">${pointings.map((d) => `<button type="button" data-field-go="${esc(d.id)}"><span>${esc(d.name)}</span><small>${esc(d.coordinates.frame)}</small></button>`).join("")}</div></div></details>`;
}

export function makeFlightNavigator(parent, { items, records, plates, targets, onSelect }) {
  const nav = document.createElement("dialog");
  nav.className = "rf-finder";
  nav.setAttribute("aria-labelledby", "rf-finder-title");
  parent.append(nav);
  const filters = [
    ["all", "All research", () => true],
    ["featured", "Featured", (d) => d.featured],
    ["image", "Images & figures", (d) => !!flightImage(d, plates)],
    ["sky", "Sky pointings", fixed],
  ].filter(([key, , test]) => key === "all" || items.some(test));
  nav.innerHTML = `<div class="rf-finder-head"><div><p class="rf-overline">Universe simulator</p><h2 id="rf-finder-title">Choose a destination</h2></div>
    <button class="rf-icon" type="button" data-dismiss aria-label="Close destinations">${icon("x")}</button></div>
    <div class="rf-finder-layout"><div class="rf-finder-chart">${skyChart(items)}</div><div class="rf-search-panel">
    <label class="rf-search">${icon("search")}<input type="search" aria-label="Search research destinations" placeholder="Search targets, instruments or findings" autocomplete="off"></label>
    <div class="rf-filters" role="group" aria-label="Destination filters">${filters.map(([key, label]) => `<button type="button" data-filter="${key}" aria-pressed="${key === "all"}">${label}</button>`).join("")}</div>
    <p class="rf-count" role="status"></p><div class="rf-results"></div></div></div>`;
  if (!items.some(fixed)) nav.classList.add("rf-no-map");
  const byId = new Map(items.map((d) => [d.id, d]));
  const ordered = [...items].sort(
    (a, b) =>
      (Date.parse(records.get(b.recordId)?.updatedAt) || 0) -
        (Date.parse(records.get(a.recordId)?.updatedAt) || 0) || a.name.localeCompare(b.name),
  );
  const index = new Map(
    ordered.map((d) => {
      const r = records.get(d.recordId),
        t = targets.get(d.id);
      return [
        d.id,
        [
          d.name,
          d.setting,
          r?.title,
          r?.topic,
          r?.summary,
          r?.status,
          t?.target,
          t?.instrument,
          d.instrument,
          ...(r?.links || []).map((l) => l.label),
        ]
          .filter(Boolean)
          .join(" ")
          .toLocaleLowerCase(),
      ];
    }),
  );
  const field = nav.querySelector("input");
  const resultList = nav.querySelector(".rf-results");
  const count = nav.querySelector(".rf-count");
  let mode = "all",
    selectedId,
    limit = 28,
    opener;
  const mapText = nav.querySelector("[data-map-label]");
  function render(more = false) {
    const previous = more ? resultList.querySelectorAll("[data-go]").length : 0;
    const words = field.value.trim().toLowerCase().split(/\s+/).filter(Boolean);
    const match = filters.find(([key]) => key === mode)[2];
    const found = ordered.filter((d) => match(d) && words.every((term) => index.get(d.id).includes(term)));
    const active = new Set(found.map((d) => d.id));
    for (const pin of nav.querySelectorAll("[data-sky-id]")) {
      pin.hidden = !active.has(pin.dataset.skyId);
      pin.classList.toggle("rf-current-pin", pin.dataset.skyId === selectedId);
    }
    for (const row of nav.querySelectorAll("[data-field-go]")) row.hidden = !active.has(row.dataset.fieldGo);
    nav
      .querySelectorAll("[data-filter]")
      .forEach((b) => b.setAttribute("aria-pressed", String(b.dataset.filter === mode)));
    count.textContent = found.length + " destination" + (found.length === 1 ? "" : "s");
    resultList.innerHTML =
      found
        .slice(0, limit)
        .map((d) => {
          const r = records.get(d.recordId),
            src = flightImage(d, plates);
          return `<button class="rf-result" data-go="${esc(d.id)}" type="button"><span class="rf-result-thumb">${src ? `<img src="${esc(src)}" alt="" loading="lazy">` : icon("book")}</span>
        <span><small>${esc(r?.status || "Research record")} · ${esc(d.setting || r?.topic || "Research")}</small><strong>${esc(d.name)}</strong></span>${icon("arrow")}</button>`;
        })
        .join("") || '<p class="rf-empty">No destinations match this search.</p>';
    if (found.length > limit)
      resultList.insertAdjacentHTML(
        "beforeend",
        '<button type="button" class="rf-more" data-more>Show more destinations</button>',
      );
    if (more) resultList.querySelectorAll("[data-go]")[previous]?.focus();
    else resultList.scrollTop = 0;
  }
  function open(trigger, current) {
    opener = trigger;
    selectedId = current;
    mode = "all";
    limit = 28;
    field.value = "";
    if (mapText) mapText.textContent = "Select a recorded pointing to fly there.";
    nav.querySelector("[data-course]")?.removeAttribute("d");
    render();
    const sky = nav.querySelector(".rf-sky-chart");
    if (sky) sky.open = matchMedia("(min-width: 761px)").matches;
    nav.showModal();
    if (matchMedia("(min-width: 741px)").matches) field.focus();
    else nav.querySelector("[data-dismiss]").focus({ preventScroll: true });
  }
  function close(restore = true) {
    if (nav.open) nav.close();
    if (restore && opener?.isConnected) opener.focus({ preventScroll: true });
  }
  nav.addEventListener("cancel", (e) => {
    e.preventDefault();
    close();
  });
  nav.addEventListener("click", (e) => {
    if (e.target === nav || e.target.closest("[data-dismiss]")) {
      close();
      return;
    }
    const btn = e.target.closest("button");
    if (!btn) return;
    if (btn.dataset.go || btn.dataset.skyId || btn.dataset.fieldGo) {
      const id = btn.dataset.go || btn.dataset.skyId || btn.dataset.fieldGo;
      close(false);
      onSelect(id);
      return;
    }
    if (btn.dataset.filter) {
      mode = btn.dataset.filter;
      limit = 28;
      render();
    } else if (btn.hasAttribute("data-more")) {
      limit += 28;
      render(true);
    }
  });
  field.addEventListener("input", () => {
    limit = 28;
    render();
  });
  function showPoint(e) {
    if (!mapText) return;
    const pin = e.target.closest("[data-sky-id],[data-field-go],[data-go]");
    if (!pin) return;
    const d = byId.get(pin.dataset.skyId || pin.dataset.fieldGo || pin.dataset.go);
    const path = nav.querySelector("[data-course]");
    if (!path) return;
    path.removeAttribute("d");
    nav
      .querySelectorAll("[data-sky-id]")
      .forEach((p) => p.classList.toggle("rf-pointed", p.dataset.skyId === d.id));
    if (!fixed(d)) {
      mapText.textContent = "Select a recorded pointing to fly there.";
      return;
    }
    const origin = byId.get(selectedId);
    const course =
      origin && fixed(origin) && origin.id !== d.id ? pointingArc(origin.coordinates, d.coordinates) : null;
    if (course) path.setAttribute("d", course.path);
    mapText.textContent =
      d.name +
      " · " +
      d.coordinates.raDeg.toFixed(3) +
      "°, " +
      d.coordinates.decDeg.toFixed(3) +
      "° · " +
      d.coordinates.frame +
      (course ? " · " + course.degrees.toFixed(1) + "° from " + origin.name : "");
  }
  nav.addEventListener("focusin", showPoint);
  nav.addEventListener("pointerover", showPoint);
  return {
    open,
    dispose() {
      close(false);
      nav.remove();
    },
  };
}

import { badge, esc } from "./evidence-view.js";
import { icon } from "./icons.js";

// Short index labels retain the full catalog topic in each button's title.
const topicLabels = new Map([
  ["The Sun, space weather and Earth's sky", "Sun, space weather & sky"],
  ["Rivers and wetlands on Earth", "Rivers & wetlands"],
  ["The Moon and Mars", "Moon & Mars"],
  ["Asteroids and comets", "Asteroids & comets"],
  ["Planets, moons and rings in our Solar System", "Solar system planets, moons & rings"],
  ["Planets around other stars", "Exoplanets"],
  ["Stars, binaries and multiple-star systems", "Stars & stellar systems"],
  ["White dwarfs", "White dwarfs"],
  ["Neutron stars, black holes and cosmic explosions", "Compact objects & explosions"],
  ["Galaxies and the distant universe", "Galaxies & the distant universe"],
  ["The project's first campaigns: galaxy spectra and survey methods", "Galaxy spectra & survey methods"],
  ["Research archive", "Additional research"],
]);

const dispositions = ["Measured", "Candidate", "Prediction", "Correction", "Null result", "Research record"];

function topicName(topic) {
  return topicLabels.get(topic) || topic;
}

export function archiveEntry(record, url, date) {
  const subject = topicName(record.topic);
  const summary = record.summary.trim();
  return `<a class="result-row index-entry" href="${esc(url)}">
    <div class="index-meta">${badge(record.status)}<span class="index-subject" title="${esc(record.topic)}">${esc(subject)}</span></div>
    <div class="row-content"><h3>${esc(record.title)}</h3>${summary && summary !== record.title ? `<p>${esc(summary)}</p>` : ""}</div>
    ${date ? `<time class="index-date" datetime="${esc(record.updatedAt)}">Updated ${esc(date)}</time>` : ""}
    <span class="index-open" aria-hidden="true">${icon("arrow")}</span></a>`;
}

// Each record occupies one position. Sorting positions by angle gives each
// disposition a continuous sector without assigning celestial coordinates.
function censusPositions(total) {
  const points = [];
  for (let i = 0; i < total; i++) {
    const fraction = (i + 0.5) / total;
    const radius = Math.sqrt(43 * 43 + fraction * (91 * 91 - 43 * 43));
    const angle = (i * Math.PI * (3 - Math.sqrt(5))) % (Math.PI * 2);
    points.push({
      angle,
      x: (100 + Math.sin(angle) * radius).toFixed(2),
      y: (100 - Math.cos(angle) * radius).toFixed(2),
    });
  }
  return points.sort((a, b) => a.angle - b.angle);
}

function censusChart(records, statuses) {
  const ordered = [...records].sort(
    (a, b) => statuses.indexOf(a.status) - statuses.indexOf(b.status) || a.id.localeCompare(b.id),
  );
  const positions = censusPositions(ordered.length);
  const markRadius = Math.min(2.35, 37 / Math.sqrt(Math.max(ordered.length, 1))).toFixed(2);
  const dots = ordered
    .map(
      (r, i) =>
        `<circle cx="${positions[i].x}" cy="${positions[i].y}" r="${markRadius}" data-record="${esc(r.id)}" data-disposition="${esc(r.status)}"/>`,
    )
    .join("");
  return `<svg viewBox="0 0 200 200" aria-hidden="true" focusable="false" class="census-plot">
    <circle class="census-orbit" cx="100" cy="100" r="97"/>
    <circle class="census-orbit" cx="100" cy="100" r="37"/>
    <path class="census-orbit" d="M100 1v5m0 188v5M1 100h5m188 0h5"/>
    ${dots}</svg>`;
}

function countBy(rows, key) {
  const counts = new Map();
  for (const row of rows) counts.set(row[key], (counts.get(row[key]) || 0) + 1);
  return counts;
}

export function prepareArchiveIndex(records, onChoice) {
  const statuses = [...new Set([...dispositions, ...records.map((r) => r.status)])].filter((status) =>
    records.some((r) => r.status === status),
  );
  const topics = [...new Set([...topicLabels.keys(), ...records.map((r) => r.topic)])].filter((topic) =>
    records.some((r) => r.topic === topic),
  );
  const census = document.querySelector("#archive-census");
  const topicsNode = document.querySelector("#archive-topics");
  census.querySelector(".census-map").insertAdjacentHTML("afterbegin", censusChart(records, statuses));
  const countValue = census.querySelector(".census-total strong");
  const countLabel = census.querySelector(".census-total small");
  const dots = Array.from(census.querySelectorAll("circle[data-record]"));
  const key = census.querySelector(".census-key");
  key.innerHTML = statuses
    .map(
      (status) =>
        `<button type="button" data-disposition="${esc(status)}" aria-pressed="false" title="Filter by ${esc(status)}"><span class="census-swatch"></span><span>${esc(status)}</span><b></b></button>`,
    )
    .join("");
  topicsNode.innerHTML = ["", ...topics]
    .map(
      (topic) =>
        `<button type="button" data-topic="${esc(topic)}" title="${esc(topic || "All topics")}" aria-pressed="false"><span>${esc(topic ? topicName(topic) : "All topics")}</span><b></b></button>`,
    )
    .join("");
  const statusButtons = [...key.querySelectorAll("button")];
  const topicButtons = [...topicsNode.querySelectorAll("button")];
  const plot = census.querySelector(".census-plot");
  function preview(status) {
    for (const dot of dots)
      dot.classList.toggle("is-muted", Boolean(status && dot.dataset.disposition !== status));
    for (const button of statusButtons)
      button.classList.toggle("is-previewed", button.dataset.disposition === status);
  }
  function chooseStatus(status) {
    const selected = statusButtons.find((button) => button.dataset.disposition === status);
    if (selected) onChoice("status", selected.getAttribute("aria-pressed") === "true" ? "" : status);
  }
  key.addEventListener("focusin", (event) => preview(event.target.closest("button")?.dataset.disposition));
  key.addEventListener("focusout", () => preview(null));
  key.addEventListener("pointerover", (event) =>
    preview(event.target.closest("button")?.dataset.disposition),
  );
  key.addEventListener("pointerleave", () => preview(null));
  key.addEventListener("click", (event) => {
    const button = event.target.closest("button");
    if (button) chooseStatus(button.dataset.disposition);
  });

  // Sector boundaries bisect adjacent marks at each change of disposition.
  // The entire sector responds to touch. Keyboard visitors can use the legend.
  const sectors = [];
  const turn = 2 * Math.PI;
  let lastAngle = 0;
  for (const dot of dots) {
    const angle = (Math.atan2(dot.cx.baseVal.value - 100, 100 - dot.cy.baseVal.value) + turn) % turn;
    const previous = sectors.at(-1);
    if (previous?.status !== dot.dataset.disposition) {
      if (previous) previous.end = (lastAngle + angle) / 2;
      sectors.push({ status: dot.dataset.disposition, end: turn });
    }
    lastAngle = angle;
  }
  function dispositionAt(event) {
    const bounds = plot.getBoundingClientRect();
    const x = (200 * (event.clientX - bounds.left)) / bounds.width - 100;
    const y = 100 - (200 * (event.clientY - bounds.top)) / bounds.height;
    const radius = Math.hypot(x, y);
    if (radius < 39 || radius > 96) return null;
    const angle = (Math.atan2(x, y) + turn) % turn;
    return event.target.dataset.disposition || sectors.find((sector) => angle < sector.end)?.status;
  }
  let hovered = null;
  plot.addEventListener("pointermove", (event) => {
    if (event.pointerType === "touch") return;
    const next = dispositionAt(event);
    if (next !== hovered) {
      hovered = next;
      preview(next);
      plot.classList.toggle("can-filter", Boolean(next));
    }
  });
  plot.addEventListener("pointerleave", () => {
    hovered = null;
    plot.classList.remove("can-filter");
    preview(null);
  });
  plot.addEventListener("click", (event) => chooseStatus(dispositionAt(event)));
  topicsNode.addEventListener("click", (event) => {
    const button = event.target.closest("button");
    if (button) onChoice("topic", button.dataset.topic);
  });
  return (matchingText, found, filters) => {
    const matchingTopic = filters.topic
      ? matchingText.filter((r) => r.topic === filters.topic)
      : matchingText;
    const matchingStatus = filters.status
      ? matchingText.filter((r) => r.status === filters.status)
      : matchingText;
    const statusCounts = countBy(matchingTopic, "status");
    const topicCounts = countBy(matchingStatus, "topic");
    const included = new Set(found.map((r) => r.id));
    for (const dot of dots) dot.classList.toggle("is-filtered", !included.has(dot.dataset.record));
    countValue.textContent = found.length.toLocaleString("en-US");
    countLabel.textContent = filters.status || filters.topic || filters.q ? "in view" : "records";
    for (const button of statusButtons) {
      const status = button.dataset.disposition;
      button.querySelector("b").textContent = statusCounts.get(status) || 0;
      button.setAttribute("aria-pressed", String(status === filters.status));
    }
    for (const button of topicButtons) {
      const topic = button.dataset.topic;
      button.querySelector("b").textContent = topic ? topicCounts.get(topic) || 0 : matchingStatus.length;
      button.setAttribute("aria-pressed", String(topic === filters.topic));
    }
  };
}

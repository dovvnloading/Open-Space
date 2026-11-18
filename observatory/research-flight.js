import { esc } from "./html.js";
import { icon } from "./icons.js";
import { destinationHash, flightImage, makeFlightNavigator } from "./flight-navigator.js";
import { mountResearchScene } from "./research-scene.js";
import { mountResearchInstruments } from "./research-instruments.js";

function flightCard(d, plates) {
  const src = flightImage(d, plates);
  return `<button type="button" class="rf-card" data-flight="${esc(d.id)}"><span class="rf-card-image">${src ? `<img src="${esc(src)}" alt="" loading="lazy">` : icon("book")}</span>
    <span><small>${esc(d.setting || d.kind)}</small><strong>${esc(d.name)}</strong></span></button>`;
}

export function mountResearchFlight(container, { catalog, plates = {}, requested, onClose, onGenerated }) {
  const previousFocus = document.activeElement;
  const dialog = document.createElement("dialog");
  dialog.className = "os-universe rf-universe";
  dialog.setAttribute("aria-label", "Universe simulator");
  const items = catalog?.destinations || [];
  const targets = new Map((catalog?.spotlights || []).map((t) => [t.id, t]));
  const records = new Map((catalog?.records || []).map((r) => [r.id, r]));
  const byId = new Map(items.map((d) => [d.id, d]));
  const preferred = ["mars", "sun", "quasar", "dwarf"];
  const featured = [
    ...preferred.map((id) => byId.get(id)).filter(Boolean),
    ...items.filter((d) => d.featured && !preferred.includes(d.id)),
  ];
  dialog.innerHTML = `<header class="rf-bar"><div class="rf-brand">${icon("orbit")}<span>Universe simulator</span></div>
    <div class="rf-modes" role="group" aria-label="Simulator mode"><span aria-current="true">Research<em> flights</em></span><button type="button" aria-label="Generated universe" data-generated>Generated<em> universe</em></button></div>
    <div class="rf-bar-actions"><button type="button" class="rf-route" data-finder aria-haspopup="dialog">${icon("search")}<span>Destinations</span></button>
    <button class="rf-icon" data-close type="button" aria-label="Close universe simulator">${icon("x")}</button></div></header>
    <div class="rf-mobile-switch" role="group" aria-label="Flight view"><button data-view="field" type="button" aria-pressed="true">${icon("focus")}Field</button><button data-view="instruments" aria-pressed="false" type="button">${icon("pulse")}Instruments</button></div>
    <div class="rf-layout"><div class="rf-port"><section class="rf-stage" aria-label="Destination field">
    <div class="rf-scene"></div><div class="rf-intro"><p class="rf-overline rf-location"></p><h1 class="rf-destination-name" tabindex="-1"></h1><p class="rf-ident"></p></div>
    <div class="rf-message" hidden></div></section>
    <nav class="rf-itinerary" aria-label="Featured flights"><div class="rf-tray-head"><h2>Featured flights</h2><div><button type="button" class="rf-share" data-share>Copy flight link</button><input type="text" class="rf-share-url" aria-label="Flight link" readonly hidden><button type="button" class="rf-icon" data-prev aria-label="Scroll featured flights back">${icon("left")}</button><button type="button" class="rf-icon" data-next aria-label="Scroll featured flights forward">${icon("right")}</button></div></div><div class="rf-strip">${featured.map((d) => flightCard(d, plates)).join("")}</div></nav>
    </div><aside class="rf-rail" aria-label="Research instruments"></aside></div><div class="visually-hidden" aria-live="polite" data-arrival></div>`;
  container.append(dialog);
  const one = (q) => dialog.querySelector(q);
  const rail = one(".rf-rail"),
    scene = one(".rf-scene"),
    strip = one(".rf-strip");
  let current,
    sequence = 0,
    clearScene,
    clearInstruments;
  const share = one("[data-share]"),
    shareUrl = one(".rf-share-url");
  function showView(which) {
    if (which === "field" && matchMedia("(max-width: 760px)").matches) rail.querySelector("image-desk")?.stop();
    dialog.dataset.view = which;
    dialog
      .querySelectorAll(".rf-mobile-switch button")
      .forEach((button) => button.setAttribute("aria-pressed", String(button.dataset.view === which)));
  }
  function stop() {
    sequence++;
    clearScene?.();
    clearInstruments?.();
    clearScene = clearInstruments = undefined;
    scene.replaceChildren();
    rail.replaceChildren();
  }
  function travel(id, focus = true) {
    if (current && current.id === id) {
      showView("field");
      if (focus) one("h1").focus({ preventScroll: true });
      return;
    }
    stop();
    share.textContent = "Copy flight link";
    shareUrl.hidden = true;
    const thisFlight = sequence;
    current = byId.get(id);
    one(".rf-message").hidden = !!current;
    one(".rf-intro").hidden = !current;
    one('.rf-mobile-switch [data-view="instruments"]').disabled = !current;
    share.disabled = !current;
    showView("field");
    strip
      .querySelectorAll("[data-flight]")
      .forEach((b) => b.setAttribute("aria-current", String(b.dataset.flight === id)));
    if (!current) {
      const available = items.length > 0;
      one(".rf-message").innerHTML =
        `<h1>${available ? "Destination unavailable" : "Flight catalog unavailable"}</h1>
        <p>${available ? "This research address does not appear in the current repository snapshot." : "Reload the catalog to visit an investigation."}</p>
        <button class="rf-route" type="button" ${available ? "data-finder" : "data-reload"}>${available ? "Choose a destination" : "Reload catalog"}</button>`;
      return;
    }
    const d = current,
      target = targets.get(id),
      record = records.get(d.recordId),
      plate = plates[id];
    const hash = destinationHash(id);
    if (location.hash !== hash) history.replaceState({ ...history.state, researchTarget: id }, "", hash);
    document.title = d.name + " · Research flight · Open Space Observatory";
    one(".rf-location").textContent = d.setting || record?.topic || "Research flight";
    one(".rf-destination-name").textContent = d.name;
    one(".rf-ident").textContent = target
      ? [target.instrument, target.epoch].filter(Boolean).join(" · ")
      : d.image?.fieldLabel || d.coordinates?.label || "";
    one("[data-arrival]").textContent = "";
    function inspectPlate({ marker } = {}) {
      if (sequence !== thisFlight) return;
      showView("instruments");
      const tab = rail.querySelector('[data-rfi-tab="images"]') || rail.querySelector('[aria-selected="true"]');
      tab?.click();
      tab?.focus({ preventScroll: true });
      if (marker) rail.querySelector("image-desk")?.focusMarker(marker);
      rail.scrollTop = 0;
    }
    clearScene = mountResearchScene(scene, {
      destination: d,
      plate,
      animate: true,
      onInspect: inspectPlate,
      onArrive: () => {
        if (sequence === thisFlight) one("[data-arrival]").textContent = "Arrived at " + d.name;
      },
    });
    clearInstruments = mountResearchInstruments(rail, { destination: d, record, target, plate, catalog });
    rail.scrollTop = 0;
    if (focus) one(".rf-destination-name").focus({ preventScroll: true });
  }
  const finder = makeFlightNavigator(dialog, { items, records, targets, plates, onSelect: travel });
  function scrollCards(direction) {
    strip.scrollBy({
      left: direction * strip.clientWidth * 0.82,
      behavior: matchMedia("(prefers-reduced-motion: reduce)").matches ? "instant" : "smooth",
    });
  }
  function syncArrows() {
    one("[data-prev]").disabled = strip.scrollLeft < 2;
    one("[data-next]").disabled = strip.scrollLeft + strip.clientWidth >= strip.scrollWidth - 2;
  }
  strip.addEventListener("scroll", syncArrows, { passive: true });
  const resize = new ResizeObserver(syncArrows);
  resize.observe(strip);
  dialog.addEventListener("click", (e) => {
    const b = e.target.closest("button");
    if (!b || b.closest(".rf-finder")) return;
    if (b.hasAttribute("data-finder")) finder.open(b, current?.id);
    else if (b.hasAttribute("data-close")) onClose();
    else if (b.hasAttribute("data-generated")) onGenerated(current?.id);
    else if (b.hasAttribute("data-share")) {
      const url = location.href,
        sharedSequence = sequence;
      const fallback = () => {
        if (sequence !== sharedSequence) return;
        shareUrl.value = url;
        shareUrl.hidden = false;
        shareUrl.focus();
        shareUrl.select();
      };
      if (navigator.clipboard?.writeText)
        navigator.clipboard.writeText(url).then(() => {
          if (sequence === sharedSequence) share.textContent = "Link copied";
        }, fallback);
      else fallback();
    } else if (b.hasAttribute("data-reload")) location.reload();
    else if (b.hasAttribute("data-prev")) scrollCards(-1);
    else if (b.hasAttribute("data-next")) scrollCards(1);
    else if (b.hasAttribute("data-flight")) travel(b.dataset.flight);
    else if (b.matches(".rf-mobile-switch button")) showView(b.dataset.view);
  });
  dialog.addEventListener("cancel", (e) => {
    if (e.target === dialog) {
      e.preventDefault();
      onClose();
    }
  });
  one("[data-finder]").disabled = !items.length;
  one(".rf-itinerary").hidden = !featured.length;
  dialog.showModal();
  travel(requested ?? (byId.has("mars") ? "mars" : items[0]?.id), false);
  if (current) one(".rf-destination-name").focus({ preventScroll: true });
  return () => {
    stop();
    resize.disconnect();
    finder.dispose();
    if (dialog.open) dialog.close();
    dialog.remove();
    if (previousFocus?.isConnected) previousFocus.focus({ preventScroll: true });
  };
}

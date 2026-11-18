const shapes = {
  deck: '<rect x="3" y="3" width="7" height="7" rx="1.4"/><rect x="14" y="3" width="7" height="7" rx="1.4"/><rect x="3" y="14" width="7" height="7" rx="1.4"/><rect x="14" y="14" width="7" height="7" rx="1.4"/>',
  orbit:
    '<circle cx="12" cy="12" r="5"/><ellipse cx="12" cy="12" rx="11" ry="3.6" transform="rotate(-35 12 12)"/>',
  archive: '<rect x="3" y="4" width="18" height="5" rx="1.4"/><path d="M5 9v11h14V9m-10 4h6"/>',
  external: '<path d="M7 17 17 7M7 7h10v10"/>',
  github:
    '<path fill="currentColor" stroke="none" d="M12 .9a11.1 11.1 0 0 0-3.51 21.63c.55.1.76-.24.76-.53l-.02-2.07c-3.08.67-3.73-1.31-3.73-1.31-.5-1.28-1.23-1.62-1.23-1.62-1-.69.08-.68.08-.68 1.1.08 1.7 1.14 1.7 1.14.99 1.69 2.59 1.2 3.22.92.1-.72.39-1.21.7-1.49-2.46-.28-5.05-1.23-5.05-5.48 0-1.21.43-2.2 1.14-2.98-.12-.28-.49-1.4.11-2.94 0 0 .94-.3 3.06 1.14a10.6 10.6 0 0 1 5.56 0c2.12-1.44 3.05-1.14 3.05-1.14.61 1.54.24 2.66.12 2.94.71.78 1.13 1.77 1.13 2.98 0 4.26-2.6 5.2-5.08 5.48.4.35.75 1.02.75 2.06l-.01 3.05c0 .3.2.64.76.53A11.11 11.11 0 0 0 12 .9Z"/>',
  clock: '<circle cx="12" cy="12" r="8.8"/><path d="M12 6.5V12l3.4 2"/>',
  focus: '<path d="M8 3H3v5m13-5h5v5M3 16v5h5m13-5v5h-5"/><circle cx="12" cy="12" r="3.4"/>',
  book: '<path d="M4 4h12a3 3 0 0 1 3 3v13H6a3 3 0 0 1-3-3V5a1 1 0 0 1 1-1Zm-1 13a3 3 0 0 1 3-3h13M7 8h6"/>',
  antenna:
    '<path d="m5 5 14 14a9.8 9.8 0 0 1-14-14Zm7 7 6-6M9 16l-2 5h10l-2-3M17 3a4.5 4.5 0 0 1 4 4"/><circle cx="18" cy="6" r="1"/>',
  calendar: '<rect x="3.5" y="5.5" width="17" height="15" rx="2"/><path d="M7 3v5m10-5v5M4 11h16"/>',
  pulse: '<path d="M2 13h5l3-8 4 14 3-8 2 2h3"/>',
  left: '<path d="m14.5 5.5-6.5 6.5 6.5 6.5"/>',
  right: '<path d="m9.5 5.5 6.5 6.5-6.5 6.5"/>',
  arrow: '<path d="M3 12h17m-6-6 6 6-6 6"/>',
  search: '<circle cx="10.5" cy="10.5" r="6.7"/><path d="m16 16 5 5"/>',
  x: '<path d="m6 6 12 12M18 6 6 18"/>',
};
export function icon(name) {
  return (
    '<svg class="icon" aria-hidden="true" xmlns="http://www.w3.org/2000/svg" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="1.6" stroke-linecap="round" stroke-linejoin="round">' +
    shapes[name] +
    "</svg>"
  );
}
export function hydrateIcons(root = document) {
  for (const placeholder of root.querySelectorAll("[data-icon]")) {
    const template = document.createElement("template");
    template.innerHTML = icon(placeholder.dataset.icon);
    const svg = template.content.firstElementChild;
    svg.classList.add(...placeholder.classList);
    placeholder.replaceWith(svg);
  }
}

"""Exercise a built observatory at a project subpath, or check a deployed URL."""

import argparse
import re
from contextlib import contextmanager
from functools import partial
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from threading import Thread

from playwright.sync_api import expect, sync_playwright


@contextmanager
def serve(directory):
    class QuietHandler(SimpleHTTPRequestHandler):
        def log_message(self, *args):
            pass

    directory = directory.resolve()
    with ThreadingHTTPServer(("127.0.0.1", 0), partial(QuietHandler, directory=directory.parent)) as server:
        thread = Thread(target=server.serve_forever, daemon=True)
        thread.start()
        try:
            yield f"http://127.0.0.1:{server.server_port}/{directory.name}/"
        finally:
            server.shutdown()
            thread.join()


def check_desktop(page, url):
    page.goto(url)
    expect(page.locator(".field-card")).to_have_count(3)
    expect(page.locator("#hero-name")).to_have_text("Mars")
    expect(page.locator("#recent-list .result-row")).to_have_count(4)
    recent = page.locator("#recent-list").bounding_box()
    assert recent["width"] > 250 and recent["height"] > 100, "Recent research is clipped"
    page.get_by_role("button", name="Measured data").click()
    expect(page.locator("#plate-figure")).to_be_visible()
    page.wait_for_function("document.querySelector('#plate-figure-image').naturalWidth > 0")
    page.get_by_role("button", name="Next featured target").click()
    expect(page.locator("#hero-name")).to_have_text("The Sun")
    expect(page.locator("#plate-figure-image")).to_have_attribute("src", "./data/figures/sun.svg")
    page.get_by_role("button", name="Context image").click()
    expect(page.locator("#plate-credit")).to_contain_text("Coronagraph")
    page.get_by_role("link", name="View image sequence").click()
    clock = page.get_by_role("slider", name="Solar frame time")
    clock.focus()
    clock.press("Home")
    expect(clock).to_have_attribute("aria-valuetext", "15:08:45 UTC")
    clock.press("ArrowRight")
    expect(clock).to_have_attribute("aria-valuetext", "16:08:45 UTC")
    page.get_by_role("button", name="Play sequence").click()
    expect(clock).not_to_have_value("1")
    page.get_by_role("button", name="Pause sequence").click()
    # An unreadable frame must never put its UTC stamp on the previous exposure.
    blocked_frame = "**/cor2/20210713_220845_n4c2A.jpg"
    previous_time = page.locator("[data-clock]").inner_text()
    page.route(blocked_frame, lambda route: route.abort())
    clock.focus()
    clock.press("End")
    retry = page.get_by_role("button", name="Retry image")
    expect(retry).to_be_visible()
    expect(page.locator("[data-clock]")).to_have_text(previous_time)
    page.unroute(blocked_frame)
    retry.click()
    expect(page.locator("[data-clock]")).to_have_text("22:08:45 UTC")
    assert "220845" in page.locator(".desk-sheet img").get_attribute("src")
    page.locator('.desk-plates a[href="#plate/quasar"]').click()
    expect(page.locator(".desk-pin")).to_have_count(2)
    zoom = page.get_by_role("button",name="Magnify 4 times")
    zoom.click()
    expect(zoom).to_have_attribute("aria-pressed","true")
    page.get_by_role("link",name="Read investigation",exact=True).click()
    expect(page.locator("#dossier-title")).to_be_focused()
    page.keyboard.press("Escape")
    expect(page.locator(".desk-pin")).to_have_count(2)
    expect(page.locator("#dossier-title")).to_be_focused()
    page.keyboard.press("Escape")
    expect(page.locator("#inspect-plate")).to_be_focused()
    page.get_by_role("button", name="Previous featured target").click()
    page.locator("#hero-detail").click()
    expect(page.get_by_role("dialog")).to_be_visible()
    expect(page.locator(".measure")).to_have_count(3)
    page.get_by_text("How these values were checked", exact=True).click()
    expect(page.locator("#dossier .methods[open] ol li")).to_have_count(3)
    for row in page.locator("#dossier .methods li").all():
        assert "/blob/" in row.locator("a").first.get_attribute("href")
    page.get_by_text("Inspect observation rows", exact=False).click()
    table = page.locator(".table-reader")
    expect(table.locator("tbody tr:visible")).to_have_count(8)
    previous = table.get_by_role("button", name="Previous observation rows")
    following = table.get_by_role("button", name="Next observation rows")
    expect(previous).to_be_disabled()
    following.click()
    expect(table.get_by_role("status")).to_have_text("Rows 9–16 of 31")
    following.click()
    following.click()
    expect(table.locator("tbody tr:visible")).to_have_count(7)
    expect(following).to_be_disabled()
    previous.click()
    expect(table.get_by_role("status")).to_have_text("Rows 17–24 of 31")
    cursor = page.get_by_role("slider",name="Inspect observation")
    cursor.focus()
    cursor.press("End")
    expect(cursor).to_have_value("31")
    page.get_by_role("button",name="Published figure",exact=True).click()
    page.locator("figure img").scroll_into_view_if_needed()
    page.wait_for_function("document.querySelector('figure img')?.naturalWidth > 0")
    assert "/blob/" in page.get_by_role("link", name="Investigation report").get_attribute("href")
    page.keyboard.press("Escape")
    expect(page.locator("#dossier")).not_to_be_visible()
    page.locator(".snapshot-link").click()
    expect(page.locator("#dossier-title")).to_have_text("How the archive was checked")
    expect(page.locator(".audit-ledger")).to_contain_text("Declared calculations")
    page.get_by_role("button", name="Close detail").click()

    page.get_by_role("link", name="Research archive", exact=True).click()
    total = page.locator(".census-plot [data-record]").count()
    expect(page.locator(".census-total strong")).to_have_text(str(total))
    prediction = page.locator('.census-key [data-disposition="Prediction"]')
    prediction.click()
    expect(page.locator("#status-filter")).to_have_value("Prediction")
    assert page.locator(".census-plot [data-record]:not(.is-filtered)").count() == int(prediction.locator("b").inner_text())
    page.locator("#clear-archive").click()
    page.get_by_role("button", name="Next", exact=True).click()
    expect(page.locator("#pagination")).to_contain_text("Page 2")
    page.locator("#query").fill("no-record-will-match-this")
    expect(page.locator("#results-count")).to_contain_text("0 records")
    page.get_by_role("button", name="Clear filters").click()
    page.locator("#status-filter").select_option("Null result")
    first = page.locator("#result-list .result-row").first
    title = first.locator("h3").inner_text()
    first.click()
    expect(page.locator("#dossier-title")).to_have_text(title)
    page.get_by_role("button", name="Close detail").click()
    expect(page.locator("#status-filter")).to_have_value("Null result")
    expect(first).to_be_focused()
    page.locator(".skip").focus()
    page.keyboard.press("Enter")
    expect(page.locator("#archive")).to_be_visible()
    expect(page.locator("#main")).to_be_focused()

    page.get_by_role("link", name="Universe simulator", exact=True).click()
    page.get_by_role("button", name="Generated universe", exact=True).click()
    expect(page.locator(".osu-levels")).to_be_visible()
    page.get_by_role("button", name="Atlas", exact=True).click()
    expect(page.locator(".osu-place")).to_have_count(4)
    page.locator('[data-atlas-index="3"]').click()
    expect(page.locator('[data-ref="name"]')).to_have_text("Planet V")
    page.locator('[data-world="0"]').click()
    expect(page.locator('[data-ref="name"]')).to_have_text("Planet I")
    for level in ("sector", "galaxy", "system", "world"):
        button = page.locator(f'[data-level="{level}"]')
        button.click()
        expect(button).to_have_attribute("aria-current", "step")
        if level == "system":
            page.get_by_role("button", name="AU scale").click()
            expect(page.locator('[data-ref="scale"]')).to_contain_text("AU")
            page.get_by_role("button", name="Advance one quarter of the selected orbit").click()
    page.get_by_role("button", name="Close universe simulator").click()
    expect(page.locator("#archive")).to_be_visible()
    page.goto("about:blank")
    page.goto(url + "#target/dwarf")
    expect(page.locator("#dossier-title")).to_be_visible()
    page.get_by_role("button", name="Close detail").click()
    expect(page.locator("#deck")).to_be_visible()
    print("Desktop: carousel, measurements, figure bytes, archive filters, pagination, keyboard focus, "
          "four simulation scales and direct links passed.")


def check_published_records(page):
    """Open every record route and compare its title, disposition and source links."""
    counts = page.evaluate("""
        async () => {
          const catalog = await (await fetch('./data/catalog.json')).json();
          async function navigate(hash) {
            // Browsers rate-limit rapid navigations within the same document.
            await new Promise(resolve => setTimeout(resolve, 150));
            const changed = new Promise(resolve => addEventListener('hashchange', resolve, {once:true}));
            location.hash = hash;
            await changed;
          }
          const targets = new Map(catalog.spotlights.map(t => [t.recordId,t]));
          let observationRows = 0;
          for (const row of catalog.records) {
            const target = targets.get(row.id);
            await navigate(target ? '#target/' + target.id : '#record/' + row.id);
            const dialog = document.querySelector('#dossier[open]');
            if (dialog?.querySelector('h2')?.textContent !== (target || row).title ||
                dialog.querySelector('.status')?.textContent !== row.status)
              throw Error('Displayed record differs from catalog: ' + row.id);
            const links = new Set([...dialog.querySelectorAll('a')].map(a => a.href));
            for (const item of [...row.links, ...(target?.links || []), row.provenance]) {
              if (!links.has(item.url) || !item.url.includes('/' + catalog.revision + '/'))
                throw Error('Missing pinned evidence: ' + row.id + ' / ' + item.url);
            }
            if (target && dialog.querySelectorAll('.methods ol li a').length < target.measurements.length)
              throw Error('Missing readout methods: ' + target.id);
            if (target) {
              const table = dialog.querySelector('.table-reader table');
              const values = [...table.tBodies[0].rows].map(tr => [...tr.cells].map(cell => cell.textContent));
              const headers = [...table.tHead.rows[0].cells].map(cell => cell.textContent);
              if (table.caption.textContent !== target.table.caption ||
                  JSON.stringify(headers) !== JSON.stringify(target.table.columns) ||
                  JSON.stringify(values) !== JSON.stringify(target.table.rows.map(row => row.map(String))) ||
                  !links.has(target.table.source.url) || !target.table.source.url.includes('/' + catalog.revision + '/'))
                throw Error('Displayed observation table differs from its source: ' + target.id);
              observationRows += values.length;
            }
          }
          return [catalog.records.length, catalog.spotlights.length, observationRows];
        }
    """)
    print(f"Compared all {counts[0]} rendered records, {counts[1]} featured dossiers and {counts[2]} observation rows with the snapshot.")


def check_universe_links(page, url, engine):
    ring = url + "#universe?view=world&at=G001-S00229&planet=5"
    archive = url + "#archive?q=MAVEN"
    name = page.locator('.os-universe [data-ref="name"]')
    atlas = page.get_by_role("button", name="Atlas", exact=True)
    page.goto(archive)
    page.get_by_role("link", name="Universe simulator", exact=True).click()
    page.get_by_role("button", name="Generated universe", exact=True).click()
    atlas.click()
    page.locator('[data-atlas-index="3"]').click()
    expect(name).to_have_text("Planet V")
    expect(page).to_have_url(ring)
    atlas.click()
    if engine == "chromium":
        page.context.grant_permissions(["clipboard-read", "clipboard-write"])
        page.get_by_role("button", name="Copy location link").click()
        expect(page.locator('[data-ref="shareStatus"]')).to_have_text("Location link copied.")
        assert page.evaluate("navigator.clipboard.readText()") == ring
    # Exercise the manual fallback in browsers which block clipboard access.
    page.evaluate("Object.defineProperty(navigator, 'clipboard', {configurable: true, value: undefined})")
    page.get_by_role("button", name="Copy location link").click()
    field = page.get_by_role("textbox", name="Location link to copy")
    expect(field).to_have_value(ring)
    expect(field).to_be_focused()
    assert field.evaluate("e => e.selectionStart === 0 && e.selectionEnd === e.value.length")
    page.go_back()
    expect(page).to_have_url(archive)
    expect(page.locator(".os-universe")).not_to_be_visible()
    page.go_forward()
    expect(page).to_have_url(ring)
    expect(name).to_have_text("Planet V")
    page.get_by_role("button", name="Close universe simulator").click()
    expect(page).to_have_url(archive)

    page.goto(url + "#target/quasar")
    expect(page.locator("#dossier")).to_be_visible()
    page.goto(ring)
    expect(name).to_have_text("Planet V")
    expect(page.locator('.osu-worlds [data-world="4"]')).to_have_attribute("aria-pressed", "true")
    page.locator('.osu-worlds [data-world="0"]').click()
    expect(name).to_have_text("Planet I")
    page.go_back()
    expect(page.locator("#dossier")).to_be_visible()
    expect(page).to_have_url(url + "#target/quasar")
    page.goto(ring)
    expect(name).to_have_text("Planet V")
    page.reload()
    expect(name).to_have_text("Planet V")
    page.goto(url + "#universe?view=world&at=G128-S32768&planet=8")
    expect(page.locator('[data-ref="locationError"]')).to_be_visible()
    expect(name).to_have_text("Vesper Spiral")
    print(f"{engine}: simulated destinations, link copy, direct/reload and history passed.")


def check_research_flights(page, url):
    page.goto(url + "#target/quasar")
    page.get_by_role("link", name="Fly to investigation").click()
    flight = page.locator(".rf-universe")
    expect(page).to_have_url(url + "#universe?target=quasar")
    expect(flight.locator("h1:visible")).to_have_text("J0826+7002")
    expect(flight.locator('[data-rfi-tab="data"]')).to_have_attribute("aria-selected", "true")
    flight.get_by_role("button", name="Inspect lensed image B", exact=True).click()
    expect(flight.locator(".desk-pin-selected")).to_have_attribute("data-pin", "B")
    expect(flight.locator('[data-zoom="2"]')).to_have_attribute("aria-pressed", "true")
    expect(flight.locator(".desk-pin")).to_have_count(2)
    flight.get_by_role("button", name="Destinations", exact=True).click()
    finder = page.get_by_role("dialog", name="Choose a destination")
    count = page.evaluate("async () => (await (await fetch('./data/catalog.json')).json()).records.length")
    expect(finder.locator(".rf-count")).to_have_text(f"{count} destinations")
    finder.get_by_role("button", name="Sky pointings", exact=True).click()
    pin = finder.locator('[data-sky-id="results-wr-117-2"]')
    pin.focus()
    expect(finder.locator("[data-map-label]")).to_contain_text("from J0826+7002")
    assert finder.locator("[data-course]").get_attribute("d").startswith("M")
    pin.click()
    expect(flight.locator("h1:visible")).to_have_text("WR 117-2")
    expect(page.locator(".rf-finder[open]")).to_have_count(0)
    flight.get_by_role("tab", name="Report", exact=True).click()
    expect(flight.get_by_text("Published WR 117-2 position", exact=True).first).to_be_visible()
    expect(flight.locator(".rfi-reading a").first).to_have_attribute("href", re.compile(r"/blob/"))
    flight.locator('[data-flight="sun"]').click()
    flight.get_by_role("tab", name="Imagery").click()
    slider = flight.get_by_role("slider", name="Solar frame time")
    slider.focus()
    slider.press("End")
    expect(flight.locator("[data-clock]")).to_have_text("22:08:45 UTC")
    flight.get_by_role("button", name="Generated universe", exact=True).click()
    page.get_by_role("button", name="Research flights", exact=True).click()
    expect(flight.locator("h1:visible")).to_have_text("The Sun")
    expect(page).to_have_url(url + "#universe?target=sun")
    # A native reload preserves the flight’s history entry in all three engines.
    with page.expect_navigation():
        page.evaluate("setTimeout(() => location.reload(), 0)")
    expect(flight.locator("h1:visible")).to_have_text("The Sun")
    page.go_back()
    expect(page.locator("#dossier-title")).to_be_visible()
    expect(page).to_have_url(url + "#target/quasar")
    page.goto(url + "#universe?target=unlisted-research-address")
    expect(flight.get_by_role("heading", name="Destination unavailable")).to_be_visible()
    flight.get_by_role("button", name="Choose a destination", exact=True).click()
    search = finder.get_by_role("searchbox")
    search.fill("nothingmatches-93852")
    expect(finder.locator(".rf-count")).to_have_text("0 destinations")
    search.fill("Enceladus")
    finder.locator('[data-go="results-cassini-uvis-enceladus-clock"]').click()
    expect(flight.locator("h1:visible")).to_have_text("Enceladus")
    expect(flight.locator("[data-arrival]")).to_have_text("Arrived at Enceladus")
    page.emulate_media(reduced_motion="no-preference")
    # Burst selections use the actual launch buttons and cancel older cameras.
    flight.evaluate("""f => {
        for (const id of ['mars','dwarf','sun']) f.querySelector('[data-flight="'+id+'"]').click();
    }""")
    expect(flight.locator(".research-scene")).to_have_count(1)
    expect(flight.locator("[data-arrival]")).to_have_text("Arrived at The Sun", timeout=8000)
    expect(flight.locator(".rs-source")).to_have_attribute("src", "./assets/sun-cor2.jpg")
    page.emulate_media(reduced_motion="reduce")
    print("Research flights: search, sky navigation, original imagery, source receipts, mode switching, deep links and history passed.")


def check_mobile(page, url):
    page.goto(url)
    expect(page.locator("#cards a")).to_have_count(3)
    assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
    page.get_by_role("button", name="Measured data").click()
    expect(page.locator("#plate-figure-link")).to_be_visible()
    page.get_by_role("button", name="Context image").click()
    page.get_by_role("link", name="Launch simulator").click()
    flight = page.locator(".rf-universe")
    expect(flight.locator(".rf-stage")).to_be_visible()
    assert flight.locator("[data-close]").evaluate("b => b.getBoundingClientRect().right <= innerWidth")
    flight.locator('.rf-mobile-switch button[data-view="instruments"]').click()
    expect(flight.get_by_role("slider", name="Inspect observation")).to_be_visible()
    assert flight.get_by_role("button", name="Interactive plot").get_attribute("aria-pressed") == "true"
    flight.get_by_role("tab", name="Imagery").click()
    expect(flight.get_by_role("button", name="Magnify 4 times")).to_be_visible()
    flight.locator('.rf-mobile-switch button[data-view="field"]').click()
    expect(flight.locator(".rf-stage")).to_be_visible()
    page.get_by_role("button", name="Generated universe", exact=True).click()
    page.get_by_role("button", name="System", exact=True).click()
    expect(page.get_by_role("button", name="Play orbital motion")).to_be_visible()
    page.get_by_role("button", name="World", exact=True).click()
    page.get_by_role("button", name="Close universe simulator").click()
    page.get_by_role("link", name="Research archive", exact=True).click()
    page.locator("#query").fill("MAVEN")
    expect(page.locator("#result-list .result-row").first).to_be_visible()
    assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
    page.goto(url + "#target/quasar")
    page.locator(".table-reader summary").click()
    cells = page.locator(".table-reader tbody tr:visible td")
    assert cells.evaluate_all("""cells => cells.every(cell => {
        const range = document.createRange();
        range.selectNodeContents(cell);
        return range.getClientRects().length === 1;
    })"""), "Mobile observation values wrap across lines"
    assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
    print("Mobile: navigation, archive search, observation table, reduced motion and layout passed.")


def run(url, engine, channel):
    failures = []
    with sync_playwright() as playwright:
        browser = getattr(playwright, engine).launch(channel=channel)
        for width, height, check in [(1440, 1000, check_desktop), (390, 844, check_mobile)]:
            page = browser.new_page(viewport={"width": width, "height": height}, reduced_motion="reduce")
            page.on("pageerror", lambda error: failures.append(str(error)))
            page.on("response", lambda response: failures.append(f"{response.status} {response.url}")
                    if response.status >= 400 and response.url.startswith(url) else None)
            check(page, url)
            if check == check_desktop:
                check_published_records(page)
                check_universe_links(page, url, engine)
                check_research_flights(page, url)
            page.close()
        browser.close()
    if failures:
        raise AssertionError("\n".join(failures))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--directory", type=Path)
    source.add_argument("--url")
    parser.add_argument("--engine", choices=("chromium", "firefox", "webkit"), default="chromium")
    parser.add_argument("--channel", help="Optional installed Chromium channel, e.g. chrome")
    args = parser.parse_args()
    if args.channel and args.engine != "chromium":
        parser.error("--channel is supported for Chromium only")
    if args.url:
        run(args.url.rstrip("/") + "/", args.engine, args.channel)
    else:
        with serve(args.directory) as local_url:
            run(local_url, args.engine, args.channel)

"""Exercise a built observatory at a project subpath, or check a deployed URL."""

import argparse
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
    page.get_by_role("button", name="Measured data").click()
    expect(page.locator("#plate-figure")).to_be_visible()
    page.wait_for_function("document.querySelector('#plate-figure-image').naturalWidth > 0")
    page.get_by_role("button", name="Next featured target").click()
    expect(page.locator("#hero-name")).to_have_text("The Sun")
    expect(page.locator("#plate-figure-image")).to_have_attribute("src", "./data/figures/sun.svg")
    page.get_by_role("button", name="Context image").click()
    expect(page.locator("#plate-credit")).to_contain_text("Coronagraph")
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
    expect(page.locator(".osu-levels")).to_be_visible()
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


def check_mobile(page, url):
    page.goto(url)
    expect(page.locator("#cards a")).to_have_count(3)
    assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
    page.get_by_role("button", name="Measured data").click()
    expect(page.locator("#plate-figure-link")).to_be_visible()
    page.get_by_role("button", name="Context image").click()
    page.get_by_role("link", name="Launch simulator").click()
    page.get_by_role("button", name="System", exact=True).click()
    expect(page.get_by_role("button", name="Play orbital motion")).to_be_visible()
    page.get_by_role("button", name="World", exact=True).click()
    page.get_by_role("button", name="Close universe simulator").click()
    page.get_by_role("link", name="Research archive", exact=True).click()
    page.locator("#query").fill("MAVEN")
    expect(page.locator("#result-list .result-row").first).to_be_visible()
    assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
    print("Mobile: navigation, archive search, reduced motion and layout passed.")


def run(url, channel):
    failures = []
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(channel=channel)
        for width, height, check in [(1440, 1000, check_desktop), (390, 844, check_mobile)]:
            page = browser.new_page(viewport={"width": width, "height": height}, reduced_motion="reduce")
            page.on("pageerror", lambda error: failures.append(str(error)))
            page.on("response", lambda response: failures.append(f"{response.status} {response.url}")
                    if response.status >= 400 and response.url.startswith(url) else None)
            check(page, url)
            if check == check_desktop:
                check_published_records(page)
            page.close()
        browser.close()
    if failures:
        raise AssertionError("\n".join(failures))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--directory", type=Path)
    source.add_argument("--url")
    parser.add_argument("--channel", help="Optional installed Chromium channel, e.g. chrome")
    args = parser.parse_args()
    if args.url:
        run(args.url.rstrip("/") + "/", args.channel)
    else:
        with serve(args.directory) as local_url:
            run(local_url, args.channel)

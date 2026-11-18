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
    page.get_by_role("button", name="Next featured target").click()
    expect(page.locator("#hero-name")).to_have_text("The Sun")
    page.get_by_role("button", name="Previous featured target").click()
    page.locator("#hero-detail").click()
    expect(page.get_by_role("dialog")).to_be_visible()
    expect(page.locator(".measure")).to_have_count(3)
    page.locator("figure img").scroll_into_view_if_needed()
    page.wait_for_function("document.querySelector('figure img')?.naturalWidth > 0")
    assert "/blob/" in page.get_by_role("link", name="Investigation report").get_attribute("href")
    page.keyboard.press("Escape")
    expect(page.locator("#dossier")).not_to_be_visible()

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
    page.get_by_role("button", name="Close universe simulator").click()
    expect(page.locator("#archive")).to_be_visible()
    page.goto("about:blank")
    page.goto(url + "#target/dwarf")
    expect(page.locator("#dossier-title")).to_be_visible()
    page.get_by_role("button", name="Close detail").click()
    expect(page.locator("#deck")).to_be_visible()
    print("Desktop: carousel, measurements, figure bytes, archive filters, pagination, keyboard focus, "
          "four simulation scales and direct links passed.")


def check_mobile(page, url):
    page.goto(url)
    expect(page.locator("#cards a")).to_have_count(3)
    assert page.evaluate("document.documentElement.scrollWidth <= innerWidth")
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

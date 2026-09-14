"""Render and check the general-public articles without rerunning numerical studies.

Run from the repository root:
  code/.venv/Scripts/python.exe code/render_readers.py
  code/.venv/Scripts/python.exe code/render_readers.py --study d4_mean_stress
"""
import argparse
import base64
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
from urllib.parse import unquote, urlsplit
from bs4 import BeautifulSoup

ROOT = Path(__file__).resolve().parents[1]
SERIES = ROOT / "code/reader_series"
WRITING = ROOT / "writing"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def quarto_path():
    candidates = [shutil.which("quarto"), ROOT / "code/.tools/bin/quarto.cmd",
                  ROOT / "code/.tools/bin/quarto"]
    return next(str(p) for p in candidates if p and Path(p).is_file())


def render(entry, quarto):
    article = WRITING / entry["folder"] / "article.qmd"
    for fmt in ("html", "latex"):
        result = subprocess.run([quarto, "render", str(article), "--to", fmt],
            cwd=ROOT, capture_output=True, text=True, encoding="utf-8")
        cache = ROOT / "code/.cache/reader-build"
        cache.mkdir(parents=True, exist_ok=True)
        (cache / f'{entry["folder"]}-{fmt}.log').write_text(
            result.stdout + result.stderr, encoding="utf-8")
        if result.returncode:
            raise RuntimeError(f'{entry["folder"]} {fmt}: {result.stderr[-2500:]}')
    print(f'Rendered {entry["folder"]}', flush=True)


def inspect(entry):
    folder = WRITING / entry["folder"]
    source = (folder / "article.qmd").read_text(encoding="utf-8")
    html = (folder / "article.html").read_text(encoding="utf-8")
    tex = (folder / "article.tex").read_text(encoding="utf-8")
    document = BeautifulSoup(html, "html.parser")
    main = document.find("main")
    if main is None:
        raise RuntimeError(f"Missing article body: {folder}")
    images = main.find_all("img")
    assert len(images) == 1 and images[0].get("src", "").startswith("data:image/")
    assert images[0].get("alt", "").strip()
    embedded = base64.b64decode(images[0]["src"].split(",", 1)[1])
    assert hashlib.sha256(embedded).hexdigest() == digest(folder/"figures/reader.png"), (
        f"Embedded reader figure is stale: {entry['folder']}")
    assert "\ufffd" not in source + main.get_text() + tex
    assert "{python}" not in html
    ids = {node["id"] for node in document.find_all(id=True)}
    links = []
    for anchor in document.find_all("a", href=True):
        href = urlsplit(anchor["href"])
        if href.scheme or href.netloc:
            continue
        if href.path:
            target = (folder / unquote(href.path)).resolve()
            assert target.exists(), f"Missing local link: {entry['folder']}: {href.path}"
            if href.fragment and target.suffix == ".html":
                destination = BeautifulSoup(target.read_text(encoding="utf-8"), "html.parser")
                assert destination.find(id=unquote(href.fragment)), (target, href.fragment)
            links.append(str(target.relative_to(ROOT)).replace("\\", "/"))
        elif href.fragment:
            assert unquote(href.fragment) in ids, href.fragment
    graphics = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", tex)
    assert len(graphics) == 1 and all((folder / name).exists() for name in graphics)
    # These are editorial guardrails, not a reading-age or scientific-validity test.
    for command in (r"\nabla", r"\partial", r"\otimes", r"\Lambda", r"\Pi", r"\asymp", r"\int", r"\sum"):
        assert command not in source, f"Advanced notation in public article: {entry['folder']}: {command}"
    headings = [h.get_text(" ", strip=True) for h in main.find_all("h2")]
    if entry["folder"] != "paper-map":
        assert "Where this fits in the bigger picture" in headings
        assert "../paper-map/article.html" in source
    if entry["technical"]:
        assert "technical-article.html" in source
        assert (folder / "technical-article.qmd").is_file()
    words = len(re.findall(r"\b[\w'-]+\b", main.get_text(" ", strip=True)))
    assert words >= 600, f"Article was reduced to an outline: {entry['folder']}"
    report = {
        "status":"passed", "audience":"General public; mostly high-school mathematics",
        "article_source_file":"article.qmd", "rendered_formats":["html","latex"],
        "source_sha256":digest(folder/"article.qmd"),
        "html_sha256":digest(folder/"article.html"), "latex_sha256":digest(folder/"article.tex"),
        "reader_png_sha256":digest(folder/"figures/reader.png"),
        "reader_svg_sha256":digest(folder/"figures/reader.svg"),
        "word_count_approximate":words, "embedded_figures":1,
        "figure_alt_text_present":True, "embedded_figure_matches_current_png":True,
        "local_links_checked":len(links),
        "latex_figure_path_resolves":True, "context_section_present":True,
        "advanced_notation_guardrail":"passed; not a readability score",
        "new_numerical_research":False, "browser_visual_review":False, "pdf_compiled":False,
        "technical_material_preserved":entry["technical"]}
    output = ROOT / "code" / entry["folder"]
    output.mkdir(exist_ok=True)
    (output / "reader-validation.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    return report


def main():
    manifest = json.loads((SERIES/"manifest.json").read_text(encoding="utf-8"))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--study", nargs="+", choices=[a["folder"] for a in manifest["articles"]])
    parser.add_argument("--skip-figures", action="store_true")
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    entries = [a for a in manifest["articles"] if not args.study or a["folder"] in args.study]
    start = time.perf_counter()
    if not args.skip_figures and not args.check_only:
        spec = importlib.util.spec_from_file_location("reader_figures", SERIES/"figures.py")
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        module.generate()
    if not args.check_only:
        quarto = quarto_path()
        with ThreadPoolExecutor(max_workers=3) as executor:
            futures = [executor.submit(render, entry, quarto) for entry in entries]
            for task in as_completed(futures):
                task.result()
    # Inspect links only after all selected destinations have been rendered.
    results = {entry["folder"]:inspect(entry) for entry in entries}
    report = {"status":"passed", "article_count":len(results),
              "all_series_articles":not bool(args.study),
              "elapsed_seconds":round(time.perf_counter()-start,3), "articles":results}
    target = SERIES / ("validation.json" if not args.study else "validation-last-selection.json")
    target.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(f"Passed: {len(results)} articles, HTML/LaTeX, figures, links, and editorial structure.", flush=True)


if __name__ == "__main__":
    main()

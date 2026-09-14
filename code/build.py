"""Validate and regenerate a deliverable in the isolated Python environment.

Usage (repository root): code/.venv/Scripts/python.exe code/build.py --study d2
Use --compute-only for the notebook/data/figures without Quarto publishing.
"""

import argparse
import hashlib
import importlib.metadata
import importlib
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
STUDIES = {
    "d1": ("d1_core_size", 3, "test_scaling.py"),
    "d2": ("d2_finite_energy", 4, "test_energy.py"),
    "d3": ("d3_residual", 3, "test_residual.py"),
    "d4": ("d4_mean_stress", 3, "test_stress.py"),
    "d5": ("d5_shearing_wave", 3, "test_shear.py"),
    "d6": ("d6_moment_correction", 3, "test_moments.py"),
    "d6b": ("d6_inner_profile", 3, "test_inner.py"),
    "d6c": ("d6_axis_pressure", 3, "test_pressure.py"),
    "d6d": ("d6_global_axis", 3, "test_axis_bounds.py"),
    "d7": ("d7_conditioning", 3, "test_conditioning.py"),
    "d9": ("d9_regime_map", 4, "test_regimes.py"),
    "d10": ("d10_compressible_vortex", 4, "test_compressible.py"),
}
sys.path.insert(0, str(CODE / "src"))


def run(command, **kwargs):
    return subprocess.run([str(x) for x in command], check=True, cwd=ROOT, **kwargs)


def find_quarto():
    on_path = shutil.which("quarto")
    candidates = [on_path, CODE / ".tools/bin/quarto.cmd", CODE / ".tools/bin/quarto"]
    for candidate in candidates:
        if candidate and Path(candidate).is_file():
            return str(candidate)
    raise RuntimeError("Quarto is needed for publishing. Install it from https://quarto.org/docs/download/ or use --compute-only.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compute-only", action="store_true")
    parser.add_argument("--study", choices=STUDIES, default="d1")
    args = parser.parse_args()
    import nbformat
    from nbclient import NotebookClient
    folder, expected_figures, test_pattern = STUDIES[args.study]
    STUDY = CODE / folder
    ARTICLE = ROOT / "writing" / folder
    # Public prose now has its own lightweight build; preserve the executable
    # computational article and its original scientific figure checks.
    article_stem = "technical-article" if (ARTICLE / "technical-article.qmd").exists() else "article"
    generate = importlib.import_module(f"nscomp.{args.study}").generate

    start = time.perf_counter()
    report_path = STUDY / "validation.json"
    report = {"status": "in_progress", "study": args.study, "python": platform.python_version(),
              "platform": platform.system(), "rendered_formats": []}
    report_path.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    try:
        tests = run([sys.executable, "-m", "unittest", "discover", "-s", CODE / "tests",
                     "-p", test_pattern, "-v"],
                    capture_output=True, text=True, encoding="utf-8")
        test_log = tests.stdout + tests.stderr
        (STUDY / "validation-tests.txt").write_text(test_log, encoding="utf-8")
        print(test_log)
        report["scientific_checks"] = "passed; see validation-tests.txt"
        if args.study == "d2":
            generate(include_animation=True)
            node = shutil.which("node")
            if node:
                result = run([node, STUDY / "check-player.cjs"], capture_output=True, text=True)
                print(result.stdout)
                report["animation_player_checks"] = "passed in Node DOM stub; not browser visual review"
            else:
                report["animation_player_checks"] = "not run: Node unavailable"
        else:
            generate()
        run([sys.executable, "-m", "ipykernel", "install", "--sys-prefix", "--name", "python3",
             "--display-name", "Python (Navier-Stokes companion)"], capture_output=True)
        # Use only this environment's kernelspec during fresh execution.
        os.environ["JUPYTER_PATH"] = str(Path(sys.prefix) / "share/jupyter")
        notebook_path = STUDY / "study.ipynb"
        notebook = nbformat.read(notebook_path, as_version=4)
        for cell in notebook.cells:
            if cell.cell_type == "code":
                cell.outputs = []
                cell.execution_count = None
        NotebookClient(notebook, timeout=120, kernel_name="python3",
                       resources={"metadata": {"path": str(STUDY)}}).execute()
        nbformat.validate(notebook)
        errors = [output for cell in notebook.cells if cell.cell_type == "code"
                  for output in cell.outputs if output.output_type == "error"]
        if errors:
            raise RuntimeError("Notebook contains error outputs")
        nbformat.write(notebook, notebook_path)
        report["notebook"] = {"fresh_kernel": True,
                              "executed_code_cells": sum(c.cell_type == "code" for c in notebook.cells),
                              "error_outputs": len(errors)}
        print("Notebook executed successfully from cleared outputs.")
        if not args.compute_only:
            quarto = find_quarto()
            report["quarto"] = subprocess.check_output([quarto, "--version"], text=True).strip()
            env = os.environ.copy()
            env["QUARTO_PYTHON"] = sys.executable
            env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")
            for output_format in ["html", "latex"]:
                completed = run([quarto, "render", ARTICLE / f"{article_stem}.qmd", "--to", output_format],
                                env=env, capture_output=True, text=True, encoding="utf-8")
                (CODE / ".cache" / f"{args.study}-quarto-{output_format}.log").parent.mkdir(parents=True, exist_ok=True)
                (CODE / ".cache" / f"{args.study}-quarto-{output_format}.log").write_text(
                    completed.stdout+completed.stderr, encoding="utf-8")
                report["rendered_formats"].append(output_format)
                print(f"Rendered {output_format}.")
            html = (ARTICLE / f"{article_stem}.html").read_text(encoding="utf-8")
            if "data:image/png;base64," not in html or "{python}" in html:
                raise RuntimeError("HTML is missing embedded figures or has unresolved inline code")
            report["html_checks"] = "figures embedded; inline Python evaluated"
            from bs4 import BeautifulSoup
            document = BeautifulSoup(html, "html.parser")
            images = document.find_all("img")
            if len(images) != expected_figures or any(not img.get("src", "").startswith("data:image/") for img in images):
                raise RuntimeError(f"Expected {expected_figures} embedded scientific figures")
            missing_links = []
            ids = {element["id"] for element in document.find_all(id=True)}
            for anchor in document.find_all("a", href=True):
                href = urlsplit(anchor["href"])
                if href.scheme or href.netloc:
                    continue
                if href.path and not (ARTICLE/unquote(href.path)).exists():
                    missing_links.append(href)
                elif not href.path and href.fragment and unquote(href.fragment) not in ids:
                    missing_links.append(href)
            if missing_links:
                raise RuntimeError(f"Broken local article links: {missing_links}")
            tex = (ARTICLE / f"{article_stem}.tex").read_text(encoding="utf-8")
            graphics = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", tex)
            if len(graphics) != expected_figures or any(not (ARTICLE/g).exists() for g in graphics):
                raise RuntimeError("LaTeX figure references did not resolve")
            report["artifact_checks"] = {"embedded_figures": len(images), "local_links_resolve": True,
                                         "latex_figure_paths_resolve": True, "pdf_compiled": False}
            report["article_source_file"] = f"{article_stem}.qmd"
            report["article_reading_copy"] = f"{article_stem}.html"
            report["article_source_sha256"] = hashlib.sha256((ARTICLE/f"{article_stem}.qmd").read_bytes()).hexdigest()
            if article_stem == "technical-article":
                run([sys.executable, CODE/"render_readers.py", "--study", folder])
                report["public_article"] = "rendered and checked separately; see reader-validation.json"
        report["dependencies"] = {
            name: importlib.metadata.version(name)
            for name in ["numpy", "matplotlib", "sympy", "scipy", "mpmath", "Pillow", "nbformat", "nbclient", "nbconvert", "ipykernel", "pypdf", "PyYAML"]
        }
        report["parameters_sha256"] = hashlib.sha256((STUDY/"parameters.json").read_bytes()).hexdigest()
        # Freeze exact installed dependencies, excluding this editable local package.
        freeze = subprocess.check_output([sys.executable, "-m", "pip", "freeze", "--exclude-editable"], text=True)
        (CODE / "requirements-lock.txt").write_text(
            "# Captured on Python " + platform.python_version() + "; install local package with -e ./code --no-deps\n" + freeze,
            encoding="utf-8")
        report["status"] = "passed"
    except Exception as error:
        report["status"] = "failed"
        report["error"] = str(error)
        if isinstance(error, subprocess.CalledProcessError):
            print(error.stdout or "")
            print(error.stderr or "")
        raise
    finally:
        report["elapsed_seconds"] = round(time.perf_counter()-start, 3)
        report_path.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    main()

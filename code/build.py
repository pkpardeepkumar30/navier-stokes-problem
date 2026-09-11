"""Validate and regenerate D1 from the repository's isolated Python environment.

Usage (repository root): code/.venv/Scripts/python.exe code/build.py
Use --compute-only for the notebook/data/figures without Quarto publishing.
"""

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import re
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
CODE = ROOT / "code"
STUDY = CODE / "d1_core_size"
ARTICLE = ROOT / "writing" / "d1_core_size"
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
    args = parser.parse_args()
    import nbformat
    from nbclient import NotebookClient
    from nscomp.d1 import generate

    start = time.perf_counter()
    report_path = STUDY / "validation.json"
    report = {"status": "in_progress", "python": platform.python_version(),
              "platform": platform.system(), "rendered_formats": []}
    report_path.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    try:
        tests = run([sys.executable, "-m", "unittest", "discover", "-s", CODE / "tests", "-v"],
                    capture_output=True, text=True, encoding="utf-8")
        test_log = tests.stdout + tests.stderr
        (STUDY / "validation-tests.txt").write_text(test_log, encoding="utf-8")
        print(test_log)
        report["scientific_checks"] = "passed; see validation-tests.txt"
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
                completed = run([quarto, "render", ARTICLE / "article.qmd", "--to", output_format],
                                env=env, capture_output=True, text=True, encoding="utf-8")
                (CODE / ".cache" / f"d1-quarto-{output_format}.log").parent.mkdir(parents=True, exist_ok=True)
                (CODE / ".cache" / f"d1-quarto-{output_format}.log").write_text(
                    completed.stdout+completed.stderr, encoding="utf-8")
                report["rendered_formats"].append(output_format)
                print(f"Rendered {output_format}.")
            html = (ARTICLE / "article.html").read_text(encoding="utf-8")
            if "data:image/png;base64," not in html or "{python}" in html:
                raise RuntimeError("HTML is missing embedded figures or has unresolved inline code")
            report["html_checks"] = "figures embedded; inline Python evaluated"
            from bs4 import BeautifulSoup
            document = BeautifulSoup(html, "html.parser")
            images = document.find_all("img")
            if len(images) != 3 or any(not img.get("src", "").startswith("data:image/") for img in images):
                raise RuntimeError("Expected three embedded scientific figures")
            missing_links = []
            for anchor in document.find_all("a", href=True):
                href = anchor["href"]
                if href.startswith((".", "source-notes")) and not (ARTICLE/href.split("#")[0]).exists():
                    missing_links.append(href)
            if missing_links:
                raise RuntimeError(f"Broken local article links: {missing_links}")
            tex = (ARTICLE / "article.tex").read_text(encoding="utf-8")
            graphics = re.findall(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", tex)
            if len(graphics) != 3 or any(not (ARTICLE/g).exists() for g in graphics):
                raise RuntimeError("LaTeX figure references did not resolve")
            report["artifact_checks"] = {"embedded_figures": len(images), "local_links_resolve": True,
                                         "latex_figure_paths_resolve": True, "pdf_compiled": False}
        report["dependencies"] = {
            name: importlib.metadata.version(name)
            for name in ["numpy", "matplotlib", "sympy", "nbformat", "nbclient", "nbconvert", "ipykernel", "pypdf", "PyYAML"]
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

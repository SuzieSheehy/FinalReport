"""Compile the draft from the supplied LaTeX, bibliography and figure PDFs."""
from pathlib import Path
import shutil
import subprocess

root = Path(__file__).resolve().parent
latex = shutil.which('pdflatex')
bibtex = shutil.which('bibtex') or shutil.which('bibtex.original')
if not latex or not bibtex:
    raise SystemExit('Install a LaTeX distribution containing pdflatex and BibTeX.')
commands = [
    [latex, '-interaction=nonstopmode', '-halt-on-error', 'main.tex'],
    [bibtex, 'main'],
    [latex, '-interaction=nonstopmode', '-halt-on-error', 'main.tex'],
    [latex, '-interaction=nonstopmode', '-halt-on-error', 'main.tex'],
]
for command in commands:
    subprocess.run(command, cwd=root, check=True)
shutil.copy2(root / 'main.pdf', root / 'Sandpile_Draft_Report.pdf')
print('Built Sandpile_Draft_Report.pdf')

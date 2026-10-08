# Two-dimensional sandpile dynamics

COMP90072 draft report by **Zhang Wenbo**, student ID **1873749**.

Date: 8 October 2026.

Read **[Sandpile_Draft_Report.pdf](Sandpile_Draft_Report.pdf)** for the compiled 16-page draft report, *Two-dimensional sandpile dynamics: Avalanche statistics and random redistribution*.

The report covers the deterministic BTW model, stabilisation theory, implementation checks, avalanche statistics on finite open grids, and a comparison with random redistribution. It includes discussion and questions for feedback before the final report.

## Report files

- `main.tex`: complete report source and compilation entry point.
- `references.bib`: bibliography database; `main.bbl` is the compiled bibliography.
- `figures/`: vector PDF figures used by the report.
- `flowchart.tex`: retained source for the algorithm diagram.
- `build_draft.py`: rebuilds the report using pdfLaTeX and BibTeX.

`main.tex` is self-contained and does not load the original template's separate chapter, cover or preamble files. The existing template files and licence remain in the repository.

## Compile

With a TeX Live installation containing the packages used in `main.tex`, run:

```text
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex
pdflatex -no-shell-escape -interaction=nonstopmode -halt-on-error main.tex
```

This generates `main.pdf`. Alternatively, `python build_draft.py` runs the compilation sequence and copies the result to `Sandpile_Draft_Report.pdf`.

On Overleaf, import the project files and select `main.tex` as the main document. The simulation software and saved trajectories belong to the separate project software package; this submission contains the report and the assets needed to compile it.

The acknowledgements and use-of-assistance statement are included in the report.

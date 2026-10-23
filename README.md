# COMP90072 final report

**Estimating Respiratory Frequency from ECG**  
Liang Zhang — student ID 1721010 — 23 October 2026

Open **main.pdf** to read the report. Compile **main.tex** to regenerate it.

## Files

| File | Purpose |
| --- | --- |
| `main.tex` | Main document and chapter order |
| `preamble.tex` | Packages, fonts, page layout and numbering |
| `cover.tex` | Title, student details and report date |
| `abstract.tex` | Abstract |
| `chapter1.tex`–`chapter10.tex` | Editable report text, equations and tables |
| `biblio.bib` | All nine references, used by BibTeX |
| `figures/` | Six report figures and the university banner from the template |
| `main.pdf` | Compiled report for assessment without recompilation |
| `SUBMISSION_GUIDE_ZH.md` | Chinese instructions for compiling and submitting on GitHub |
| `LICENSE` | Original course template's BSD 3-Clause license |

## Build

Use a LaTeX installation with **pdfLaTeX, BibTeX and latexmk**, then run this command in the directory containing `main.tex`:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Alternatively:

```sh
pdflatex main.tex
bibtex main
pdflatex main.tex
pdflatex main.tex
```

All figures and bibliography data needed to build the report are included. The PDF and numerical analysis do not need to be rerun to view this report. After editing the LaTeX source, rebuild and upload the updated `main.pdf` together with the changed sources.

## Scope of this package

The report was converted from the supplied **Report 90072 .pdf**. Its substantive text, student details, six figures, six tables, three equations and nine references were retained. The figures were extracted from that PDF; no plots or numerical results were regenerated. The title page, table of contents, typography and bibliography were reformatted for LaTeX.

This is the **report-source package**. The analysis scripts, notebook, data and CSV results mentioned in Chapter 10 are separate project files and are not included in this package. Add the existing analysis project to your submission if the course requires it. The analysis commands in Chapter 10 are for that separate project; the LaTeX build commands above compile this report.

## Template attribution

Adapted from the COMP90072 [FinalReport template](https://github.com/SuzieSheehy/FinalReport), commit `671f572c1d5da799ae8b5064f71243c9c2df61e7`.

The `report` document class, A4 format, 2 cm margins, `main.tex` / `preamble.tex` / `cover.tex` structure, BibTeX bibliography and university banner follow the template. Example physics text and figures were replaced by this report. Unused packages were removed. The template's unsupported standard-class `10.5pt` option was changed to `11pt`. Figure, table and equation numbers retain the supplied report's sequence.

The original template license is retained in `LICENSE`. The banner remains the university's original asset from the template.

The supplied submission screenshot names `SuzieSheehy/FinalReport` in its written link and shows `clreichardt/FinalReport` in the example image. Use the repository linked by your current course instructions as the submission target.

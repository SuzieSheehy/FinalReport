# FinalReport
Final Report template

The PDF for the report should be generated from main.tex:



Note for Latex or Overleaf users: (Credit to Antonio Spataro)
You may need to move a line between the preamble.tex file to the beginning of main.tex (for Overleaf) or vice versa:
\documentclass[a4paper, 10.5pt, oneside]{report}


## Project report

ECG-Based Respiratory Rate Estimation and Apnea Record Classification by Kaibo Qiao (student ID 1652063).

Compile from the repository root:

```text
pdflatex -interaction=nonstopmode -halt-on-error main.tex
bibtex main
pdflatex -interaction=nonstopmode -halt-on-error main.tex
pdflatex -interaction=nonstopmode -halt-on-error main.tex
```

The compiled report is main.pdf. Required PDF figures are in figures/. Bibliography entries are in biblio.bib. AI use is disclosed in the report appendix and AI_USAGE_DISCLOSURE.md.

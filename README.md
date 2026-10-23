# Respiratory Rates: Expanded report source

Upload this ZIP as a new Overleaf project and choose `main.tex` as the main document. Compile with pdfLaTeX. All required figures and seven table fragments are included.

For local compilation, run `pdflatex main.tex` three times. Bibliographic entries are in `references.tex`; `biblio.bib` retains the same source entries for reference. Biber is not required.

The report text is in `chapter1.tex`, and formatting is in `cover.tex` and `preamble.tex`. All seven figures and seven tables are copied directly from the verified outputs of the matching Expanded Synced software package.

To regenerate them, put the separately supplied `data.zip` beside the software's `run_all.py` and run `python run_all.py --data-zip data.zip`. This runs 35 automated tests, all four analysis stages and `make_report_tables.py`, which generates all seven tables. Copy the regenerated table files and the seven figure PDFs used here into this report project's `tables/` and `figures/` directories before recompiling.

The software README maps each report figure and table to its script and saved CSV. The former separate `build_additional_tables.py` step has been integrated into the software's table generator.

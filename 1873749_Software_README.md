# COMP90072 Project software — Zhang Wenbo (1873749)

`1873749_Software.zip` contains the sandpile project software prepared on
18 October 2026. This pull request makes that software package available to
Suzie alongside the earlier draft-report submission (PR #5).

## Download and run

Download `1873749_Software.zip` using GitHub's **Download raw file** button,
extract it, and open a terminal in `Sandpile_Project_Software`.

```text
python -m venv .venv
```

Activate the environment with `.venv\Scripts\Activate.ps1` in PowerShell,
`.venv\Scripts\activate.bat` in Command Prompt, or
`source .venv/bin/activate` on macOS/Linux. Then run:

```text
python -m pip install -r software/requirements.txt
python software/sandpile.py test --out checks
python software/sandpile.py demo --out demo_results
```

The reference environment uses Python 3.13.5, NumPy 2.3.5, SciPy 1.17.0 and
Matplotlib 3.10.8. The archive's root `README.md` documents full reproduction,
analysis and optional report-compilation commands. Use a new simulation output
directory and do not run the checks with Python's `-O` flag.

## Package contents and integrity

The archive includes the simulation and analysis code, pinned dependencies,
README files, 38 reference NPZ records, 12 compressed readable CSV data files,
summary tables, figures, implementation and analysis notes, verification scripts,
and the final report PDF with its LaTeX source and figure assets.

- File size: 10,284,331 bytes.
- SHA-256: `11edb535d5a0e1b8049ff583a39a5aa84facda0da9a2fac71a300720f19f712f`.
- 185 files, including `MANIFEST.sha256`; all 184 listed file hashes were checked.

This is the same ZIP supplied for the Canvas **Project software** upload.
The final-report PDF is also supplied separately for the Canvas **Final report**
assignment. This GitHub pull request does not establish either Canvas submission
status.

The report and archive README include the acknowledgement of AI assistance used
in implementation, checking, analysis and writing.

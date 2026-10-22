NV-CENTRE OPTICAL READOUT SIMULATION
====================================

Project title:
From Deterministic Rate Equations to Stochastic Optical Readout of
Nitrogen-Vacancy Centres

Author: Gu Gu
Subject: COMP90072 - The Art of Scientific Computing


1. PROJECT OVERVIEW
-------------------

This submission implements a seven-state rate-equation model for the optical
readout of nitrogen-vacancy (NV) centres. The code progresses through four
computational investigations:

1. deterministic state-population and photoluminescence dynamics;
2. optimisation of the optical readout duration;
3. stochastic Poisson photon-count simulation and threshold classification;
4. recovery of selected transition parameters from synthetic time-binned data.

All transition rates are stored in MHz (equivalently, inverse microseconds),
and all times are stored in microseconds. No external data files or internet
connection are required to reproduce the numerical results.


2. REQUIRED FILENAMES AND DIRECTORY LAYOUT
------------------------------------------

Place the following six Python files in the same directory:

    CONSTANTS.py
    NV_base_model.py
    Population_states.py
    Readout_optimisation.py
    Stochastic_readout.py
    Parameter_inference.py

The filenames and capitalisation are important. NV_base_model.py imports
CONSTANTS, and the four result scripts import NV_base_model. On case-sensitive
systems, changing the capitalisation will cause an import error.

The scripts write figures into the current working directory. Running a script
again will replace an existing output image with the same name.


3. SOFTWARE REQUIREMENTS
------------------------

The submission was syntax-checked and all four result scripts were run
successfully with the following exact environment:

    Python       3.13.5
    NumPy        2.3.3
    SciPy        1.16.2
    Matplotlib   3.10.7

To check the version of your current library, add:
    --version
behind the library and run as a command in the terminal. 

If you want to run using the exact versions of libraries, run these commands:

    pip install NumPy==2.3.3
    pip install SciPy==1.16.2
    pip install Matplotlib==3.10.7

The code uses established NumPy, SciPy and Matplotlib interfaces and may work
with other recent versions, but only the versions listed above were tested for
this submission.


4. HOW TO RUN THE CODE
----------------------

Open a terminal in the directory containing the six canonical filenames. The
scripts do not accept command-line arguments; the reproducibility settings are
defined explicitly near the top of each script.

Run the result scripts in this order:

    python Population_states.py
    python Readout_optimisation.py
    python Stochastic_readout.py
    python Parameter_inference.py

CONSTANTS.py and NV_base_model.py are support modules and do not need to be run
separately. Each result script starts a new Python process and imports them.


The scripts save figures but do not call plt.show(), so no plot window is
required.


5. EXPECTED CONSOLE OUTPUT
--------------------------

The fixed settings and random seed should reproduce the following values.
Minor changes in the last printed digit may occur with different library
versions or hardware.

Population_states.py:

    largest normalisation errors: approximately 1.02e-14 and 8.66e-15
    expected photons, initial m_s=0:   13.005622
    expected photons, initial m_s=-1:  11.973622
    photon-count difference:            1.032000

    Figure 3 -> State_populations_initial_spin0.png
		State_populations_initial_spin1.png
    Figure 4 -> PL_comparison.png

Readout_optimisation.py:

    optimal readout time: 196.0 ns
    expected m_s=0 count:   2.791968
    expected m_s=-1 count:  1.936564
    single-shot SNR:         0.393376

    Figure 5 -> Cumulative_photon_counts.png
		Readout_window_optimisation.png

Stochastic_readout.py:

    m_s=0 deterministic mean:  2.7920
    m_s=0 sample mean:         2.7689
    m_s=0 sample variance:     2.7608
    m_s=-1 deterministic mean: 1.9366
    m_s=-1 sample mean:        1.9402
    m_s=-1 sample variance:    1.9144
    best threshold:            3 photons
    classification accuracy:  0.6053

    Figure 6 -> Poisson_photon_count_probabilities.png

Parameter_inference.py:

    beta:  generating 0.7000, fitted 0.6886
    k_35:  generating 10.0000 MHz, fitted 9.6535 MHz
    k_45:  generating 55.0000 MHz, fitted 53.2065 MHz
    optimiser success: True
    objective-function evaluations: 7

    Figure 7 -> Parameter_inference_fit.png




# Project: Finite graphene flakes using Hückel

A project generating a finite honeycomb lattice and solving a nearest neighbour π electron Hückel Hamiltonian, calculating orbital energies and HOMO-LUMO separations and visualising molecular orbitals. 
The code investigates how the electronic structure changes between small molecules and larger graphene flakes. It permits an incomplete final ring, which is convenient mathematically, but isn't necessarily representative of a physical system. 

Note: The code takes around 2 minutes for me to run, when plotting everything in the default settings


# Project Files:

    Project/
    ├── Final_code_AoSC.py          # Main program (functions and main settings) 
    ├── README.md                   # This guide
    └── Figures/                    # Output figures (created automatically)


# Dependencies 

The project was developed on Windows using python 3.13.12 and the following packages:


| Package     |  Version  |  Purpose                                        |
|:--------|:--------|:--------|
| NumPy       |  2.5.4    |  Arrays, Hamiltonian matrices and eigenvalues   |
| Matplotlib  |  3.11.2   |  2D geometry, spectra and contour plots         |
| PyVista     |  0.49.1   |  3D orbital isosurfaces                         |
| SciPy       |  1.18.1   |  Efficient calculations of central eigenvalues  |


# Default Parameters and Units

Each carbon atom contributes one π orbital. 


| Parameter    | Default value   |  Meaning                                                                |
|:--------|:--------|:--------|
| alpha        |  0.0            |  On site energy, in model units                                         |
| beta         |  -1.0           |  Nearest neighbour hopping, in model units                              |
| a_cc         |  1.42           |  Carbon-Carbon separation, Å                                            |
| tolerance    |  0.1            |  Bond distance tolerance, Å                                             |
| l            |  1              |  Angular momentum used for visualising Carbon 2p orbitals               |
| beta_eV      |  2.4            |  Assumed physical magnitude of hopping parameter for energy conversion  |
| lam          |  1.5            |  Slater radial decay parameter in orbitals()                            |
| grid_points  |  40             |  Grid points per axis in mesh()                                         |

Huckel eigenvalues are in model energy units, not eV. The HOMO-LUMO separation is in |beta| units. 


# Choosing Parameters in main()

|Parameter             |  Example Value    |  Data type    |  Meaning                                                 |
|:--------|:--------|:--------|:--------|
|n_sites               |  24               |  Int          |  Number of atoms in lattice                              |
|n_values              |  [6, 24, 54, 96]  |  List of Int  |  Flakes compared in multi-spectrum figure                |
|max_sites             |  294              |  Int          |  Largest number of sites in intermediate gap plot        |
|orb_index             |  20               |  Int          |  Individual orbital to display, valid from 1 to n_sites  |

# Plot switches
|Switch             |  Example Value    |  Data type    |  Meaning                                                 |
|:--------|:--------|:--------|:--------|
|PLOT_POSITIONS        |  True             |  Boolean      |  Lattice diagram                                         |
|PLOT_CONTOUR_MO       |  True             |  Boolean      |  2D cross sectional plots of molecular orbitals          |
|PLOT_ENERGY_SPECTRUM  |  True             |  Boolean      |  One flake orbital energy spectrum                       |
|PLOT_ENERGIES         |  True             |  Boolean      |  Multiple flake energy spectra                           |
|PLOT_3D_ORBITAL       |  True             |  Boolean      |  One 3D molecular orbital                                |
|PLOT_3D_ORBITALS_ALL  |  True             |  Boolean      |  All 3D moelcular orbitals, three views each (slow)      |
|PLOT_HOMO_GAP         |  True             |  Boolean      |  Gap vs size comparison                                  |
|PLOT_INCOMPLETE       |  True             |  Boolean      |  Show all even size series in gap comparison             |
|PLOT_COMPLETE         |  True             |  Boolean      |  Highlight completed shell gaps in comparison            |
|PLOT_LITERATURE       |  True             |  Boolean      |  Add reference gap markers                               |
|PLOT_SCALING_TEST     |  True             |  Boolean      |  Log–log fit for completed shell gaps                    |

To produce only one selected 3D orbital, set PLOT_3D_ORBITAL = True and PLOT_3D_ORBITALS_ALL = False. PLOT_INCOMPLETE plots every even N from 6 through to max_sites, including completed shells.


# Outputs

To reproduce a result, use the script version, parameter values, plot flags and filenames. The following configurations correspond to existing plotting functions:


|Result                             |  Required settings/call                                 |  Saved output                                                |
|:--------|:--------|:--------|
|Benzene energy spectrum            |  n_sites=6, PLOT_ENERGY_SPECTRUM=True                   |  6_energy_spectrum.png                                       |
|Geometry of N sized flake          |  n_sites=N, PLOT_POSITIONS=True                         |  N_positions.png                                             |
|Energy spectra comparison          |  n_values=[6,24,54,96], PLOT_ENERGIES=True              |  4_plots_energies.png                                        |
|Gap comparison, any combination    |  PLOT_HOMO_GAP=True, choose the three gap series flags  |  294_gap_comparisons_gap_vs_size_plot.png for max_sites=294  |
|Completed shell scaling fit        |  PLOT_SCALING_TEST=True, max_sites=294                  |  294_complete_shell_scaling_test.png                         |
|Single 3D orbitals                 |  n_sites=N, orb_index=k, PLOT_3D_ORBITAL=True           |  k_MO_individual.svg                                         |
|All 3D orbitals, three views each  |  n_sites=N, PLOT_3D_ORBITALS_ALL=True                   |  N_MOs.svg                                                   |
|2D MO contours                     |  n_sites=N, PLOT_CONTOUR_MO=True                        |  N_MOs_1-4.png, N_MOs_5-8.png, etc.                          |


# Key Functions

|Function                 |  Input                                                          |  Output/Purpose                                          |
|:--------|:--------|:--------|
|generate_graphene_shell  |  n_sites, a_cc                                                  |  Array of planar atomic coordinates (positions = (N,2))  |
|find_bonds               |  positions, a_cc, tolerance                                     |  List of nearest neighbour index pairs                   |
|slater, mesh, orbitals   |  l, lam, r, y, grid_points, n_sites, positions, e_vec           |  Recontructs 3D wavefunction and orbital isosurfaces     |
|H_graphene               |  n_sites, positions, alpha, beta, a_cc, tolerance               |  Symmetric Hamiltonian shape (N,N)                       |
|solve_huckel             |  n_sites, alpha, beta, a_cc, tolerance                          |  Positions, eigenvalues and eigenvector matrix           |
|add_lattice_to_plotter   |  plotter, positions, a_cc, tolerance, atom_radius, bond_radius  |  Add bonds and atoms to 3D orbitals                      |
|homo_lumo_gap            |  e_val                                                          |  HOMO energy, LUMO energy, separation                    |
|calculate_gaps           |  n_values, alpha, beta, a_cc, tolerance, beta_eV                |  Compute the central gaps in eV for even sizes           |
|plot_*                   |  results and output_dir                                         |  Saves 2D or 3D graphics                                 |


# Model Limitations

- Incomplete final rings may create physically unrealistic molecules.
- There are no electron-electron interactions, spin polarisation or longer range hopping terms.
- The 3D Slater basis is used to visualise the Hückel coefficients, not to derive or refine the Hamiltonian.
- MO wavefunction colours represent sign/phase, not positive and negative charge density.


# Code References

Jeong & Lee (2024), Calculation of Electronic Structures of Graphene Nanoflakes by Hückel Molecular Orbital Method, Korean Journal of Materials Research 34. DOI: 10.3740/MRSK.2024.34.12.620.

Heiskanen, Manninen & Akola (2008), Electronic structure of triangular, hexagonal and round graphene flakes near the Fermi level, New Journal of Physics 10, 103015. DOI: 10.1088/1367-2630/10/10/103015.

Castro Neto et al. (2009), The electronic properties of graphene, Reviews of Modern Physics 81, 109. DOI: 10.1103/RevModPhys.81.109.



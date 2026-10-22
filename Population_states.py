import CONSTANTS as c
import numpy as np
import matplotlib.pyplot as plt

from NV_base_model import (constant_laser, 
    STATE_LABELS, 
    photoluminescence, 
    Initial_ms0, 
    Initial_ms1, 
    integrated_counts)


################################
# SIMULATION SETTINGS

# For these results, using 0 to 1001 microsecond scale. 
t_us = np.linspace(0.0, 1.0, 1001)

P_ms0 = constant_laser(t_us, Initial_ms0)
P_ms1 = constant_laser(t_us, Initial_ms1)

PL_ms0 = photoluminescence(P_ms0)
PL_ms1 = photoluminescence(P_ms1)

assert np.allclose(P_ms0.sum(axis=1), 1.0, atol=1e-10)
assert np.allclose(P_ms1.sum(axis=1), 1.0, atol=1e-10)

# Check for numerical rounging errors
print(f"Largest normalization error: {np.max(P_ms0.sum(axis=1) - 1)}")
print(f"Largest normalization error: {np.max(P_ms1.sum(axis=1) - 1)}")

################################
plt.figure()
for i, labels in enumerate(STATE_LABELS):
    plt.plot(t_us, P_ms0[:, i], label=labels)


plt.xlabel(r"Laser time $t$ ($\mu$s)")
plt.ylabel("Population")
plt.title("State populations: initial NV$^-$ in $m_s=0$")
plt.grid()
plt.legend(fontsize = "small")
plt.savefig("State_populations_initial_spin0.png")

################################
plt.figure()
for i, labels in enumerate(STATE_LABELS):
    plt.plot(t_us, P_ms1[:, i], label=labels)


plt.xlabel(r"Laser time $t$ ($\mu$s)")
plt.ylabel("Population")
plt.title("State populations: initial NV$^-$ in $m_s=1$")
plt.grid()
plt.legend(fontsize = "small")
plt.savefig("State_populations_initial_spin1.png")

################################
plt.figure()
normalisation = max(PL_ms0.max(), PL_ms1.max())
plt.plot(t_us, PL_ms0 / normalisation, label="initial NV$^-$: $m_s=0$")
plt.plot(t_us, PL_ms1 / normalisation, label="initial NV$^-$: $m_s=1$")
plt.title("Photoluminescence comparison of different initialisations") # Compare with Di Wang's plot, looks very similar
plt.legend()
plt.xlim(0, 0.6)
plt.xlabel(r"Laser time $t$ ($\mu$s)")
plt.ylabel("Normalised Photoluminescence")
plt.savefig("PL_comparison.png", dpi=300)

################################
# Considering the total expected photon counts of both initial states.
counts_ms0 = integrated_counts(t_us=t_us, populations=P_ms0)

counts_ms1 = integrated_counts(t_us=t_us, populations=P_ms1)

print(f"Expected photons for initial ms=0: {counts_ms0:.6f}")
print(f"Expected photons for initial ms=1: {counts_ms1:.6f}")
print(f"Photon-count difference: {counts_ms0 - counts_ms1:.6f}")
# Note: photon counts are experimentally integer values, however here they display as fractions. This is 
# allowed because we are explicitly looking at expected photons, providing an estimate rather than discrete
# capturing of photons. 

################################



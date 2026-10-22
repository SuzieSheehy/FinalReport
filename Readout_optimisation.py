import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

from NV_base_model import (constant_laser,
    photoluminescence, 
    Initial_ms0, 
    Initial_ms1)


################################
# SIMULATION SETTINGS

# Simulate one microsecond with one-nanosecond spacing.
t_us = np.linspace(0.0, 1.0, 1001)

# Illustrative detector parameters.
# Still using default photoluminescence detection conditions. Conditions are subject to experimental 
# calibration, and thus we ignore this for now. 

P_ms0 = constant_laser(t_us, Initial_ms0)
P_ms1 = constant_laser(t_us, Initial_ms1)

PL_ms0 = photoluminescence(P_ms0)
PL_ms1 = photoluminescence(P_ms1)


################################
# Cumulative expected photon counts

counts_ms0 = cumulative_trapezoid(
    PL_ms0,
    x=t_us,
    initial=0.0
)

counts_ms1 = cumulative_trapezoid(
    PL_ms1,
    x=t_us,
    initial=0.0
)


################################
# Shot-noise-limited signal-to-noise ratio

count_difference = np.abs(counts_ms0 - counts_ms1)

# For independent Poisson photon counts:
# variance(N0 - N1) = variance(N0) + variance(N1) = N0 + N1
shot_noise = np.sqrt(counts_ms0 + counts_ms1)

# Avoid division by zero at t = 0.
SNR = np.divide(
    count_difference,
    shot_noise,
    out=np.zeros_like(count_difference),
    where=shot_noise > 0
)


################################
# Locate the optimum readout duration

# Return the argument at which the maximum SNR was detected. 
optimal_index = np.argmax(SNR)

optimal_time_us = t_us[optimal_index]
optimal_SNR = SNR[optimal_index]
optimal_counts_ms0 = counts_ms0[optimal_index]
optimal_counts_ms1 = counts_ms1[optimal_index]

print(f"Optimal readout time: {optimal_time_us * 1000:.1f} nanoseconds")
print(f"Expected ms=0 count: {optimal_counts_ms0:.6f}")
print(f"Expected ms=1 count: {optimal_counts_ms1:.6f}")
print(f"Single-shot SNR: {optimal_SNR:.6f}")


################################
# Plot the result

# Cumulative photon counts
plt.figure()
plt.plot(
    t_us,
    counts_ms0,
    label=r"Initial $m_s=0$"
)

plt.plot(
    t_us,
    counts_ms1,
    label=r"Initial $m_s=-1$"
)

plt.axvline(
    optimal_time_us,
    color="black",
    linestyle="--",
    label="Optimal readout time"
)

plt.ylabel("Expected detected photons")
plt.xlabel(r"Readout duration ($\mu$s)")
plt.title("Cumulative photon counts")
plt.legend()
plt.grid()
plt.savefig("Cumulative_photon_counts.png")


# Signal-to-noise ratio
plt.figure()
plt.plot(
    t_us,
    SNR,
    color="tab:purple"
)

plt.axvline(
    optimal_time_us,
    color="black",
    linestyle="--"
)

plt.scatter(
    optimal_time_us,
    optimal_SNR,
    color="red",
    zorder=3,
    label=(rf"$t_\mathrm{{opt}}={optimal_time_us * 1000:.0f}$ ns")
)

plt.xlabel(r"Readout duration ($\mu$s)")
plt.ylabel("Shot-noise-limited SNR")
plt.title("Readout-window optimisation")
plt.legend()
plt.grid()
plt.savefig("Readout_window_optimisation.png", dpi=300)



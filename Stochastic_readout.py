import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import cumulative_trapezoid

from NV_base_model import (
    constant_laser,
    photoluminescence,
    Initial_ms0,
    Initial_ms1)


################################
# DETERMINISTIC READOUT MEANS -> Similar to Readout_optimisation set up

t_us = np.linspace(0.0, 1.0, 1001)

P_ms0 = constant_laser(t_us, Initial_ms0)
P_ms1 = constant_laser(t_us, Initial_ms1)

PL_ms0 = photoluminescence(
    P_ms0,
    detection_scale=1.0,
    background=0.0,
    alpha=0.3)

PL_ms1 = photoluminescence(
    P_ms1,
    detection_scale=1.0,
    background=0.0,
    alpha=0.3)

counts_ms0 = cumulative_trapezoid(
    PL_ms0,
    x=t_us,
    initial=0.0)

counts_ms1 = cumulative_trapezoid(
    PL_ms1,
    x=t_us,
    initial=0.0)

count_difference = np.abs(counts_ms0 - counts_ms1)
shot_noise = np.sqrt(counts_ms0 + counts_ms1)

# Only considering the Singal-to-Noise ratio with a positive noice value.
SNR = np.divide(
    count_difference,
    shot_noise,
    out=np.zeros_like(count_difference),
    where=shot_noise > 0)

optimal_index = np.argmax(SNR)
optimal_time_us = t_us[optimal_index]
optimal_counts_ms0 = counts_ms0[optimal_index]
optimal_counts_ms1 = counts_ms1[optimal_index]


################################
# STOCHASTIC PHOTON COUNTS

# Number of repeated readout experiments.
number_of_shots = 10_000

# Fixed seed makes the result reproducible.
# Use a random number generator to simulate poisson type distribution. 
random_generator = np.random.default_rng(seed=676767)

# The deterministic integrated counts become the Poisson means.
# Different means for different initial states require use to center the distribution for 
# each separately. 
photon_counts_ms0 = random_generator.poisson(
    lam=optimal_counts_ms0,
    size=number_of_shots)

photon_counts_ms1 = random_generator.poisson(
    lam=optimal_counts_ms1,
    size=number_of_shots)


################################
# Compare stochastic and deterministic results

print("\nPoisson photon-count simulation")

print(f"ms=0 deterministic mean: {optimal_counts_ms0:.4f}")

print(f"ms=0 stochastic mean:    {np.mean(photon_counts_ms0):.4f}")

print(f"ms=0 stochastic variance:" f"{np.var(photon_counts_ms0, ddof=1):.4f}")

print(f"ms=-1 deterministic mean: {optimal_counts_ms1:.4f}")

print(f"ms=-1 stochastic mean:    {np.mean(photon_counts_ms1):.4f}")

print(f"ms=-1 stochastic variance:" f"{np.var(photon_counts_ms1, ddof=1):.4f}")


################################
# Find the best photon count threshold

maximum_count = max(photon_counts_ms0.max(), photon_counts_ms1.max())

possible_thresholds = np.arange(maximum_count + 2)
threshold_accuracies = []

for threshold in possible_thresholds:

    # Higher photon counts are identified as ms=0.
    accuracy_ms0 = np.mean(photon_counts_ms0 >= threshold)

    # Lower photon counts are identified as ms=-1.
    accuracy_ms1 = np.mean(photon_counts_ms1 < threshold)

    # Assume both initial spin states are equally likely.
    total_accuracy = 0.5 * (accuracy_ms0 + accuracy_ms1)

    threshold_accuracies.append(total_accuracy)

threshold_accuracies = np.asarray(threshold_accuracies)

best_threshold_index = np.argmax(threshold_accuracies)
best_threshold = possible_thresholds[best_threshold_index]
best_accuracy = threshold_accuracies[best_threshold_index]

print(f"\nBest photon-count threshold: {best_threshold}")
print(f"Estimated classification accuracy: {best_accuracy:.4f}")


################################
# Discrete photon-count plot

# When choosing our horizontal range, ignore extreme tail of data. 
maximum_display_count = int(np.ceil(max(
    np.quantile(photon_counts_ms0, 0.999),
    np.quantile(photon_counts_ms1, 0.999),)))

# +1 is included because counting from 0. 
photon_numbers = np.arange(maximum_display_count + 1)

# np.bincount counts how many times each integer occurs.
# Probability of k photons = number of experiments containing k photons/total number of experiments.

probability_ms0 = (
    np.bincount(photon_counts_ms0, 
    minlength=maximum_display_count + 1)[:maximum_display_count + 1]/ number_of_shots)

probability_ms1 = (
    np.bincount(photon_counts_ms1, 
    minlength=maximum_display_count + 1)[:maximum_display_count + 1]/ number_of_shots)

# Used to serperate the two different bar type data.
bar_width = 0.38

plt.bar(
    photon_numbers - bar_width / 2,
    probability_ms0,
    width=bar_width,
    color="cornflowerblue",
    edgecolor="black",
    linewidth=0.5,
    label=r"Initial $m_s=0$")

plt.bar(
    photon_numbers + bar_width / 2,
    probability_ms1,
    width=bar_width,
    color="coral",
    edgecolor="black",
    linewidth=0.5,
    label=r"Initial $m_s=-1$")

# A threshold of 3 means:
# 0, 1, 2 -> classify as ms=-1
# 3, 4,...-> classify as ms=0
# The boundary is therefore drawn at 2.5 to seperate these classifciations. 
# If below 3 -> ms=-1
# If above or equal to 3 -> ms=0

plt.axvline(
    best_threshold - 0.5,
    color="black",
    linestyle="--",
    linewidth=1.5,
    label=f"Decision threshold = {best_threshold}")

plt.xticks(photon_numbers)
plt.xlim(-0.7, maximum_display_count + 0.7)

plt.xlabel("Photons observed in one readout")
plt.ylabel("Empirical probability")
plt.title("Simulated photon-count probability distributions")
plt.legend(fontsize="small")

plt.savefig("Poisson_photon_count_probabilities.png", dpi=300, bbox_inches="tight")


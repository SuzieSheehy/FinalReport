import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

from NV_base_model import (
    constant_laser,
    photoluminescence,
    Initial_ms0,
    Initial_ms1)


#-------------------------------------------------------------#
# 1. Values used to generate the artificial experiment.
#-------------------------------------------------------------#

TRUE_PARAMETERS = np.array([
    0.7,      # beta
    10.0,     # k_35, MHz
    55.0,])   # k_45, MHz

PARAMETER_NAMES = ["beta", "k_35", "k_45"]

# Allowing bounds for parameters to sweep through
LOWER_BOUNDS = np.array([
    0.03,     # beta
    1.0,      # k_35
    5.0])     # k_45

UPPER_BOUNDS = np.array([
    3.0,      # beta
    50.0,     # k_35
    150.0])   # k_45



#-------------------------------------------------------------#
# 2. Predict the mean counts in every time bin.
#-------------------------------------------------------------#

def expected_binned_counts(
    parameters,
    bin_centres_us,
    bin_width_us,
    repetitions):

    # Bin the expected photon detection together, where spin 0 is appended first. 
    # This gives 1 big array, and can be separated later. 

    beta, k_35, k_45 = parameters

    predicted_counts = []

    for initial_state in (Initial_ms0, Initial_ms1):

        populations = constant_laser(
            bin_centres_us,
            initial_state,
            beta=beta,
            k_35=k_35,
            k_45=k_45)

        fluorescence_rate = photoluminescence(
            populations,
            detection_scale=1.0,
            background=0.0,
            alpha=0.3)
        

        # fluorescence_rate has units photons per microsecond.
        # Multiplying by the bin width gives photons per repetition.
        # Multiplying by repetitions gives total expected counts.

        mean_counts = (repetitions * bin_width_us * fluorescence_rate)

        predicted_counts.append(mean_counts)

    return np.concatenate(predicted_counts)


#-------------------------------------------------------------#
# 3. Residual appropriate for Poisson-distributed data.
#-------------------------------------------------------------#
# For every time bin, we want separate data:
# observed -> measured/simulated photon count
# expected -> prediction made by the mdoel
# residual -> observed - expected

# A property of a poisson distributed observation is that its mean and variance are the same, given by
# a probability mass function of:
# P(X = K) = lambda^K * exp{-lambda} / K! 
# where mean = lambda = variance. 

def poisson_deviance_residual(observed, expected):

    observed = np.asarray(observed, dtype=float)
    expected = np.maximum(np.asarray(expected, dtype=float), 1.0e-12)

    log_term = np.zeros_like(observed)

    nonzero = observed > 0  # prevent dividing with zero

    log_term[nonzero] = (observed[nonzero] * np.log(observed[nonzero] / expected[nonzero]))

    # Formula for deviance of Poisson distribution
    deviance = 2.0 * (expected - observed + log_term)

    # Ensuring we eliminate artificial negative values caused by rounding
    deviance = np.maximum(deviance, 0.0)

    # np.sign recovers the direction of the residual. 
    return (np.sign(observed - expected) * np.sqrt(deviance))

# If the residual is positive: model underpredicts the observed counts. 
# If the residual is negative: model overpredicts the observed counts. 

#-------------------------------------------------------------#
# 4. Set up the artificial experiment
#-------------------------------------------------------------#

# Setting up the bins
bin_edges_us = np.linspace(0.0, 0.6, 121)

bin_width_us = bin_edges_us[1] - bin_edges_us[0]

bin_centres_us = (bin_edges_us[:-1] + 0.5 * bin_width_us)

repetitions = 5000

random_generator = np.random.default_rng(seed=676767)


# Calculate the exact mean values for the true parameters.

true_mean_counts = expected_binned_counts(
    TRUE_PARAMETERS,
    bin_centres_us,
    bin_width_us,
    repetitions)


# Generating one noisy artificial experiment. This is fixed when fitting the data. 
observed_counts = random_generator.poisson(true_mean_counts)


#-------------------------------------------------------------#
# 5. Define the objective function for the optimiser.
#-------------------------------------------------------------#


def objective_function(log_parameters):

    # Using a log to exponential transform, we can effectively recover the same information, with the
    # added benefit that the transition rates remain positive and thus physical. Optimising logarithms
    # also allows us to deal with numerical scales that differ considerably, such as the beta scale and
    # the transition rate scale. 

    parameters = np.exp(log_parameters)

    model_counts = expected_binned_counts(
        parameters,
        bin_centres_us,
        bin_width_us,
        repetitions)

    # These are expected values, and we therefore allow for fractional photons. 
    return poisson_deviance_residual(observed_counts, model_counts)


# Deliberately start away from the true parameters.
initial_guess = np.array([
    0.5,     # beta
    7.0,     # k_35
    40.0])   # k_45



#-------------------------------------------------------------#
# 6. Perform the parameter fit.
#-------------------------------------------------------------#

# Trial log parameters -> NV population dynamics -> Photoluminescence -> Predicted photon counts
# -> Poisson Residuals -> Updates trial parameters 

# Setting x_scale to 'jac' helps with hugely differing scales. 
fit_result = least_squares(
    objective_function,
    x0 = np.log(initial_guess),
    bounds =(
        np.log(LOWER_BOUNDS),
        np.log(UPPER_BOUNDS)),
    x_scale = "jac") 

fitted_parameters = np.exp(fit_result.x)

if not fit_result.success:
    raise RuntimeError("Parameter fit failed")


#-------------------------------------------------------------#
# 7. Display the fitted values.
#-------------------------------------------------------------#-

print("\nParameter inference results")
print("---------------------------")

print(f"beta: true = {TRUE_PARAMETERS[0]:.4f}," f"fitted = {fitted_parameters[0]:.4f}")

print(f"k_35: true = {TRUE_PARAMETERS[1]:.4f} MHz," f"fitted = {fitted_parameters[1]:.4f} MHz")

print(f"k_45: true = {TRUE_PARAMETERS[2]:.4f} MHz," f"fitted = {fitted_parameters[2]:.4f} MHz")

print(f"\nOptimiser success: {fit_result.success}")

# How many times our residual function was called in the process of fitting the data.
print(f"Function evaluations: {fit_result.nfev}")


#-------------------------------------------------------------#
# 8. Compare fitted curves with the artificial data.
#-------------------------------------------------------------#

fitted_mean_counts = expected_binned_counts(
    fitted_parameters,
    bin_centres_us,
    bin_width_us,
    repetitions)

number_of_bins = bin_centres_us.size

# Splicing results here.
observed_ms0 = observed_counts[:number_of_bins]
observed_ms1 = observed_counts[number_of_bins:]

fitted_ms0 = fitted_mean_counts[:number_of_bins]
fitted_ms1 = fitted_mean_counts[number_of_bins:]


# Convert total counts back into an average detected rate.
normalisation = repetitions * bin_width_us

observed_rate_ms0 = observed_ms0 / normalisation
observed_rate_ms1 = observed_ms1 / normalisation

fitted_rate_ms0 = fitted_ms0 / normalisation
fitted_rate_ms1 = fitted_ms1 / normalisation


# Calculate the final residuals.
final_residuals = poisson_deviance_residual(observed_counts, fitted_mean_counts)

residuals_ms0 = final_residuals[:number_of_bins]
residuals_ms1 = final_residuals[number_of_bins:]

#-------------------------------------------------------------#
# 9. Validate with deterministic model. 
#-------------------------------------------------------------#
true_rate_ms0 = (true_mean_counts[:number_of_bins] / normalisation)

true_rate_ms1 = (true_mean_counts[number_of_bins:] / normalisation)


#-------------------------------------------------------------#
# Plot the result
#-------------------------------------------------------------#

# Plotting the simulated data and the residuals on two separate plots.
figure, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

axes[0].scatter(
    bin_centres_us,
    observed_rate_ms0,
    s=14,
    alpha=0.65,
    label=r"Simulated data: $m_s=0$")

axes[0].scatter(
    bin_centres_us,
    observed_rate_ms1,
    s=14,
    alpha=0.65,
    label=r"Simulated data: $m_s=-1$")

axes[0].plot(
    bin_centres_us,
    fitted_rate_ms0,
    linewidth=2,
    label=r"Fitted model: $m_s=0$")

axes[0].plot(
    bin_centres_us,
    fitted_rate_ms1,
    linewidth=2,
    label=r"Fitted model: $m_s=-1$")

axes[0].plot(
    bin_centres_us,
    true_rate_ms0,
    color="darkorchid",
    linestyle="--",
    linewidth=1.5,
    label=r"True deterministic model: $m_s=0$")

axes[0].plot(
    bin_centres_us,
    true_rate_ms1,
    color="orangered",
    linestyle="--",
    linewidth=1.5,
    label=r"True deterministic model: $m_s=-1$")

axes[0].set_ylabel(r"Detected rate ($\mu\mathrm{s}^{-1}$)")

axes[0].set_title("Recovery of NV-centre transition rates")

axes[0].legend()
axes[0].grid(alpha=0.3)


axes[1].scatter(
    bin_centres_us,
    residuals_ms0,
    s=14,
    alpha=0.7,
    label=r"$m_s=0$")

axes[1].scatter(
    bin_centres_us,
    residuals_ms1,
    s=14,
    alpha=0.7,
    label=r"$m_s=-1$")

axes[1].axhline(
    0.0,
    color="black",
    linestyle="--")

axes[1].set_xlabel(r"Time after laser turn-on ($\mu$s)")

axes[1].set_ylabel("Poisson deviance residual")

axes[1].legend()
axes[1].grid(alpha=0.3)


plt.tight_layout()

plt.savefig("Parameter_inference_fit.png", dpi=300, bbox_inches="tight")



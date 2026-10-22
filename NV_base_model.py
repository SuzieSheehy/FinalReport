import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp, trapezoid
from scipy.linalg import expm
import CONSTANTS as c

# Transition rates measured in MHz
# Beta terms are associated with laser intensity

# Transition matrix, governs the population state. Some variables are redefined locally because 
# we need to look at different rates without changing them gloabally later. 
def rate_matrix(beta = c.beta,
    k_35 = c.k_35,
    k_45 = c.k_45):

    if beta < 0:
        raise ValueError("beta must be non-negative")

    A = np.array([
        [-beta*c.k_13, 0, c.k_31, 0, c.k_51, 0, beta*c.k_71] , 
        [0, -beta*c.k_24, 0, c.k_42, c.k_52, 0, beta*c.k_72] , 
        [beta*c.k_13, 0, -c.k_31-beta*c.k_36-k_35, 0, 0, 0, 0] , 
        [0, beta*c.k_24, 0, -c.k_42-k_45-beta*c.k_46, 0, 0, 0] , 
        [0, 0, k_35, k_45, -c.k_52-c.k_51, 0, 0] , 
        [0, 0, beta*c.k_36, beta*c.k_46, 0, -beta*c.k_67, c.k_76] , 
        [0, 0, 0, 0, 0, beta*c.k_67, -c.k_76-beta*c.k_71-beta*c.k_72]])
    return A

A = rate_matrix()
# print(A)

# check that all the columns does add to 0 for conserved transition rates. 
assert np.allclose(A.sum(axis=0), 0.0, atol=1e-12)

# Intialising
# Assume entire population resides in states 1, 2, 6. These are the ground states for NV0 and NV minus including both spin states for NV minus. 
# Construct and arbitrary conserved initial state.
def Initial_state(NV_minus_m0_i : float, 
    NV_minus_m1_i : float,
    NV_minus_excited_m0_i : float, 
    NV_minus_excited_m1_i : float, 
    NV_minus_meta_i : float, 
    NV_0_i : float, 
    NV_0_excited_i : float):
    
    # to keep NV populations conserved:

    NV_minus_i = NV_minus_m0_i + NV_minus_m1_i + NV_minus_excited_m0_i + NV_minus_excited_m1_i + NV_minus_meta_i
    
    if NV_minus_i +  NV_0_i + NV_0_excited_i != 1:
        raise ValueError("Population fraction must add up to 1")
        
    if not 0.0 <= NV_minus_m0_i <= 1.0:
        raise ValueError("Initial NV minus spin 0 population must lie within [0, 1]")

    if not 0.0 <= NV_minus_m1_i <= 1.0:
        raise ValueError("Initial NV minus spin 1 population must lie within [0, 1]")

    if not 0.0 <= NV_minus_excited_m0_i <= 1.0:
        raise ValueError("Initial excited NV minus spin 0 population  must lie within [0, 1]")

    if not 0.0 <= NV_minus_excited_m1_i <= 1.0:
        raise ValueError("Initial NV minus metastable population must lie within [0, 1]")

    if not 0.0 <= NV_minus_meta_i <= 1.0:
        raise ValueError("Initial NV minus population must lie within [0, 1]")

    if not 0.0 <= NV_0_i <= 1.0:
        raise ValueError("Initial NV 0 population must lie within [0, 1]")

    if not 0.0 <= NV_0_excited_i <= 1.0:
        raise ValueError("Initial NV 0 excited population must lie within [0, 1]")
    
        
    Initial = np.array([NV_minus_m0_i,
        NV_minus_m1_i,
        NV_minus_excited_m0_i,
        NV_minus_excited_m1_i,
        NV_minus_meta_i,
        NV_0_i,
        NV_0_excited_i])

    return Initial


# Considering a 7:3 ratio of NV minus to NV0 initial population
Initial_ms0 = Initial_state(0.7,0,0,0,0,0.3,0)
Initial_ms1 = Initial_state(0,0.7,0,0,0,0.3,0)
# print(Initial_ms0)

# For a constant laser, can find time evolution by taking the exponential of the transition matrix. 
def constant_laser(t_us, Initial, beta = c.beta, k_35=c.k_35, k_45=c.k_45,):
    
    # Returns array with shape (number of times, 7)
    t_us = np.array(t_us, dtype=float)
    Initial = np.array(Initial, dtype=float)

    # Stacking the matrices to store as data
    A_constant = rate_matrix(beta, k_35=k_35, k_45=k_45)
    return np.vstack([expm(A_constant * t) @ Initial for t in t_us])

# The fastest tabulated rates are about 100 MHz (a 10 ns timescale).
# Set t_us to be 1 micro second with 1000 intervals such that each interval represents 1 nanoseconds. 

# t_us = np.linspace(0.0, 1.0, 1001)
# B = constant_laser(t_us = np.linspace(0.0, 1.0, 1001), Initial = Initial_ms0)
# print(B)

# Detection is via photoluminescence.
# I need to create some way such that the function converts the state populations into the expected rates at which photons are detected. 
# These should in theory, match actual real experimental results. 
def photoluminescence(
    populations: np.ndarray,
    detection_scale: float = 1,
    background: float = 0,
    alpha: float = 0.3): 
    
    P = np.array(populations)

    # The transitions that emit photons are:
    # |3> -> |1>
    # |4> -> |2>
    # |7> -> |6>

    # According to the Di's paper, NV0 fluorescence is mostly filtered out and only left with 30%
    pl_model =  c.k_31 * P[:, 2] + c.k_42 * P[:, 3] + alpha * c.k_76 * P[:, 6]

    # Also, included potential background noise. 
    
    return detection_scale * pl_model + background

# Photoluminescence looks at counts per unit time. To get an understanding of the total counts, we
# can apply an integration. Integrate PL over time. 
# Returns the expected number of photons over the supplied time interval. 

def integrated_counts(
    t_us: np.ndarray,
    populations: np.ndarray,
    detection_scale: float = 1,
    background: float = 0,
    alpha: float = 0.3,
):

    t_us = np.asarray(t_us, dtype=float)
    P = np.asarray(populations, dtype=float)

    # Some easy self-consistency checks: 
    # Check the time array.
    if t_us.ndim != 1:
        raise ValueError("t_us must be a one-dimensional array")

    if len(t_us) < 2:
        raise ValueError("At least two time points are required")

    if np.any(np.diff(t_us) <= 0):
        raise ValueError("t_us must be strictly increasing")

    # Check the population array.
    if P.ndim != 2 or P.shape[1] != 7:
        raise ValueError(
            "populations must have shape (number_of_times, 7)"
        )

    if P.shape[0] != len(t_us):
        raise ValueError(
            "The number of population rows must equal the number of time points"
        )

    # Photon rate at every supplied time, just used PL described previously. 
    photon_rate = photoluminescence(
        populations=P,
        detection_scale=detection_scale,
        background=background,
        alpha=alpha,
    )

    # Numerically integrate photon rate with respect to time.
    expected_counts = trapezoid(photon_rate, x=t_us)

    return float(expected_counts)


STATE_LABELS = [
    r"$|1\rangle$: NV$^-$ ground, $m_s=0$",
    r"$|2\rangle$: NV$^-$ ground, $m_s=-1$",
    r"$|3\rangle$: NV$^-$ excited, $m_s=0$",
    r"$|4\rangle$: NV$^-$ excited, $m_s=-1$",
    r"$|5\rangle$: metastable state",
    r"$|6\rangle$: NV$^0$ ground",
    r"$|7\rangle$: NV$^0$ excited",
]

def numerical_solver(
    t_us: np.ndarray,
    Initial: np.ndarray,
    beta_t,
    step_size: float = 0.001):

    # Define a single generation of the new state.
    # t is the current time, and p is the current state vector. 
    # Placeholders such that we can numerically solve this. 
    # This allows us to implement sequences of laser pulses!! Much more interesting. 
    def single(t, p):
        return rate_matrix(float(beta_t(t))) @ p    
        
    # Integrating differential equation, however now the laser is also time dependent -> beta(t)
    solution = solve_ivp(
        fun=single,
        t_span=(float(t_us[0]), float(t_us[-1])),
        y0=np.array(Initial),
        t_eval=t_us,
        method="RK45", # might potentially look at other integration techniques if time allows. 
        rtol=1e-9,
        atol=1e-12,
        max_step=step_size)

    # Convenient to transpose our data such that P[:,0] represents the population of state 1 
    # at every time, consistent with what we already defined. 
    return solution.y.T





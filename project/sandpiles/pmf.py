import numpy as np
from scipy.special import logsumexp
# ---
# Distribution pmf/pdfs
# ---

# 1. Discrete power law distribution
# p(x; s, x_max) = 1/C x^(-s)
#   x = 1, 2, ... x_max
#   s in (0, inf)

def power_law_pmf(
        params: tuple[float],   # parameter s
        x: np.ndarray[int],     # data
        x_max: int,             # hard cutoff (inclusive)
        loglik: bool = False    # log-pmf if true
) -> np.ndarray[float]:         # outputs vector of (log-)pmf values computed
    s, = params

    # compute the constant using log-sum-exp to avoid numeric over/underflow
    k = np.arange(1, x_max + 1)         # all possible values for x
    log_C = logsumexp(-s * np.log(k))   # sum over all values

    if loglik:
        return -log_C - s * np.log(x)
    else:
        return 1 / np.exp(log_C) * np.power(x, -s)


# 2. Exponentially Truncated Discrete Power Law Distribution
# p(x; tau, c, x_max) = 1 / C x^(-tau) e^(- x/c)
#   x = 1, 2, ... x_max
#   tau in (0, inf)
#   c in (0, inf)

def exp_power_law_pmf(
        params: tuple[float, float],    # parameters tau, c
        x: np.ndarray[int],             # data
        x_max: int,                     # hard cutoff (inclusive)
        loglik: bool = False            # log-pmf if true
) -> np.ndarray[float]:                 # outputs vector of (log-)pmf values computed
    tau, c = params

    # compute the constant using log-sum-exp to avoid numeric over/underflow
    k = np.arange(1, x_max + 1)
    log_C = logsumexp(-tau * np.log(k) - k/c)

    if loglik:
        return -log_C - tau * np.log(x) - x/c
    
    else: 
        return 1/np.exp(log_C) * np.power(x, -tau) * np.exp(-x/c)


# 3. Pooled Exp-Truncated Discrete Power Law Distribution
# p(x, L; tau, gamma, D, x_max) = 1 / C x^(-tau) e^(-x / (gamma L^D))
#   x = 1, 2, ..., x_max
#   tau in (0, inf)
#   gamma in (0, inf)
#   D in (0, inf)
#   L in 1, 2, ..., 

def pooled_exp_power_law_pmf(
        params: tuple[float, float],        # parameters tau, gamma, D
        data: dict[int, np.ndarray[int]],   # dictionary of (L, observations) to associate correct L with data
        x_max: dict[int, int],              # hard cutoff (depends on L)
        loglik: bool = False                # log-pmf if true
) -> np.ndarray[float]:                     # outputs MERGED vector of (log-)pmf values computed, length is sum of lengths of observations

    tau, gamma, D = params

    # use previous likelihood function and evaluate for all provided (L, x) pairs
    if loglik:
        lik = [exp_power_law_pmf(
                params=(tau, gamma * np.power(L, D)), 
                x=x, 
                x_max=x_max[L], 
                loglik=True) 
                for L, x in data.items()]

    else:
        lik = [exp_power_law_pmf(
                params=(tau, gamma * np.power(L, D)), 
                x=x, 
                x_max=x_max[L], 
                loglik=False) 
                for L, x in data.items()]

    return np.concatenate(lik)


# 4. Stretched Exponentially Truncated Discrete Power Law Distribution
# p(x; tau, c, x_max) = 1 / C x^(-tau) e^(- (x/c)^b)
#   x = 1, 2, ... x_max
#   tau in (0, inf)
#   c in (0, inf)
#   b in (0, inf)

def stretched_exp_power_law_pmf(
        params: tuple[float, float],    # parameters tau, c, b
        x: np.ndarray[int],             # data
        x_max: int,                     # hard cutoff (inclusive)
        loglik: bool = False            # log-pmf if true
) -> np.ndarray[float]:                 # outputs vector of (log-)pmf values computed
    tau, c, b = params

    # compute the constant using log-sum-exp to avoid numeric over/underflow
    k = np.arange(1, x_max + 1)
    log_C = logsumexp(-tau * np.log(k) - np.power(k/c, b))

    if loglik:
        return -log_C - tau * np.log(x) - np.power(x/c, b)
    
    else: 
        return 1/np.exp(log_C) * np.power(x, -tau) * np.exp(-np.power(x/c, b))



# 5. Pooled Exp-Truncated Discrete Power Law Distribution
# p(x, L; tau, gamma, D, beta, x_max) = 1 / C x^(-tau) e^(- (x / (gamma L^D))^beta)
#   x = 1, 2, ..., x_max
#   tau in (0, inf)
#   gamma in (0, inf)
#   D in (0, inf)
#   beta in (0, inf)
#   L in 1, 2, ..., 

def pooled_stretched_exp_power_law_pmf(
        params: tuple[float, float],        # parameters tau, gamma, D
        data: dict[int, np.ndarray[int]],   # dictionary of (L, observations) to associate correct L with data
        x_max: dict[int, int],              # hard cutoff (depends on L)
        loglik: bool = False                # log-pmf if true
) -> np.ndarray[float]:                     # outputs MERGED vector of (log-)pmf values computed, length is sum of lengths of observations

    tau, gamma, D, beta = params

    # use previous likelihood function and evaluate for all provided (L, x) pairs
    if loglik:
        lik = [stretched_exp_power_law_pmf(
                params=(tau, gamma * np.power(L, D), beta), 
                x=x, 
                x_max=x_max[L], 
                loglik=True) 
                for L, x in data.items()]

    else:
        lik = [stretched_exp_power_law_pmf(
                params=(tau, gamma * np.power(L, D), beta), 
                x=x, 
                x_max=x_max[L], 
                loglik=False) 
                for L, x in data.items()]

    return np.concatenate(lik)
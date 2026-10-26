import numpy as np
from collections.abc import Callable
import scipy.optimize

# ---
# Optimisation Procedure to compute the MLE
# The procedure used is the L-BFGS-B method implemented in scipy.optimise.minimise
# link: https://docs.scipy.org/doc/scipy/reference/optimize.minimize-lbfgsb.html#optimize-minimize-lbfgsb
#
# This method was originally introduced by Byrd et. al in 1995:
# Byrd, R. H., Lu, P., Nocedal, J., & Zhu, C. (1995). A limited memory algorithm for bound constrained optimization. 
# SIAM Journal on scientific computing, 16(5), 1190-1208.
# ---

def compute_MLE(
    pmf: Callable,                          # pmf function f(params, x, log=True)
    p0: tuple[float, ],                     # parameter initial guess
    bounds: tuple[tuple[float, float], ],   # parameter space
    scale: tuple[float] = None,             # rescale parameters to help optimisation
    pmf_kwargs: dict = {},                  # keyword arguments to pass to pmf
    **kwargs                                # other kwargs to pass to optimisation (e.g. method, gradient, etc)
) -> np.ndarray[float]:                     # returns a numpy array of estimated parameters

    # if no scale provided, don't modify it, just set the scale to 1
    if scale is None:
        scale = np.ones(len(p0))
    else:
        scale = np.array(scale)

    # minimise negative log-likelihood = maximise likelihood
    def neg_loglik(params: tuple[float, ]) -> float:
        params = np.array(params)
        scaled_params = params*scale
        return - np.sum(pmf(params=scaled_params, loglik=True, **pmf_kwargs))

    # run the optimisation procedure (default L-BFGS-B)
    optim = scipy.optimize.minimize(
        fun = neg_loglik,
        x0 = p0,
        bounds=bounds,
        **kwargs
    )
    
    # test for convergence
    assert optim.success, optim.message
    
    return optim.x * scale

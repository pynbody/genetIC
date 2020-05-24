import numpy as np

def powerlaw_covariance(k, plaw = -0.0, k0=0.5, k1=np.inf):
    """Calculate the Fourier-space covariance corresponding to a power-law power spectrum P(k) propto k
    
    Cut-offs are imposed at k<k0 and k>k1 where required. The cut-off at k0 is also smoothed."""

    try:
        cv = k ** plaw
    except ZeroDivisionError:
        return 0.0
    except ValueError:
        return 0.0

    if plaw>-1.0:
        k0 = 0

    if isinstance(cv, np.ndarray):
        cv[cv==np.inf]=0
        cv[cv!=cv]=0
        cv[k > k1] = 0
    else:
        if k>k1 : cv = 0.0

    cv*=1./(1. + k0 * k ** (plaw - 0.5))

    return cv
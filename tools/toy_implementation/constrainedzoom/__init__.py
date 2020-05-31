from importlib import reload

import pylab as p
import scipy.fftpack
import scipy.integrate
import scipy.interpolate

from constrainedzoom.methods import ZoomConstrained
from constrainedzoom.power_spectrum import powerlaw_covariance
from . import fft_wrapper
from . import methods

# auto-reload to make development easier
reload(fft_wrapper)
reload(methods)
#reload(plot)

from .fft_wrapper import *


def constraint_vector(scale=100,length=768,position=None) :
    """Generate a constraint vector corresponding to the Gaussian-filtered
    density at the given position."""
    if position is None :
        position = length//2

    pixel_vals = np.arange(0.,length)
    constraint = np.exp(-(pixel_vals-position)**2/(2*scale))
    constraint/=constraint.sum()
    return constraint.view(FFTArray)

def deriv_constraint_vector(smooth=None,length=768,position=None) :
    """Constraint derivative at given position. First arg is ignored
    but could later define smoothing scale"""

    """
    if position is None :
        position = length//2


    constraint = np.zeros(length)
    if smooth is None:
        smooth = 1

    constraint[position+smooth]=0.5
    constraint[position-smooth]=-0.5

    constraint/=np.sqrt(np.dot(constraint,constraint))
    """

    constraint = np.gradient(constraint_vector(smooth,length,position))
    return constraint


def cov2cor(matr):
    return matr/np.sqrt(np.outer(matr.diagonal(),matr.diagonal()))


def demo(val=2.0,seed=1,plaw=-1.5, deriv=False, nP=1024, nW=256, k_cut=0.2, scale=4, smooth=10):
    cv_gen = deriv_constraint_vector if deriv else constraint_vector
    cov_this = functools.partial(powerlaw_covariance, plaw=plaw)

    # set up zoom solution
    if seed is not None:
        np.random.seed(seed)
    G = ZoomConstrained(cov_this, k_cut=k_cut, nW=nW, nP=nP, hires_window_scale=scale, offset=(nP * (scale - 1)) // (2 * scale))

    G.add_constraint(val, cv_gen(smooth, nW))

    # set up ideal (uniform) solution
    nP_eff = nW*scale
    print("nP_eff=",nP_eff)

    Gs = ZoomConstrained(cov_this, k_cut=10000, nW=nP_eff, nP=nP_eff, hires_window_scale=1, offset=0)
    Gs.add_constraint(val, cv_gen(smooth, nP_eff))
    _, r_ideal = Gs.realization(no_random=True)
    x_ideal = (np.arange(nP_eff)+0.5)/(nP_eff/nP)

    # make plots
    p.clf()
    p.subplot(211)
    p.plot(x_ideal,r_ideal,linewidth=4,color="#ffdddd")
    x0, x1 = G.xs()
    r0, r1 = G.realization(no_random=True)
    p.plot(x0,r0,'k:')
    p.plot(x1,r1,'k-')
    p.ylabel("Solution")

    p.subplot(212)
    ideal_sln = scipy.interpolate.interp1d(x_ideal,r_ideal)
    x0_max = max(abs(ideal_sln(x0)))
    x1_max = max(abs(ideal_sln(x1)))
    p.plot(x0,(r0-ideal_sln(x0))/x0_max,'k:')
    p.plot(x1,(r1-ideal_sln(x1))/x1_max,'k')
    p.ylim(-0.05,0.05)
    p.ylabel("$\Delta$ solution / max |solution|")
    return G, Gs

"""
def cov_constraint_demo(downgrade_view=False,plaw=-1.5):
    cov_this = functools.partial(globals()['cov'],plaw=plaw)
    G = ZoomConstrained(cov_this,nW=256)
    #G.add_constraint(0.0,constraint_vector())
    #G.constraints_real.append(np.ones(768))
    cov, means, stds = G.estimate_cov(with_means=True)
    print("Mean of constraint:",means)
    print("Std-dev of constraints:",stds)
    plot_real_space(G, cov, downgrade_view)
"""


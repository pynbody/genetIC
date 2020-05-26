import pynbody
import numpy as np
import pylab as p
import copy
import math
import warnings

class GaussianLowPass(object):
    def __init__(self, k0):
        self.k0 = k0

    def __call__(self, k):
        return np.exp(-k**2/self.k0**2)


class Complementary(object):
    def __init__(self, underlying):
        self.underlying = underlying

    def __call__(self, k):
        return 1.-self.underlying(k)

class SharpHighPass(object):
    def __init__(self, k0):
        self.k0 = k0

    def __call__(self, k):
        return k>self.k0

class SharpLowPass(object):
    def __init__(self, k0):
        self.k0 = k0

    def __call__(self, k):
        return k<self.k0

class FermiLowPass(object):
    def __init__(self, k0, sharpness=10):
        self.k0 = k0
        self.temp = k0/sharpness

    def __call__(self, k):
        return 1. / (1. + np.exp((k - self.k0) / self.temp));

class VelocityFromDensity(object):
    def __init__(self, OmegaM0=0.3, z=15.7):
        self._OmegaM0 = OmegaM0
        self._z = z

    def __call__(self, k):
        E_z2 = self._OmegaM0 * (1+self._z)**3 + (1.-self._OmegaM0)# H(z)^2/H_0^2
        result = 1e4 * E_z2/(k**2 * (1+self._z)**2) # H(z)^2/((1+z)^2 k^2)
        result[np.isinf(result)] = 0
        return result

def subsampling_plot(x,y,*args,**kwargs):
    length = len(x)
    if length>2**10:
        subsampler = np.array(2**np.arange(0,math.log(length)/math.log(2),0.025),dtype=int)
        x = x[subsampler]
        y = y[subsampler]
    p.plot(x,y,*args,**kwargs)


class PowerSpectrum(object):
    def __init__(self, pspec=None, spacing_Mpc=0.01, npix=2 ** 19, z=0.0):
        self.spacing_Mpc = spacing_Mpc
        self.npix = npix
        self.k = 2 * np.pi * np.fft.fftfreq(npix, spacing_Mpc)
        self.z = z

        k_sanitized = copy.copy(self.k)
        k_sanitized[k_sanitized == 0] = 1
        if pspec is None:
            f = pynbody.new()
            f.properties['z'] = self.z
            pspec = pynbody.analysis.hmf.PowerSpectrumCAMBLive(f)
            print("sigma8 = %.2f"%pspec.get_sigma8())
        if not hasattr(pspec, "__call__"):
            self.Pk = 1e2*k_sanitized**pspec*(k_sanitized>1e-2)*(k_sanitized<1e2)
        else:
            self.Pk = pspec(abs(k_sanitized))
            self.Pk[self.k == 0] = 0



    @property
    def correlation(self):
        return CorrelationFunction(self)

    def filter(self, filt):
        return FilteredPowerSpectrum(self, filt)

    def boxsize(self, boxsize):
        return self.filter(SharpHighPass(2*np.pi/boxsize))

    def plot(self, scale=False):
        subsampling_plot(self.k[:self.npix//2], self.Pk[:self.npix // 2])
        p.loglog()
        if scale:
            p.ylim(2e2, 2e5)
            p.xlim(1e-3,1)

class FilteredPowerSpectrum(PowerSpectrum):
    def __init__(self, underlying, filt):
        self.spacing_Mpc = underlying.spacing_Mpc
        self.npix = underlying.npix
        self.k = underlying.k
        self.Pk = underlying.Pk*filt(abs(self.k))


class CorrelationFunction(object):
    def __init__(self, underlying):
        self.spacing_Mpc = underlying.spacing_Mpc
        self.npix = underlying.npix
        self.underlying = underlying
        self.k=underlying.k
        r = np.linspace(0, self.spacing_Mpc * self.npix // 2, self.npix // 2)
        self.r = np.concatenate([r, -r[::-1]])
        self.xi_filt = 1.0

    def __call__(self, value):
        offset_index = np.argmin(abs(self.r - value))
        offset = self.xi[offset_index]
        return offset

    @property
    def xi(self):
        return self._get_xi()

    @property
    def r_xi(self):
        result = np.fft.ifft(self.underlying.Pk * self.k, norm='ortho').imag
        delta_k = self.underlying.k[1] - self.underlying.k[0]
        result *= delta_k * np.sqrt(self.underlying.npix) / (2 * np.pi ** 2)
        return result

    def _get_xi(self):
        xi_r = self.r_xi / self.r
        xi_r[xi_r != xi_r] = 0
        xi_r[xi_r == np.inf] = 0
        return xi_r

    def plot(self, remove_dc=False, r2_weighted=True, r_offset=0, scale=False):
        xi_r = self.xi[:self.npix // 2]
        r = self.r[:self.npix // 2]
        if remove_dc:
            xi_r = xi_r-np.average(xi_r,weights=r**2)

        if r2_weighted:
            xi_r = xi_r*r**2
        subsampling_plot(r-r_offset, xi_r)

        if scale:
            p.xlim(0, 200)
            p.ylim(-100,100)

    @property
    def spectrum(self):
        return PowerSpectrumFromCorrelation(self)

    def filter(self, filt):
        return FilteredCorrelationFunction(self, filt)

    def boxsize(self, boxsize, dc_mode=False):
        filtered = self.filter(SharpLowPass(boxsize/2))

        if dc_mode:
            offset = OffsetCorrelationFunction(filtered, self(boxsize/2))
        else:
            offset = filtered

        return offset

class FilteredCorrelationFunction(CorrelationFunction):
    def __init__(self, underlying, filt):
        self._setup_from_underlying(underlying)
        self._r_xi_filtered = filt(abs(self.r))*underlying.r_xi
        self._xi_filtered = self._r_xi_filtered / self.r

    def _setup_from_underlying(self, underlying):
        self.spacing_Mpc = underlying.spacing_Mpc
        self.npix = underlying.npix
        self.r = underlying.r
        self.k = underlying.k

    def _get_xi(self):
        return self._xi_filtered

    @property
    def r_xi(self):
        return self._r_xi_filtered


class OffsetCorrelationFunction(FilteredCorrelationFunction):
    def __init__(self, underlying, offset):
        self._setup_from_underlying(underlying)
        self._xi_filtered = underlying.xi - offset




class PowerSpectrumFromCorrelation(PowerSpectrum):
    def __init__(self, underlying):
        self.spacing_Mpc = underlying.spacing_Mpc
        self.npix = underlying.npix
        self.underlying = underlying
        self.k = underlying.k
        self.Pk = self._get_Pk_from_r_xi(underlying.r_xi)
        self.Pk[self.Pk!=self.Pk]=0

    def _get_Pk_from_r_xi(self, r_xi):
        P_k = -np.fft.fft(r_xi, norm='ortho').imag / self.k
        delta_k = self.k[1] - self.k[0]
        P_k /= delta_k * np.sqrt(self.npix) / (2 * np.pi ** 2)
        return P_k


class FilterExplorer(object):
    def __init__(self, pspec=None, spacing_Mpc=0.1, npix=2**16):
        if pspec is None:
            pspec = pynbody.analysis.hmf.PowerSpectrumCAMB(pynbody.new())
        self.npix = npix
        self.k = 2*np.pi*np.fft.fftfreq(npix, spacing_Mpc)
        r = np.linspace(0,spacing_Mpc*npix/2,npix//2)
        self.r = np.concatenate([r,-r[::-1]])
        k_sanitized = copy.copy(self.k)
        k_sanitized[k_sanitized==0] = 1
        self._Pk = pspec(abs(k_sanitized))
        self._Pk[self.k == 0]=0
        self.filt = 1.0
        self.r_filt=1.0

    def plot_spectrum(self):
        k = self.k[:self.npix//2]
        p.plot(k, self.Pk[:self.npix // 2])
        p.loglog()

    def plot_spectral_index(self):
        k = self.k[:self.npix // 2]
        Pk = self.Pk[:self.npix // 2]
        log_gradient = np.diff(np.log(Pk))/np.diff(np.log(k))
        p.plot((k[:-1]+k[1:])/2, log_gradient)
        p.semilogx()

    def plot_roundtrip_spectrum(self):
        k = self.k[:self.npix // 2]
        Pk = self.get_Pk_from_xi(self.get_xi()*self.r_filt)


    def set_k_filter(self, filt):
        self.filt=filt(abs(self.k))

    def set_r_filter(self, filt):
        self.r_filt = filt(abs(self.r))

    @property
    def Pk(self):
        return self._Pk*self.filt

    def plot_xi(self):
        xi_r = (self.r_filt*self.get_xi())[:self.npix//2]
        r = self.r[:self.npix//2]
        p.plot(r, xi_r*r**2)
        p.xlim(0,200)

    def get_xi(self):
        xi_r = np.fft.ifft(self.Pk * self.k, norm='ortho').imag / self.r
        delta_k = self.k[1]-self.k[0]
        xi_r *= delta_k * np.sqrt(self.npix) / (2 * np.pi**2)
        xi_r[xi_r!=xi_r]=0
        xi_r[xi_r==np.inf]=0
        return xi_r

    def get_Pk_from_xi(self, xi):
        P_k = -np.fft.fft(xi*self.r, norm='ortho').imag / self.k
        delta_k = self.k[1] - self.k[0]
        P_k /= delta_k * np.sqrt(self.npix) / (2 * np.pi ** 2)
        return P_k



def plot_filter_effect(filt, filter_real_space=False,
                       plot_real_space=True, ns=None, coarse_grid_space=0.1):
    FRACTIONAL_K_SPLIT = 0.3



    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        if ns:
            ps = PowerSpectrum(ns,spacing_Mpc=5e-4,npix=2**21)
        else:
            ps = PowerSpectrum(spacing_Mpc=5e-4, npix=2 ** 21)

        if filter_real_space:
            filt = filt(coarse_grid_space / FRACTIONAL_K_SPLIT)
            f1 = ps.correlation.filter(filt)
            f2 = ps.correlation.filter(Complementary(filt))
        else:
            filt = filt(FRACTIONAL_K_SPLIT * 2.0 * np.pi / coarse_grid_space)
            f1 = ps.filter(filt).correlation
            f2 = ps.filter(Complementary(filt)).correlation


        subsampler = np.array(2**np.arange(0,21,0.025),dtype=int)

        p.clf()
        p.subplot(211)
        weight = 1./ps.correlation.xi #f1.r**2
        p.plot(f1.r[subsampler],(abs(f1.xi)*weight)[subsampler],"k--")
        p.plot(f1.r[subsampler],(abs(f2.xi)*weight)[subsampler],"k:")
        p.plot(f1.r[subsampler],abs(ps.correlation.xi*weight)[subsampler],"r")
        p.xlim(0.001,1000.0)
        p.loglog()

        p.xlim(1e-3 * coarse_grid_space / 0.2, 10 * coarse_grid_space / 0.2)

        p.ylim(1e-3,4.0)

        p.gca().xaxis.tick_top()
        p.gca().xaxis.set_label_position('top')
        p.xlabel("$r/Mpc\, h^{-1}$")
        p.ylabel(r"$\xi_i(r) / \xi(r)$")

        p.subplot(212, sharex=p.gca())
        f1 = f1.spectrum
        f2 = f2.spectrum
        k_recip = 2.0*np.pi/f1.k
        p.plot(k_recip[subsampler], (f1.Pk/ps.Pk)[subsampler], "k--")
        p.plot(k_recip[subsampler], (f2.Pk/ps.Pk)[subsampler], "k:")
        p.plot(k_recip[subsampler], (ps.Pk/ps.Pk)[subsampler], "r")

        p.xlim(1e-3 * coarse_grid_space / 0.2, 10 * coarse_grid_space / 0.2)
        p.ylim(1e-3,4)
        p.loglog()

        p.axvline(coarse_grid_space, color=(0.8, 0.8, 0.8), linewidth=2)
        p.axvline(2 * coarse_grid_space, color=(0.85, 0.85, 0.85), linewidth=2)
        p.axvline(3 * coarse_grid_space, color=(0.9, 0.9, 0.9), linewidth=2)
        p.axvline(4 * coarse_grid_space, color=(0.95, 0.95, 0.95), linewidth=2)

        p.xlabel(r"$2\pi/k/Mpc\, h^{-1}$")
        p.ylabel("$P_i(k)/P(k)$")



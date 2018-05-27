"""Module for exploring the power spectrum covariance of fixed fields (as opposed to Gaussian fields)"""

import fft_wrapper
import numpy as np
import numpy.fft
import pylab as p

class EitherField(object):
    def __init__(self, pspec=-2, spacing_Mpc=0.01, npix=1024, npix_subbox=512):
        self.spacing_Mpc = spacing_Mpc
        self.npix = npix
        self.npix_subbox = npix_subbox
        self.k = 2.0*np.pi*np.fft.fftfreq(npix, spacing_Mpc)
        self.Pk = self.get_power_spectrum(pspec, self.k)
        self.k_subbox = 2.0*np.pi*np.fft.fftfreq(npix_subbox, spacing_Mpc)
        self.Pk_subbox = self.get_power_spectrum(pspec, self.k_subbox)

    def get_power_spectrum(self, pspec, k):
        k_sanitized = np.abs(k)
        k_sanitized[k_sanitized == 0] = 1
        Pk = 1e2 * k_sanitized ** pspec
        Pk[k == 0] = 0

        return Pk

    def apply_window(self, vec: fft_wrapper.FFTArray):
        return vec.in_real_space()[:self.npix_subbox].in_fourier_space()

    def subbox_covariance(self):
        """Returns <delta_i delta_j> for the subbox"""
        covariance = np.zeros((self.npix_subbox, self.npix_subbox))
        for i in range(self.npix):
            cov_element = np.zeros(self.npix).view(fft_wrapper.FFTArray)
            cov_element.fourier = True
            cov_element[i]=np.sqrt(self.Pk[i])
            cov_element_subbox = self.apply_window(cov_element)
            covariance+=np.outer(cov_element_subbox, cov_element_subbox.conj())
        return covariance

    def expectation_field4(self):
        """Returns <delta_i delta_i delta_j delta_j> for the subbox"""
        raise NotImplementedError("To use this method, you need to instantiate either a GaussianField or a FixedField")

    def plot_subbox_power(self):
        """Plot the ensemble average power in the subbox.

        Overplot (as a dashed line) what it would be if the subbox contained its own realization of a 
        random field with the same theory power spectrum
        """
        restrict = slice(1, self.npix_subbox//2)
        Pk_empirical = self.subbox_covariance().diagonal()

        p.plot(self.k_subbox[restrict], Pk_empirical[restrict])
        p.plot(self.k_subbox[restrict], self.Pk_subbox[restrict], '--')
        p.loglog()

    def plot_subbox_power_std(self):
        """Plot the standard deviation of the power in the subbox.

        Overplot (as a dashed line) what it would be if the subbox contained its own realization of a Gaussian
        random field with the same theory power spectrum (NB this idealised comparison line is Gaussian even if 
        the parent field is Fixed)."""
        power = self.subbox_covariance().diagonal()
        power_var = np.sqrt(self.expectation_field4().diagonal() - power**2)
        p.plot(self.k_subbox[1:self.npix_subbox//2], power_var[1:self.npix_subbox//2])
        p.plot(self.k_subbox[1:self.npix_subbox//2], np.sqrt(2)*self.Pk_subbox[1:self.npix_subbox//2], '--')
        p.loglog()

    def plot_subbox_correlation(self):
        p.clf()
        cov = self.subbox_covariance()[:self.npix_subbox//2, :self.npix_subbox//2]
        pow = cov.diagonal()
        cor = cov/np.sqrt(np.outer(pow,pow))
        p.imshow(cor)


class GaussianField(EitherField):
    def expectation_field4(self):
        """Returns <delta_i delta_i delta_j delta_j> for the subbox"""
        subbox_cov = self.subbox_covariance()
        subbox_pow = subbox_cov.diagonal()
        field4 = 2*subbox_cov**2 + np.outer(subbox_pow, subbox_pow)
        return field4

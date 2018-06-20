"""Module for exploring the power spectrum covariance of fixed fields (as opposed to Gaussian fields)"""

import fft_wrapper
import numpy as np
import scipy.linalg
import pylab as p

from functools import lru_cache

class EitherField(object):
    def __init__(self, pspec=-2, spacing_Mpc=0.01, npix=1024, npix_subbox=512, 
                 subbox_policy="small", bandwidth=2):
        """

        *params*
         pspec - the power law of the power spectrum
         spacing_Mpc - a notional spacing for the pixels; normalises the k values
         npix - the number of pixels in the big box
         npix_subbox - the number of pixels in the subbox
         subbox_policy - either 'small' or 'zeroed'
                         small implies that subbox values are stored in a npix_subbox vector
                         zero implies that subbox values are setored in a npix vector with zero padding
         bandwidth - the number of Fourier big-box pixels in a band for band-power estimates
        """
        self.spacing_Mpc = spacing_Mpc
        self.npix = npix
        if subbox_policy=="small":
            self.npix_subbox_vector = npix_subbox
        elif subbox_policy=="zeroed":
            self.npix_subbox_vector = npix
        self.npix_subbox_mask = npix_subbox
        self.bandwidth = bandwidth
        self.nbands = npix//(self.bandwidth*2) # number of bandpowers to estimate
        self.k = 2.0*np.pi*np.fft.fftfreq(npix, spacing_Mpc)
        self.Pk = self.get_power_spectrum(pspec, self.k)

        self.k_subbox = 2.0*np.pi*np.fft.fftfreq(self.npix_subbox_vector, spacing_Mpc)
 
        self.k_bandpower = self.k[self.bandwidth//2:self.npix//2:self.bandwidth]
        self.Pk_bandpower = self.Pk[self.bandwidth//2:self.npix//2:self.bandwidth]
        assert len(self.k_bandpower)==self.nbands

    def get_power_spectrum(self, pspec, k):
        k_sanitized = np.abs(k)
        k_sanitized[k_sanitized == 0] = 1
        Pk = 1e2 * k_sanitized ** pspec
        Pk[k == 0] = 0

        return Pk

    def apply_window(self, vec: fft_wrapper.FFTArray):
        vec.in_real_space()
        vec[self.npix_subbox_mask:]=0
        vec = vec[:self.npix_subbox_vector]
        vec = vec.in_fourier_space()
        assert len(vec)==self.npix_subbox_vector
        return vec

    @lru_cache()
    def subbox_covariance_convolution_matrix(self):
        matrix = np.zeros((self.npix_subbox_vector, self.nbands))
        for band in range(self.nbands):
            
            i_big_box_min = band*self.bandwidth
            i_big_box_max = (band+1)*self.bandwidth
            
            WXW = np.zeros((self.npix_subbox_vector ,self.npix_subbox_vector ))
            for i_big_box in range(i_big_box_min, i_big_box_max):
                cov_element : fft_wrapper.FFTArray = np.zeros(self.npix).view(fft_wrapper.FFTArray) 
                cov_element.fourier = True
                cov_element[i_big_box] = 1.0
                cov_element = self.apply_window(cov_element)
                WXW += np.outer(cov_element, cov_element.conj())*self.Pk[i_big_box]
            if self.Pk_bandpower[band]>0:
                WXW/=self.Pk_bandpower[band]
            matrix[:,band] += WXW.diagonal()
        return matrix

    @lru_cache()
    def subbox_covariance_deconvolution_matrix(self):
        return np.linalg.pinv(self.subbox_covariance_convolution_matrix(),rcond=1e-5)

    @lru_cache()
    def subbox_covariance(self):
        """Returns <delta_i delta_j> for the subbox"""
        covariance = np.zeros((self.npix_subbox_vector, self.npix_subbox_vector))
        for i in range(self.npix):
            cov_element = np.zeros(self.npix).view(fft_wrapper.FFTArray)
            cov_element.fourier = True
            cov_element[i]=np.sqrt(self.Pk[i])
            cov_element_subbox = self.apply_window(cov_element)
            covariance+=np.outer(cov_element_subbox, cov_element_subbox.conj())
        return covariance

    @lru_cache()
    def expectation_field4(self):
        """Returns <delta_i delta_i delta_j delta_j> for the subbox"""
        raise NotImplementedError("To use this method, you need to instantiate either a GaussianField or a FixedField")

    def plot_subbox_power(self, deconvolved=True):
        """Plot the ensemble average power in the subbox.

        Overplot (as a dashed line) what it would be if the subbox contained its own realization of a 
        random field with the same theory power spectrum
        """
        restrict = slice(1, self.npix_subbox_vector // 2)

        if deconvolved:
            matr = self.subbox_covariance_deconvolution_matrix()
            Pk_empirical = np.dot(matr, self.subbox_covariance().diagonal())
            Pk_expected = self.Pk_bandpower
            k_expected = self.k_bandpower
        else:
            Pk_empirical = self.subbox_covariance().diagonal()[:self.npix_subbox_vector//2]
            matr = self.subbox_covariance_convolution_matrix()
            Pk_expected = np.dot(matr, self.Pk_bandpower)[:self.npix_subbox_vector//2]
            k_expected = self.k_subbox[:self.npix_subbox_vector//2]

        p.plot(k_expected, Pk_empirical)
        p.plot(k_expected, Pk_expected, '--')
        p.loglog()

    def subbox_power_covar(self, deconvolved=True):
        cov = self.subbox_covariance()
        field4 = self.expectation_field4()
        power_covar = field4 - cov**2
        if deconvolved:
            decon = self.subbox_covariance_deconvolution_matrix()
            return np.dot(decon, np.dot(power_covar, decon.T))
        else:
            return power_covar

    def subbox_power_var(self, deconvolved=True):
        return self.subbox_power_covar(deconvolved).diagonal()

    def plot_subbox_power_std(self, deconvolved=True):
        """Plot the standard deviation of the power in the subbox.

        Overplot (as a dashed line) what it would be if the subbox contained its own realization of a Gaussian
        random field with the same theory power spectrum (NB this idealised comparison line is Gaussian even if 
        the parent field is Fixed)."""
        power_var = self.subbox_power_var(deconvolved)
        if deconvolved:
            k = self.k_bandpower
            Pk_expected = self.Pk_bandpower
        else:
            k = self.k_subbox[:self.npix_subbox_vector//2]
            matr = self.subbox_covariance_convolution_matrix()
            Pk_expected = np.dot(matr, self.Pk_bandpower)[:self.npix_subbox_vector//2]

        p.plot(k, np.sqrt(power_var))
        p.plot(k, np.sqrt(2) * Pk_expected, '--')
        p.loglog()

    def plot_subbox_correlation(self):
        p.clf()
        cov = self.subbox_covariance()[:self.npix_subbox_vector // 2, :self.npix_subbox_vector // 2]
        pow = cov.diagonal()
        cor = cov/np.sqrt(np.outer(pow,pow))
        p.imshow(cor)

    def plot_subbox_power_correlation(self):
        power = self.subbox_covariance().diagonal()
        cov = self.expectation_field4() - np.outer(power,power)
        std = power**2
        cor = cov/np.sqrt(np.outer(std,std))
        p.imshow(cor,vmin=-1.0,vmax=1)
        p.set_cmap('RdBu')
        p.colorbar()


class GaussianField(EitherField):
    def expectation_field4(self):
        """Returns <delta_i delta_i delta_j delta_j> for the subbox"""
        subbox_cov = self.subbox_covariance()
        subbox_pow = subbox_cov.diagonal()
        field4 = 2*subbox_cov**2 + np.outer(subbox_pow, subbox_pow)
        return field4

class FixedField(EitherField):
    def expectation_field4_diagonal_correction(self):
        """In the full box:

        <delta_i^4> = 3<delta_i^2>^2 for Gaussian
        <delta_i^4> = <delta_i>^2 for fixed

        Correction: -2<delta_i>^2

        Correction to subbox <P_i P_j>

        sum_k <Tr P_i W delta_k delta_k W P_j W delta_k delta_k W>
         -> sum_k Tr P_i W C_kk W P_j W C_kk W

        """
        covariance = np.zeros((self.npix_subbox_vector, self.npix_subbox_vector))
        for k in range(self.npix):
            cov_element = np.zeros(self.npix).view(fft_wrapper.FFTArray)
            cov_element.fourier = True
            cov_element[k] = np.sqrt(self.Pk[k])
            cov_element_subbox = self.apply_window(cov_element)
            WC_kkW = np.outer(cov_element_subbox, cov_element_subbox.conj())
            covariance+=WC_kkW**2
        return -2*covariance

    def expectation_field4(self):
        gaussian = GaussianField.expectation_field4(self)
        return gaussian+self.expectation_field4_diagonal_correction()



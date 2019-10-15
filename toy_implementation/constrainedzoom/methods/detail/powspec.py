import numpy as np
from functools import partial
import scipy.integrate, scipy.fftpack
from ...fft_wrapper import unitary_fft, unitary_inverse_fft, FFTArray, in_fourier_space

class Powspec:
    def __init__(self, cov_fn, nP, nW, hires_window_scale, offset):
        if cov_fn is None:
            from ... import powerlaw_covariance
            cov_fn = partial(powerlaw_covariance, plaw=-1.0)
        self.cov_fn = cov_fn

        self.k_low = scipy.fftpack.rfftfreq(self.nP, d=self.delta_low)
        self.k_high = scipy.fftpack.rfftfreq(self.nW, d=self.delta_high)

        self.C_low = self._get_variance_k(self.k_low)
        self.C_high = self._get_variance_k(self.k_high)

        # covariance matrix for potential
        with np.errstate(divide='ignore'):
            self.C_low_potential = self._C_delta_to_potential(self.C_low, self.k_low)
            self.C_high_potential = self.pixel_size_ratio * self._C_delta_to_potential(self.C_high, self.k_high)

    def _C_delta_to_potential(self, C, k):
        with np.errstate(divide='ignore'):
            Cpot = C/k**4
            Cpot[k==0] = 0
        return Cpot

    def _get_variance_k(self, k, k_filter=lambda x:1):
        """This routine converts a continuous power spectrum [cov_fn] into a discrete version.

        Note that this is only technically possible if the power is band-limited, i.e. there cannot be any
        power on or below the scale of an individual cell. In practice, that requirement is routinely violated
        but then there is no principled way to make the mapping.

        To understand how the mapping is made in the case where the power is band-limited, consider the variance of
        the field at a single point. From the continuous power spectrum (adopting unitary FT conventions), one gets

        sigma^2 = (2 pi)^(-1) \int dk P(k).

        From a discrete power spectrum, P_i, one extracts the diagonal part in real space by writing e.g.

        sigma^2 = (N_cells)^(-1) Tr C

        where C is the covariance matrix in real space. Again using unitary FFTs, this is just

        sigma^2 = (N_cells)^(-1) Sum_i P_i

        Thus (ignoring proportionalities which we don't care about here) we come to the mapping

        P_i \propto N_cells \int_{k_i - delta_k/2}^{k_i + delta_k/2} P(k)
        """
        integrand = lambda k1: self.cov_fn(k1)*k_filter(k1)
        delta_k = k[1]-k[0]
        Cv = np.zeros_like(k)
        for i,ki in enumerate(k):
            if ki==0:
                lower_k = 0.0
            else:
                lower_k = ki-delta_k/2
            upper_k = ki + delta_k / 2
            Cv[i] = scipy.integrate.quad(integrand, lower_k, upper_k)[0]

        return Cv * len(k)

    def set_Chigh_realspace(self):
        self.C_high = self._calc_transfer_fn_realspace(self.nW)
        self.C_high_potential = self._calc_transfer_fn_realspace(self.nW, potential=True)

    def _calc_transfer_fn_realspace(self, num_pixels, potential=False):
        fullbox_nW = self.nP * self.pixel_size_ratio
        k_high_full = scipy.fftpack.rfftfreq(fullbox_nW, d=1. / fullbox_nW)

        # pretend high-res is across the full box temporarily; get the TF

        C_high_full = self._get_variance_k(k_high_full)
        if potential:
            C_high_full = self._C_delta_to_potential(C_high_full, k_high_full)
        transfer = self._cov_to_transfer(C_high_full)
        # now truncate the TF to the correct size
        transfer_hi = np.concatenate((transfer[:num_pixels // 2], transfer[-num_pixels // 2:]))

        return self._transfer_to_cov(transfer_hi)

    def _cov_to_transfer(self, cov):
        # take real part of C_high (i.e. zero imaginary components)
        sqrt_C_high_re = np.sqrt(cov)
        sqrt_C_high_re[2::2] = 0

        # use this to calculate the real-space transfer function
        T_high = unitary_inverse_fft(sqrt_C_high_re)/np.sqrt(len(cov))

        return T_high

    def _transfer_to_cov(self, transfer):
        sqrt_C_high_apo_re = unitary_fft(transfer)* np.sqrt(len(transfer))

        # copy back imaginary parts ready for convolution
        if len(sqrt_C_high_apo_re) % 2 == 0:
            sqrt_C_high_apo_re[2::2] = sqrt_C_high_apo_re[1:-1:2]
        else:
            sqrt_C_high_apo_re[2::2] = sqrt_C_high_apo_re[1::2]

        return sqrt_C_high_apo_re ** 2

    @in_fourier_space
    def _apply_transfer_function(self, noiseP, noiseW=None):
        """Apply the transfer function to the pixelized and windowed noise

        :param noiseP: the pixelized whitenoise
        :param noiseW: the windowed whitenoise
        :return: deltaP, deltaW: the convolved fields
        """
        deltaP = FFTArray(noiseP * np.sqrt(self.C_low))
        deltaP.fourier = True
        if noiseW is not None:
            deltaW = FFTArray(noiseW * np.sqrt(self.C_high))
            deltaW.fourier = True
            return deltaP, deltaW
        else:
            return deltaP
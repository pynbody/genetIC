import numpy as np
from functools import partial
import scipy.integrate
from ...fft_wrapper import unitary_fft, unitary_inverse_fft

class VarianceCalculationTools:
    def __init__(self, cov_fn, n1, n2, hires_window_scale, offset):
        if cov_fn is None:
            from ... import powerlaw_covariance
            cov_fn = partial(powerlaw_covariance, plaw=-1.0)
        self.cov_fn = cov_fn

    def _C_delta_to_potential(self, C, k):
        with np.errstate(divide='ignore'):
            Cpot = C/k**4
            Cpot[k==0] = 0
        return Cpot

    def _get_variance_k(self, k, k_filter=lambda x:1):
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

        return Cv

    def set_Chigh_realspace(self):
        self.C_high = self._calc_transfer_fn_realspace(self.n2)
        self.C_high_potential = self._calc_transfer_fn_realspace(self.n2, potential=True)

    def _calc_transfer_fn_realspace(self, num_pixels, potential=False):
        fullbox_n2 = self.n1 * self.pixel_size_ratio
        k_high_full = scipy.fftpack.rfftfreq(fullbox_n2, d=1. / fullbox_n2)

        # pretend high-res is across the full box temporarily; get the TF

        C_high_full = self._get_variance_k(k_high_full) * self.n1
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
import copy

import scipy.fftpack
from Cython.Includes.numpy import __init__

from constrainedzoom import UnfilteredZoomConstrained, in_real_space, FFTArray


class MLZoomConstrained(UnfilteredZoomConstrained):
    def __init__(self, *args, **kwargs):
        self._approximations = kwargs.pop('approximations',[])
        super().__init__(*args, **kwargs)
        self.n_underlying = self.n2 * self.window_size_ratio
        self.k_underlying = scipy.fftpack.rfftfreq(self.n_underlying, d=self.delta_high)
        self.C_underlying = self._get_variance_k(self.k_underlying) * float(self.n2)

        self.k_high2 = scipy.fftpack.rfftfreq(self.n2*2, d=self.delta_high)
        self.C_high2 = self._calc_transfer_fn_realspace(self.n2*2)/4

    def realization(self, verbose=False, no_random=False, white_noise_lo=None, white_noise_hi=None):
        white_noise_lo, white_noise_hi = self._get_whitenoise(white_noise_lo, white_noise_hi)
        white_noise_lo, white_noise_hi = self._separate_whitenoise(white_noise_lo, white_noise_hi)

        delta_low, delta_high = self._apply_transfer_function(white_noise_lo, white_noise_hi)

        return delta_low.in_real_space(), delta_high.in_real_space()

    @in_real_space
    def _separate_whitenoise(self, delta_low, delta_high):
        assert not delta_high.fourier
        assert not delta_low.fourier

        delta_high -= self.upsample_zeroorder(self.downsample(delta_high))

        lo_modes_for_hi_window = self.upsample_zeroorder(delta_low)

        delta_high += lo_modes_for_hi_window

        return delta_low, delta_high

    def add_constraint(self, val=0.0, hr_vec=None):
        raise RuntimeError("MLZoomConstrained does not support constraints")

    @in_real_space
    def _apply_transfer_function(self, white_noise_lo, white_noise_hi):

        C = self.C_underlying

        approx_i = 'i' in self._approximations
        approx_ii = 'ii' in self._approximations
        approx_iii = 'iii' in self._approximations
        approx_iv = 'iv' in self._approximations

        if approx_i:
            term_i : FFTArray = np.zeros(self.n2*2).view(FFTArray)
            term_i.fourier=False
            term_i[self.n2//2:3*self.n2//2] = white_noise_hi
            term_i.in_fourier_space()
            term_i*=self.C_high2**0.5
            term_i.in_real_space()
            term_i = term_i[self.n2//2:3*self.n2//2]


        else:
            term_i = self.place_window(white_noise_hi).in_fourier_space()
            term_i*=C**0.5
            term_i = self.extract_window(term_i.in_real_space())

        if approx_ii:
            term_ii = copy.copy(white_noise_lo).in_fourier_space()
            term_ii*=self.C_low**0.5
            term_ii.in_real_space()
        else:
            term_ii = self.upsample_zeroorder(white_noise_lo, in_window=False).in_fourier_space()
            term_ii*=C**0.5
            term_ii = self.downsample(term_ii,in_window=False)

        if approx_iii:
            term_iii = copy.copy(white_noise_lo)
            term_iii[self.offset:self.offset + self.n1 // self.window_size_ratio] = 0
            term_iii.in_fourier_space()
            term_iii*=self.C_low**0.5
            term_iii=self.upsample_cubic(term_iii.in_real_space())

        else:
            term_iii = copy.copy(white_noise_lo)
            term_iii[self.offset:self.offset + self.n1 // self.window_size_ratio] = 0
            term_iii = self.upsample_zeroorder(term_iii, in_window=False).in_fourier_space()
            term_iii *= C ** 0.5
            term_iii = self.extract_window(term_iii.in_real_space())


        if approx_iv:
            term_iv = 0
        else:
            term_iv = copy.copy(white_noise_hi)
            term_iv -= self.upsample_zeroorder(self.downsample(term_iv))
            term_iv = self.place_window(term_iv).in_fourier_space()
            term_iv*=C**0.5
            term_iv = self.downsample(term_iv.in_real_space(), in_window=False)

        return term_ii+term_iv, term_i+term_iii
import copy

import scipy.fftpack
import numpy as np

from ..fft_wrapper import in_real_space, FFTArray
from ..methods import UnfilteredZoomConstrained


class MLZoomConstrained(UnfilteredZoomConstrained):
    description = "Maximum Likelihood"

    def __init__(self, *args, **kwargs):
        self._approximations = kwargs.pop('approximations',[])
        super().__init__(*args, **kwargs)
        self.n_underlying = self.nW * self.window_size_ratio
        self.k_underlying = scipy.fftpack.rfftfreq(self.n_underlying, d=self.delta_high)
        self.C_underlying = self._get_variance_k(self.k_underlying)

    def get_default_plot_padding(self):
        return 0

    def realization(self, verbose=False, no_random=False, noiseP=None, noiseW=None):
        noiseP, noiseW = self._get_whitenoise(noiseP, noiseW)
        noiseP, noiseW = self._modify_whitenoise(noiseP, noiseW)

        delta_low, delta_high = self._apply_transfer_function(noiseP, noiseW)

        return delta_low.in_real_space(), delta_high.in_real_space()

    @in_real_space
    def _modify_whitenoise(self, delta_low, delta_high):
        assert not delta_high.fourier
        assert not delta_low.fourier

        delta_high -= self.upsample_zeroorder(self.downsample(delta_high))

        lo_modes_for_hi_window = self.upsample_zeroorder(delta_low) / np.sqrt(self.pixel_size_ratio)

        delta_high += lo_modes_for_hi_window

        return delta_low, delta_high

    def add_constraint(self, val=0.0, hr_vec=None, potential=False):
        raise RuntimeError("MLZoomConstrained does not support constraints")

    @in_real_space
    def _apply_transfer_function(self, noiseP, noiseW=None):

        C = self.C_underlying

        approx_i = 'i' in self._approximations
        approx_ii = 'ii' in self._approximations
        approx_iii = 'iii' in self._approximations
        approx_iv = 'iv' in self._approximations

        if approx_i:
            term_i : FFTArray = np.zeros(self.nW*2).view(FFTArray)
            term_i.fourier=False
            term_i[self.nW//2:3*self.nW//2] = noiseW
            term_i.in_fourier_space()
            term_i*=self.C_high2**0.5
            term_i.in_real_space()
            term_i = term_i[self.nW//2:3*self.nW//2]


        else:
            term_i = self.place_window(noiseW).in_fourier_space()
            term_i*=C**0.5
            term_i = self.extract_window(term_i.in_real_space())

        if approx_ii:
            term_ii = copy.copy(noiseP).in_fourier_space()
            term_ii*=self.C_low**0.5
            term_ii.in_real_space()
        else:
            term_ii = self.upsample_zeroorder(noiseP, in_window=False).in_fourier_space()
            term_ii*=C**0.5
            term_ii = self.downsample(term_ii, input_unpadded=False)
            term_ii/=self.pixel_size_ratio**0.5

        if approx_iii:
            term_iii = copy.copy(noiseP)
            term_iii[self.offset:self.offset + self.nP // self.window_size_ratio] = 0
            term_iii.in_fourier_space()
            term_iii*=self.C_low**0.5
            term_iii=self.upsample_cubic(term_iii.in_real_space())

        else:
            term_iii = copy.copy(noiseP)
            term_iii[self.offset:self.offset + self.nP // self.window_size_ratio] = 0
            term_iii = self.upsample_zeroorder(term_iii, in_window=False).in_fourier_space()
            term_iii *= C ** 0.5
            term_iii = self.extract_window(term_iii.in_real_space())
            term_iii /= self.pixel_size_ratio ** 0.5


        if approx_iv:
            term_iv = 0
        else:
            term_iv = copy.copy(noiseW)
            term_iv -= self.upsample_zeroorder(self.downsample(term_iv))
            term_iv = self.place_window(term_iv).in_fourier_space()
            term_iv*=C**0.5
            term_iv = self.downsample(term_iv.in_real_space(), input_unpadded=False)

        return term_ii+term_iv, term_i+term_iii
import copy
import numpy as np

from . import ZoomConstrained, UnfilteredZoomConstrained
from ..fft_wrapper import FFTArray
from ..power_spectrum import powerlaw_covariance

class IdealizedZoomConstrained(ZoomConstrained):
    """Calculate the low-res/high-res split by making a full box at the high resolution,
    then downgrading the resolution of the low-res region"""

    def __init__(self, cov_fn=powerlaw_covariance, n1=256, n2=768, hires_window_scale=4, offset=10):
        super().__init__(cov_fn, n1, n2, hires_window_scale, offset)
        self._underlying = UnfilteredZoomConstrained(cov_fn, n2*hires_window_scale, n2, hires_window_scale, offset)

    def _iter_cov_elements(self):
        test_field_lo = np.zeros(self._underlying.n1)
        test_field_hi = np.zeros(self._underlying.n2)

        for i in range(self._underlying.n1):
            test_field_lo[i] = 1.0
            yield test_field_lo, test_field_hi
            test_field_lo[i] = 0.0

    def _iter_one_cov_element(self, offset):
        test_field_lo = np.zeros(self._underlying.n1)
        test_field_hi = np.zeros(self._underlying.n2)

        if offset>=0:
            # place test delta function in high-res region
            test_field_lo[self.offset*self.pixel_size_ratio + self.n2 // 2 + self.pixel_size_ratio//2 + offset] = 1.0
        else:
            # in low-res region
            test_field_lo[(self.offset+offset)*self.pixel_size_ratio:(self.offset+offset+1)*self.pixel_size_ratio] = 1.0

        yield FFTArray(test_field_lo).in_fourier_space(), test_field_hi


    def realization(self, *args, **kwargs):
        if 'underlying_realization' in kwargs:
            underlying = kwargs.pop('underlying_realization')
        else:
            underlying,_ = self._underlying.realization(*args, **kwargs)

        underlying_lores = copy.copy(underlying)

        underlying_lores.in_fourier_space()
        underlying_lores[self.n1:]=0
        underlying_lores.in_real_space()

        lores = underlying_lores.reshape((self.n1,self.pixel_size_ratio)).mean(axis=1)
        hires = underlying[self.offset*self.pixel_size_ratio:self.offset*self.pixel_size_ratio+self.n2]
        return lores, hires

    def _apply_constraints(self, delta_low_k, delta_high_k, verbose):
        return delta_low_k, delta_high_k # constraints are applied in underlying object

    def _modify_whitenoise(self, wn_lo, wn_hi):
        return wn_lo, wn_hi

    def _separate_fields(self, delta_low_k, delta_high_k):
        return delta_low_k, delta_high_k, None

    def _recombine_fields(self, delta_low, delta_high, memos):
        return delta_low, delta_high





class FastIdealizedZoomConstrained(IdealizedZoomConstrained):

     def get_cov(self, one_element=False):
            # the base class implementation does work, but this is optimized for the idealized case

            if one_element:
                return super().get_cov(one_element)

            test_field_lo : FFTArray = np.zeros(self._underlying.n1).view(FFTArray)
            test_field_hi = np.zeros(self._underlying.n2)

            # in the underlying white noise, just put a single spike at coordinate zero. We will then cycle the result
            test_field_lo[0] = np.sqrt(2.0)
            test_field_lo.in_fourier_space()
            test_field_lo[0]/=np.sqrt(2.0) # constant mode
            if len(test_field_lo)%2==0:
                test_field_lo[-1]/=np.sqrt(2.0) # nyquist mode

            underlying_out_0, _ = self._underlying.realization(white_noise_lo=test_field_lo, white_noise_hi = test_field_hi)

            cov = np.zeros((self.n1 + self.n2, self.n1 + self.n2))

            for i in range(self._underlying.n1):
                underlying_out_i = np.roll(underlying_out_0,i)

                out_lo, out_hi = self.realization(underlying_realization=underlying_out_i)

                cov[:self.n1, :self.n1] += np.outer(out_lo, out_lo)
                cov[self.n1:, self.n1:] += np.outer(out_hi, out_hi)
                cov[:self.n1, self.n1:] += np.outer(out_lo, out_hi)

            cov[self.n1:, :self.n1] = cov[:self.n1, self.n1:].T
            return cov
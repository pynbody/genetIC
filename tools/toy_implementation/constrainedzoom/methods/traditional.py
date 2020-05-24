import copy
import numpy as np
import numpy.testing
from ..fft_wrapper import in_real_space, FFTArray, complex_dot
from .geometric import ZoomConstrainedWithGeometricConstraints


class TraditionalZoomConstrained(ZoomConstrainedWithGeometricConstraints):
    description = "Traditional"
    realspace_convolution = True
    constrain_noise_directly = True # i.e. constraints are applied directly to independent noise vectors
    convolve_lores_at_hires = True # i.e. bring the lores noise in the centre of W into nW BEFORE convolving; otherwise interpolation errors will be greater

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.realspace_convolution:
            self.set_Chigh_realspace()

    def get_default_plot_padding(self):
        """The default plot padding (in coarse pixels) to hide from the high-res region"""
        return self.nW // (self.pixel_size_ratio) // 4




    @in_real_space
    def _modify_whitenoise(self, delta_low: FFTArray, delta_high: FFTArray):
        assert not delta_high.fourier
        assert not delta_low.fourier

        delta_high = self._remove_lr_pixel_info(delta_high)

        delta_low_zeroed = copy.copy(delta_low)

        if self.convolve_lores_at_hires:
            delta_low_zeroed[self._B_window_slice] = 0


        self._delta_low_residual = (delta_low-delta_low_zeroed).view(FFTArray)
        self._delta_low_residual.fourier = False

        lo_modes_for_hi_window = self.upsample_zeroorder(self._delta_low_residual) / np.sqrt(self.pixel_size_ratio)

        delta_high += lo_modes_for_hi_window

        return delta_low_zeroed, delta_high

    def _remove_lr_pixel_info(self, hr_vector):
        return hr_vector - self._only_lr_pixel_info(hr_vector)

    def _only_lr_pixel_info(self, hr_vector):
        return self.upsample_zeroorder(self.downsample(hr_vector))

    @in_real_space
    def _recombine_fields(self, deltaP, deltaW):
        deltaW += self.upsample_cubic(deltaP)
        deltaP += self._apply_transfer_function(self._delta_low_residual.in_fourier_space()).in_real_space()

        return deltaP, deltaW

    @in_real_space
    def covector_to_vector(self, lr_covec, hr_covec):
        hr_vec = hr_covec #self._remove_lr_pixel_info(hr_covec)
        # remove bits outside W but inside the 'pad' region
        hr_vec[:self.nW//4]=0
        hr_vec[3*self.nW//4]=0
        return lr_covec,hr_vec

    def zoom_covec_from_uniform_covec_in_window(self, hr_covec, potential):
        if potential:
            C_high = self.C_high_potential
            C_low = self.C_low_potential
        else:
            C_high = self.C_high
            C_low = self.C_low
        if hr_covec is None:
            hr_covec = self._default_constraint_hr_vec()
        # Unable to constrain in the 'pad' region
        if not (np.allclose(hr_covec[:self.nW // 4], 0, atol=hr_covec.max() * 1e-3)
                and np.allclose(hr_covec[(self.nW * 3) // 4:], 0, atol=hr_covec.max() * 1e-3)):
            raise ValueError(
                "The constraint extends outside the valid part of the high resolution window (which is only half of nW for traditional algorithms)")
        # filter is ( (I-W^+ W) C_P^(1/2) W  \hat{P} m   + P m W^+ XC^{1/2}X^{+}  )
        lr_covec = FFTArray(self.downsample_cubic(hr_covec, pad_around_window=True))
        lr_covec.in_fourier_space()
        lr_covec *= C_low ** 0.5 * self.pixel_size_ratio
        lr_covec.in_real_space()
        lr_covec[self._B_window_slice] = 0
        lr_covec2 = FFTArray(copy.copy(hr_covec))
        lr_covec2.in_fourier_space()
        lr_covec2 *= C_high ** 0.5 * self.pixel_size_ratio ** 0.5
        lr_covec2.in_real_space()
        lr_covec2[:self.nW // 4] = 0
        lr_covec2[3 * self.nW // 4:] = 0
        lr_covec += self.downsample(lr_covec2, input_unpadded=True)
        # filter is (I-WP^+PW^+) XC^{1/2}X+
        hr_covec = FFTArray(copy.copy(hr_covec))
        hr_covec.in_fourier_space()
        hr_covec *= C_high ** 0.5
        hr_covec.in_real_space()
        hr_covec = self._remove_lr_pixel_info(hr_covec)
        hr_covec[:self.nW // 4] = 0
        hr_covec[3 * self.nW // 4] = 0

        lr_covec.in_fourier_space()
        hr_covec.in_fourier_space()
        return lr_covec, hr_covec


class TraditionalZoomSingleLevelConstraint(TraditionalZoomConstrained):
    def zoom_covec_from_uniform_covec_in_window(self, hr_covec, potential):
        if potential:
            C_high = self.C_high_potential
            C_low = self.C_low_potential
        else:
            C_high = self.C_high
            C_low = self.C_low
        if hr_covec is None:
            hr_covec = self._default_constraint_hr_vec()

        # C_P^(1/2) W  \hat{P} m
        lr_covec = FFTArray(self.downsample_cubic(hr_covec))
        lr_covec.in_fourier_space()
        lr_covec *= C_low ** 0.5 * self.pixel_size_ratio

        hr_covec = FFTArray(np.zeros_like(hr_covec))
        hr_covec.in_fourier_space()

        return lr_covec, hr_covec

class BertschingerZoomConstrained(TraditionalZoomConstrained):
    description = "Naive Traditional"
    realspace_convolution = False
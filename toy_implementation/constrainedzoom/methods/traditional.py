import copy
import numpy as np
from ..fft_wrapper import in_real_space, FFTArray
from . import UnfilteredZoomConstrained


class TraditionalZoomConstrained(UnfilteredZoomConstrained):
    description = "Hahn/Abel"
    realspace_convolution = True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.realspace_convolution:
            self.set_Chigh_realspace()

    def get_default_plot_padding(self):
        """The default plot padding (in coarse pixels) to hide from the high-res region"""
        return self.n2 // (self.pixel_size_ratio) // 4

    @property
    def _B_window_slice(self):
        B_window_size = self.n2 // (self.pixel_size_ratio) // 2
        offset_B = self.offset + self.n2 // self.pixel_size_ratio // 2 - B_window_size // 2
        return slice(offset_B,offset_B+B_window_size)


    @in_real_space
    def _modify_whitenoise(self, delta_low: FFTArray, delta_high: FFTArray):
        assert not delta_high.fourier
        assert not delta_low.fourier

        delta_high -= self.upsample_zeroorder(self.downsample(delta_high))

        delta_low_zeroed = copy.copy(delta_low)

        delta_low_zeroed[self._B_window_slice] = 0


        self._delta_low_residual = (delta_low-delta_low_zeroed).view(FFTArray)
        self._delta_low_residual.fourier = False

        lo_modes_for_hi_window = self.upsample_zeroorder(self._delta_low_residual)

        delta_high += lo_modes_for_hi_window

        return delta_low_zeroed, delta_high

    def _separate_fields(self, delta_low: FFTArray, delta_high: FFTArray):
        return delta_low, delta_high, None

    @in_real_space
    def _recombine_fields(self, delta_low, delta_high, _):
        delta_high += self.upsample_cubic(delta_low)
        delta_low += self._apply_transfer_function(self._delta_low_residual.in_fourier_space()).in_real_space()

        return delta_low, delta_high


class BertschingerZoomConstrained(TraditionalZoomConstrained):
    description = "Bertschinger"
    realspace_convolution = False
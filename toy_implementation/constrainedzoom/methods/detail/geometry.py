from ...fft_wrapper import FFTArray, unitary_fft, unitary_inverse_fft, in_fourier_space, in_real_space, complex_dot
from functools import lru_cache, partial
import numpy as np
import scipy.fftpack

class GeometryAndPixelization:
    def __init__(self, cov_fn, n1, n2, hires_window_scale, offset):
        super().__init__(cov_fn, n1, n2, hires_window_scale, offset)
        assert n1 % hires_window_scale == 0, "Scale must divide n1 to fit pixels exactly"
        assert n2 % hires_window_scale == 0, "Scale must divide n2 to fit pixels exactly"
        self.n1 = n1
        self.n2 = n2
        self.window_size_ratio = hires_window_scale
        self.pixel_size_ratio = (self.window_size_ratio * self.n2) // self.n1
        self.offset = offset
        self.delta_low = 1. / self.n1
        self.delta_high = 1. / (self.n2 * self.window_size_ratio)

    def xs(self):
        """Return the real-space coordinates of the two outputs"""
        return (np.arange(self.n1) + 0.5) / self.n1, \
               (self.offset + (np.arange(self.n2) + 0.5) / self.pixel_size_ratio) / self.n1

    def boundary_xs(self):
        """Return the real-space coordinates of the pixel boundaries for the outputs"""
        centre1, centre2 = self.xs()
        return centre1 - 0.5 / self.n1, centre2 - 0.5 / (self.n1 * self.pixel_size_ratio)

    @in_real_space
    def extract_window(self, delta_highres_unwindowed):
        return delta_highres_unwindowed[
               self.offset * self.pixel_size_ratio:self.offset * self.pixel_size_ratio + self.n2]

    @in_real_space
    def place_window(self, delta_highres_windowed: FFTArray) -> FFTArray:
        delta_highres_unwindowed = np.zeros(self.n1 * self.pixel_size_ratio).view(FFTArray)
        assert delta_highres_unwindowed.fourier is False
        delta_highres_unwindowed[
        self.offset * self.pixel_size_ratio:self.offset * self.pixel_size_ratio + self.n2] = delta_highres_windowed
        return delta_highres_unwindowed

    @in_real_space
    def upsample_zeroorder(self, delta_low: FFTArray, in_window=True) -> FFTArray:
        """Take a low-res vector and put it in the high-res region without interpolating"""

        delta_highres = np.zeros(self.n1 * self.pixel_size_ratio)
        delta_highres = delta_highres.view(type=FFTArray)
        delta_highres.fourier = False

        for i in range(self.pixel_size_ratio):
            delta_highres[i::self.pixel_size_ratio] = delta_low

        if in_window:
            return self.extract_window(delta_highres)
        else:
            return delta_highres

    @in_real_space
    def upsample_linear(self, delta_low) -> FFTArray:
        """Take a low-res vector and interpolate it into the high-res region"""

        delta_highres = np.zeros(self.n1 * self.pixel_size_ratio)
        delta_low_left = np.roll(delta_low, 1)
        delta_low_right = np.roll(delta_low, -1)

        for i in range(self.pixel_size_ratio):
            sub_offset = (float(i) + 0.5) / self.pixel_size_ratio
            weight_cen = 1 - abs(sub_offset - 0.5)
            weight_left = 0.5 - sub_offset
            if weight_left < 0: weight_left = 0
            weight_right = sub_offset - 0.5
            if weight_right < 0: weight_right = 0

            delta_highres[i::self.pixel_size_ratio] = delta_low * weight_cen
            delta_highres[i::self.pixel_size_ratio] += delta_low_left * weight_left
            delta_highres[i::self.pixel_size_ratio] += delta_low_right * weight_right

        result = delta_highres[self.offset * self.pixel_size_ratio:self.offset * self.pixel_size_ratio + self.n2]
        return result.view(type=FFTArray)

    @in_real_space
    def upsample_cubic(self, delta_low) -> FFTArray:
        "Take a low-res vector and interpolate it into the high-res region - cubic interpolation"

        x_vals_low, x_vals_high = self.xs()
        delta_highres = scipy.interpolate.interp1d(x_vals_low, delta_low, kind='cubic')(x_vals_high)

        return delta_highres.view(type=FFTArray)

    @in_real_space
    def downsample(self, hires_vector: FFTArray, input_unpadded=True, output_padded=True) -> FFTArray:
        """Take a high-res region vector and downsample it onto the low-res grid"""
        vec_lr = np.zeros(self.n1).view(type=FFTArray)
        if input_unpadded:
            vec_lr[self.offset:self.offset + self.n2 // self.pixel_size_ratio] = \
                hires_vector.reshape((self.n2 // self.pixel_size_ratio, self.pixel_size_ratio)).mean(axis=1)
        else:
            vec_lr = hires_vector.reshape((self.n1, self.pixel_size_ratio)).mean(axis=1)

        if output_padded:
            return vec_lr
        else:
            return vec_lr[self.offset:self.offset + self.n2 // self.pixel_size_ratio]

    @property
    @lru_cache()
    def _upsample_cubic_matrix(self):
        matr = np.zeros((self.n2, self.n1))
        for i in range(self.n1):
            test = np.zeros(self.n1)
            test[i] = 1.0
            matr[:, i] = self.upsample_cubic(test)

        test = np.random.uniform(0, 1, self.n1)
        result = self.upsample_cubic(test)
        np.testing.assert_allclose(np.dot(matr, test), result, atol=1e-5)

        return matr

    @in_real_space
    def downsample_cubic(self, hr_vec, pad_around_window=True):
        """Take a high-res region vector and cubic downsample onto the low-res grid as an inverse to cubic upsampling"""
        matr = self._upsample_cubic_matrix.T / self.pixel_size_ratio
        lr_vec = np.dot(matr, hr_vec)
        if pad_around_window:
            return lr_vec.view(FFTArray)
        else:
            return lr_vec[self._B_window_slice].view(FFTArray)

    @property
    def _B_window_slice(self):
        B_window_size = self.n2 // (self.pixel_size_ratio) // 2
        offset_B = self.offset + self.n2 // self.pixel_size_ratio // 2 - B_window_size // 2
        return slice(offset_B, offset_B + B_window_size)

    def _get_ks(self):
        pixel_dx_low = 1. / self.n1
        pixel_dx_high = 1. / (self.n2 * self.window_size_ratio)
        k_low = scipy.fftpack.rfftfreq(self.n1, d=pixel_dx_low)
        k_high = scipy.fftpack.rfftfreq(self.n2, d=pixel_dx_high)
        return k_high, k_low

from .geometric import ZoomConstrainedWithGeometricConstraints
from ..fft_wrapper import in_fourier_space, in_real_space, FFTArray
import numpy as np
import copy


class FilteredZoomConstrained(ZoomConstrainedWithGeometricConstraints):
    """The reference implementation for fast-filter zooms, which should match exactly the description in the notes"""
    description = "Fast Filter"
    constrain_noise_directly = True

    def __init__(self, *args, **kwargs):
        """In addition to the standard arguments presented by ZoomConstrained, this accepts
        :param k_cut (0.5): the fraction of the low-res grid nyquist frequency used as the filter frequency
        :param T (0.1): the fractional temperature of the Fermi filter function, relative to the cut-off frequency
        """
        k_cut_fractional = kwargs.pop('k_cut', 0.5)
        T = kwargs.pop('T', 0.1)
        super().__init__(*args, **kwargs)
        self.k_cut = k_cut_fractional * self.k_low[-1]
        self.T = self.k_cut * T

    def filter_low(self, k):
        return 1./(1.+np.exp((k-self.k_cut)/self.T))

    def filter_high(self, k):
        return (1.-self.filter_low(k)**2)**0.5

    def get_default_plot_padding(self):
        return 4

    @in_real_space
    def zoom_covec_from_uniform_covec_in_window(self, hr_covec, potential):
        # C_high = self.C_high_potential if potential else self.C_high
        # C_low = self.C_low_potential if potential else self.C_low

        hr_covec.in_real_space()

        high = copy.copy(hr_covec)
        low = self.downsample_cubic(hr_covec) * self.pixel_size_ratio
        low.in_fourier_space()
        high.in_fourier_space()

        low *= self.filter_low(self.k_low)
        high *= self.filter_high(self.k_high)

        if potential:
            low *= self.C_low_potential ** 0.5
            high *= self.C_high_potential ** 0.5
        else:
            low *= self.C_low ** 0.5
            high *= self.C_high ** 0.5

        return low, high

    def _filter_in_window(self, lr: FFTArray, window_first) -> FFTArray:
        """This implements low-pass filtering only in the windowed region.

        This becomes important when converting vectors back to covectors (basically for getting chi^2). See notes.

        Depending on whether the filter appears as a Hermitian conjugate (or not), the windowing either occurs
        after the filtering (or before)."""
        lr.in_fourier_space()
        lr_filtered = lr.copy()

        if window_first:
            lr_filtered*=self.filter_low(self.k_low)
        lr_filtered.in_real_space()

        lr_no_window = self.zero_window(lr)
        lr.in_fourier_space()
        lr_filtered_in_window = lr_filtered - self.zero_window(lr_filtered)

        lr_filtered_in_window.in_fourier_space()
        if not window_first:
            lr_filtered_in_window*=self.filter_low(self.k_low)

        lr_no_window.in_fourier_space()

        return lr_filtered_in_window + lr_no_window

    @in_fourier_space
    def _apply_metric(self, low, high, use_windowing = False):
        """Applies the metric. In its simplified form the inverse metric is identical to the metric itself.

        Strictly, when applying the low-pass filter, it should only apply inside the high resolution window, since the
        low-res field does have some high-frequency information which is used outside that window.

        However when converting covectors to vectors one gets away with this because they are forced to be localised
        inside the high-res window and therefore there will not be any high frequency information outside that region.

        The approximation becomes very poor if one tries to convert the overdensity field to a covector (in order to
        calculate a chi^2), and so there is an option to include the correct windowing (set use_windowing=True) for
        this case. Note that the windowing order would be different if one wanted to include it for the
        covector -> vector case.
        """
        if use_windowing:
            # The following is only the correct form for vector->covector transforms
            # Otherwise the window ordering needs to be reversed
            low_from_low = self._filter_in_window(self._filter_in_window(low, window_first=True), window_first=False)
        else:
            low_from_low = low * self.filter_low(self.k_low) ** 2

        high_from_high = high * self.filter_high(self.k_high) ** 2

        low_from_high = high * self.filter_low(self.k_high) * self.filter_high(self.k_high)
        assert low_from_high.fourier
        low_from_high = self.downsample_cubic(
            low_from_high.in_real_space()).in_fourier_space() * self.pixel_size_ratio ** 0.5

        high_from_low = low * self.filter_low(self.k_low) * self.filter_high(
            self.k_low) / self.pixel_size_ratio ** 0.5
        assert high_from_low.fourier
        high_from_low = self.upsample_cubic(high_from_low.in_real_space()).in_fourier_space()
        assert high_from_low.fourier

        return low_from_low + low_from_high, \
               high_from_low + high_from_high

    def covector_to_vector(self, low, high):
        return self._apply_metric(low, high)

    def vector_to_covector(self, low, high):
        return self._apply_metric(low, high, use_windowing=True)

    @in_fourier_space
    def estimate_chi2_from_whitenoise(self, low, high):
        """A test to check we can get a sensible chi^2 for the fields we generate"""
        lowv, highv = self.vector_to_covector(low, high)
        return self.covector_vector_inner_product(low, high, lowv, highv)/2 # white noise norm: sigma^2 = 2

    @in_fourier_space
    def _recombine_fields(self, delta_low_k, delta_high_k):
        low_k_for_window = delta_low_k*self.filter_low(self.k_low)
        low_k_in_window = self.upsample_cubic(low_k_for_window.in_real_space())
        low_k_in_window.in_fourier_space()

        delta_high_k = delta_high_k*self.filter_high(self.k_high) + low_k_in_window

        return delta_low_k, delta_high_k


    def _modify_whitenoise(self, noiseP, noiseW):
        return noiseP, noiseW # filtering method regards initial white-noise fields as independent; no modification required


class FilteredZoomConstrainedInDeltaBasis(FilteredZoomConstrained):
    """A version of FilteredZoomConstrained that works in the delta_Z basis instead of the noise_Z basis"""
    description = "Fast Filter (delta basis)"
    constrain_noise_directly = False

    @in_real_space
    def zoom_covec_from_uniform_covec_in_window(self, hr_covec, potential):
        # C_high = self.C_high_potential if potential else self.C_high
        # C_low = self.C_low_potential if potential else self.C_low

        hr_covec.in_real_space()

        high = copy.copy(hr_covec)
        low = self.downsample_cubic(hr_covec) * self.pixel_size_ratio
        low.in_fourier_space()
        high.in_fourier_space()

        low*=self.filter_low(self.k_low)
        high*=self.filter_high(self.k_high)

        if potential:
            low *= (self.C_low_potential / self.C_low) ** 0.5
            high *= (self.C_high_potential / self.C_high) ** 0.5

        return low, high

    @in_fourier_space
    def covector_to_vector(self, low, high):

        low_from_low = low * self.C_low * self.filter_low(self.k_low) ** 2 / self.pixel_size_ratio
        high_from_high = high * self.C_high * self.filter_high(self.k_high) ** 2

        low_from_high = high * self.filter_low(self.k_high) * self.filter_high(self.k_high) * self.C_high
        assert low_from_high.fourier
        low_from_high = self.downsample_cubic(low_from_high.in_real_space()).in_fourier_space()

        high_from_low = low * self.C_low  * self.filter_low(self.k_low) * self.filter_high(self.k_low) / self.pixel_size_ratio
        assert high_from_low.fourier
        high_from_low = self.upsample_cubic(high_from_low.in_real_space()).in_fourier_space()
        assert high_from_low.fourier

        return low_from_low + low_from_high, \
               high_from_low + high_from_high


class HybridZoomConstrained(FilteredZoomConstrained):
    """This is a trial implementation to see if a hybrid approach can improve accuracy further (it doesn't seem to help)"""
    description = "Hybrid trad+filter"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_Chigh_realspace()


    @in_real_space
    def _modify_whitenoise(self, delta_low: FFTArray, delta_high: FFTArray):
        assert not delta_high.fourier
        assert not delta_low.fourier

        delta_high -= self.upsample_zeroorder(self.downsample(delta_high))
        delta_high += self.upsample_zeroorder(delta_low)

        delta_high.in_fourier_space()
        delta_high*=self.filter_high(self.k_high)

        delta_low.in_fourier_space()
        delta_low*=self.filter_low(self.k_low)

        return delta_low, delta_high


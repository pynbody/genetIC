from .geometric import ZoomConstrainedWithGeometricConstraints
from ..fft_wrapper import in_fourier_space, in_real_space, unitary_inverse_fft, unitary_fft, complex_dot, FFTArray
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

    @in_fourier_space
    def covector_to_vector(self, low, high):
        low_from_low = low * self.filter_low(self.k_low) ** 2 / self.pixel_size_ratio
        high_from_high = high * self.filter_high(self.k_high) ** 2

        low_from_high = high * self.filter_low(self.k_high) * self.filter_high(self.k_high)
        assert low_from_high.fourier
        low_from_high = self.downsample_cubic(low_from_high.in_real_space()).in_fourier_space()

        high_from_low = low * self.filter_low(self.k_low) * self.filter_high(
            self.k_low) / self.pixel_size_ratio
        assert high_from_low.fourier
        high_from_low = self.upsample_cubic(high_from_low.in_real_space()).in_fourier_space()
        assert high_from_low.fourier

        return low_from_low + low_from_high, \
               high_from_low + high_from_high

    @in_fourier_space
    def _recombine_fields(self, delta_low_k, delta_high_k):
        low_k_for_window = delta_low_k*self.filter_low(self.k_low)
        low_k_in_window = self.upsample_cubic(low_k_for_window.in_real_space())
        low_k_in_window.in_fourier_space()

        delta_high_k = delta_high_k*self.filter_high(self.k_high) + low_k_in_window

        return delta_low_k, delta_high_k


    def _modify_whitenoise(self, wn_lo, wn_hi):
        return wn_lo, wn_hi # filtering method regards initial white-noise fields as independent; no modification required


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






class FilteredZoomConstrainedOriginal(FilteredZoomConstrained):
    """The original implementation of fast-filter zooms, which corresponds more closely to the genetIC code but
    shifts various factors around and therefore looks quite different from the FilteredZoomConstrained implementation
    above."""
    description = "Fast Filter (original implementation)"
    constrain_noise_directly = False

    @in_fourier_space
    def _separate_fields(self, delta_low_k, delta_high_k):
        k_high, k_low = self._get_ks()
        delta_high_k *= np.sqrt(1. - self.filter_low(k_high) ** 2)  # keep original power spectrum
        self.delta_low_supplement = delta_low_k * (1.-self.filter_low(
            k_low))  # store the filtered-out components of the low-res field
        delta_low_k *= self.filter_low(k_low)
        return delta_low_k, delta_high_k

    @in_real_space
    def _recombine_fields(self, delta_low, delta_high):
        delta_high += self.upsample_cubic(delta_low)
        delta_low += self.delta_low_supplement.in_real_space()
        return delta_low, delta_high

    @in_fourier_space
    def high_k_vector_from_low_k_vector(self, low_harmonics):
        pixelized_highres = self.upsample_cubic(unitary_inverse_fft(low_harmonics))
        return unitary_fft(pixelized_highres)

    @in_fourier_space
    def low_k_vector_from_high_k_vector(self, high_harmonics):
        pixelized_lowres = self.downsample(unitary_inverse_fft(high_harmonics))
        return unitary_fft(pixelized_lowres)

    @in_fourier_space
    def covector_vector_inner_product(self, low1, high1, low2, high2, more_accurate=True):

        # The low.low part picks up a pixel_size_ratio factor. We can see this as follows.
        # Take the ideal case where low1 = f_low v1, high1 = f_high v1, f_low+f_high=1 and f_low f_high = 0.
        # Without varying pixel sizes, v1 = (low1+high1), and v1.v2 = (low1.low1)+(high1.high1), exactly.
        # Now, let's downsample the low1 pixel scale. The convention we've adopted is that the big
        # low1 pixels take on the mean value of the original, finer low1 pixels. So, the new dot product
        # in real space has been multiplied by 1/pixel_size_ratio. (Because FFTs are unitary, this
        # applies also to the harmonic space dot product). We need to multiply by pixel_size_ratio to
        # cancel out this change.
        product = complex_dot(low1,low2)* self.pixel_size_ratio+complex_dot(high1,high2)
        if more_accurate:
            # add in the low1^dagger C high2 + high1 C low2^dagger terms
            low1_as_high = self.high_k_vector_from_low_k_vector(low1)
            low2_as_high = self.high_k_vector_from_low_k_vector(low2)
            product+=(complex_dot(low1_as_high,high2) + complex_dot(high1, low2_as_high))
        return product

    @in_fourier_space
    def covector_to_vector(self, low, high):
        return self.C_low * low, self.C_high * high

    @in_fourier_space
    def zoom_covec_from_uniform_covec_in_window(self, hr_covec, potential):

        high = hr_covec * (1. - self.filter_low(self.k_high) )
        low = self.downsample(hr_covec - high)

        high.in_fourier_space()
        low.in_fourier_space()

        if potential:
            high *= (self.C_high_potential / self.C_high) ** 0.5
            low *= (self.C_low_potential / self.C_low) ** 0.5

        return low, high

    @in_fourier_space
    def _apply_constraints(self, noise_or_delta_low_k, noise_or_delta_high_k):
        """This is a slightly ugly re-intepretation of _apply_constraints from the base class, but is required
        to handle the delta_low_supplement part"""
        for (al_low_k, al_high_k), d in zip(self.constraints, self.constraints_val):
            scale = d - self.covector_vector_inner_product(al_low_k, al_high_k, noise_or_delta_low_k,
                                                           noise_or_delta_high_k)
            vec_low_k, vec_high_k = self.covector_to_vector(al_low_k, al_high_k)
            noise_or_delta_low_k += vec_low_k * scale
            noise_or_delta_high_k += vec_high_k * scale

            self.delta_low_supplement += self.low_k_vector_from_high_k_vector(vec_high_k * scale)

        return noise_or_delta_low_k, noise_or_delta_high_k


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


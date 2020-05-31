import numpy as np

from constrainedzoom.fft_wrapper import in_fourier_space, in_real_space, unitary_inverse_fft, unitary_fft, complex_dot
from constrainedzoom.methods import FilteredZoomConstrained


class FilteredZoomConstrainedOriginal(FilteredZoomConstrained):
    """The original implementation of fast-filter zooms, which corresponds more closely to the genetIC code but
    shifts various factors around and therefore looks quite different from the FilteredZoomConstrained implementation
    above."""
    description = "Fast Filter (original implementation)"
    constrain_noise_directly = False

    @in_fourier_space
    def _separate_fields(self, deltaP, deltaW):
        k_high, k_low = self._get_ks()
        deltaW *= np.sqrt(1. - self.filter_low(k_high) ** 2)  # keep original power spectrum
        self.delta_low_supplement = deltaP * (1. - self.filter_low(
            k_low))  # store the filtered-out components of the low-res field
        deltaP *= self.filter_low(k_low)
        return deltaP, deltaW

    @in_real_space
    def _recombine_fields(self, deltaP, deltaW):
        deltaW += self.upsample_cubic(deltaP)
        deltaP += self.delta_low_supplement.in_real_space()
        return deltaP, deltaW

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
from . import ZoomConstrained
from ..fft_wrapper import in_fourier_space, in_real_space, unitary_inverse_fft, unitary_fft, complex_dot, FFTArray
import numpy as np
import copy
import math

class FilteredZoomConstrained(ZoomConstrained):
    description = "Fast Filter"

    def get_default_plot_padding(self):
        return 4

    def __init__(self, *args, **kwargs):
        """Initialise a ZoomConstrained instance that uses filtering to combine fields on different grids.
        
        In addition to the standard arguments presented by ZoomConstrained, this accepts
        :param k_cut (0.5): the fraction of the low-res grid nyquist frequency used as the filter frequency"""

        k_cut_fractional = kwargs.pop('k_cut', 0.3)
        super().__init__(*args, **kwargs)
        self.k_cut = k_cut_fractional * self.k_low[-1]

    def _modify_whitenoise(self, wn_lo, wn_hi):
        return wn_lo, wn_hi # filtering method regards initial white-noise fields as independent

    def filter_low(self, k):
        T = self.k_cut/10
        return 1./(1.+np.exp((k-self.k_cut)/T))

    def filter_high(self, k):
        return 1.-self.filter_low(k)

    @in_fourier_space
    def _separate_fields(self, delta_low_k, delta_high_k):
        k_high, k_low = self._get_ks()
        delta_high_k *= np.sqrt(1. - self.filter_low(k_high) ** 2)  # keep original power spectrum
        delta_low_k_plus = delta_low_k * self.filter_high(
            k_low)  # store the filtered-out components of the low-res field
        delta_low_k *= self.filter_low(k_low)
        return delta_low_k, delta_high_k, delta_low_k_plus

    @in_real_space
    def _recombine_fields(self, delta_low, delta_high, delta_low_k_plus):
        delta_high += self.upsample_cubic(delta_low)
        delta_low += delta_low_k_plus.in_real_space()
        return delta_low, delta_high

    @in_fourier_space
    def _apply_constraints(self, delta_low_k, delta_high_k, verbose):
        for (al_low_k, al_high_k), d in zip(self.constraints, self.constraints_val):
            scale = d - self.covector_vector_inner_product(al_low_k, al_high_k, delta_low_k, delta_high_k)
            vec_low_k, vec_high_k = self.covector_to_vector(al_low_k, al_high_k)
            delta_low_k += vec_low_k * scale
            delta_high_k += vec_high_k * scale
            self.delta_low_supplement+=self.low_k_vector_from_high_k_vector(vec_high_k*scale)

        return delta_low_k, delta_high_k


    @in_fourier_space
    def high_k_vector_from_low_k_vector(self, low_harmonics):
        pixelized_highres = self.upsample_cubic(unitary_inverse_fft(low_harmonics))
        return unitary_fft(pixelized_highres)

    @in_fourier_space
    def low_k_vector_from_high_k_vector(self, high_harmonics):
        pixelized_lowres = self.downsample(unitary_inverse_fft(high_harmonics))
        return unitary_fft(pixelized_lowres)

    @in_fourier_space
    def covector_norm(self, low, high):
        return self.covector_covector_inner_product(low, high, low, high)

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
    def covector_covector_inner_product(self, low1, high1, low2, high2):
        return self.covector_vector_inner_product(low1, high1, *self.covector_to_vector(low2, high2))

    @in_fourier_space
    def covector_to_vector(self, low, high):
        return self.C_low*low, self.C_high*high

    def add_constraint(self, val=0.0, hr_covec=None, potential=False):
        """Add a constraint, specifying the constraint covector in the high-res region and the value it should take
        """

        if hr_covec is None:
            hr_covec = self._default_constraint_hr_vec()

        self.constraints_real.append(hr_covec) # stored only for information - not part of the algorithm


        high = unitary_fft(hr_covec) * (1. - self.filter_low(self.k_high)) ** 0.5
        low = unitary_fft(self.downsample(hr_covec - unitary_inverse_fft(high)))

        if potential:
            high*=(self.C_high_potential/self.C_high)**0.5
            low*=(self.C_low_potential/self.C_low)**0.5

        # perform Gram-Schmidt orthogonalization
        for (la, ha),va in zip(self.constraints,self.constraints_val):
            dotprod = self.covector_covector_inner_product(la,ha,low,high)
            low-=dotprod*la
            high-=dotprod*ha
            val-=dotprod*va

        norm = self.covector_norm(low, high)
        low/=math.sqrt(norm)
        high/=math.sqrt(norm)
        val/=math.sqrt(norm)

        self.constraints.append((low,high))
        self.constraints_val.append(val)



class HybridZoomConstrained(FilteredZoomConstrained):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.set_Chigh_realspace()
        #self._apodize_Chigh()

    @property
    def _B_window_slice(self):
        B_window_size = self.n2 // (self.pixel_size_ratio) // 2
        offset_B = self.offset + self.n2 // self.pixel_size_ratio // 2 - B_window_size // 2
        return slice(offset_B, offset_B + B_window_size)

    @in_real_space
    def _modify_whitenoise(self, delta_low: FFTArray, delta_high: FFTArray):
        assert not delta_high.fourier
        assert not delta_low.fourier

        delta_high -= self.upsample_zeroorder(self.downsample(delta_high))
        delta_high += self.upsample_zeroorder(delta_low)

        return delta_low, delta_high

    @in_fourier_space
    def _separate_fields(self, delta_low_k, delta_high_k):
        return delta_low_k, delta_high_k, None

    def _recombine_fields(self, delta_low: FFTArray, delta_high: FFTArray, _):
        delta_high.in_fourier_space()
        delta_high *= self.filter_high(self.k_high)

        delta_low_upsampled = copy.copy(delta_low)
        delta_low_upsampled.in_fourier_space()

        delta_low_upsampled = self.upsample_cubic(delta_low_upsampled)
        delta_low_upsampled.in_fourier_space()
        delta_low_upsampled*=self.filter_low(self.k_high)
        
        delta_high += delta_low_upsampled
        return delta_low, delta_high


class StepFilterZoomConstrained(FilteredZoomConstrained):
    @in_fourier_space
    def _separate_fields(self, delta_low_k, delta_high_k):
        k_high, k_low = self._get_ks()
        delta_high_k[k_high<=k_low.max()] = 0
        return delta_low_k, delta_high_k, None

    def _recombine_fields(self, delta_low: FFTArray, delta_high: FFTArray, _):
        delta_low_window = copy.copy(delta_low.in_real_space()[self.offset:self.offset+self.n1//self.window_size_ratio])
        delta_low_window.in_fourier_space()

        delta_high.in_fourier_space()
        delta_high[:len(delta_low_window)] = delta_low_window*np.sqrt(len(delta_high)/len(delta_low_window))

        return delta_low, delta_high
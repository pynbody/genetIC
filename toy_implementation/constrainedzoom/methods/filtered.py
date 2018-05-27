from . import ZoomConstrained
from fft_wrapper import in_fourier_space, in_real_space, unitary_inverse_fft, complex_dot, FFTArray
import numpy as np
import copy

class FilteredZoomConstrained(ZoomConstrained):
    description = "Fast Filter"

    def get_default_plot_padding(self):
        return 4

    def __init__(self, *args, **kwargs):
        """Initialise a ZoomConstrained instance that uses filtering to combine fields on different grids.
        
        In addition to the standard arguments presented by ZoomConstrained, this accepts
        :param k_cut (0.5): the fraction of the low-res grid nyquist frequency used as the filter frequency"""

        k_cut_fractional = kwargs.pop('k_cut', 0.5)
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

            if verbose:
                print("lowdot=", complex_dot(al_low_k, delta_low_k) * self.pixel_size_ratio)
                print("highdot=", complex_dot(al_high_k, delta_high_k))
                print("sum=", complex_dot(al_low_k, delta_low_k) * self.pixel_size_ratio + complex_dot(al_high_k,
                                                                                                       delta_high_k))
                al_low, al_high = self._harmonic_to_combined_pixel(al_low_k, al_high_k)
                print("RS simple dot=", np.dot(unitary_inverse_fft(al_high_k), unitary_inverse_fft(delta_high_k)))
                print("RS dot=", np.dot(al_high, delta_high_k))

            scale = d - complex_dot(al_low_k, delta_low_k) * self.pixel_size_ratio - complex_dot(al_high_k,
                                                                                                 delta_high_k)

            if verbose:
                print("scale=", scale)

            delta_low_k += self.C_low * al_low_k * scale
            delta_high_k += self.C_high * al_high_k * scale
        return delta_low_k, delta_high_k


    def _harmonic_to_combined_pixel(self, f_low_k, f_high_k):
        delta_low, delta_high = self.harmonic_to_pixel(f_low_k, f_high_k)# delta_high does not contain low-frequency contributions.

        # interpolate the low frequency contribution into the high-res window
        # prepare the data range just around the window ...
        delta_lowhigh = self.upsample_linear(delta_low)

        delta_high+=delta_lowhigh

        return delta_low, delta_high

    def harmonic_to_pixel(self, f_low_k, f_high_k):
        delta_low = unitary_inverse_fft(f_low_k)  # *self.filter_low(self.k_low))
        if f_high_k is not None:
            delta_high = unitary_inverse_fft(f_high_k)
        else:
            delta_high = np.zeros(self.n2)

        return delta_low, delta_high



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
        #delta_low_upsampled*=self.filter_low(self.k_low)**0.1

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
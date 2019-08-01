import copy
import abc
import math
from ..fft_wrapper import FFTArray, complex_dot, in_fourier_space
from . import ZoomConstrained


class ZoomConstrainedWithGeometricConstraints(ZoomConstrained):
    """A base class for those zoom implementations which adopt the covector/vector geometry described in the notes."""

    @abc.abstractmethod
    def zoom_covec_from_uniform_covec_in_window(self, hr_covec, potential):
        """Create a zoom covector from the uniform-basis covector within the zoom window.

        Should apply T^\dagger W^+ in terms of the generalised form in the notes.

        Returns u_P, u_W: the P and W parts of the covector in Fourier space"""
        pass

    @abc.abstractmethod
    def covector_to_vector(self, lr_covec, hr_covec):
        """Turn a covector into a vector.

        Should multiply by the metric g in terms of the generalised form in the notes.

        Returns u^P, u^W: the P and W parts of the vector in Fourier space"""
        pass

    @in_fourier_space
    def covector_vector_inner_product(self, low1, high1, low2, high2):
        """Perform the covector vector inner product"""
        return complex_dot(low1, low2) + complex_dot(high1, high2)

    @in_fourier_space
    def covector_norm(self, low, high):
        """Return the covector norm"""
        return self.covector_covector_inner_product(low, high, low, high)

    @in_fourier_space
    def covector_covector_inner_product(self, low1, high1, low2, high2):
        """Perform a covector/covector inner product (by first converting one of the covectors to a vector)"""
        low2v, high2v = self.covector_to_vector(low2, high2)
        return self.covector_vector_inner_product(low1, high1, low2v.in_fourier_space(), high2v.in_fourier_space())

    @in_fourier_space
    def _apply_constraints(self, noise_or_delta_low_k, noise_or_delta_high_k):
        for (al_low_k, al_high_k), d in zip(self.constraints, self.constraints_val):
            scale = d - self.covector_vector_inner_product(al_low_k, al_high_k, noise_or_delta_low_k,
                                                           noise_or_delta_high_k)
            vec_low_k, vec_high_k = self.covector_to_vector(al_low_k, al_high_k)
            noise_or_delta_low_k += vec_low_k.in_fourier_space() * scale
            noise_or_delta_high_k += vec_high_k.in_fourier_space() * scale

        return noise_or_delta_low_k, noise_or_delta_high_k

    def add_constraint(self, val=0.0, hr_covec=None, potential=False):

        if hr_covec is None:
            hr_covec = self._default_constraint_hr_vec()

        hr_covec: FFTArray = copy.copy(hr_covec).view(FFTArray)
        self.constraints_real.append(hr_covec)  # stored only for information - not part of the algorithm

        low, high = self.zoom_covec_from_uniform_covec_in_window(hr_covec, potential)

        # perform Gram-Schmidt orthogonalization
        for (la, ha), va in zip(self.constraints, self.constraints_val):
            dotprod = self.covector_covector_inner_product(la, ha, low, high)
            low -= dotprod * la
            high -= dotprod * ha
            val -= dotprod * va

        norm = self.covector_norm(low, high)

        low /= math.sqrt(norm)
        high /= math.sqrt(norm)
        val /= math.sqrt(norm)

        self.constraints.append((low, high))
        self.constraints_val.append(val)
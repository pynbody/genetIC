import numpy as np
from ...fft_wrapper import FFTArray
import abc

class CovarianceCalculation(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def __init__(self):
        # indicate the fields expected to be present
        self.nP = 0
        self.nW = 0
        self.constraints_real = []
        self.pixel_size_ratio = 0

    @abc.abstractmethod
    def realization(self, no_random=False, noiseP: FFTArray = None, noiseW: FFTArray = None) -> [FFTArray, FFTArray]:
        return FFTArray(), FFTArray()

    def get_cov(self, one_element=None) -> np.ndarray:
        """Get the exact covariance matrix using a deterministic algorithm.

        Returns a (nP+nW) x (nP+nW) matrix, with the upper left bit being the deltaP covariance, bottom right being
        the deltaW covariance, and the upper-right/bottom-left being the deltaP deltaW cross-covariance.

        See estimate_cov for an estimation algorithm approach to the same calculation.

        :param one_element: if set, inject noise only into a single pixel (if None, get the entire covariance matrix).
        """
        cov = np.zeros((self.nP + self.nW, self.nP + self.nW))

        element_iterator = self._iter_one_cov_element(
            one_element) if one_element is not None else self._iter_cov_elements()

        for test_field_lo, test_field_hi in element_iterator:
            out_lo, out_hi = self.realization(noiseP=test_field_lo, noiseW=test_field_hi)
            cov[:self.nP, :self.nP] += np.outer(out_lo, out_lo)
            cov[self.nP:, self.nP:] += np.outer(out_hi, out_hi)
            cov[:self.nP, self.nP:] += np.outer(out_lo, out_hi)

        cov[self.nP:, :self.nP] = cov[:self.nP, self.nP:].T
        return cov

    def estimate_cov(self, Ntrials=2000, with_means=False):
        """Estimate the covariance matrix using the specified number of random trials

        :param with_means: if True, return also the mean field and the variance in the constrained values
        (which should be near-zero)"""

        cov = np.zeros((self.nP + self.nW, self.nP + self.nW))
        cstr_means = np.zeros(len(self.constraints_real))
        constraint_variance = np.zeros(len(self.constraints_real))
        means = np.zeros(self.nW)
        for i in range(Ntrials):
            r1, r2 = self.realization()
            cov[:self.nP, :self.nP] += np.outer(r1, r1)
            cov[:self.nP, self.nP:] += np.outer(r1, r2)
            cov[self.nP:, self.nP:] += np.outer(r2, r2)
            means += r2
            cstr_means += [np.dot(cstr, r2) for cstr in self.constraints_real]
            constraint_variance += [np.dot(cstr, r2) ** 2 for cstr in self.constraints_real]

        cstr_means /= Ntrials
        constraint_variance /= Ntrials
        cov /= Ntrials
        means /= Ntrials

        constraint_variance -= cstr_means ** 2

        cov[self.nP:, :self.nP] = cov[:self.nP, self.nP:].T

        if with_means:
            return cov, means, np.sqrt(constraint_variance)
        else:
            return cov

    def _iter_cov_elements(self):
        """Helper for get_cov: iterates over sample 'white noise' fields such that the covariance is exactly the sum
        of the resulting output fields"""
        test_field_lo = np.zeros(self.nP)
        test_field_hi = np.zeros(self.nW)

        for i in range(self.nP):
            test_field_lo[i] = 1.0
            yield test_field_lo, test_field_hi
            test_field_lo[i] = 0.0

        for i in range(self.nW):
            test_field_hi[i] = 1.0
            yield test_field_lo, test_field_hi
            test_field_hi[i] = 0.0

    def _iter_one_cov_element(self, offset):
        """Helper for get_cov testing: acts like _iter_cov_elements but only generates a single pixel. If offset is
        positive, this is a delta-function at the centre of the high-res box plus offset pixels. If offset is negative,
        it is a delta-function to the left of the high-res box by the specified amount."""
        test_field_lo = np.zeros(self.nP)
        test_field_hi = np.zeros(self.nW)

        if offset >= 0:
            # place test delta-function in high-res region
            test_field_hi[self.nW // 2 + self.pixel_size_ratio // 2 + offset] = self.pixel_size_ratio ** 0.5
            test_field_lo = self.downsample(test_field_hi)
        else:
            # place in low-res region just next to high-res region
            test_field_lo[self.offset + offset] = 1.0

        yield FFTArray(test_field_lo).in_fourier_space(), FFTArray(test_field_hi).in_fourier_space()
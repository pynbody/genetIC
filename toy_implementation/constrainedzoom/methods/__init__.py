import abc
import math
import numpy as np
import scipy.fftpack
import scipy.interpolate

from ..fft_wrapper import FFTArray, in_fourier_space, in_real_space, complex_dot
from .detail.powspec import Powspec
from .detail.geometry import GeometryAndPixelization
from .detail.covcalc import CovarianceCalculation


class ZoomConstrained(GeometryAndPixelization, Powspec, CovarianceCalculation,
                      metaclass=abc.ABCMeta):

    description = "Unknown Method"
    constrain_noise_directly = False

    def __init__(self, cov_fn = None, nP=256, nW=256, hires_window_scale=4, offset = 10):
        super().__init__(cov_fn, nP, nW, hires_window_scale, offset)

        self.constraints =[]
        self.constraints_val = []
        self.constraints_real = []

    def get_default_plot_padding(self):
        """The default plot padding (in coarse pixels) to hide from the high-res region"""
        return 0

    def realization(self, no_random=False, noiseP: FFTArray=None, noiseW: FFTArray=None) -> [FFTArray, FFTArray] :
        """Generate a realization of the Gaussian random field with constraints, if present.
        
        This implementation calls a generic sequence of operations which are overriden in base class implementations to
        achieve different zoom techniques.

        :param no_random: if True, create a zeroed random field to see effect of constraints in isolation
        :param noiseP: if present, the array of white noise for the low-frequency (nP) box
        :param noiseW:  if present, the array of white noise for the high-frequency (nW) sub-box
        :return: delta_low, delta_high: the fields for the low-frequency (nP) box and high-frequency (nW) sub-box respectively
        """

        k_high, k_low = self._get_ks()

        if no_random:
            noiseP = np.zeros_like(k_low).view(FFTArray)
            noiseW = np.zeros_like(k_high).view(FFTArray)
        else:
            # Generate white-noise, or correctly encapsulate the white-noise samples passed in
            noiseP, noiseW = self._get_whitenoise(noiseP, noiseW)

        # In trad approach, apply constraints directly to the noise; covariance lives in the constraint covectors
        if self.constrain_noise_directly and len(self.constraints)>0:
            self._apply_constraints(noiseP, noiseW)

        # Traditional approaches make the low-frequency and high-frequency white noise compatible
        noiseP, noiseW = self._modify_whitenoise(noiseP, noiseW)

        # Apply the transfer function as stored in Fourier space as self.C_high and self.C_low
        delta_low, delta_high = self._apply_transfer_function(noiseP, noiseW)

        # In the original implementation of the filtering approach, at this point the relevant filters were applied
        # to delta_high and delta_low. This original implementation is in the `old' module. In the new implementation,
        # as well as in the traditional approach, nothing actually happens here.
        delta_low, delta_high = self._separate_fields(delta_low, delta_high)

        # Now we can apply the constraints in the filter approach
        if (not self.constrain_noise_directly) and len(self.constraints)>0:
            delta_low, delta_high = self._apply_constraints(delta_low, delta_high)

        # As a final step, relevant information must be copied from the low-resolution to the high-resolution regions.
        delta_low, delta_high = self._recombine_fields(delta_low, delta_high)

        return delta_low.in_real_space(), delta_high.in_real_space()


    @in_fourier_space
    def _get_whitenoise(self, noiseP: FFTArray=None, noiseW: FFTArray=None):
        """Generate white noise fields nP and nW, with unit variance per mode.

        :param noiseP: set to provide a specified pixelized noise vector
        :param noiseW: set to provide a specified windowed noise vector
        :return: nP, nW: the generated noise fields
        """

        if noiseP is None:
            noiseP, noiseW = np.random.normal(0.0, 1.0, size=self.nP), np.random.normal(0.0, 1.0, size=self.nW)
            # only one nyquist mode is included in DFT but two are physically present
            if len(noiseP) % 2 == 0:
                noiseP[-1] *= np.sqrt(2)
            if len(noiseW) % 2 == 0:
                noiseW[-1] *= np.sqrt(2)

        noiseP = np.copy(noiseP).view(FFTArray)
        noiseW = np.copy(noiseW).view(FFTArray)
        noiseP.fourier=noiseW.fourier=True

        return noiseP, noiseW
    
    @abc.abstractmethod
    def _modify_whitenoise(self, noiseP: FFTArray, noiseW: FFTArray) -> [FFTArray, FFTArray]:
        """Give subclass an opportunity to modify the whitenoise fields before the transfer function is applied
 
        :param noiseP: the low-frequency pixelised whitenoise
        :param noiseW: the high-frequency windowed whitenoise
        :return: noiseP, noiseW: the same fields, modified
        """
        pass

    def _separate_fields(self, deltaP: FFTArray, deltaW: FFTArray):
        """Give subclass an opportunity to apply low and high-pass filters to the fields
        before constraints are applied. Default implementation does nothing.
        
        :param deltaP: the field in the pixelised basis, to be low-pass filtered
        :param deltaW: the field in the windowed region, to be high-pass filtered
        :return: deltaP, deltaW - the two fields post-separation
        """
        return deltaP, deltaW

    def _recombine_fields(self, deltaP: FFTArray, deltaW: FFTArray):
        """Give subclass an opportunity to combine information from the low-frequency into the high-frequency field.

        Default implementation does nothing.
        
        :param deltaP: the field in the pixelised basis, without any high-frequency information
        :param deltaW: the field in the window region, without any low-frequency information
        :return: deltaP, deltaW: the final version of the fields
        """
        return deltaP, deltaW

    def _default_constraint_hr_vec(self):
        """Generate a simple constraint covector (in the high resolution window) for quick tests"""
        hr_vec = np.zeros(self.nW)
        cval = self.nW // 2
        hr_vec[cval + 2] = 1.0
        return hr_vec

    @abc.abstractmethod
    def add_constraint(self, val=0.0, hr_covec=None, potential=False):
        """Add a constraint, specifying the constraint covector in the high-res region and the value it should take.

        If hr_covec is None, use a default constraint covector consisting of a delta function in the centre of the window
        """
        pass

    @abc.abstractmethod
    def _apply_constraints(self, noise_or_delta_low_k, noise_or_delta_high_k):
        """Apply the constraints either to the noise or delta zoom vectors.

        Returns the modified vectors.

        If constrain_noise_directly is True in the class, the noise vectors n_W and n_P are passed; otherwise
        delta_W and delta_P"""
        pass



class UnfilteredZoomConstrained(ZoomConstrained):
    @in_real_space
    def _separate_fields(self, delta_low, delta_high):
        delta_high -= self.upsample_zeroorder(self.downsample(delta_high))
        return delta_low, delta_high

    @in_real_space
    def _recombine_fields(self, deltaP, deltaW):
        deltaW += self.upsample_zeroorder(deltaP)
        return deltaP, deltaW

    def _modify_whitenoise(self, noiseP, noiseW):
        return noiseP, noiseW

    def add_constraint(self, val=0.0, lr_covec=None, potential=False):
        """Add a constraint to the LR field (note this differs from other classes which apply to the HR field)
        :param potential:
        """

        lr_covec = lr_covec.view(FFTArray).in_fourier_space()

        if potential:
            lr_covec*=(self.C_low_potential/self.C_low)**0.5

        # perform Gram-Schmidt orthogonalization
        if len(self.constraints)!=0:
            raise RuntimeError("Orthogonalization not implemented yet for UnfilteredZoomConstrained!")

        norm = complex_dot(lr_covec, self.C_low*lr_covec)
        lr_covec /= math.sqrt(norm)
        val /= math.sqrt(norm)

        self.constraints.append(lr_covec)
        self.constraints_val.append(val)

    @in_fourier_space
    def _apply_constraints(self, delta_low_k, delta_high_k):
        # implementation works purely on low-res part; ignores window
        for al_low_k, d in zip(self.constraints, self.constraints_val):
            scale = d - complex_dot(al_low_k, delta_low_k)
            delta_low_k += self.C_low * al_low_k * scale

        return delta_low_k, delta_high_k



from . import idealized, ml, traditional, filtered

# auto-reload submodules to make development easier
from importlib import reload
reload(idealized)
reload(ml)
reload(traditional)
reload(filtered)

from .idealized import FastIdealizedZoomConstrained, IdealizedZoomConstrained
from .ml import MLZoomConstrained
from .traditional import TraditionalZoomConstrained
from .filtered import FilteredZoomConstrained

__all__ = ['FastIdealizedZoomConstrained', 'IdealizedZoomConstrained',
           'MLZoomConstrained', 'TraditionalZoomConstrained', 'FilteredZoomConstrained']

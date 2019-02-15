import importlib
import math

import abc

import scipy.fftpack
import scipy.integrate
import scipy.interpolate
import numpy as np

from ..fft_wrapper import FFTArray, unitary_fft, unitary_inverse_fft, in_fourier_space, in_real_space, complex_dot


class ZoomConstrained(metaclass=abc.ABCMeta):
    description = "Unknown Method"
    constrain_noise_directly = False

    def __init__(self, cov_fn = None, n1=256, n2=256, hires_window_scale=4, offset = 10):

        self.cov_fn = cov_fn
        assert n1%hires_window_scale==0, "Scale must divide n1 to fit pixels exactly"
        assert n2%hires_window_scale==0, "Scale must divide n2 to fit pixels exactly"
        self.n1 = n1
        self.n2 = n2
        self.window_size_ratio = hires_window_scale
        self.pixel_size_ratio = (self.window_size_ratio * self.n2) // self.n1
        self.offset = offset
        self.delta_low = 1./self.n1
        self.delta_high = 1./(self.n2 * self.window_size_ratio)
        self.k_low = scipy.fftpack.rfftfreq(self.n1,d=self.delta_low)
        self.k_high = scipy.fftpack.rfftfreq(self.n2,d=self.delta_high)

        self.C_low = self._get_variance_k(self.k_low) * self.n1
        self.C_high = self._get_variance_k(self.k_high) * self.n1 / self.window_size_ratio

        self.constraints =[]
        self.constraints_val = []
        self.constraints_real = []

    def get_default_plot_padding(self):
        """The default plot padding (in coarse pixels) to hide from the high-res region"""
        return 0

    # Covariance calculation utilities
    
    def set_Chigh_realspace(self):
        self.C_high = self._calc_transfer_fn_realspace(self.n2)

    def _calc_transfer_fn_realspace(self, num_pixels):
        fullbox_n2 = self.n1 * self.pixel_size_ratio
        k_high_full = scipy.fftpack.rfftfreq(fullbox_n2, d=1. / fullbox_n2)

        # pretend high-res is across the full box temporarily; get the TF

        transfer = self._cov_to_transfer(self._get_variance_k(k_high_full) * self.n1)
        # now truncate the TF to the correct size
        transfer_hi = np.concatenate((transfer[:num_pixels // 2], transfer[-num_pixels // 2:]))

        return self._transfer_to_cov(transfer_hi)

    def _cov_to_transfer(self, cov):
        # take real part of C_high (i.e. zero imaginary components)
        sqrt_C_high_re = np.sqrt(cov)
        sqrt_C_high_re[2::2] = 0

        # use this to calculate the real-space transfer function
        T_high = unitary_inverse_fft(sqrt_C_high_re)/np.sqrt(len(cov))

        return T_high

    def _transfer_to_cov(self, transfer):
        sqrt_C_high_apo_re = unitary_fft(transfer)* np.sqrt(len(transfer))

        # copy back imaginary parts ready for convolution
        if len(sqrt_C_high_apo_re) % 2 == 0:
            sqrt_C_high_apo_re[2::2] = sqrt_C_high_apo_re[1:-1:2]
        else:
            sqrt_C_high_apo_re[2::2] = sqrt_C_high_apo_re[1::2]

        return sqrt_C_high_apo_re ** 2

    def _get_variance_k(self, k, k_filter=lambda x:1):
        integrand = lambda k1: self.cov_fn(k1)*k_filter(k1)
        delta_k = k[1]-k[0]
        Cv = np.zeros_like(k)
        for i,ki in enumerate(k):
            if ki==0:
                lower_k = 0.0
            else:
                lower_k = ki-delta_k/2
            upper_k = ki + delta_k / 2
            Cv[i] = scipy.integrate.quad(integrand, lower_k, upper_k)[0]

        #if len(Cv)%2==0:
        #    Cv[-1]*=np.sqrt(2)

        return Cv


    def xs(self):
        """Return the real-space coordinates of the two outputs"""
        return (np.arange(self.n1)+0.5)/self.n1, \
               (self.offset + (np.arange(self.n2)+0.5)/self.pixel_size_ratio)/self.n1

    def boundary_xs(self):
        """Return the real-space coordinates of the pixel boundaries for the outputs"""
        centre1, centre2 = self.xs()
        return centre1-0.5/self.n1, centre2-0.5/(self.n1*self.pixel_size_ratio)


    def realization(self, verbose=False, no_random=False, white_noise_lo=None, white_noise_hi=None) -> [FFTArray, FFTArray] :
        """Generate a realization of the Gaussian random field with constraints, if present.
        
        This implementation calls a generic sequence of operations which are overriden in base class implementations to
        achieve different zoom techniques.
        
        :param verbose: if True, print more output
        :param no_random: if True, create a zeroed random field to see effect of constraints in isolation
        :param white_noise_lo: if present, the array of white noise for the low-frequency (n1) box
        :param white_noise_hi:  if present, the array of white noise for the high-frequency (n2) sub-box
        :return: delta_low, delta_high: the fields for the low-frequency (n1) box and high-frequency (n2) sub-box respectively
        """

        k_high, k_low = self._get_ks()

        if no_random:
            white_noise_lo = np.zeros_like(k_low).view(FFTArray)
            white_noise_hi = np.zeros_like(k_high).view(FFTArray)
        else:
            # Generate white-noise, or correctly encapsulate the white-noise samples passed in
            white_noise_lo, white_noise_hi = self._get_whitenoise(white_noise_lo, white_noise_hi)

        # In trad approach, apply constraints directly to the noise; covariance lives in the constraint covectors
        if self.constrain_noise_directly and len(self.constraints)>0:
            if verbose:
                print("Constraint NOISE")
            self._apply_constraints(white_noise_lo, white_noise_hi, verbose)

        # Traditional approaches make the low-frequency and high-frequency white noise compatible
        white_noise_lo, white_noise_hi = self._modify_whitenoise(white_noise_lo, white_noise_hi)

        # Apply the transfer function as stored in Fourier space as self.C_high and self.C_low
        delta_low, delta_high = self._apply_transfer_function(white_noise_lo,white_noise_hi)

        # Filtering approach now pulls only the high-frequency part from delta_high and low-f part from delta_low
        delta_low, delta_high, memos = self._separate_fields(delta_low, delta_high)

        # Now we can apply the constraints in the filter approach
        if (not self.constrain_noise_directly) and len(self.constraints)>0:
            if verbose:
                print("Constraint DELTA")
            delta_low, delta_high = self._apply_constraints(delta_low, delta_high, verbose)

        # Filtering approach now needs to recombine the low-frequency and high-frequency information for the final output
        delta_low, delta_high = self._recombine_fields(delta_low, delta_high, memos)

        return delta_low.in_real_space(), delta_high.in_real_space()


    def _get_whitenoise(self, white_lo=None, white_hi=None):
        if white_lo is None:
            white_lo, white_hi = np.random.normal(0.0,1.0,size=self.n1), np.random.normal(0.0,self.pixel_size_ratio**0.5,size=self.n2)
        # only one nyquist mode is included in DFT but two are physically present
        if len(white_lo) % 2 == 0:
            white_lo[-1] *= np.sqrt(2)
        if len(white_hi) % 2 == 0:
            white_hi[-1] *= np.sqrt(2)

        white_lo = FFTArray(np.copy(white_lo))
        white_hi = FFTArray(np.copy(white_hi))
        white_lo.fourier=white_hi.fourier=True

        return white_lo, white_hi
    
    @abc.abstractmethod
    def _modify_whitenoise(self, wn_lo, wn_hi):
        """Give subclass an opportunity to modify the whitenoise fields before the transfer function is applied
 
        :param wn_lo: the low-frequency pixelised whitenoise
        :param wn_hi: the high-frequency windowed whitenoise
        :return: wn_lo, wn_hi: the same fields, modified
        """
        pass

    @in_fourier_space
    def _apply_transfer_function(self, white_noise_lo, white_noise_hi=None):
        result_lo = FFTArray(white_noise_lo * np.sqrt(self.C_low))
        result_lo.fourier = True
        if white_noise_hi is not None:
            result_hi = FFTArray(white_noise_hi * np.sqrt(self.C_high) )
            result_hi.fourier = True
            return result_lo, result_hi
        else:
            return result_lo

    @abc.abstractmethod
    def _separate_fields(self, delta_low_k, delta_high_k):
        """Give subclass an opportunity to apply low and high-pass filters to the specified fields
        
        :param delta_low_k: the field in the pixelised basis, to be low-pass filtered
        :param delta_high_k: the field in the windowed region, to be high-pass filtered
        :return: delta_low_k, delta_high_k, memos - the two fields plus any additional information 
                 required for recombining the fields later
        """
        pass

    @abc.abstractmethod
    def _recombine_fields(self, delta_low: FFTArray, delta_high: FFTArray, memos):
        """Give subclass an opportunity to combine information from the low-frequency into the high-frequency field
        
        :param delta_low: the field in the pixelised basis, without any high-frequency information
        :param delta_high: the field in the window region, without any low-frequency information
        :param memos: any memos returned from _separate_fields at an earlier stage
        :return: delta_low, delta_high: the final version of the fields
        """
        pass

    @abc.abstractmethod
    def _apply_constraints(self, delta_low_k, delta_high_k, verbose):
        """Apply the constraints to the low-frequency and high-frequency fields"""

    def _get_ks(self):
        pixel_dx_low = 1. / self.n1
        pixel_dx_high = 1. / (self.n2 * self.window_size_ratio)
        k_low = scipy.fftpack.rfftfreq(self.n1, d=pixel_dx_low)
        k_high = scipy.fftpack.rfftfreq(self.n2, d=pixel_dx_high)
        return k_high, k_low

    def get_cov(self, one_element=None):
        """Get the exact covariance matrix using a deterministic algorithm rather than sampling (see estimate_cov for the latter)"""
        cov = np.zeros((self.n1 + self.n2, self.n1 + self.n2))

        element_iterator = self._iter_one_cov_element(one_element) if one_element is not None else self._iter_cov_elements()

        for test_field_lo, test_field_hi in element_iterator:
            out_lo, out_hi = self.realization(white_noise_lo=test_field_lo, white_noise_hi=test_field_hi)
            cov[:self.n1, :self.n1]+=np.outer(out_lo, out_lo)
            cov[self.n1:, self.n1:] += np.outer(out_hi, out_hi)
            cov[:self.n1, self.n1:] += np.outer(out_lo, out_hi)

        cov[self.n1:,:self.n1]=cov[:self.n1,self.n1:].T
        return cov

    def estimate_cov(self, Ntrials=2000, with_means=False):
        """Estimate the covariance matrix using the specified number of random trials
        
        :param with_means: if True, return also the mean field and the variance in the constrained values 
        (which should be near-zero)"""

        cov = np.zeros((self.n1+self.n2,self.n1+self.n2))
        cstr_means = np.zeros(len(self.constraints_real))
        constraint_variance = np.zeros(len(self.constraints_real))
        means = np.zeros(self.n2)
        for i in range(Ntrials):
            r1,r2 = self.realization()
            cov[:self.n1,:self.n1]+=np.outer(r1,r1)
            cov[:self.n1,self.n1:]+=np.outer(r1,r2)
            cov[self.n1:,self.n1:]+=np.outer(r2,r2)
            means+=r2
            cstr_means+=[np.dot(cstr,r2) for cstr in self.constraints_real]
            constraint_variance+=[np.dot(cstr,r2)**2 for cstr in self.constraints_real]

        cstr_means/=Ntrials
        constraint_variance/=Ntrials
        cov/=Ntrials
        means/=Ntrials

        constraint_variance-=cstr_means**2

        cov[self.n1:,:self.n1]=cov[:self.n1,self.n1:].T

        if with_means:
            return cov, means, np.sqrt(constraint_variance)
        else:
            return cov

    def _iter_cov_elements(self):
        """Helper for get_cov: iterates over sample 'white noise' fields such that the covariance is exactly the sum
        of the resulting output fields"""
        test_field_lo = np.zeros(self.n1)
        test_field_hi = np.zeros(self.n2)

        for i in range(self.n1):
            test_field_lo[i] = 1.0
            yield test_field_lo, test_field_hi
            test_field_lo[i] = 0.0

        for i in range(self.n2):
            test_field_hi[i] = self.pixel_size_ratio**(0.5)
            yield test_field_lo, test_field_hi
            test_field_hi[i] = 0.0

    def _iter_one_cov_element(self, offset):
        """Helper for get_cov testing: acts like _iter_cov_elements but only generates a single pixel. If offset is
        positive, this is a delta-function at the centre of the high-res box plus offset pixels. If offset is negative,
        it is a delta-function to the left of the high-res box by the specified amount."""
        test_field_lo = np.zeros(self.n1)
        test_field_hi = np.zeros(self.n2)

        if offset>=0:
            # place test delta-function in high-res region
            test_field_hi[self.n2//2+self.pixel_size_ratio//2+offset]=self.pixel_size_ratio**0.5
            test_field_lo = self.downsample(test_field_hi)
        else:
            # place in low-res region just next to high-res region
            test_field_lo[self.offset+offset]=1.0

        yield FFTArray(test_field_lo).in_fourier_space(), FFTArray(test_field_hi).in_fourier_space()







    @in_real_space
    def extract_window(self, delta_highres_unwindowed):
        return delta_highres_unwindowed[self.offset * self.pixel_size_ratio:self.offset * self.pixel_size_ratio + self.n2]

    @in_real_space
    def place_window(self, delta_highres_windowed: FFTArray) -> FFTArray:
        delta_highres_unwindowed = np.zeros(self.n1*self.pixel_size_ratio).view(FFTArray)
        assert delta_highres_unwindowed.fourier is False
        delta_highres_unwindowed[self.offset * self.pixel_size_ratio:self.offset * self.pixel_size_ratio + self.n2] = delta_highres_windowed
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

        delta_highres = np.zeros(self.n1*self.pixel_size_ratio)
        delta_low_left = np.roll(delta_low,1)
        delta_low_right = np.roll(delta_low,-1)

        for i in range(self.pixel_size_ratio):
            sub_offset = (float(i)+0.5)/self.pixel_size_ratio
            weight_cen = 1-abs(sub_offset-0.5)
            weight_left = 0.5-sub_offset
            if weight_left<0 : weight_left = 0
            weight_right = sub_offset-0.5
            if weight_right<0: weight_right = 0

            delta_highres[i::self.pixel_size_ratio] = delta_low * weight_cen
            delta_highres[i::self.pixel_size_ratio] += delta_low_left * weight_left
            delta_highres[i::self.pixel_size_ratio] += delta_low_right * weight_right



        result = delta_highres[self.offset*self.pixel_size_ratio:self.offset*self.pixel_size_ratio+self.n2]
        return result.view(type=FFTArray)

    @in_real_space
    def upsample_cubic(self, delta_low) -> FFTArray:
        "Take a low-res vector and interpolate it into the high-res region - cubic interpolation"

        x_vals_low, x_vals_high = self.xs()
        delta_highres = scipy.interpolate.interp1d(x_vals_low, delta_low, kind='cubic')(x_vals_high)

        return delta_highres.view(type=FFTArray)

    @in_real_space
    def downsample(self, hires_vector: FFTArray, in_window=True) -> FFTArray:
        """Take a high-res region vector and downsample it onto the low-res grid"""
        vec_lr = np.zeros(self.n1).view(type=FFTArray)
        if in_window:
            vec_lr[self.offset:self.offset+self.n2//self.pixel_size_ratio] = \
                  hires_vector.reshape((self.n2//self.pixel_size_ratio,self.pixel_size_ratio)).mean(axis=1)
        else:
            vec_lr = hires_vector.reshape((self.n1,self.pixel_size_ratio)).mean(axis=1)
        return vec_lr


    @abc.abstractmethod
    def add_constraint(self, val=0.0, hr_vec=None):
        """Add a constraint, specifying the constraint covector in the high-res region and the value it should take.

        If hr_vec is None, use a default constraint covector consisting of a delta function in the centre of the window"""

    def _default_constraint_hr_vec(self):
        hr_vec = np.zeros(self.n2)
        cval = self.n2 // 2
        hr_vec[cval+2] = 1.0
        print(hr_vec.sum())
        #hr_vec[cval+1] = -1.0
        return hr_vec


class UnfilteredZoomConstrained(ZoomConstrained):
    @in_real_space
    def _separate_fields(self, delta_low, delta_high):
        delta_high -= self.upsample_zeroorder(self.downsample(delta_high))
        return delta_low, delta_high, None

    @in_real_space
    def _recombine_fields(self, delta_low, delta_high, _):
        delta_high += self.upsample_zeroorder(delta_low)
        return delta_low, delta_high

    def _modify_whitenoise(self, wn_lo, wn_hi):
        return wn_lo, wn_hi

    def add_constraint(self, val=0.0, lr_covec=None):
        """Add a constraint to the LR field (note this differs from other classes which apply to the HR field)"""

        lr_covec = lr_covec.view(FFTArray).in_fourier_space()

        # perform Gram-Schmidt orthogonalization
        if len(self.constraints)!=0:
            raise RuntimeError("Orthogonalization not implemented yet for UnfilteredZoomConstrained!")

        norm = complex_dot(lr_covec, self.C_low*lr_covec)
        lr_covec /= math.sqrt(norm)
        val /= math.sqrt(norm)

        self.constraints.append(lr_covec)
        self.constraints_val.append(val)

    @in_fourier_space
    def _apply_constraints(self, delta_low_k, delta_high_k, verbose):
        # implementation works purely on low-res part; ignores window
        for al_low_k, d in zip(self.constraints, self.constraints_val):
            scale = d - complex_dot(al_low_k, delta_low_k)
            if verbose:
                print("scale=", scale)

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

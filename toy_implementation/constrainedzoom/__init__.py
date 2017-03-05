import numpy as np
import math
import scipy.fftpack
import scipy.integrate
import scipy.interpolate
import pylab as p
import functools
from . import fft_wrapper
import copy

# for simple reloading
reload(fft_wrapper)

from .fft_wrapper import *

"""
# old self-integrating version (deprecated because integral of filters is also important)
def cov(x,plaw = -1.5):
    if plaw ==0.0:
        return 1.0
    delta = x[1]-x[0]
    if plaw!=-1 :
        cv = abs(x**(plaw+1)-(x-delta)**(plaw+1))
    else:
        cv = np.log(x)-np.log(x-delta)
    cv[cv==np.inf]=0
    cv[cv!=cv]=0
    return cv
"""


def cov(x,plaw = -0.0,k0=0.5,k1=np.inf):
    try:
        cv = x**plaw
    except ZeroDivisionError:
        return 0.0
    except ValueError:
        return 0.0

    if plaw>-1.0:
        k0 = 0

    if isinstance(cv, np.ndarray):
        cv[cv==np.inf]=0
        cv[cv!=cv]=0
        cv[x>k1] = 0
    else:
        if x>k1 : cv = 0.0

    cv*=1./(1.+k0*x**(plaw-0.5))

    return cv



def complex_dot(x,y):
    """Dot product for packed FFT complex coeffs"""
    return np.dot(x,y) + np.dot(x[1:],y[1:]) # counts +ve and -ve modes

class ZoomConstrained(object):
    def __init__(self, cov_fn = cov, k_cut=0.3, n1=256, n2=768, hires_window_scale=4, offset = 10):

        self.cov_fn = cov_fn
        assert n1%hires_window_scale==0, "Scale must divide n1 to fit pixels exactly"
        assert n2%hires_window_scale==0, "Scale must divide n2 to fit pixels exactly"
        self.n1 = n1
        self.n2 = n2
        self.window_size_ratio = hires_window_scale
        self.offset = offset
        self.k_cut = self.n1*k_cut
        self.delta_low = 1./self.n1
        self.delta_high = 1./(self.n2 * self.window_size_ratio)
        self.k_low = scipy.fftpack.rfftfreq(self.n1,d=self.delta_low)
        self.k_high = scipy.fftpack.rfftfreq(self.n2,d=self.delta_high)

        # Naive approach:
        self.C_low = self._get_variance_k(self.k_low) * float(self.n1)
        self.C_high = self._get_variance_k(self.k_high) * float(self.n2)


        self.pixel_size_ratio = (self.window_size_ratio * self.n2) / self.n1
        self.constraints =[]
        self.constraints_val = []
        self.constraints_real = []


    def set_Chigh_realspace(self):
        fullbox_n2 = self.n1*self.pixel_size_ratio
        k_high_full = scipy.fftpack.rfftfreq(fullbox_n2,d=1./fullbox_n2)

        # pretend high-res is across the full box temporarily; get the TF

        transfer = self._cov_to_transfer(self._get_variance_k(k_high_full)*float(fullbox_n2))
        # now truncate the TF to the correct size
        transfer_hi = np.concatenate((transfer[:self.n2/2],transfer[-self.n2/2:]))

        self.C_high = self._transfer_to_cov(transfer_hi)


    def _cov_to_transfer(self, cov):
        # take real part of C_high (i.e. zero imaginary components)
        sqrt_C_high_re = np.sqrt(cov)
        sqrt_C_high_re[2::2] = 0

        # use this to calculate the real-space transfer function
        T_high = unitary_inverse_fft(sqrt_C_high_re)/np.sqrt(len(cov))

        return T_high

    def _transfer_to_cov(self, transfer):
        print len(transfer)
        sqrt_C_high_apo_re = unitary_fft(transfer)* np.sqrt(len(transfer))

        # copy back imaginary parts ready for convolution
        if len(sqrt_C_high_apo_re) % 2 == 0:
            sqrt_C_high_apo_re[2::2] = sqrt_C_high_apo_re[1:-1:2]
        else:
            sqrt_C_high_apo_re[2::2] = sqrt_C_high_apo_re[1::2]

        return sqrt_C_high_apo_re ** 2



    def filter_low(self, k):
        #return k<self.k_cut
        T = self.k_cut/10
        return 1./(1.+np.exp((k-self.k_cut)/T))

    def filter_high(self, k):
        return 1.-self.filter_low(k)

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

    _wn_psr_power = 0.5

    def _get_whitenoise(self, white_lo=None, white_hi=None):
        if white_lo is None:
            white_lo, white_hi = np.random.normal(0.0,1.0,size=self.n1), np.random.normal(0.0,self.pixel_size_ratio**self._wn_psr_power,size=self.n2)

        white_lo = FFTArray(np.copy(white_lo))
        white_hi = FFTArray(np.copy(white_hi))
        white_lo.fourier=white_hi.fourier=True

        # only one nyquist mode is included in DFT but two are physically present
        if len(white_lo) % 2 == 0:
            white_lo[-1] *= np.sqrt(2)
        if len(white_hi) % 2 == 0:
            white_hi[-1] *= np.sqrt(2)

        return white_lo, white_hi

    def _separate_whitenoise(self, wn_lo, wn_hi):
        return wn_lo, wn_hi

    @in_fourier_space
    def _apply_transfer_function(self, white_noise_lo, white_noise_hi=None):
        result_lo = FFTArray(white_noise_lo * np.sqrt(self.C_low))
        result_lo.fourier = True
        if white_noise_hi is not None:
            result_hi = FFTArray(white_noise_hi * np.sqrt(self.C_high)/self.pixel_size_ratio**self._wn_psr_power)
            result_hi.fourier = True
            return result_lo, result_hi
        else:
            return result_lo



    def xs(self):
        return np.arange(self.n1)+0.5, \
               self.offset + (np.arange(self.n2)+0.5)/self.pixel_size_ratio




    def realization(self,test=False,no_random=False,white_noise_lo=None,white_noise_hi=None):

        k_high, k_low = self._get_ks()

        if no_random:
            white_noise_lo = np.zeros_like(k_low)
            white_noise_hi = np.zeros_like(k_high)
        else:
            white_noise_lo, white_noise_hi = self._get_whitenoise(white_noise_lo, white_noise_hi)

        white_noise_lo, white_noise_hi = self._separate_whitenoise(white_noise_lo, white_noise_hi)

        delta_low, delta_high = self._apply_transfer_function(white_noise_lo,white_noise_hi)


        delta_low, delta_high, memos = self._separate_fields(delta_low, delta_high)

        delta_low, delta_high = self._apply_constraints(delta_low, delta_high, test)

        delta_low, delta_high = self._recombine_fields(delta_low, delta_high, memos)

        return delta_low.in_real_space(), delta_high.in_real_space()

    @in_real_space
    def _recombine_fields(self, delta_low, delta_high, delta_low_k_plus):
        delta_high += self.upsample_cubic(delta_low)
        delta_low += delta_low_k_plus.in_real_space()
        return delta_low, delta_high

    @in_fourier_space
    def _apply_constraints(self, delta_low_k, delta_high_k, verbose):
        for (al_low_k, al_high_k), d in zip(self.constraints, self.constraints_val):
            if self._real_space_filter:
                raise NotImplementedError, "Can't apply constraints when the filter is in real space"

            if verbose:
                print "lowdot=", complex_dot(al_low_k, delta_low_k) * self.pixel_size_ratio
                print "highdot=", complex_dot(al_high_k, delta_high_k)
                print "sum=", complex_dot(al_low_k, delta_low_k) * self.pixel_size_ratio + complex_dot(al_high_k,
                                                                                                       delta_high_k)
                al_low, al_high = self.harmonic_to_combined_pixel(al_low_k, al_high_k)
                print "RS simple dot=", np.dot(unitary_inverse_fft(al_high_k), unitary_inverse_fft(delta_high_k))
                print "RS dot=", np.dot(al_high, delta_high)

            scale = d - complex_dot(al_low_k, delta_low_k) * self.pixel_size_ratio - complex_dot(al_high_k,
                                                                                                 delta_high_k)

            if verbose:
                print "scale=", scale

            delta_low_k += self.C_low * al_low_k * scale
            delta_high_k += self.C_high * al_high_k * scale
        return delta_low_k, delta_high_k

    @in_fourier_space
    def _separate_fields(self, delta_low_k, delta_high_k):
        k_high, k_low = self._get_ks()
        delta_high_k *= np.sqrt(1. - self.filter_low(k_high) ** 2)  # keep original power spectrum
        delta_low_k_plus = delta_low_k * self.filter_high(
            k_low)  # store the filtered-out components of the low-res field
        delta_low_k *= self.filter_low(k_low)
        return delta_low_k, delta_high_k, delta_low_k_plus

    def _filter_fields_real(self, delta_low, delta_high):
        raise RuntimeError, "Incorrect call to _filter_fields_real - should be filtering in harmonic space"

    def _get_ks(self):
        pixel_dx_low = 1. / self.n1
        pixel_dx_high = 1. / (self.n2 * self.window_size_ratio)
        k_low = scipy.fftpack.rfftfreq(self.n1, d=pixel_dx_low)
        k_high = scipy.fftpack.rfftfreq(self.n2, d=pixel_dx_high)
        return k_high, k_low

    def estimate_cov(self, Ntrials=2000, with_means=False):
        cov = np.zeros((self.n1+self.n2,self.n1+self.n2))
        cstr_means = np.zeros(len(self.constraints_real))
        cstr_vars = np.zeros(len(self.constraints_real))
        means = np.zeros(self.n2)
        for i in xrange(Ntrials):
            r1,r2 = self.realization()
            cov[:self.n1,:self.n1]+=np.outer(r1,r1)
            cov[:self.n1,self.n1:]+=np.outer(r1,r2)
            cov[self.n1:,self.n1:]+=np.outer(r2,r2)
            means+=r2
            cstr_means+=[np.dot(cstr,r2) for cstr in self.constraints_real]
            cstr_vars+=[np.dot(cstr,r2)**2 for cstr in self.constraints_real]

        cstr_means/=Ntrials
        cstr_vars/=Ntrials
        cov/=Ntrials
        means/=Ntrials

        cstr_vars-=cstr_means**2

        cov[self.n1:,:self.n1]=cov[:self.n1,self.n1:].T

        if with_means:
            return cov, means, np.sqrt(cstr_vars)
        else:
            return cov

    def _iter_one_cov_element(self, offset):
        test_field_lo = np.zeros(self.n1)
        test_field_hi = np.zeros(self.n2)
        print self.pixel_size_ratio
        test_field_hi[self.n2/2+self.pixel_size_ratio/2+offset]=self.pixel_size_ratio**self._wn_psr_power
        test_field_lo = self.downsample(test_field_hi)

        yield FFTArray(test_field_lo).in_fourier_space(), FFTArray(test_field_hi).in_fourier_space()

    def _iter_cov_elements(self):
        test_field_lo = np.zeros(self.n1)
        test_field_hi = np.zeros(self.n2)

        for i in xrange(self.n1):
            test_field_lo[i] = 1.0
            yield test_field_lo, test_field_hi
            test_field_lo[i] = 0.0

        for i in xrange(self.n2):
            test_field_hi[i] = self.pixel_size_ratio**(self._wn_psr_power)
            yield test_field_lo, test_field_hi
            test_field_hi[i] = 0.0

    def get_cov(self, one_element=False):
        cov = np.zeros((self.n1 + self.n2, self.n1 + self.n2))

        element_iterator = self._iter_one_cov_element(one_element) if one_element else self._iter_cov_elements()

        for test_field_lo, test_field_hi in element_iterator:
            out_lo, out_hi = self.realization(white_noise_lo=test_field_lo, white_noise_hi=test_field_hi)
            cov[:self.n1, :self.n1]+=np.outer(out_lo, out_lo)
            cov[self.n1:, self.n1:] += np.outer(out_hi, out_hi)
            cov[:self.n1, self.n1:] += np.outer(out_lo, out_hi)

        cov[self.n1:,:self.n1]=cov[:self.n1,self.n1:].T
        return cov


    def hi_U(self):
        """Return the unitary matrix mapping pixel to harmonic space"""

        U = np.zeros((self.n2,self.n2))

        U_cplx = np.exp(- 1.j * 2 * np.outer(np.arange(self.n2/2+1),np.arange(self.n2)) * math.pi/self.n2)/np.sqrt(self.n2/2+1)
        U[2::2] = np.imag(U_cplx[1:-1])
        U[1::2] = np.real(U_cplx[1:])
        U[0] = np.real(U_cplx[0])

        return U


    def harmonic_to_combined_pixel(self, f_low_k, f_high_k):
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

    @in_real_space
    def upsample_zeroorder(self, delta_low):
        """Take a low-res vector and put it in the high-res region without interpolating"""


        delta_highres = np.zeros(self.n1 * self.pixel_size_ratio)
        delta_highres = delta_highres.view(type=FFTArray)
        delta_highres.fourier = False

        delta_low_left = np.roll(delta_low, 1)
        delta_low_right = np.roll(delta_low, -1)

        for i in range(self.pixel_size_ratio):
            delta_highres[i::self.pixel_size_ratio] = delta_low

        return delta_highres[self.offset * self.pixel_size_ratio:self.offset * self.pixel_size_ratio + self.n2]

    @in_real_space
    def upsample_linear(self, delta_low):
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
    def upsample_cubic(self, delta_low):
        "Take a low-res vector and interpolate it into the high-res region - cubic interpolation"

        x_vals_low, x_vals_high = self.xs()
        delta_highres = scipy.interpolate.interp1d(x_vals_low, delta_low, kind='cubic')(x_vals_high)

        return delta_highres.view(type=FFTArray)

    @in_real_space
    def downsample(self, hires_vector):
        """Take a high-res region vector and downsample it onto the low-res grid"""
        vec_lr = np.zeros(self.n1).view(type=FFTArray)
        vec_lr[self.offset:self.offset+self.n2/self.pixel_size_ratio] = \
              hires_vector.reshape((self.n2/self.pixel_size_ratio,self.pixel_size_ratio)).mean(axis=1)
        return vec_lr

    def high_k_vector_from_low_k_vector(self, low_harmonics):
        pixelized_highres = self.harmonic_to_combined_pixel(low_harmonics, None)[1]
        return unitary_fft(pixelized_highres)

    def hr_pixel_to_harmonic(self, vec=None):
        if vec is None:
            vec = np.zeros(self.n2)
            vec[self.n2/2]=1.0

        # FA . vec
        vec_k_high = unitary_fft(vec)*self.filter_high(self.k_high)

        vec_lr = self.downsample(vec-unitary_inverse_fft(vec_k_high))

        vec_k_low = unitary_fft(vec_lr)

        return vec_k_low, vec_k_high


    def norm(self, low, high):
        return self.C_weighted_inner_product(low,high,low,high)

    def inner_product(self, low1, high1, low2, high2, more_accurate=True):

        # The low.low part picks up a pixel_size_ratio factor. We can see this as follows.
        # Take the ideal case where low1 = f_low v1, high1 = f_high v1, f_low+f_high=1 and f_low f_high = 0.
        # Without varying pixel sizes, v1 = (low1+high1), and v1.v2 = (low1.low1)+(high1.high1), exactly.
        # Now, let's downsample the low1 pixel scale. The convention we've adopted is that the big
        # low1 pixels take on the mean value of the original, finer low1 pixels. So, the new dot product
        # in real space has been multiplied by 1/pixel_size_ratio. (Because FFTs are unitary, this
        # applies also to the harmonic space dot product). We need to multiply by pixel_size_ratio to
        # cancel out this change.
        product = self.pixel_size_ratio*complex_dot(low1,low2)+complex_dot(high1,high2)
        if more_accurate:
            # add in the low1^dagger C high2 + high1 C low2^dagger terms
            low1_as_high = self.high_k_vector_from_low_k_vector(low1)
            low2_as_high = self.high_k_vector_from_low_k_vector(low2)
            product+=complex_dot(low1_as_high,high2) + complex_dot(high1, low2_as_high)
        return product


    def C_weighted_inner_product(self, low1, high1, low2, high2):
        return self.inner_product(low1,high1,low2*self.C_low,high2*self.C_high)

    def add_constraint(self, val=0.0, hr_vec=None):


        self.constraints_real.append(hr_vec) # stored only for information - not part of the algorithm


        low, high = self.hr_pixel_to_harmonic(hr_vec)


        # perform Gram-Schmidt orthogonalization
        for (la, ha),va in zip(self.constraints,self.constraints_val):
            dotprod = self.xCy(la,ha,low,high)
            low-=dotprod*la
            high-=dotprod*ha
            val-=dotprod*va
            print self.xCy(la,ha,low,high)


        norm = self.norm(low,high)
        low/=math.sqrt(norm)
        high/=math.sqrt(norm)
        val/=math.sqrt(norm)

        self.constraints.append((low,high))
        self.constraints_val.append(val)





class IdealizedZoomConstrained(ZoomConstrained):
    """Calculate the low-res/high-res split by making a full box at the high resolution,
    then downgrading the resolution of the low-res region"""

    def __init__(self, cov_fn=cov, k_cut=0.2, n1=256, n2=768, hires_window_scale=4, offset=10):
        super(IdealizedZoomConstrained, self).__init__(cov_fn, k_cut, n1, n2, hires_window_scale, offset)
        self._underlying = ZoomConstrained(cov_fn, k_cut, n2*hires_window_scale, n2, hires_window_scale, offset)

    def _iter_cov_elements(self):
        test_field_lo = np.zeros(self._underlying.n1)
        test_field_hi = np.zeros(self._underlying.n2)

        for i in xrange(self._underlying.n1):
            test_field_lo[i] = 1.0
            yield test_field_lo, test_field_hi
            test_field_lo[i] = 0.0

    def _iter_one_cov_element(self, offset):
        test_field_lo = np.zeros(self._underlying.n1)
        test_field_hi = np.zeros(self._underlying.n2)

        test_field_lo[self.offset*self.pixel_size_ratio + self.n2 / 2 + self.pixel_size_ratio/2 + offset] = 1.0

        yield FFTArray(test_field_lo).in_fourier_space(), test_field_hi

    def realization(self, *args, **kwargs):
        underlying,_ = self._underlying.realization(*args, **kwargs)


        lores = underlying.reshape((self.n1,self.pixel_size_ratio)).mean(axis=1)
        hires = underlying[self.offset*self.pixel_size_ratio:self.offset*self.pixel_size_ratio+self.n2]
        return lores, hires


class UnfilteredZoomConstrained(ZoomConstrained):
    @in_real_space
    def _separate_fields(self, delta_low, delta_high):
        delta_high -= self.upsample_zeroorder(self.downsample(delta_high))
        return delta_low, delta_high, None

    @in_real_space
    def _recombine_fields(self, delta_low, delta_high, _):
        delta_high += self.upsample_zeroorder(delta_low)
        return delta_low, delta_high



class HybridZoomConstrained(ZoomConstrained):
    def __init__(self, *args, **kwargs):
        super(HybridZoomConstrained, self).__init__(*args, **kwargs)
        self.set_Chigh_realspace()


    @property
    def _B_window_slice(self):
        B_window_size = self.n2 / (self.pixel_size_ratio)
        offset_B = self.offset + self.n2 / self.pixel_size_ratio / 2 - B_window_size / 2
        return slice(offset_B,offset_B+B_window_size)

    @in_real_space
    def _separate_whitenoise(self, delta_low, delta_high):
        assert not delta_high.fourier
        assert not delta_low.fourier

        delta_high -= self.upsample_zeroorder(self.downsample(delta_high))


        lo_modes_for_hi_window = self.upsample_zeroorder(delta_low)

        delta_high += lo_modes_for_hi_window

        return delta_low, delta_high

    @in_fourier_space
    def _separate_fields(self, delta_low_k, delta_high_k):
        return delta_low_k, delta_high_k, 0

    def _recombine_fields(self, delta_low, delta_high, delta_low_plus):
        """
        delta_high -= self.upsample_zeroorder(self.downsample(delta_high))


        lo_modes_for_hi_window = self.upsample_zeroorder(delta_low).in_fourier_space()

        delta_high += lo_modes_for_hi_window

        """


        remove_modes = delta_high*self.filter_low(self.k_high)**2
        add_modes = self.upsample_cubic(delta_low.in_fourier_space()*self.filter_low(self.k_low)).in_fourier_space()*self.filter_low(self.k_high)
        delta_high.in_fourier_space()

        #delta_high = 0
        delta_high-=remove_modes
        delta_high+=add_modes


        return delta_low, delta_high



class HahnAbelZoomConstrained(UnfilteredZoomConstrained):
    def __init__(self, *args, **kwargs):
        super(HahnAbelZoomConstrained, self).__init__(*args, **kwargs)
        self.set_Chigh_realspace()
        #self._apodize_Chigh()

    @property
    def _B_window_slice(self):
        B_window_size = self.n2 / (self.pixel_size_ratio) / 2
        offset_B = self.offset + self.n2 / self.pixel_size_ratio / 2 - B_window_size / 2
        return slice(offset_B,offset_B+B_window_size)

    @in_real_space
    def _separate_whitenoise(self, delta_low, delta_high):
        assert not delta_high.fourier
        assert not delta_low.fourier


        delta_high -= self.upsample_zeroorder(self.downsample(delta_high))

        delta_low_zeroed = copy.copy(delta_low)

        delta_low_zeroed[self._B_window_slice] = 0


        self._delta_low_residual = delta_low-delta_low_zeroed
        lo_modes_for_hi_window = self.upsample_zeroorder(self._delta_low_residual)

        delta_high += lo_modes_for_hi_window

        return delta_low_zeroed, delta_high

    def _separate_fields(self, delta_low, delta_high):
        return delta_low, delta_high, None

    @in_real_space
    def _recombine_fields(self, delta_low, delta_high, _):
        delta_high += self.upsample_cubic(delta_low)
        delta_low += self._apply_transfer_function(self._delta_low_residual.in_fourier_space()).in_real_space()
        return delta_low, delta_high




def constraint_vector(scale=100,length=768,position=None) :
    """Generate a constraint vector corresponding to the Gaussian-filtered
    density at the given position."""
    if position is None :
        position = length/2

    pixel_vals = np.arange(0.,length)
    constraint = np.exp(-(pixel_vals-position)**2/(2*scale))
    constraint/=constraint.sum()
    return constraint

def deriv_constraint_vector(smooth=None,length=768,position=None) :
    """Constraint derivative at given position. First arg is ignored
    but could later define smoothing scale"""
    if position is None :
        position = length/2


    constraint = np.zeros(length)
    """
    constraint[position+1]=0.5
    constraint[position-1]=-0.5
    if smooth is not None:
        X = Ufft(constraint)
        k = scipy.fftpack.rfftfreq(length,d=1.0)
        X*=np.exp(-k**2*smooth)
        constraint = Uifft(X)
    """

    if smooth is None:
        smooth = 1

    constraint[position+smooth]=0.5
    constraint[position-smooth]=-0.5

    constraint/=np.sqrt(np.dot(constraint,constraint))
    return constraint

def display_cov(G, cov, downgrade=False, vmin=None, vmax=None, pad=0, show_hh=True):
    p.set_cmap('PuOr_r')
    vmin = vmin if vmin is not None else np.min(cov)
    vmax = vmax if vmax is not None else np.max(cov)

    C11 = cov[:G.n1,:G.n1]
    C22 = cov[G.n1:,G.n1:]
    C12 = cov[:G.n1,G.n1:]

    zoom_width = G.n1/G.window_size_ratio
    offset = G.offset
    n1 = G.n1

    if pad>0:
        pixel_scale = (G.n2*G.window_size_ratio)/G.n1
        pad_fine = pad*pixel_scale
        C22 = C22[pad_fine:-pad_fine,pad_fine:-pad_fine]
        C12 = C12[:,pad_fine:-pad_fine]
        zoom_width-=pad*2
        offset+=pad

    p.imshow(C11,extent=[0,G.n1,G.n1,0],vmin=vmin,vmax=vmax,interpolation='nearest')
    if downgrade:
        zoom_fac = G.window_size_ratio*(G.n2/G.n1)
        print "zoom_fac=",zoom_fac
        C22new=0
        for i in range(zoom_fac):
            for j in range(zoom_fac):
                C22new += C22[i::zoom_fac,j::zoom_fac]
        C22 = C22new/zoom_fac**2

        C12new=0
        for i in range(zoom_fac):
            C12new+=C12[:,i::zoom_fac]
        C12 = C12new/zoom_fac

    p.imshow(C12.T,extent=[0,n1,offset+zoom_width,offset],vmin=vmin,vmax=vmax,interpolation='nearest')
    if show_hh:
        p.imshow(C22,extent=[offset,offset+zoom_width,offset+zoom_width,offset],vmin=vmin,vmax=vmax,interpolation='nearest')

    p.plot([0,G.n1],[offset,offset],'k:')
    p.plot([0,G.n1],[offset+zoom_width,offset+zoom_width],'k:')
    if show_hh:
        p.plot([offset,offset],[0,G.n1],'k:')
        p.plot([offset+zoom_width,offset+zoom_width],[0,G.n1],'k:')

    if show_hh:
        p.text(offset+zoom_width,offset,'h-h',horizontalalignment='right',verticalalignment='top',color='black')
    p.text(G.n1,offset,'h-l',horizontalalignment='right',verticalalignment='top',color='black')
    p.text(G.n1,offset+zoom_width,'l-l',horizontalalignment='right',verticalalignment='top',color='black')
    if show_hh:
        p.text(offset+zoom_width,offset+zoom_width,'l-l',horizontalalignment='right',verticalalignment='top',color='black')
    p.xlim(0,G.n1)
    p.ylim(G.n1,0)


def cov2cor(matr):
    return matr/np.sqrt(np.outer(matr.diagonal(),matr.diagonal()))

def WC_vs_CW(plaw=0.0, k_cut = 0.2, part='lo', log=False):

    cov_this = functools.partial(cov,plaw=plaw)

    G = ZoomConstrained(cov_this,n2=256,k_cut=k_cut)


    U = G.hi_U()

    if part=='lo':
        cv_noW = G.get_lo_cov_in_subbox_ideal_limit()
        cv_W = G.get_lo_cov_in_subbox()
    elif part=='hi':
        cv_noW = G.get_hi_cov_ideal_limit()
        cv_W = G.get_hi_cov()
    elif part=='sum':
        cv_noW = G.get_hi_cov_ideal_limit() + G.get_lo_cov_in_subbox_ideal_limit()
        cv_W = G.get_hi_cov() + G.get_lo_cov_in_subbox()

    cv_noW_orig = cv_noW
    cv_W_orig = cv_W

    if log:
        cv_noW = np.log10(abs(cv_noW))
        cv_W = np.log10(abs(cv_W))
        vmax = cv_W.max()
        vmin = vmax-3

        vmax_diff = vmax
        vmin_diff = vmin

    else:
        vmin = cv_noW.min()
        vmax = cv_noW.max()

        vmin_diff = vmin/10
        vmax_diff = vmax/10
    #vmin = vmax = None

    p.subplot(231)

    p.imshow(cv_noW,vmin=vmin,vmax=vmax)
    p.ylabel("Real space")
    p.title("Ideal covariance")
    p.subplot(232)
    p.imshow(cv_W,vmin=vmin,vmax=vmax)
    p.title("Actual covariance")
    p.subplot(233)
    p.imshow(cv_W-cv_noW,vmin=vmin_diff,vmax=vmax_diff)
    p.title("Residuals x 10")
    p.draw()


    cv_noW = np.dot(U,np.dot(cv_noW_orig,U.T))
    cv_W = np.dot(U,np.dot(cv_W_orig,U.T))



    if log:
        cv_noW = np.log10(abs(cv_noW))
        cv_W = np.log10(abs(cv_W))
        vmax = cv_noW.max()
        vmin = vmax-3
        vmin_diff = vmin
        vmax_diff = vmax
    else:
        vmin = cv_noW.min()
        vmax = cv_noW.max()
        vmin_diff = vmin/10
        vmax_diff = vmax/10

    p.subplot(234)
    p.ylabel("Harmonic space")
    p.imshow(cv_noW,vmin=vmin,vmax=vmax)
    p.subplot(235)
    p.imshow(cv_W,vmin=vmin,vmax=vmax)
    p.subplot(236)
    p.imshow(cv_W-cv_noW,vmin=vmin_diff,vmax=vmax_diff)


def demo(val=2.0,seed=1,plaw=-1.5, deriv=False, n1=1024, n2=256, k_cut=0.2, scale=4, smooth=10):
    cv_gen = deriv_constraint_vector if deriv else constraint_vector
    cov_this = functools.partial(cov,plaw=plaw)

    # set up zoom solution
    if seed is not None:
        np.random.seed(seed)
    G = ZoomConstrained(cov_this, k_cut=k_cut, n2=n2, n1=n1, hires_window_scale=scale, offset=(n1*(scale-1))/(2*scale))

    G.add_constraint(val,cv_gen(smooth,n2))

    # set up ideal (uniform) solution
    n1_eff = n2*scale
    print "n1_eff=",n1_eff

    Gs = ZoomConstrained(cov_this, k_cut=10000,n2=n1_eff,n1=n1_eff, hires_window_scale=1,offset=0)
    Gs.add_constraint(val,cv_gen(smooth,n1_eff))
    _, r_ideal = Gs.realization(no_random=True)
    x_ideal = (np.arange(n1_eff)+0.5)/(n1_eff/n1)

    # make plots
    p.clf()
    p.subplot(211)
    p.plot(x_ideal,r_ideal,linewidth=4,color="#ffdddd")
    x0, x1 = G.xs()
    r0, r1 = G.realization(no_random=True)
    p.plot(x0,r0,'k:')
    p.plot(x1,r1,'k-')
    p.ylabel("Solution")

    p.subplot(212)
    ideal_sln = scipy.interpolate.interp1d(x_ideal,r_ideal)
    x0_max = max(abs(ideal_sln(x0)))
    x1_max = max(abs(ideal_sln(x1)))
    p.plot(x0,(r0-ideal_sln(x0))/x0_max,'k:')
    p.plot(x1,(r1-ideal_sln(x1))/x1_max,'k')
    p.ylim(-0.05,0.05)
    p.ylabel("$\Delta$ solution / max |solution|")
    return G, Gs

def cov_constraint_demo(downgrade_view=False,plaw=-1.5):
    cov_this = functools.partial(globals()['cov'],plaw=plaw)
    G = ZoomConstrained(cov_this,n2=256)
    #G.add_constraint(0.0,constraint_vector())
    #G.constraints_real.append(np.ones(768))
    cov, means, stds = G.estimate_cov(with_means=True)
    print "Mean of constraint:",means
    print "Std-dev of constraints:",stds
    display_cov(G, cov, downgrade_view)


def cov_zoom_demo(n1=64, n2=32, hires_window_scale=4, estimate=False, Ntrials=2000,
                  plaw=-1.5, cl=ZoomConstrained,pad=0,vmin=None,vmax=None,errors=False,
                  show_hh=True,k_cut=0.3,one_element=False):
    p.clf()
    cov_this = functools.partial(globals()['cov'], plaw=plaw)
    X = cl(cov_this, n1=n1, n2=n2, hires_window_scale=hires_window_scale, k_cut=k_cut)
    if estimate:
        cv_est = X.estimate_cov(Ntrials)
    else:
        cv_est = X.get_cov(one_element=one_element)

    if errors:
        Y = IdealizedZoomConstrained(cov_this, n1=n1, n2=n2, hires_window_scale=hires_window_scale)
        true_cov = Y.get_cov(one_element=one_element)
        cv_est-=true_cov
        cv_est/=true_cov.max() # variance in hi-res region

    display_cov(X, cv_est,pad=pad,vmin=vmin,vmax=vmax,show_hh=show_hh)
    p.colorbar()

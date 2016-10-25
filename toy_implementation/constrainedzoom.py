import numpy as np
import math
import scipy.fftpack
import scipy.integrate
import scipy.interpolate
import pylab as p
import functools

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


def cov(x,plaw = -0.0,k0=0.0):
    try:
        cv = x**plaw
    except ZeroDivisionError:
        return 0.0
    except ValueError:
        return 0.0

    if isinstance(cv, np.ndarray):
        cv[cv==np.inf]=0
        cv[cv!=cv]=0
        cv[x<k0] = 0.0
    else:
        if x<k0 : cv = 0.0

    return cv

def unitary_fft(x):
    return scipy.fftpack.rfft(x)/math.sqrt(float(len(x)))

def unitary_inverse_fft(x):
    return scipy.fftpack.irfft(x)*math.sqrt(float(len(x)))

def complex_dot(x,y):
    """Dot product for packed FFT complex coeffs"""
    return np.dot(x,y) + np.dot(x[1:],y[1:]) # counts +ve and -ve modes

class ZoomConstrained(object):
    _real_space_filter = False

    def __init__(self, cov_fn = cov, k_cut=0.2, n1=256, n2=768, hires_window_scale=4, offset = 10):

        self.cov_fn = cov_fn
        assert n1%hires_window_scale==0, "Scale must divide n1 to fit pixels exactly"
        assert n2%hires_window_scale==0, "Scale must divide n2 to fit pixels exactly"
        self.n1 = n1
        self.n2 = n2
        self.scale = hires_window_scale
        self.offset = offset
        self.k_cut = self.n1*k_cut
        self.weights = self.interpolation_coeffs()
        self.delta_low = 1./self.n1
        self.delta_high = 1./(self.n2*self.scale)
        self.k_low = scipy.fftpack.rfftfreq(self.n1,d=self.delta_low)
        self.k_high = scipy.fftpack.rfftfreq(self.n2,d=self.delta_high)

        # The covariances have a pixel size dependence that is automatically
        # included in get_cov, which integrates over the relevant part of the
        # power spectrum. They also have a window size dependence, because
        # <x_W x_W^dagger> = <x x^dagger> * W, where W is the fractional length
        # of the window. This can be most easily seen by thinking in real space.
        # So, only C_high picks up any explicit dependence here.
        self.C_low = self.get_variance_k(self.k_low)
        self.C_high = self.get_variance_k(self.k_high) / hires_window_scale

        self.pixel_size_ratio = (self.scale*self.n2)/self.n1
        self.constraints =[]
        self.constraints_val = []
        self.constraints_real = []

    def filter_low(self, k):
        #return k<self.k_cut
        T = self.k_cut/10
        return 1./(1.+np.exp((k-self.k_cut)/T))

    def filter_high(self, k):
        return 1.-self.filter_low(k)

    def C_filter_low(self,k):
        return self.filter_low(k)**2

    def C_filter_high(self,k):
        return 1.-self.C_filter_low(k)

    def get_variance_k(self, k, k_filter=lambda x:1):
        integrand = lambda k1: self.cov_fn(k1)*k_filter(k1)
        delta_k = k[1]-k[0]
        Cv = np.zeros_like(k)
        for i,ki in enumerate(k):
            if ki==0:
                Cv[i]=0
            else:
                upper_k = ki+delta_k/2
                Cv[i] = scipy.integrate.quad(integrand,ki-delta_k/2,upper_k)[0]
        return Cv


    def randoms(self, C):
        size = len(C)
        r = np.random.normal(0.0,1.0,size=size)
        r*=np.sqrt(C)
        return r

    def xs(self):
        return np.arange(self.n1)+0.5, \
               self.offset + (np.arange(self.n2)+0.5)/self.pixel_size_ratio



    def lowres_covariance(self):
        # calculate the pixel-space diagonal variance and just-off-diagonal covariance
        # for the lo-res grid
        k = scipy.fftpack.rfftfreq(self.n1,d=1./self.n1)
        C = self.get_variance_k(k)

        pa = np.zeros(self.n1)
        pb = np.zeros(self.n1)
        pa[0]=1
        pb[1]=1

        pa_k = unitary_fft(pa)
        pb_k = unitary_fft(pb)

        covm = np.zeros((2,2))

        covm[0,0] = complex_dot(pa_k,C*pa_k)
        covm[0,1] = complex_dot(pa_k,C*pb_k)
        covm[1,0] = covm[0,1]
        covm[1,1] = covm[0,0]

        return covm

    def interpolation_coeffs(self):
        grid_scalefactor = (self.scale*self.n2)/self.n1

        C = self.lowres_covariance()
        X = C[1,0]/C[0,0]

        # linear interpolation version:
        weights = np.arange(grid_scalefactor+1,dtype=float)/grid_scalefactor
        return weights
        """
        # attempt to keep power the same (but goes screwy immediately off-diagonal):

        # following is exact for X=0
        weights = (np.arange(grid_scalefactor+1,dtype=float)/grid_scalefactor)**0.5

        for i in range(10):
            xweights = weights
            weights = -weights[::-1] * X + np.sqrt(weights[::-1]**2*(X**2-1)+1)
            weights = (weights+xweights)/2


        return weights
        """


    def realization(self,test=False,no_random=False):

        k_high, k_low = self._get_ks()

        if no_random:
            delta_low_k = np.zeros_like(k_low)
            delta_high_k = np.zeros_like(k_high)
        else:
            delta_low_k = self.randoms(self.C_low)
            delta_high_k = self.randoms(self.C_high)

        if not self._real_space_filter:
            delta_high_k, delta_low_k, delta_low_k_plus = self._filter_fields(delta_high_k, delta_low_k)

        for (al_low_k,al_high_k), d in zip(self.constraints, self.constraints_val):
            if self._real_space_filter:
                raise NotImplementedError, "Can't apply constraints when the filter is in real space"

            if test:
                print "lowdot=",complex_dot(al_low_k,delta_low_k)*self.pixel_size_ratio
                print "highdot=", complex_dot(al_high_k, delta_high_k)
                print "sum=",complex_dot(al_low_k,delta_low_k)*self.pixel_size_ratio+complex_dot(al_high_k, delta_high_k)
                al_low, al_high = self.harmonic_to_combined_pixel(al_low_k, al_high_k)
                print "RS simple dot=",np.dot(unitary_inverse_fft(al_high_k),unitary_inverse_fft(delta_high_k))
                print "RS dot=",np.dot(al_high,delta_high)

            scale = d - complex_dot(al_low_k,delta_low_k)*self.pixel_size_ratio - complex_dot(al_high_k, delta_high_k)

            if test:
                print "scale=",scale

            delta_low_k += self.C_low*al_low_k * scale
            delta_high_k += self.C_high*al_high_k  * scale

        delta_low, delta_high = self.harmonic_to_pixel(delta_low_k,
                                                                delta_high_k)

        if self._real_space_filter:
            pass
            #delta_low, delta_high = self._filter_fields_real(delta_low, delta_high)
            #delta_high += self.upsample(delta_low)
        else:
            delta_high += self.upsample(delta_low)
            delta_low+=unitary_inverse_fft(delta_low_k_plus)



        return delta_low, delta_high

    def _filter_fields(self, delta_high_k, delta_low_k):
        k_high, k_low = self._get_ks()
        delta_high_k *= np.sqrt(1. - self.filter_low(k_high) ** 2)  # keep original power spectrum
        delta_low_k_plus = delta_low_k * self.filter_high(
            k_low)  # store the filtered-out components of the low-res field
        delta_low_k *= self.filter_low(k_low)
        return delta_high_k, delta_low_k, delta_low_k_plus

    def _filter_fields_real(self, delta_low, delta_high):
        raise RuntimeError, "Incorrect call to _filter_fields_real - should be filtering in harmonic space"

    def _get_ks(self):
        pixel_dx_low = 1. / self.n1
        pixel_dx_high = 1. / (self.n2 * self.scale)
        k_low = scipy.fftpack.rfftfreq(self.n1, d=pixel_dx_low)
        k_high = scipy.fftpack.rfftfreq(self.n2, d=pixel_dx_high)
        return k_high, k_low

    def estimate_cov(self, Ntrials=2000, with_means=False):
        cov = np.zeros((self.n1+self.n2,self.n1+self.n2))
        cstr_means = np.zeros(len(self.constraints_real))
        cstr_vars = np.zeros(len(self.constraints_real))
        for i in xrange(Ntrials):
            r1,r2 = self.realization()
            cov[:self.n1,:self.n1]+=np.outer(r1,r1)
            cov[:self.n1,self.n1:]+=np.outer(r1,r2)
            cov[self.n1:,self.n1:]+=np.outer(r2,r2)
            cstr_means+=[np.dot(cstr,r2) for cstr in self.constraints_real]
            cstr_vars+=[np.dot(cstr,r2)**2 for cstr in self.constraints_real]

        cstr_means/=Ntrials
        cstr_vars/=Ntrials
        cov/=Ntrials

        cstr_vars-=cstr_means**2

        cov[self.n1:,:self.n1]=cov[:self.n1,self.n1:].T

        if with_means:
            return cov, cstr_means, np.sqrt(cstr_vars)
        else:
            return cov

    def get_cov(self, ideal=False):
        cov = np.zeros((self.n1 + self.n2, self.n1 + self.n2))

        if ideal:
            cov[:self.n1, :self.n1] = self.get_lo_cov_in_superbox()
            cov[self.n1:, self.n1:] = self.get_hi_cov_ideal_limit() + self.get_lo_cov_in_subbox_ideal_limit()
        else:
            cov[:self.n1, :self.n1] = self.get_lo_cov_in_superbox()
            locov_sub, locov_cross = self.get_lo_cov_in_subbox(with_cross=True)
            cov[:self.n1, self.n1:] = locov_cross
            cov[self.n1:, :self.n1] = locov_cross.T
            cov[self.n1:, self.n1:] = self.get_hi_cov() + locov_sub

        return cov

    def get_lo_cov_in_subbox(self, white=False, with_cross=False):
        """Returns the covariance matrix for the low-res part of the
        high-res window"""

        if self._real_space_filter:
            Cs = self.C_low
        else:
            Cs = self.C_low * self.filter_low(self.k_low) ** 2
        C_sub = np.zeros((self.n2,self.n2))
        C_cross = np.zeros((self.n1,self.n2))
        for i,Ci in enumerate(Cs):
            T = np.zeros(self.n1)
            T[i]=1.0
            if white:
                Ci = 1.0
            pspace_super, pspace_sub = self.harmonic_to_combined_pixel(T, None)
            if self._real_space_filter:
                pspace_super, pspace_sub = self._filter_fields_real(pspace_super, pspace_sub)
                pspace_sub+=self.upsample(pspace_super)
            C_sub+=Ci*np.outer(pspace_sub,pspace_sub)
            C_cross+=Ci*np.outer(pspace_super, pspace_sub)

        if with_cross:
            return C_sub, C_cross
        else:
            return C_sub


    def get_hi_cov(self,white=False):
        """Returns the covariance matrix for the high-res segment.

        Constraints are ignored."""

        C = np.zeros((self.n2,self.n2))
        dummy_superbox = np.zeros(self.n1)

        if self._real_space_filter:
            Cs = self.C_high
        else:
            Cs = self.C_high * (1 - self.filter_low(self.k_high) ** 2)

        for i,Ci in enumerate(Cs):
            T = np.zeros(self.n2)
            T[i]=1.0
            T = unitary_inverse_fft(T)
            if white:
                Ci = 1.0

            if self._real_space_filter:
                _, T = self._filter_fields_real(dummy_superbox,T)
            C+=Ci*np.outer(T,T)

        return C

    def get_lo_cov_in_superbox(self,white=False):
        """Returns the covariance matrix for the low-res box.

        Constraints are ignored."""

        C = np.zeros((self.n1,self.n1))
        Cs = self.C_low
        for i,Ci in enumerate(Cs):
            T = np.zeros(self.n1)
            T[i]=1.0
            T = unitary_inverse_fft(T)
            if white:
                Ci = 1.0
            C+=Ci*np.outer(T,T)

        return C


    def get_lo_cov_in_subbox_ideal_limit(self):
        """Returns the covariance matrix for the low-res part of the
        high-res window, IN THE IDEAL CASE"""

        return self.get_hi_cov_ideal_limit(lo=True)

    def get_hi_cov_ideal_limit(self, white=False, lo=False):
        """Returns the covariance matrix for the high-res segment as
        calculated for the FULL box size (i.e. the ideal case where we could
        actually calculate the whole box at the full resolution, then just
        extract a small segment)

        Constraints are ignored."""
        C = np.zeros((self.n2*self.pixel_size_ratio,self.n2*self.pixel_size_ratio))

        super_k = scipy.fftpack.rfftfreq(self.n2*self.pixel_size_ratio,d=self.delta_high)
        if lo:
            Cs = self.get_variance_k(super_k, self.C_filter_low)
        else:
            Cs = self.get_variance_k(super_k, self.C_filter_high)

        for i,Ci in enumerate(Cs):
            T = np.zeros(self.n2*self.pixel_size_ratio)
            T[i]=1.0
            T = unitary_inverse_fft(T)
            if white:
                Ci = 1.0
            else:
                Ci*=self.pixel_size_ratio
            C+=Ci*np.outer(T,T)

        return C[self.offset*self.pixel_size_ratio:self.offset*self.pixel_size_ratio+self.n2,
                 self.offset*self.pixel_size_ratio:self.offset*self.pixel_size_ratio+self.n2]

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
        delta_lowhigh = self.upsample(delta_low)

        delta_high+=delta_lowhigh

        return delta_low, delta_high

    def harmonic_to_pixel(self, f_low_k, f_high_k):
        delta_low = unitary_inverse_fft(f_low_k)  # *self.filter_low(self.k_low))
        if f_high_k is not None:
            delta_high = unitary_inverse_fft(f_high_k)
        else:
            delta_high = np.zeros(self.n2)

        return delta_low, delta_high

    def upsample(self, delta_low):
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



        return delta_highres[self.offset*self.pixel_size_ratio:self.offset*self.pixel_size_ratio+self.n2]

    def downsample(self, hires_vector):
        """Take a high-res region vector and downsample it onto the low-res grid"""
        vec_lr = np.zeros(self.n1)
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


class UnfilteredZoomConstrained(ZoomConstrained):
    _real_space_filter = True

    def filter_low(self, k):
        raise RuntimeError, "Not available - unfiltered"

    def filter_high(self,k ):
        raise RuntimeError, "Not available - unfiltered"

    def C_filter_low(self, k):
        raise RuntimeError, "Not available - unfiltered"

    def C_filter_high(self,k):
        raise RuntimeError, "Not available - unfiltered"

    def _filter_fields_real(self, delta_low, delta_high):
        delta_high -= self.upsample(self.downsample(delta_high))
        return delta_low, delta_high

    def _filter_fields(self, delta_high_k, delta_low_k):
        raise RuntimeError, "Incorrect call to _filter_fields - should be filtering in real space"

    def upsample(self, delta_low):
        """Take a low-res vector and interpolate it into the high-res region"""

        delta_highres = np.zeros(self.n1*self.pixel_size_ratio)


        for i in range(self.pixel_size_ratio):
            delta_highres[i::self.pixel_size_ratio] = delta_low



        return delta_highres[self.offset*self.pixel_size_ratio:self.offset*self.pixel_size_ratio+self.n2]



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

def display_cov(G, cov, downgrade=False, vmin=None, vmax=None):
    vmin = vmin or np.min(cov)
    vmax = vmax or np.max(cov)

    C11 = cov[:G.n1,:G.n1]
    C22 = cov[G.n1:,G.n1:]
    C12 = cov[:G.n1,G.n1:]

    zoom_width = G.n1/G.scale

    p.imshow(C11,extent=[0,G.n1,G.n1,0],vmin=vmin,vmax=vmax,interpolation='nearest')
    if downgrade:
        zoom_fac = G.scale*(G.n2/G.n1)
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

    p.imshow(C12.T,extent=[0,G.n1,G.offset+zoom_width,G.offset],vmin=vmin,vmax=vmax,interpolation='nearest')
    p.imshow(C22,extent=[G.offset,G.offset+zoom_width,G.offset+zoom_width,G.offset],vmin=vmin,vmax=vmax,interpolation='nearest')

    p.plot([0,G.n1],[G.offset,G.offset],'w:')
    p.plot([0,G.n1],[G.offset+zoom_width,G.offset+zoom_width],'w:')
    p.plot([G.offset,G.offset],[0,G.n1],'w:')
    p.plot([G.offset+zoom_width,G.offset+zoom_width],[0,G.n1],'w:')

    p.text(G.offset+zoom_width,G.offset,'h-h',horizontalalignment='right',verticalalignment='top',color='white')
    p.text(G.n1,G.offset,'h-l',horizontalalignment='right',verticalalignment='top',color='white')
    p.text(G.n1,G.offset+zoom_width,'l-l',horizontalalignment='right',verticalalignment='top',color='white')
    p.text(G.offset+zoom_width,G.offset+zoom_width,'l-l',horizontalalignment='right',verticalalignment='top',color='white')
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

def cov_demo(downgrade_view=False,plaw=-1.5):
    cov_this = functools.partial(globals()['cov'],plaw=plaw)
    G = ZoomConstrained(cov_this,n2=256)
    #G.add_constraint(0.0,constraint_vector())
    #G.constraints_real.append(np.ones(768))
    cov, means, stds = G.estimate_cov(with_means=True)
    print "Mean of constraint:",means
    print "Std-dev of constraints:",stds
    display_cov(G, cov, downgrade_view)

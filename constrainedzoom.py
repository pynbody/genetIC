import numpy as np
import numpy.linalg
import math
import scipy.fftpack
import scipy.integrate
import pylab as p
import copy
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


def cov(x,plaw = -1.5):
    try:
        cv = x**plaw
    except ZeroDivisionError:
        return 0.0
    except ValueError:
        return 0.0

    if isinstance(cv, np.ndarray):
        cv[cv==np.inf]=0
        cv[cv!=cv]=0

    return cv

def Ufft(x):
    """Unitary FFT"""
    return scipy.fftpack.rfft(x)/math.sqrt(float(len(x)))

def Uifft(x):
    """Unitary inverse FFT"""
    return scipy.fftpack.irfft(x)*math.sqrt(float(len(x)))

def complex_dot(x,y):
    """Dot product for packed FFT complex coeffs"""
    return np.dot(x,y) + np.dot(x[1:],y[1:]) # counts +ve and -ve modes

class ZoomConstrained(object):
    def __init__(self, cov_fn = cov, k_cut=0.2, n1=256, n2=768, scale=4, offset = 10):
        self.cov_fn = cov_fn
        self.n1 = n1
        self.n2 = n2
        self.scale = scale
        self.offset = offset
        self.k_cut = self.n1*k_cut
        self.weights = self.interpolation_coeffs()
        self.delta_low = 1./self.n1
        self.delta_high = 1./(self.n2*self.scale)
        self.k_low = scipy.fftpack.rfftfreq(self.n1,d=self.delta_low)
        self.k_high = scipy.fftpack.rfftfreq(self.n2,d=self.delta_high)
        self.C_low = self.get_cov(self.k_low)
        self.C_high = self.get_cov(self.k_high)
        self.zoom_fac = self.scale*(self.n2/self.n1)
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

    def get_cov(self, k, k_filter=lambda x:1):
        integrand = lambda k1: self.cov_fn(k1)*k_filter(k1)
        delta_k = k[1]-k[0]
        Cv = np.zeros_like(k)
        for i,ki in enumerate(k):
            if ki==0:
                Cv[i]=0
            else:
                Cv[i] = scipy.integrate.quad(integrand,ki-delta_k/2,ki+delta_k/2)[0]
        return Cv


    def randoms(self, k):
        size = len(k)
        r = np.random.normal(0.0,1.0,size=size)
        r*=np.sqrt(self.get_cov(k)*size)
        return r

    def xs(self):
        delta_fine = float(self.n1)/(self.scale*self.n2)
        return np.arange(self.n1)+0.5, \
               np.linspace(self.offset+delta_fine,
                           self.offset+self.n1/self.scale-delta_fine,self.n2)


    def lowres_covariance(self):
        # calculate the pixel-space diagonal variance and just-off-diagonal covariance
        # for the lo-res grid
        k = scipy.fftpack.rfftfreq(self.n1,d=1./self.n1)
        C = self.get_cov(k)

        pa = np.zeros(self.n1)
        pb = np.zeros(self.n1)
        pa[0]=1
        pb[1]=1

        pa_k = Ufft(pa)
        pb_k = Ufft(pb)

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


    def realization(self,test=False):

        grid_scalefactor = (self.scale*self.n2)/self.n1


        delta_low = 1./self.n1
        delta_high = 1./(self.n2*self.scale)
        k_low = scipy.fftpack.rfftfreq(self.n1,d=delta_low)
        k_high = scipy.fftpack.rfftfreq(self.n2,d=delta_high)


        delta_low_k = self.randoms(k_low)
        delta_high_k = self.randoms(k_high)

        delta_high_k*=np.sqrt(1.-self.filter_low(k_high)**2) # keep original power spectrum
        delta_low_k_plus = delta_low_k*self.filter_high(k_low) # store the filtered-out components of the low-res field
        delta_low_k*=self.filter_low(k_low)

        # TEST - no low-k component
        # delta_low_k[:]=0

        # THIS SHOULD GO ultimately to save time:
        if test:
            delta_low, delta_high = self.harmonic_to_pixel(delta_low_k,
                                                           delta_high_k)




        for (al_low_k,al_high_k), d in zip(self.constraints, self.constraints_val):
            #al_low_k[:]=0

            if test:
                print "lowdot=",complex_dot(al_low_k,delta_low_k)*self.zoom_fac
                print "highdot=", complex_dot(al_high_k, delta_high_k)
                print "sum=",complex_dot(al_low_k,delta_low_k)*self.zoom_fac+complex_dot(al_high_k, delta_high_k)
                al_low, al_high = self.harmonic_to_pixel(al_low_k,al_high_k)
                print "RS simple dot=",np.dot(Uifft(al_high_k),Uifft(delta_high_k))
                print "RS dot=",np.dot(al_high,delta_high)

            scale = d - complex_dot(al_low_k,delta_low_k)*grid_scalefactor - complex_dot(al_high_k, delta_high_k)

            if test:
                print "scale=",scale

            delta_low_k += self.C_low*al_low_k * scale
            delta_high_k += self.C_high*al_high_k  * scale




        delta_low, delta_high = self.harmonic_to_pixel(delta_low_k,
                                                       delta_high_k)


        # not sure how good this is - but add in the high-k power into the low-res region
        delta_low2 = Uifft(delta_low_k_plus)

        return delta_low2+delta_low, delta_high

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

        cstr_vars-=cstr_means**2

        cov[self.n1:,:self.n1]=cov[:self.n1,self.n1:].T

        if with_means:
            return cov, cstr_means, np.sqrt(cstr_vars)
        else:
            return cov


    def get_lo_cov(self,white=False):
        """Returns the covariance matrix for the low-res part of the
        high-res window"""

        Cs = self.C_low*self.filter_low(self.k_low)**2
        C = np.zeros((self.n2,self.n2))

        for i,Ci in enumerate(Cs):
            T = np.zeros(self.n1)
            T[i]=1.0
            if white:
                Ci = 1.0
            pspace = self.harmonic_to_pixel(T,None)[1]
            C+=Ci*np.outer(pspace,pspace)

        return C



    def get_lo_cov_unwindowed(self):
        """Returns the covariance matrix for the low-res part of the
        high-res window, IN THE IDEAL CASE"""

        return self.get_hi_cov_unwindowed(lo=True)

    def get_hi_cov(self,white=False):
        """Returns the covariance matrix for the high-res segment.

        Constraints are ignored."""

        C = np.zeros((self.n2,self.n2))
        Cs = self.C_high*(1-self.filter_low(self.k_high)**2)
        for i,Ci in enumerate(Cs):
            T = np.zeros(self.n2)
            T[i]=1.0
            T = Uifft(T)
            if white:
                Ci = 1.0
            C+=Ci*np.outer(T,T)

        return C

    def get_hi_cov_unwindowed(self, white=False, lo=False):
        """Returns the covariance matrix for the high-res segment as
        calculated for the FULL box size (i.e. the ideal case where we could
        actually calculate the whole box at the full resolution, then just
        extract a small segment)

        Constraints are ignored."""
        C = np.zeros((self.n2*self.zoom_fac,self.n2*self.zoom_fac))

        super_k = scipy.fftpack.rfftfreq(self.n2*self.zoom_fac,d=self.delta_high)
        if lo:
            Cs = self.get_cov(super_k,self.C_filter_low)
        else:
            Cs = self.get_cov(super_k,self.C_filter_high)

        for i,Ci in enumerate(Cs):
            T = np.zeros(self.n2*self.zoom_fac)
            T[i]=1.0
            T = Uifft(T)
            if white:
                Ci = 1.0
            else:
                Ci*=self.zoom_fac
            C+=Ci*np.outer(T,T)

        return C[self.offset*self.zoom_fac:self.offset*self.zoom_fac+self.n2,
                 self.offset*self.zoom_fac:self.offset*self.zoom_fac+self.n2]

    def hi_U(self):
        """Return the unitary matrix mapping pixel to harmonic space"""

        U = np.zeros((self.n2,self.n2))

        U_cplx = np.exp(- 1.j * 2 * np.outer(np.arange(self.n2/2+1),np.arange(self.n2)) * math.pi/self.n2)/np.sqrt(self.n2/2+1)
        print U_cplx.shape
        U[2::2] = np.imag(U_cplx[1:-1])
        U[1::2] = np.real(U_cplx[1:])
        U[0] = np.real(U_cplx[0])

        return U


    def harmonic_to_pixel(self, f_low_k, f_high_k):
        delta_low = Uifft(f_low_k*self.filter_low(self.k_low))
        if f_high_k is not None:
            delta_high = Uifft(f_high_k)
        else:
            delta_high = np.zeros(self.n2)

        # delta_high does not contain low-frequency contributions.

        # interpolate the low frequency contribution into the high-res window
        # prepare the data range just around the window ...

        delta_low_window = delta_low[self.offset-1:self.offset+self.n1/self.scale+1]
        delta_lowhigh = np.zeros(self.n2+self.zoom_fac)

        #  ... and interpolate linearly ...

        for i in range(self.zoom_fac):
            weight_hi = self.weights[i]
            weight_lo = self.weights[self.zoom_fac-i]
            delta_lowhigh[i::self.zoom_fac] = weight_lo*delta_low_window[:-1] +\
                weight_hi*delta_low_window[1:]

        # trim to size
        upper = -self.zoom_fac/2+1
        if upper==0: upper = None
        delta_lowhigh = delta_lowhigh[self.zoom_fac/2+1:upper]

        return delta_low, delta_lowhigh+delta_high

    def hr_pixel_to_harmonic(self, vec=None):
        if vec is None:
            vec = np.zeros(self.n2)
            vec[self.n2/2]=1.0

        # FA . vec
        vec_k_high = Ufft(vec)

        # down-sample
        vec_lr = np.zeros(self.n1)
        vec_lr[self.offset:self.offset+self.n2/self.zoom_fac] = \
              vec.reshape((self.n2/self.zoom_fac,self.zoom_fac)).mean(axis=1)

        vec_k_low = Ufft(vec_lr)

        return vec_k_low*self.filter_low(self.k_low), \
               vec_k_high*self.filter_high(self.k_high)

        # FB . vec

    def norm(self, low, high):
        return self.zoom_fac*complex_dot(low,low*self.C_low)+complex_dot(high,high*self.C_high)

    def xCy(self, low1, high1, low2, high2):
        return self.zoom_fac*complex_dot(low1,low2*self.C_low)+complex_dot(high1,high2*self.C_high)

    def add_constraint(self, val=0.0, hr_vec=None):


        self.constraints_real.append(hr_vec) # stored only for information - not part of the algorithm


        low, high = self.hr_pixel_to_harmonic(hr_vec)


        # perform G-S

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


def constraint_vector(scale=100,length=768,position=None) :
    """Generate a constraint vector corresponding to the Gaussian-filtered
    density at the given position."""
    if position is None :
        position = length/2

    pixel_vals = np.arange(0.,length)
    constraint = np.exp(-(pixel_vals-position)**2/(2*scale))
    constraint/=constraint.sum()
    return constraint


def display_cov(G, cov, downgrade=False):
    vmin = np.min(cov)
    vmax = np.max(cov)

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

def WC_vs_CW(plaw=0.0, k_cut = 0.2, lo=False):

    cov_this = functools.partial(cov,plaw=plaw)

    G = ZoomConstrained(cov_this,n2=256,k_cut=k_cut)


    U = G.hi_U()

    if lo:
        cv_noW = G.get_lo_cov_unwindowed()
        cv_W = G.get_lo_cov()
    else:
        cv_noW = G.get_hi_cov_unwindowed()
        cv_W = G.get_hi_cov()

    vmin = cv_noW.min()
    vmax = cv_noW.max()
    #vmin = vmax = None

    p.subplot(231)

    p.imshow(cv_noW,vmin=vmin,vmax=vmax)
    p.ylabel("Real space")
    p.title("Ideal HR covariance")
    p.subplot(232)
    p.imshow(cv_W,vmin=vmin,vmax=vmax)
    p.title("Actual HR covariance")
    p.subplot(233)
    p.imshow(cv_W-cv_noW,vmin=vmin/10,vmax=vmax/10)
    p.title("Residuals x 10")
    p.draw()
    cv_noW = np.dot(U,np.dot(cv_noW,U.T))
    cv_W = np.dot(U,np.dot(cv_W,U.T))
    vmin = cv_noW.min()
    vmax = cv_noW.max()
    p.subplot(234)
    p.ylabel("Harmonic space")
    p.imshow(cv_noW,vmin=vmin,vmax=vmax)
    p.subplot(235)
    p.imshow(cv_W,vmin=vmin,vmax=vmax)
    p.subplot(236)
    p.imshow(cv_W-cv_noW,vmin=vmin/10,vmax=vmax/10)


def demo(val=2.0,seed=1,plaw=-1.5):
    cov_this = functools.partial(cov,plaw=plaw)
    if seed is not None:
        np.random.seed(seed)
    G = ZoomConstrained(cov_this)
    #G.add_constraint(val,constraint_vector())
    #G.add_constraint(-val,constraint_vector(500))
    x0, x1 = G.xs()
    r0, r1 = G.realization()
    p.plot(x0,r0,':')
    p.plot(x1,r1,'.')
    #p.plot([42.05,42.05],[-20,20])


def cov_demo(downgrade_view=False,plaw=-1.5):
    cov_this = functools.partial(globals()['cov'],plaw=plaw)
    G = ZoomConstrained(cov_this,n2=256)
    #G.add_constraint(0.0,constraint_vector())
    #G.constraints_real.append(np.ones(768))
    cov, means, stds = G.estimate_cov(with_means=True)
    print "Mean of constraint:",means
    print "Std-dev of constraints:",stds
    display_cov(G, cov, downgrade_view)

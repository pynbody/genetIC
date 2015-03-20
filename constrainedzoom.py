import numpy as np
import numpy.linalg
import math
import scipy.fftpack
import pylab as p

def cov(x,plaw = -1.0):
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

def Ufft(x):
    """Unitary FFT"""
    return scipy.fftpack.rfft(x)/math.sqrt(float(len(x)))

def Uifft(x):
    """Unitary inverse FFT"""
    return scipy.fftpack.irfft(x)*math.sqrt(float(len(x)))

class ZoomConstrained(object):
    def __init__(self, cov_fn = cov, n1=256, n2=768, scale=4, offset = 10):
        self.cov_fn = cov_fn
        self.n1 = n1
        self.n2 = n2
        self.scale = scale
        self.offset = offset
        self.k_cut = self.n1*0.05
        self.weights = self.interpolation_coeffs()
        delta_low = 1./self.n1
        delta_high = 1./(self.n2*self.scale)
        self.k_low = scipy.fftpack.rfftfreq(self.n1,d=delta_low)
        self.k_high = scipy.fftpack.rfftfreq(self.n2,d=delta_high)
        self.C_low = self.cov_fn(self.k_low)
        self.C_high = self.cov_fn(self.k_high)
        self.zoom_fac = self.scale*(self.n2/self.n1)
        self.constraints =[]
        self.constraints_val = []

    def filter_low(self, k):
        return k<self.k_cut

    def filter_high(self, k):
        return k>=self.k_cut

    def randoms(self, k):
        size = len(k)
        r = np.random.normal(0.0,1.0,size=size)
        r*=np.sqrt(self.cov_fn(k)*size)
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
        C = self.cov_fn(k)

        pa = np.zeros(self.n1)
        pb = np.zeros(self.n1)
        pa[0]=1
        pb[1]=1

        pa_k = Ufft(pa)
        pb_k = Ufft(pb)

        covm = np.zeros((2,2))

        covm[0,0] = np.dot(pa_k,C*pa_k)
        covm[0,1] = np.dot(pa_k,C*pb_k)
        covm[1,0] = covm[0,1]
        covm[1,1] = covm[0,0]

        return covm

    def interpolation_coeffs(self):
        grid_scalefactor = (self.scale*self.n2)/self.n1

        C = self.lowres_covariance()
        X = C[1,0]/C[0,0]

        # following is exact for X=0
        weights = (np.arange(grid_scalefactor+1,dtype=float)/grid_scalefactor)**0.5

        print X
        for i in range(10):
            xweights = weights
            weights = -weights[::-1] * X + np.sqrt(weights[::-1]**2*(X**2-1)+1)
            weights = (weights+xweights)/2


        return np.arange(grid_scalefactor+1,dtype=float)/grid_scalefactor


    def realization(self):

        grid_scalefactor = (self.scale*self.n2)/self.n1


        delta_low = 1./self.n1
        delta_high = 1./(self.n2*self.scale)
        k_low = scipy.fftpack.rfftfreq(self.n1,d=delta_low)
        k_high = scipy.fftpack.rfftfreq(self.n2,d=delta_high)


        delta_low_k = self.randoms(k_low)
        delta_high_k = self.randoms(k_high)

        delta_high_k*=self.filter_high(k_high)
        delta_low_k_plus = delta_low_k+self.filter_high(k_low) # store the filtered-out components of the low-res field
        delta_low_k*=self.filter_low(k_low)


        for (al_low,al_high), d in zip(self.constraints, self.constraints_val):
            scale = d - np.dot(al_low,delta_low_k) - np.dot(al_high, delta_high_k)
            delta_low_k += self.C_low*al_low * scale
            delta_high_k += self.C_high*al_high  * scale
            print np.dot(al_low,delta_low_k)+np.dot(al_high,delta_high_k),d

        delta_low, delta_high = self.harmonic_to_pixel(delta_low_k,
                                                       delta_high_k)




        # not sure how good this is - but add in the high-k power into the low-res region
        delta_low2 = 1.0*Uifft(delta_low_k_plus)

        return delta_low2+delta_low, delta_high

    def estimate_cov(self, Ntrials=2000):
        cov = np.zeros((self.n1+self.n2,self.n1+self.n2))
        for i in xrange(Ntrials):
            r1,r2 = self.realization()
            cov[:self.n1,:self.n1]+=np.outer(r1,r1)
            cov[:self.n1,self.n1:]+=np.outer(r1,r2)
            cov[self.n1:,self.n1:]+=np.outer(r2,r2)

        cov[self.n1:,:self.n1]=cov[:self.n1,self.n1:].T
        return cov

    def harmonic_to_pixel(self, f_low_k, f_high_k):
        delta_low = Uifft(f_low_k*self.filter_low(self.k_low))
        delta_high = Uifft(f_high_k)

        #  interpolate the low frequency contribution into the high-res window
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
        return np.dot(low,low*self.C_low)+np.dot(high,high*self.C_high)

    def xCy(self, low1, high1, low2, high2):
        return np.dot(low1,low2*self.C_low)+np.dot(high1,high2*self.c_high)

    def add_constraint(self, val=0.0, hr_vec=None):
        if len(self.constraints)>0:
            print "oops - G-S not implemented yet for multiple constraints"

        low, high = self.hr_pixel_to_harmonic(hr_vec)
        norm = self.norm(low,high)
        low/=math.sqrt(norm)
        high/=math.sqrt(norm)
        val/=math.sqrt(norm)

        self.constraints.append((low,high))
        self.constraints_val.append(val)



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

    p.imshow(C12.T,extent=[0,G.n1,G.offset+zoom_width,G.offset],vmin=vmin,vmax=vmax,interpolation='nearest')
    p.imshow(C22,extent=[G.offset,G.offset+zoom_width,G.offset+zoom_width,G.offset],vmin=vmin,vmax=vmax,interpolation='nearest')

    p.plot([0,G.n1],[G.offset,G.offset],'w:')
    p.plot([0,G.n1],[G.offset+zoom_width,G.offset+zoom_width],'w:')
    p.plot([G.offset,G.offset],[0,G.n1],'w:')
    p.plot([G.offset+zoom_width,G.offset+zoom_width],[0,G.n1],'w:')

    p.xlim(0,G.n1)
    p.ylim(G.n1,0)


def demo(val=2.0):
    np.random.seed(1)
    G = ZoomConstrained()
    G.add_constraint(val)
    x0, x1 = G.xs()
    r0, r1 = G.realization()
    p.plot(x0,r0)
    p.plot(x1,r1)


def cov_demo():
    G = ZoomConstrained()

    cov = G.estimate_cov()
    display_cov(G, cov)

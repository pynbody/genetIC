import numpy as np
import numpy.linalg
import math
import scipy.fftpack
import pylab as p

def cov(x):
    #return 1.
    delta = x[1]-x[0]
    cv = x**0.5-(x-delta)**0.5
    cv[cv==np.inf]=0
    cv[cv!=cv]=0
    return cv

class ZoomConstrained(object):
    def __init__(self, cov_fn = cov, n1=256, n2=768, scale=4, offset = 10):
        self.cov_fn = cov_fn
        self.n1 = n1
        self.n2 = n2
        self.scale = scale
        self.offset = offset
        self.k_cut = self.n1*0.1

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


    def realization(self):

        grid_scalefactor = (self.scale*self.n2)/self.n1


        delta_low = 1./self.n1
        delta_high = 1./(self.n2*self.scale)
        k_low = scipy.fftpack.rfftfreq(self.n1,d=delta_low)
        k_high = scipy.fftpack.rfftfreq(self.n2,d=delta_high)


        delta_low_k = self.randoms(k_low)
        delta_high_k = self.randoms(k_high)

        delta_high_k*=self.filter_high(k_high)

        delta_low = scipy.fftpack.irfft(delta_low_k*self.filter_low(k_low))
        delta_low2 = scipy.fftpack.irfft(delta_low_k)
        delta_high = scipy.fftpack.irfft(delta_high_k)

        # now interpolate the low frequency contribution into the high-res window
        # prepare the data range just around the window ...

        delta_low_window = delta_low[self.offset-1:self.offset+self.n1/self.scale+1]
        delta_lowhigh = np.zeros(self.n2+grid_scalefactor)

        #  ... and interpolate linearly ...

        for i in range(grid_scalefactor):
            weight_hi = float(i)/grid_scalefactor
            weight_lo = 1.0-weight_hi
            delta_lowhigh[i::grid_scalefactor] = weight_lo*delta_low_window[:-1] +\
                weight_hi*delta_low_window[1:]

        # trim to size
        upper = -grid_scalefactor/2+1
        if upper==0: upper = None
        delta_lowhigh = delta_lowhigh[grid_scalefactor/2+1:upper]

        return delta_low2, delta_high+delta_lowhigh

    def estimate_cov(self, Ntrials=2000):
        cov = np.zeros((self.n1+self.n2,self.n1+self.n2))
        for i in xrange(Ntrials):
            r1,r2 = self.realization()
            cov[:self.n1,:self.n1]+=np.outer(r1,r1)
            cov[:self.n1,self.n1:]+=np.outer(r1,r2)
            cov[self.n1:,self.n1:]+=np.outer(r2,r2)

        cov[self.n1:,:self.n1]=cov[:self.n1,self.n1:].T
        return cov


def demo():
    #np.random.seed(1)
    G = ZoomConstrained()

    x0, x1 = G.xs()
    r0,r1 = G.realization()
    p.plot(x0,r0)
    p.plot(x1,r1)


def cov_demo():
    G = ZoomConstrained()

    cov = G.estimate_cov()

    vmin = np.min(cov)
    vmax = np.max(cov)

    C11 = cov[:G.n1,:G.n1]
    C22 = cov[G.n1:,G.n1:]
    C12 = cov[:G.n1,G.n1:]

    zoom_width = G.n1/G.scale

    p.imshow(C11,extent=[0,G.n1,G.n1,0],vmin=vmin,vmax=vmax,interpolation='nearest')
    p.imshow(C12.T,extent=[0,G.n1,G.offset+zoom_width,G.offset],vmin=vmin,vmax=vmax,interpolation='nearest')
    p.imshow(C22,extent=[G.offset,G.offset+zoom_width,G.offset+zoom_width,G.offset],vmin=vmin,vmax=vmax,interpolation='nearest')

    p.plot([0,G.n1],[G.offset,G.offset],'w:')
    p.plot([0,G.n1],[G.offset+zoom_width,G.offset+zoom_width],'w:')
    p.plot([G.offset,G.offset],[0,G.n1],'w:')
    p.plot([G.offset+zoom_width,G.offset+zoom_width],[0,G.n1],'w:')

    p.xlim(0,G.n1)
    p.ylim(G.n1,0)

import numpy as np
import numpy.linalg
import math
import scipy.fftpack
import pylab as p

def cov(x):
    return 1.0
    cv = x**-0.75
    cv[cv==np.inf]=0
    cv[cv!=cv]=0
    return cv

class MultiscaleGaussian(object):
    def __init__(self, cov_fn=cov, n_refine=0, size=256):
        self.cov_fn = cov
        self.n_refine = n_refine
        self.size = size


    def realization(self, level=0, means=None):
        r = np.random.normal(0.0,1.0,size=self.size)
        r*=np.sqrt(cov(scipy.fftpack.rfftfreq(self.size,d=0.5**(level)/self.size)))
        r*=np.exp(np.random.uniform(0.0,2*math.pi,size=self.size)*1.j)
        r = scipy.fftpack.irfft(r)/np.sqrt(self.size)


        """
        if means is not None:
            means/=np.sqrt(2)
            means_realized = (r[::2]+r[1::2])/2
            r[::2]-=means_realized
            r[1::2]-=means_realized
            r[::2]+=means[:self.size/2]
            r[1::2]+=means[:self.size/2]
        """

        if level==self.n_refine:
            return [r]
        else:
            return [r]+self.realization(level+1,r)

    def xs(self):
        xs = []
        for i in range(self.n_refine+1):
            xs.append(np.arange(self.size,dtype=float)/2**i)
        return xs

    def estimate_cov(self, level, Ntrials=10000):
        print scipy.fftpack.rfftfreq(self.size,d=2.0**level)
        cov = np.zeros((self.size,self.size))
        for i in xrange(Ntrials):
            r = self.realization()[level]
            cov+=np.outer(r,r)
        return cov


def demo():
    G = MultiscaleGaussian(cov, 1)

    G2 = MultiscaleGaussian(cov, 0,512)

    x0, x1 = G.xs()
    r0,r1 = G.realization()
    p.plot(x0,r0)
    p.plot(x1,r1)

    x0 = G2.xs()[0]/2
    r0, = G2.realization()
    p.plot(x0,r0)

def demo_cov():
    G = MultiscaleGaussian(cov, 1)
    cv = G.estimate_cov(1)
    p.subplot(121)
    vmin = cv.min()
    vmax = cv.max()
    p.imshow(cv,vmin=vmin,vmax=vmax)

    G = MultiscaleGaussian(cov,0, 512)
    cv = G.estimate_cov(0)

    p.subplot(122)
    p.imshow(cv,vmin=vmin,vmax=vmax)
    p.xlim(0,256)
    p.ylim(256,0)

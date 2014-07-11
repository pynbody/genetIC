import numpy as np
import numpy.linalg
import math

class Gaussian(object) :
    def __init__(self,C,x0=None) :
        """Initialize a Gaussian random field generator with mean x0
        and covariance C.

        If x0 is None, the Gaussian is taken to be zero-mean"""
        C = np.asanyarray(C)
        if x0 is not None :
            x0 = np.asanyarray(x0)
        else :
            x0 = np.zeros(len(C))
        assert C.shape[0]==C.shape[1]==len(x0), "Length of x0 and shape of C do not match"
        
     

        vals, vecs = np.linalg.eigh(C)
        # sometimes numerical errors give negative eigenvalue that should really just be zero.
        # Fix it within some tolerance
        tol = abs(vals.mean())*1.e-8
        vals[(vals>-tol) * (vals<0)] = 0
        
        assert (vals>=0).all(), "Not a covariance matrix - some eigenvalues are negative"
        
        self._vals = vals
        self._vecs = vecs
        self._C = C
        self._x0 = x0

    def realization(self) :
        realization = np.zeros_like(self._vals)
        realization[self._vals>0] = np.random.normal(0.0,np.sqrt(self._vals[self._vals>0]))

        return np.dot(realization,self._vecs.T)+self._x0

    def constrained(self,alpha, d) :
        """Return a new Gaussian object subject to the given linear constraint.

        i.e. the return Gaussian will always realize x such that
        alpha.x = d.  The constraint vector does not need to be
        normalized correctly; this is done internally.
        """

        alpha = np.array(alpha) # take copy as we'll be modifying this
        C = self._C
        x0 = self._x0

        assert len(alpha)==len(x0), "Constraint vector does not match length of data"

        # Normalize to our convention:
        norm = np.dot(alpha,np.dot(C,alpha))
        assert norm!=0, "Attempting to construct an overconstrained system"
        
        alpha/=np.sqrt(norm)
        d/=np.sqrt(norm)

        x1 = x0 - np.dot(np.dot(C,np.outer(alpha,alpha)),x0) + d * np.dot(C, alpha)
        C1 = C - np.dot(C,np.dot(np.outer(alpha,alpha),C))
        return Gaussian(C1,x1)

    def iterative_constrained(self, alpha, d, niter=4) :
        """Return a new object subject to the given linear constraint, using
        an iterative algorithm that never explicitly constructs the new
        covariance matrix."""
        return IterativeConstrainedGaussian(self, alpha, d, niter)

    
class IterativeConstrainedGaussian(object) :
    def __init__(self, underlying, constraint, value, niter) :
        """Initialize a constrained Gaussian without ever calculating
        its covariance explicitly.

        *args*

        underlying - the underlying Gaussian object from which realizations
                     will be drawn

        constraint - the constraint vector

        value - the value the vector should take (so that data.constraint = value)

        niter - the number of iterations to use when calculating"""
        self.underlying = underlying

        # normalize the constraint
        norm = np.dot(constraint,np.dot(underlying._C,constraint))
    
        self.cv = constraint/np.sqrt(norm)
        self.x1 = underlying._x0 - np.dot(np.dot(underlying._C,np.outer(self.cv,self.cv)),underlying._x0) + value * np.dot(underlying._C, self.cv)
        self.niter = niter

    def realization(self) :
        
        def calculate(i,d_vec) :
            AATvec = self.cv*np.dot(d_vec,self.cv)
            if i==0 :
                return -0.5*AATvec
                        
            C_Mi_d = np.dot(self.underlying._C,calculate(i-1,d_vec))
            Mi_C_Mi_d = calculate(i-1,C_Mi_d)

            return -0.5*(AATvec + Mi_C_Mi_d)

        data = self.underlying.realization() - self.underlying._x0
        return data + np.dot(self.underlying._C,calculate(self.niter,data)) + self.x1
        #return data - self.underlying._C

    


def pixel_to_harmonic_matrix(length=500) :
    """Generate the unitary transformation matrix between pixel and harmonic space"""
    p_vec = np.arange(0.,length)
    return np.exp(2.j*math.pi*np.outer(p_vec,p_vec)/length)/np.sqrt(length)

def powerlaw_covariance(power_index=-0.5, length=500) :
    """Generate a covariance matrix for a given power law spectrum"""
    M = pixel_to_harmonic_matrix(length)
    pspec = np.arange(0.,length)**power_index
    pspec[0]=0

    C = np.dot(M.conj().T,np.dot(np.diag(pspec),M)).real
    return C

def constraint_vector(scale=100,length=500,position=None) :
    """Generate a constraint vector corresponding to the Gaussian-filtered
    density at the given position."""
    if position is None :
        position = length/2

    pixel_vals = np.arange(0.,length)
    constraint = np.exp(-(pixel_vals-position)**2/(2*scale))
    constraint/=constraint.sum()
    return constraint


def estimate_covariance(G, nsamples=1000) :
    est_cov = 0
    for n in xrange(nsamples) :
        d = G.realization()
        est_cov+=np.outer(d,d)/nsamples
    return est_cov

                                    
def demo1() :
    import pylab as p
    cv1 = constraint_vector(150)
    cv2 = constraint_vector(20)
    
    G = Gaussian(powerlaw_covariance())
    G1 = G.constrained(cv1,-0.20)
    G2 = G1.constrained(cv2,+0.20)
    
    realization = G.realization()
    p.plot(realization,label="Unconstrained")
    print "Unconstrained - values of constraints are: ",realization.dot(cv1), realization.dot(cv2)

    realization = G1.realization()
    p.plot(realization,label="Large-scale void")
    print "One constraint - values of constraints are: ",realization.dot(cv1), realization.dot(cv2)
     

    realization = G2.realization()
    p.plot(realization,label="+Local density enhancement")
    print "Two constraints - values of constraints are: ",realization.dot(cv1), realization.dot(cv2)
    p.legend()


def demo2(nsamples=10000) :
    import pylab as p
    cv1 = constraint_vector(50)
    G = Gaussian(powerlaw_covariance(-1.0))
    G1 = G.constrained(cv1,0.0)
    

    p.subplot(241)
    vmin = G1._C.min()
    vmax = G1._C.max()
    p.imshow(G1._C,vmin=vmin,vmax=vmax)
    p.title("Exact covariance")
    p.draw()

    
    p.subplot(242)
    cov = estimate_covariance(G1,nsamples=nsamples)
    p.imshow(cov, vmin=vmin,vmax=vmax)
    p.xticks([])
    p.yticks([])
    p.title("MC exact")
    
    p.subplot(246)
    p.imshow(cov-G1._C, vmin=-vmax/10,vmax=vmax/10)
    p.xticks([])
    p.yticks([])
    p.ylabel("Differences (x10)")
    p.draw()
    
    p.subplot(243)
    G1_approx = G.iterative_constrained(cv1,0.0,niter=5)
    cov = estimate_covariance(G1_approx,nsamples=nsamples)
    p.imshow(cov, vmin=vmin,vmax=vmax)
    p.xticks([])
    p.yticks([])
    p.title("MC iterative\n iterations = 5")

    p.subplot(247)
    p.imshow(cov-G1._C, vmin=-vmax/10,vmax=vmax/10)
    p.xticks([])
    p.yticks([])
    p.draw()
    
    p.subplot(244)
    G1_approx = G.iterative_constrained(cv1,0.0,niter=0)
    cov = estimate_covariance(G1_approx,nsamples=nsamples)
    p.imshow(cov, vmin=vmin,vmax=vmax)
    p.xticks([])
    p.yticks([])
    p.title("MC iterative\n iterations = 0")
    p.colorbar().set_label("C colourscale")

    p.subplot(248)
    p.imshow(cov-G1._C, vmin=-vmax/10,vmax=vmax/10)
    p.draw()
    p.xticks([])
    p.yticks([])
    p.colorbar().set_label("diff colourscale")

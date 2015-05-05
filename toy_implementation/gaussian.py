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

    def realization(self,underlying_realization=None) :
        if underlying_realization is not None:
            print "Warning - no underlying field description, ignoring underlying_realization"
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
        covariance matrix.

        This no longer seems to be a useful technique -- superceded by
        projection_constrained -- but is kept here for reference."""
        return IterativeConstrainedGaussian(self, alpha, d, niter)

    def projection_constrained(self, alpha, d):
        """Return a new object subject to the given linear constraint, using
        the non-iterative projection method."""
        return ProjectionConstrainedGaussian(self, alpha, d)

    def cov_dot(self, vec):
        """Return the covariance matrix applied to a given vector. Note that
        in some subclasses the covariance matrix may not be known explicitly,
        so the appropriate calculation will occur."""

        return np.dot(self._C,vec)

    def get_mean(self):
        return self._x0

    def chi2(self, vec):
        return np.dot(vec,np.dot(np.linalg.inv(self._C),vec))


    def get_chi2_matrix(self, *constraints):
        n_constraints = len(constraints)
        A = np.zeros((n_constraints,n_constraints))
        for i in range(n_constraints):
            for j in range(n_constraints):
                A[i,j] = np.dot(constraints[i],self.cov_dot(constraints[j]))

        return np.linalg.inv(A)


    def build_cov(self):
        """Work out the covariance matrix, even if it's not explicitly stored"""

        N = len(self.get_mean())
        cov = np.zeros((N,N))
        v = np.zeros(N)
        for i in xrange(N):
            v[i]=1.0
            cov[i,:]=self.cov_dot(v)
            v[i]=0.0
        return cov


class ProjectionConstrainedGaussian(Gaussian) :
    def __init__(self, underlying, constraint, value) :
        """Initialize a constrained Gaussian without ever calculating
        its covariance explicitly.

        *args*

        underlying - the underlying Gaussian object from which realizations
                     will be drawn

        constraint - the constraint vector

        value - the value the vector should take (so that data.constraint = value)
        """
        self.underlying = underlying

        # normalize the constraint
        norm = np.dot(constraint,underlying.cov_dot(constraint))
        self.alpha = constraint/np.sqrt(norm)
        self.C0_dot_alpha = underlying.cov_dot(self.alpha)


        x0 = underlying.get_mean()
        self.x1 = x0 \
                  - self.C0_dot_alpha*np.dot(self.alpha,x0) \
                  + self.C0_dot_alpha*(value/np.sqrt(norm))



    def cov_dot(self, vec):
        C0_dot_vec = self.underlying.cov_dot(vec)
        return C0_dot_vec - self.C0_dot_alpha*np.dot(self.alpha,C0_dot_vec)

    def get_mean(self):
        return self.x1

    def realization(self,underlying_realization=None) :

        if underlying_realization is None:
            underlying_realization = self.underlying.realization()

        data = underlying_realization - self.underlying.get_mean()
        """
        print "alpha.data=",np.dot(self.alpha,data)
        print "C0_dot_alpha norm=",self.underlying.chi2(self.C0_dot_alpha)
        """


        """
        u = self.underlying
        CdA = np.dot(u._C,self.alpha)
        print "OR",np.dot(CdA,np.dot(np.linalg.inv(u._C),CdA))
        print "OR",np.dot(self.alpha,np.dot(u._C,np.dot(np.linalg.inv(u._C),CdA)))
        print "alpha C0 alpha = ",np.dot(self.alpha,CdA)
        """

        return data - np.dot(self.alpha,data)*self.C0_dot_alpha + self.x1


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

    def realization(self,underlying_realization=None) :

        def calculate(i,d_vec) :
            AATvec = self.cv*np.dot(d_vec,self.cv)
            if i==0 :
                return -1.0*AATvec

            C_Mi_d = np.dot(self.underlying._C,calculate(i-1,d_vec))
            Mi_C_Mi_d = calculate(i-1,C_Mi_d)

            return -0.5*(AATvec + Mi_C_Mi_d)

        if underlying_realization is None:
            underlying_realization = self.underlying.realization()

        data = underlying_realization - self.underlying._x0
        return data + np.dot(self.underlying._C,calculate(self.niter,data)) + self.x1
        #return data - self.underlying._C




def pixel_to_harmonic_matrix(length=500) :
    """Generate the unitary transformation matrix between pixel and harmonic space"""
    p_vec = np.arange(0.,length)
    return np.exp(2.j*math.pi*np.outer(p_vec,p_vec)/length)/np.sqrt(length)

def powerlaw_covariance(power_index=-0.5, length=500, cutoff=None) :
    """Generate a covariance matrix for a given power law spectrum"""
    M = pixel_to_harmonic_matrix(length)
    pspec = np.arange(0.,length)**power_index
    pspec[0]=pspec[1]

    if cutoff:
        pspec*=1./(np.exp((np.arange(0.,length)-cutoff)/(50))+1)


    C = np.dot(M.conj().T,np.dot(np.diag(pspec),M)).real


    return C

def random_covariance(length=500):
    import scipy
    A = scipy.random.rand(length,length)
    return np.dot(A,A.transpose())

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


def demo1(projection=False) :
    """Demonstration of a double-constrained problem"""

    import pylab as p
    cv1 = constraint_vector(150)
    cv2 = constraint_vector(20)

    G = Gaussian(powerlaw_covariance())
    if projection:
        G1 = G.projection_constrained(cv1,-2.0)
        G2 = G1.projection_constrained(cv2,2.0)
    else:
        G1 = G.constrained(cv1,-2.0)
        G2 = G1.constrained(cv2,0.0)

    realization = G.realization()
    p.plot(realization,label="Unconstrained")
    print "Unconstrained - values of constraints are: ",realization.dot(cv1), realization.dot(cv2)
    print "                covariance along constraints:",np.dot(cv1,G1.cov_dot(cv1)), np.dot(cv2,G1.cov_dot(cv2))


    realization = G1.realization(underlying_realization=realization)


    p.plot(realization,label="Large-scale void")
    print "One constraint - values of constraints are: ",realization.dot(cv1), realization.dot(cv2)
    print "                covariance along constraints:",np.dot(cv1,G1.cov_dot(cv1)), np.dot(cv2,G1.cov_dot(cv2))

    realization = G2.realization(underlying_realization=realization)
    p.plot(realization,label="+Local density enhancement")
    print "Two constraints - values of constraints are: ",realization.dot(cv1), realization.dot(cv2)
    print "                covariance along constraints:",np.dot(cv1,G2.cov_dot(cv1)), np.dot(cv2,G2.cov_dot(cv2))
    p.legend()


def demo2(nsamples=10000) :
    """Demonstration of the covariance matrices"""

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


def demo3(width=50,plaw=-1.0,maxiter=10,size=500):
    """Demonstration of the convergence properties"""

    import pylab as p

    if plaw is None:
        C = random_covariance(size)
    else:
        C = powerlaw_covariance(plaw,size)

    G = Gaussian(C)

    if width is None:
        cv1 = np.random.uniform(size=size)
    else:
        cv1 = constraint_vector(width)

    R0 = G.realization()
    orig = np.dot(cv1,R0)

    lc = p.plot([-1],[1.0],'x')[0].get_color()

    for niter in range(maxiter):
        GI = G.iterative_constrained(cv1,0.0,niter=niter)
        R1 = GI.realization(underlying_realization=R0)
        p.plot([niter],[np.dot(cv1,R1)/orig],lc+'x')
        p.xlim(-2,maxiter)
        p.ylim(-0.5,1.5)
        p.draw()



def demo4(width=50,plaw=-1.0,niter=4,size=500):
    """Comparison of exact vs explicit iterative covariance matrix"""

    import pylab as p

    if plaw is None:
        C = random_covariance(size)
    else:
        C = powerlaw_covariance(plaw,size)

    G = Gaussian(C)

    if width is None:
        cv1 = np.random.uniform(size=size)
    else:
        cv1 = constraint_vector(width)

    Citer = G.projection_constrained(cv1, 0).build_cov()
    Cactual = G.constrained(cv1,0)._C

    p.clf()
    p.subplot(311)
    a = Cactual.min()
    b = Cactual.max()


    p.imshow(Cactual,vmin=a,vmax=b)
    p.colorbar()
    p.subplot(312)
    p.imshow(Citer,vmin=a,vmax=b)
    p.colorbar()
    p.subplot(313)
    p.imshow((Citer-Cactual)/b)
    p.colorbar()

def demo5(width1=50,width2=100,plaw=-1.0,niter=4,size=500):
    """Comparison of exact vs explicit iterative covariance matrix"""

    import pylab as p

    if plaw is None:
        C = random_covariance(size)
    else:
        C = powerlaw_covariance(plaw,size)

    G = Gaussian(C)


    cv1 = constraint_vector(width1)
    cv2 = constraint_vector(width2)

    Citer = G.projection_constrained(cv1, 10).projection_constrained(cv2,-10).build_cov()
    Cactual = G.constrained(cv1,10).constrained(cv2,-10)._C

    p.clf()
    p.subplot(311)
    a = Cactual.min()
    b = Cactual.max()


    p.imshow(Cactual,vmin=a,vmax=b)
    p.colorbar()
    p.subplot(312)
    p.imshow(Citer,vmin=a,vmax=b)
    p.colorbar()
    p.subplot(313)
    p.imshow((Citer-Cactual)/b)
    p.colorbar()



def demo_fabio(*heights,**kwargs):
    seed = kwargs.pop('seed',0)
    import pylab as p
    maxwidth = 300

    cv1 = constraint_vector(100)
    cv2 = constraint_vector(300)



    cvs = (cv1,cv2)[:len(heights)]

    G = Gaussian(powerlaw_covariance(-1.0,cutoff=50))


    G1 = G
    for height_i, cv_i in zip(heights,cvs):
        G1 = G1.projection_constrained(cv_i,height_i)

    np.random.seed(seed)
    R = G.realization()
    np.random.seed(seed)
    R1 = G1.realization()
    p.plot(R)
    p.plot(R1)



    d0 = [np.dot(R,cv_i) for cv_i in cvs] # original (unmodified) 'heights'

    M = G.get_chi2_matrix(*cvs)
    print "chi2 actual=",G.chi2(R1)-G.chi2(R)
    print "chi2 expect=",np.dot(heights,np.dot(M,heights)) - np.dot(d0,np.dot(M,d0))

def chi2_fabio(plaw=-1.0, height2_fac=1.0, offset=0):
    import pylab as p

    cv1 = constraint_vector(50)
    cv2 = constraint_vector(5,position=250+offset)

    G = Gaussian(powerlaw_covariance(plaw))

    M = G.get_chi2_matrix(cv1,cv2)

    np.random.seed(1)

    height2 = np.sqrt(np.dot(cv2,G.cov_dot(cv2)))*height2_fac
    height1_scale = np.sqrt(np.dot(cv1,G.cov_dot(cv1)))

    heights1 = np.arange(0.1,5.0,0.01)*height1_scale

    chi2 = [np.dot([height1,height2],np.dot(M,[height1,height2])) for height1 in heights1]
    chi2 = np.array(chi2)
    chi2-=chi2.min()
    p.plot(heights1/height1_scale,chi2)

    p.ylabel(r"$\Delta \chi^2$")

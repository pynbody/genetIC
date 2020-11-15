import numpy as np
import scipy.fftpack
import functools
import math

def unitary_fft(x):
    return scipy.fftpack.rfft(x)/math.sqrt(float(len(x)))

def unitary_inverse_fft(x):
    return scipy.fftpack.irfft(x)*math.sqrt(float(len(x)))

def complex_dot(x,y):
    """Dot product for packed FFT complex coeffs"""
    if len(x)%2==0:
        return np.dot(x,y) + np.dot(x[1:-1],y[1:-1]) # counts +ve and -ve modes; NB Nyquist is only ONE mode
    else:
        return np.dot(x, y) + np.dot(x[1:], y[1:])  # counts +ve and -ve modes

class FFTArray(np.ndarray):
    def __new__(subtype, data, **kwargs):
        X = np.asarray(data,**kwargs).view(subtype)
        X.fourier = False
        return X

    def __array_finalize__(self, obj):
        self.fourier = getattr(obj, 'fourier', False)

    def __init__(self, *args, **kwargs):
        assert hasattr(self, 'fourier')

    def in_fourier_space(self):
        if not self.fourier:
            self[:] = unitary_fft(self)
            self.fourier = True
            return self

        return self

    def in_real_space(self):
        if self.fourier:
            self[:] = unitary_inverse_fft(self)
            self.fourier = False

        return self

    def norm(self):
        if self.fourier:
            return math.sqrt(complex_dot(self,self))
        else:
            return np.linalg.norm(self)

    @property
    def k(self):
        return scipy.fftpack.rfftfreq(len(self))

def unitary_fft_matrix(n):
    ki, xi = np.mgrid[:n,:n]
    return np.exp(-2.j*np.pi*ki*xi/n)/np.sqrt(n)


def _converter(fn, call='in_fourier_space'):
    @functools.wraps(fn)
    def wrapped(*args, **kwargs):
        new_args = []
        new_kwargs = {}
        for a in args:
            if hasattr(a, call):
                new_args.append(getattr(a, call)())
            else:
                new_args.append(a)

        for k, v in kwargs.items():
            if hasattr(v, call):
                new_kwargs[k] = getattr(v, call)()
            else:
                new_kwargs[k] = v

        return fn(*args, **kwargs)

    return wrapped

def in_fourier_space(fn):
    """Function decorator to ensure all FFTArray inputs are given in Fourier space"""
    return _converter(fn, 'in_fourier_space')

def in_real_space(fn):
    """Function decorator to ensure all FFTArray inputs are given in real space"""
    return _converter(fn, 'in_real_space')

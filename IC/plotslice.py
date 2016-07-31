import numpy as np
import pylab as p

def plot1dslice(prefix="output/",ps="-",slice_z=None,slice_y=None,plot_b=True,scale_b=1.0,vmin=-0.15,vmax=0.15):
    a = np.load(prefix+"grid-0.npy")
    if plot_b:
        b = np.load(prefix+"grid-1.npy")

    ax,ay,az,aL = [float(x) for x in open(prefix+"grid-info-0.txt").readline().split()]

    if slice_z is None:
        slice_z = aL/2

    if slice_y is None:
        slice_y = aL/2

    if plot_b:
        bx,by,bz,bL = [float(x) for x in open(prefix+"grid-info-1.txt").readline().split()]

    a_sl_z = int(len(a)*((slice_z-az)/aL))
    a_sl_y = int(len(a)*((slice_y-ay)/aL))

    print "a_sl_z,a_sl_y=",a_sl_z,a_sl_y

    dx = aL/len(a)

    a_vals = np.linspace(ax+dx/2,ax+aL-dx/2,len(a))



    if plot_b:
        b_sl_z = int(len(b)*((slice_z-bz)/bL))
        b_sl_y = int(len(b)*((slice_y-by)/bL))
        print "b_sl_z,b_sl_y=",b_sl_z,b_sl_y
        dx = bL/len(b)
        b_vals = np.linspace(bx+dx/2,bx+bL-dx/2,len(b))
        a[(a_vals>b_vals.min()) * (a_vals<b_vals.max())] = np.nan
        p.plot(a_vals,a[:,a_sl_y,a_sl_z].real,ps)
        p.plot(b_vals,b[:,b_sl_y,b_sl_z].real*scale_b,ps)



    else :
        p.plot(a_vals,a[:,a_sl_y,a_sl_z].real,ps)

    p.xlim(0,aL)



def plotslice(prefix="output/",slice=None,plot_b=True,vmin=-0.15,vmax=0.15):
    a = np.load(prefix+"grid-0.npy")
    if plot_b:
        b = np.load(prefix+"grid-1.npy")

    ax,ay,az,aL = [float(x) for x in open(prefix+"grid-info-0.txt").readline().split()]
    if slice is None:
        slice = aL/2

    if plot_b:
        bx,by,bz,bL = [float(x) for x in open(prefix+"grid-info-1.txt").readline().split()]

    a_sl = int(len(a)*((slice-az)/aL))

    if vmin is None:
        vmin = a[a_sl].real.min()
    if vmax is None:
        vmax = a[a_sl].real.max()

    p.imshow(a[a_sl].real,extent=(0,aL,aL,0),vmin=vmin,vmax=vmax,interpolation='nearest')

    if plot_b:
        b_sl = int(len(b)*((slice-bz)/bL))
        if b_sl<0 or b_sl>=len(b):
            print "Note that level 1 is out of range"
        else:
            p.imshow(b[b_sl].real,extent=(bx,bx+bL,by+bL,by),vmin=vmin,vmax=vmax,interpolation='nearest')
            p.plot([bx,bx+bL,bx+bL,bx,bx],[by,by,by+bL,by+bL,by],'k:')
    p.xlim(0,aL)
    p.ylim(0,aL)


def plotslice_pynbody(f, slice=0.0,vmin=-0.15,vmax=0.15):
    import pynbody
    f = pynbody.load(f)
    slice/=f.properties['boxsize'].in_units('Mpc a h^-1',**f.conversion_context())
    rho_mean = f.dm['mass'].sum()/f.properties['boxsize']**3 # should be numerically equal to omegaM0
    print "rho_mean=",rho_mean
    f.dm['delta'] = (f.dm['rho']-rho_mean)/rho_mean

    assert abs(f.dm['z'].max() - f.dm['z'].min() - 1.0)<0.03

    f.dm['z']-=slice

    pynbody.plot.sph.image(f.dm,qty='delta',width=1,log=False,vmin=vmin,vmax=vmax,denoise=True)

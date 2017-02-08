import numpy as np
import pylab as p
import glob

def plot1dslice(prefix="output/",ps="-",slice_z=None,slice_y=None,maxgrid=2,vmin=-0.15,vmax=0.15,thisgrid=0):
    a = np.load(prefix+"grid-%d.npy"%thisgrid)

    ax,ay,az,aL = [float(x) for x in open(prefix+"grid-info-%d.txt"%thisgrid).readline().split()]

    if slice_z is None:
        slice_z = aL/2

    if slice_y is None:
        slice_y = aL/2

    a_sl_z = int(len(a)*((slice_z-az)/aL))
    a_sl_y = int(len(a)*((slice_y-ay)/aL))

    print "a_sl_z,a_sl_y=",a_sl_z,a_sl_y

    dx = aL/len(a)

    a_vals = np.linspace(ax+dx/2,ax+aL-dx/2,len(a))



    if thisgrid<maxgrid:
        bx,by,bz,bL = [float(x) for x in open(prefix+"grid-info-%d.txt"%(thisgrid+1)).readline().split()]
        b = np.load(prefix+"grid-%d.npy"%(thisgrid+1))
        dx = bL/len(b)
        b_vals = np.linspace(bx+dx/2,bx+bL-dx/2,len(b))
        a[(a_vals>b_vals.min()) * (a_vals<b_vals.max())] = np.nan
        p.plot(a_vals,a[:,a_sl_y,a_sl_z].real,ps)

        plot1dslice(prefix,ps,slice_z,slice_y,maxgrid,vmin,vmax,thisgrid+1)

    else :
        p.plot(a_vals,a[:,a_sl_y,a_sl_z].real,ps)

    p.xlim(0,aL)


def plotslice_onegrid(prefix="output/",grid=0,slice=None,vmin=-0.15,vmax=0.15,padcells=0,offset=None):
    new_plot = p.gcf().axes == []
    a = np.load(prefix+"grid-%d.npy"%grid)

    ax,ay,az,aL = [float(x) for x in open(prefix+"grid-info-%d.txt"%grid).readline().split()]

    if offset is None:
        offset=(-aL/2,-aL/2)
    ax+=offset[0]
    ay+=offset[1]

    if slice is None:
        slice = az+aL/2

    a_sl = int(len(a)*((slice-az)/aL))

    if a_sl<0 or a_sl>=len(a):
        print "Grid %d is not contained in this z-slice"%grid
        return

    a = a[a_sl]


    if vmin is None:
        vmin = a.real.min()
    if vmax is None:
        vmax = a.real.max()

    dx=aL/len(a)
    aL-=dx*padcells*2
    ax+=dx*padcells
    ay+=dx*padcells

    if padcells>0:
        a = a[padcells:-padcells,padcells:-padcells]

    p.imshow(a.real,extent=(ax,ax+aL,ay+aL,ay),vmin=vmin,vmax=vmax,interpolation='nearest')
    p.plot([ax,ax+aL,ax+aL,ax,ax],[ay,ay,ay+aL,ay+aL,ay],'k:')

    if new_plot:
        p.xlim(ax,ax+aL)
        p.ylim(ay,ay+aL)
    return slice, vmin, vmax, offset


def plotslice(prefix="output/",maxgrid=10,slice=None,plot_b=True,vmin=-0.15,vmax=0.15,padcells=4, offset=None):
    maxgrid_on_disk = len(glob.glob(prefix+"grid-?.npy"))
    print maxgrid_on_disk
    if maxgrid_on_disk<maxgrid:
        maxgrid = maxgrid_on_disk

    for level in range(maxgrid):
        slice, vmin, vmax, offset = plotslice_onegrid(prefix,level,slice,vmin,vmax,padcells=0 if level==0 else padcells, offset=offset)



def plotslice_pynbody(f, slice=0.0,vmin=-0.15,vmax=0.15):
    import pynbody
    f = pynbody.load(f)
    f.physical_units("Mpc a h^-1")
    slice/=f.properties['boxsize'].in_units('Mpc a h^-1',**f.conversion_context())
    rho_mean = f.dm['mass'].sum()/f.properties['boxsize']**3 # should be numerically equal to omegaM0
    print "rho_mean=",rho_mean
    f.dm['delta'] = (f.dm['rho']-rho_mean)/rho_mean
    f.dm['delta'].convert_units("1")

    assert abs(f.dm['z'].max().in_units(f.properties['boxsize'])
               - f.dm['z'].min().in_units(f.properties['boxsize']) - 1.0)<0.03

    f.dm['z']-=slice

    pynbody.plot.sph.image(f.dm,qty='delta',width=f.properties['boxsize'],log=False,vmin=vmin,vmax=vmax,denoise=True,
                           show_cbar=False)

    return f

def plot_ps(f):
    for level in range(3):
        search = f+"/*_%d.ps"%level
        ps_fname = glob.glob(search)
        if len(ps_fname)==0:
            return
        print search,"->",ps_fname
        k, Pk = np.loadtxt(ps_fname[0],unpack=True,usecols=(0,3))
        p.plot(k,Pk)
    p.loglog()

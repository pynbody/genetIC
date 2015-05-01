import numpy as np
import pylab as p

def plotslice(prefix="output/",slice=0.0,plot_b=True,vmin=-0.15,vmax=0.15):
    a = np.load(prefix+"grid-0.npy")
    b = np.load(prefix+"grid-1.npy")

    ax,ay,az,aL = [float(x) for x in open(prefix+"grid-info-0.txt").readline().split()]
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

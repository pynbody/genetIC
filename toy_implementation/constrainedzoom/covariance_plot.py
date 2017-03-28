import functools

import numpy as np
import pylab as p

from .power_spectrum import powerlaw_covariance
from .methods import ZoomConstrained
from .methods.filtered import FilteredZoomConstrained
from .methods.idealized import FastIdealizedZoomConstrained
from . import fft_wrapper


def plot_covariance_slice(G: ZoomConstrained, cov: np.ndarray, element_offset=0):
    """Plots the 1D slice of the covariance matrix specified by element_offset from the middle of the high-res region"""
    hires_element = G.n2//2+element_offset+G.pixel_size_ratio//2
    lores_element = G.offset+hires_element//G.pixel_size_ratio
    C11 = cov[:G.n1, :G.n1]
    C22 = cov[G.n1:, G.n1:]
    C12 = cov[:G.n1, G.n1:]

    x1, x2 = G.xs()

    p.plot(x1,C11[lores_element],"k")
    p.plot(x2,C22[hires_element],"r")


def overplot_boundary_real_space(G: ZoomConstrained, cov: np.ndarray,
                                 zoomin_size=5, plot_size=0.333, pad=0,
                                 vmin=None, vmax=None, show_hh=True, left=0.05, bottom=0.95):
    """Overplots a zoomed representation of the covariance matrix at the edge of the zoom region

    :param G - the ZoomConstrained object that originally generated the covariance matrix
    :param cov - the (n1+n2) x (n1+n2) representation of the covariance matrix
    :param zoomin_size - the number of (coarse) pixels on either side of the edge to show
    :param pad - the number of (coarse) pixels ommitted on the edges of the fine region
    :param plot_size - the size of the panel to be added
    :param bottom - the location of the bottom of the panel to be added
    :param left - the location of the left of the panel to be added

    Typically you want to call this after having first called display_cov
    """

    C11 = cov[:G.n1, :G.n1]
    C22 = cov[G.n1:, G.n1:]
    C12 = cov[:G.n1, G.n1:]

    offset_end = G.offset + G.n1//G.window_size_ratio - pad
    slice_1 = slice(offset_end-zoomin_size,offset_end+zoomin_size)
    if pad>0:
        slice_2 = slice(-zoomin_size*G.pixel_size_ratio-pad*G.pixel_size_ratio,-pad*G.pixel_size_ratio)
    else:
        slice_2 = slice(-zoomin_size * G.pixel_size_ratio, None)

    xs, _ = G.xs()
    source_start = xs[slice_1.start]
    source_end = xs[slice_1.stop-1]


    p.imshow(C11[slice_1,slice_1], extent=[left, left + plot_size,
                                           bottom , bottom - plot_size], vmin=vmin, vmax=vmax)
    p.imshow(C12[slice_1,slice_2].T, extent=[left, left + plot_size,
                                             bottom-plot_size/2 , bottom- plot_size ], vmin=vmin, vmax=vmax)
    p.imshow(C22[slice_2,slice_2], extent=[left, left + plot_size / 2,
                                           bottom-plot_size/2 , bottom - plot_size ], vmin=vmin, vmax=vmax)


    p.plot(left + np.array([0, plot_size, plot_size, 0, 0]),
           bottom - np.array([plot_size, plot_size, 0, 0, plot_size]), "k")

    p.plot([left,left+plot_size],[bottom-plot_size/2,bottom-plot_size/2],"k:")
    p.plot([left+plot_size/2,left+plot_size/2],[bottom-plot_size/2,bottom-plot_size],'k:')


    p.plot([source_start,source_end,source_end,source_start,source_start],
           [source_start,source_start,source_end,source_end,source_start],
           "k")

    p.plot([source_start,left],[source_end,bottom-plot_size],"k-")
    p.plot([source_end,left+plot_size],[source_end,bottom-plot_size],'k-')

    p.xlim(0,1.0)
    p.ylim(1.0,0)


def plot_power_spec(G: ZoomConstrained, cov:np.ndarray, pad=0):
    C11 = cov[:G.n1, :G.n1]
    C12 = cov[:G.n1, G.n1:]
    C22 = cov[G.n1:, G.n1:]

    if pad > 0:
        C22 = C22[pad:-pad, pad:-pad]
        C12 = C12[:, pad:-pad]

    k_nyq_1 = G.n1/2  # N pi / L, but L=1
    k_nyq_2 = G.n2 * G.window_size_ratio/2

    U1 = fft_wrapper.unitary_fft_matrix(len(C11))
    U2 = fft_wrapper.unitary_fft_matrix(len(C22))

    # keep only positive frequency components
    U1 = U1[:len(U1) // 2]
    U2 = U2[:len(U2) // 2]

    C11 = U1 @ C11 @ U1.conj().T
    C22 = U2 @ C22 @ U2.conj().T
    C12 = U1 @ C12 @ U2.conj().T

    k1 = np.linspace(0,k_nyq_1-1,len(C11))
    k2 = np.linspace(0,k_nyq_2-1,len(C22))

    p.plot(k1, C11.diagonal()/k1[1]/2,"r")
    p.plot(k2, C22.diagonal()/k2[1]/2,"g")

    p.plot(G.k_low, G.C_low/G.k_low[1], "m:")
    p.plot(G.k_high, G.C_high/G.k_high[1], "y:")
    p.loglog()


def plot_fourier_space(G: ZoomConstrained, cov:np.ndarray, pad=0, vmin=None, vmax=None,
                       hilo_axis=None, hihi_axis=None) :
    """Display the covariance matrix in fourier space

    :param G: the ZoomConstrained object that originally generated the covariance matrix
    :param cov: the covariance matrix in real space
    :param pad: the number of (coarse) pixels on either edge of the zoom region to ignore in real space (before FTing)
    :param vmin: passed to pyplot's imshow
    :param vmax: passed to pyplot's imshow
    :param hilo_axis: the axes to use for the fine x coarse covariance matrix, or None to create fresh ones
    :param hihi_axis: the axes to use for the fine x fine covariance matrix (if hilo_axis is not None)
    """
    C11 = cov[:G.n1, :G.n1]
    C12 = cov[:G.n1, G.n1:]
    C22 = cov[G.n1:, G.n1:]

    if pad>0:
        pad2 = pad*G.pixel_size_ratio
        C22  = C22[pad2:-pad2,pad2:-pad2]
        C12 = C12[:,pad2:-pad2]

    k_nyq_1 = G.n1 # N pi / L, but L=1
    k_nyq_2 = G.n2*G.window_size_ratio

    U1 = fft_wrapper.unitary_fft_matrix(len(C11))
    U2 = fft_wrapper.unitary_fft_matrix(len(C22))

    # keep only positive frequency components
    U1 = U1[:len(U1)/2]
    U2 = U2[:len(U2)/2]

    C11 = U1@C11@U1.conj().T
    C22 = U2@C22@U2.conj().T
    C12 = U1@C12@U2.conj().T

    k_pixel_size_1 = k_nyq_1 / len(C11)
    k_pixel_size_2 = k_nyq_2 / len(C22)

    if not hilo_axis:
        p.clf()
        fig, (hilo_axis, hihi_axis) = p.subplots(1, 2, num = p.gcf().number,
                                                 gridspec_kw={'width_ratios': [k_nyq_1, k_nyq_2]},
                                                 squeeze=True, sharey=True)

    p.sca(hihi_axis)

    p.imshow(C22.real,extent=[0,k_nyq_2+k_pixel_size_2,k_nyq_2+k_pixel_size_2,0],interpolation='none',vmin=vmin,vmax=vmax)
    p.plot([k_nyq_1,k_nyq_1,0],[0,k_nyq_1,k_nyq_1],'k:')

    p.xlabel("Wavenumber/$\pi$")

    p.xlim(0, k_nyq_2)
    p.ylim(k_nyq_2, 0)

    p.text(k_nyq_2 - k_nyq_2 * 0.02, k_nyq_2 * 0.02, 'Fine', horizontalalignment='right',
           verticalalignment='top', color='black', fontsize=15)

    p.sca(hilo_axis)
    text_offset = k_nyq_2*0.02

    p.imshow(abs(C12).T, extent=[0,k_nyq_1+k_pixel_size_1,k_nyq_2+k_pixel_size_2,0], interpolation='none', vmin=vmin, vmax=vmax)
    #p.imshow(C11.real, extent=[0,k_nyq_1+k_pixel_size_1,k_nyq_1+k_pixel_size_1,0], interpolation='none', vmin=vmin,vmax=vmax)

    p.text(k_nyq_1-text_offset , k_nyq_1+text_offset, 'Fine $\\times$ \n Coarse', horizontalalignment='right',
           verticalalignment='top', color='black', fontsize=15)
    p.text(k_nyq_1-text_offset, text_offset, 'Coarse', horizontalalignment='right',
           verticalalignment='top', color='black', fontsize=15)

    p.plot([0, k_nyq_1], [k_nyq_1, k_nyq_1], 'k:')


    p.xlim(0,k_nyq_1)
    p.ylim(k_nyq_2,0 )

    p.xlabel("Wavenumber/$\pi$")
    p.ylabel("Wavenumber/$\pi$")


def plot_real_space(G: ZoomConstrained, cov: np.ndarray,
                    downgrade=False, vmin=None, vmax=None, pad=0, show_hh=True):
    """Display the covariance matrix in real space

    :param G: the ZoomConstrained object that originally generated the covariance matrix
    :param cov: the covariance matrix
    :param pad: the number of (coarse) pixels on either edge of the zoom region to ignore
    :param vmin: passed to pyplot's imshow
    :param vmax: passed to pyplot's imshow
    :param show_hh: whether to show the fine x fine component (otherwise show only coarse and fine x coarse)
    :param downgrade: if True, downgrade the fine region to the same resolution as the coarse region
    :return:
    """

    p.set_cmap('PuOr_r')
    vmin = vmin if vmin is not None else np.min(cov)
    vmax = vmax if vmax is not None else np.max(cov)

    C11 = cov[:G.n1,:G.n1]
    C22 = cov[G.n1:,G.n1:]
    C12 = cov[:G.n1,G.n1:]

    zoom_width = G.n1//G.window_size_ratio
    offset = G.offset
    n1 = G.n1

    if pad>0:
        pixel_scale = (G.n2*G.window_size_ratio)//G.n1
        pad_fine = pad*pixel_scale
        C22 = C22[pad_fine:-pad_fine,pad_fine:-pad_fine]
        C12 = C12[:,pad_fine:-pad_fine]
        zoom_width-=pad*2
        offset+=pad

    p.imshow(C11,extent=[0,1.0,1.0,0],vmin=vmin,vmax=vmax,interpolation='none')
    if downgrade:
        zoom_fac = G.window_size_ratio*(G.n2//G.n1)
        print("zoom_fac=",zoom_fac)
        C22new=0
        for i in range(zoom_fac):
            for j in range(zoom_fac):
                C22new += C22[i::zoom_fac,j::zoom_fac]
        C22 = C22new/zoom_fac**2

        C12new=0
        for i in range(zoom_fac):
            C12new+=C12[:,i::zoom_fac]
        C12 = C12new/zoom_fac

    zoom_window_start_coordinate = offset/n1
    zoom_window_end_coordinate = (offset+zoom_width)/n1

    p.imshow(C12.T,extent=[0,1.0,zoom_window_end_coordinate,zoom_window_start_coordinate],vmin=vmin,vmax=vmax,interpolation='none')
    if show_hh:
        p.imshow(C22,extent=[zoom_window_start_coordinate,zoom_window_end_coordinate,
                             zoom_window_end_coordinate,zoom_window_start_coordinate],
                 vmin=vmin,vmax=vmax,interpolation='none')

    p.plot([0,1.0],[zoom_window_start_coordinate]*2,'k:')
    p.plot([0,1.0],[zoom_window_end_coordinate]*2,'k:')
    if show_hh:
        p.plot([zoom_window_start_coordinate]*2,
               [zoom_window_start_coordinate,zoom_window_end_coordinate],'k:')
        p.plot([zoom_window_end_coordinate]*2,
               [zoom_window_start_coordinate,zoom_window_end_coordinate],'k:')

    if show_hh:
        p.text(zoom_window_end_coordinate-0.01,zoom_window_start_coordinate+0.01,
               'Fine',horizontalalignment='right',verticalalignment='top',color='black',fontsize=15)

    p.text(0.99,zoom_window_start_coordinate+0.01,
           r'Fine $\times$ Coarse',horizontalalignment='right',verticalalignment='top',color='black',fontsize=15)
    p.text(0.99,zoom_window_end_coordinate+0.01,
           'Coarse',horizontalalignment='right',verticalalignment='top',color='black',fontsize=15)

    p.xlim(0,1.0)
    p.ylim(1.0,0)
    p.xlabel("Position")
    #p.ylabel("Position")




def combined_plots(G: ZoomConstrained, cov:np.ndarray, pad=0, vmin=None, vmax=None):
    """Create a three-panel plot of the provided covariance in harmonic and real space

    :param G: the ZoomConstrained object that originally generated the covariance matrix
    :param cov: the covariance matrix
    :param pad: the number of (coarse) pixels on either edge of the zoom region to ignore
    :param vmin: passed to pyplot's imshow
    :param vmax: passed to pyplot's imshow
    """
    k_nyq_1 = G.n1  # N pi / L, but L=1
    k_nyq_2 = G.n2 * G.window_size_ratio

    # ideal - but does not seem to work on mac os:
    p.gcf().set_size_inches(11.275, 4.925, forward=True)

    gs_fourier = p.GridSpec(2, 1, hspace=0)


    fig, (hilo_axis, hihi_axis, real_axis, colbar_axis) = p.subplots(1, 4, num=p.gcf().number,
                                             gridspec_kw={'width_ratios': [k_nyq_1, k_nyq_2, k_nyq_2, k_nyq_2/15]},
                                             squeeze=True, sharey=False)

    #real_axis.yaxis.set_label_position("right")
    #real_axis.yaxis.set_ticks_position("right")
    hihi_axis.yaxis.set_ticklabels([])

    pos = hihi_axis.get_position()
    pos.x0-=0.02
    pos.x1-=0.02
    hihi_axis.set_position(pos)

    plot_fourier_space(G, cov, pad, vmin, vmax, hilo_axis, hihi_axis)
    p.sca(hihi_axis)
    p.title("Fourier space")
    p.sca(real_axis)
    plot_real_space(G, cov, vmin=vmin, vmax=vmax, pad=pad)
    overplot_boundary_real_space(G, cov, pad=pad, vmin=vmin, vmax=vmax)
    p.title("Real space")
    p.xticks([0,0.5,1.0])
    p.yticks([0,0.5,1.0])
    cbar = p.colorbar(cax=colbar_axis)
    cbar.set_label("Fractional error")
    cbar.set_ticks([vmin,0,vmax])
    cbar.set_ticklabels(["$%.0f$%%" % (vmin*100),
                          "0",
                          "$%.0f$%%" % (vmax * 100)])


def cov_zoom_demo(n1=256, n2=256,
                  hires_window_scale=4,
                  estimate=False, Ntrials=2000,
                  plaw=-1.5, cl=FilteredZoomConstrained,pad=0,vmin=None,vmax=None,errors=False,
                  show_hh=True,one_element=None,subplot=False,
                  display_fourier=False,
                  initialization_kwargs={}):
    if not subplot:
        p.clf()

    cov_this = functools.partial(powerlaw_covariance, plaw=plaw)
    X = cl(cov_this, n1=n1, n2=n2, hires_window_scale=hires_window_scale, **initialization_kwargs)
    if estimate:
        cv_est = X.estimate_cov(Ntrials)
    else:
        cv_est = X.get_cov(one_element=one_element)

    if errors:
        Y = FastIdealizedZoomConstrained(cov_this, n1=n1, n2=n2, hires_window_scale=hires_window_scale)
        true_cov = Y.get_cov(one_element=one_element)
        cv_est-=true_cov
        cv_est/=true_cov.max() # variance in hi-res region
        if vmin is None:
            vmin = -max(abs(cv_est[1, 1]), abs(cv_est[-n2 / 2, -n2 / 2]))
        if vmax is None:
            vmax = max(abs(cv_est[1, 1]), abs(cv_est[-n2 / 2, -n2 / 2]))

    else:
        cv_est/=cv_est.max()
        if vmin is None:
            vmin = -1.0
        if vmax is None:
            vmax = 1.0



    if display_fourier:
        combined_plots(X, cv_est, pad, vmin=vmin, vmax=vmax)
    else:
        if one_element is None:
            plot_real_space(X, cv_est, pad=pad, vmin=vmin, vmax=vmax, show_hh=show_hh)
            overplot_boundary_real_space(X, cv_est, pad=pad, vmin=vmin, vmax=vmax, show_hh=show_hh)
            if not subplot:
                cbar = p.colorbar()
                if errors:
                    cbar.set_label("Fractional error in covariance")
                else:
                    cbar.set_label("Covariance (arbitrary units)")
            p.xlabel("Real-space pixel 1")
            p.ylabel("Real-space pixel 2")
            p.subplots_adjust(0.04, 0.14, 0.93, 0.95)
        else:
            plot_covariance_slice(X, cv_est, one_element)

    return X
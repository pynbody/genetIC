import functools

import numpy as np
import pylab as p

from .power_spectrum import powerlaw_covariance
from .methods import ZoomConstrained
from .methods.filtered import FilteredZoomConstrained
from .methods.idealized import FastIdealizedZoomConstrained
from . import fft_wrapper
from . import methods


def plot_covariance_slice(G: ZoomConstrained, cov: np.ndarray, element_offset=0):
    """Plots the 1D slice of the covariance matrix specified by element_offset from the middle of the high-res region"""
    hires_element = G.nW//2+element_offset+G.pixel_size_ratio//2
    lores_element = G.offset+hires_element//G.pixel_size_ratio
    C11 = cov[:G.nP, :G.nP]
    C22 = cov[G.nP:, G.nP:]
    C12 = cov[:G.nP, G.nP:]

    x1, x2 = G.xs()

    p.plot(x1,C11[lores_element],"k")
    p.plot(x2,C22[hires_element],"r")


def overplot_boundary_real_space(G: ZoomConstrained, cov: np.ndarray,
                                 zoomin_size=5, plot_size=0.333, pad=0,
                                 vmin=None, vmax=None, show_hh=True, left=0.05, bottom=0.95):
    """Overplots a zoomed representation of the covariance matrix at the edge of the zoom region

    :param G - the ZoomConstrained object that originally generated the covariance matrix
    :param cov - the (nP+nW) x (nP+nW) representation of the covariance matrix
    :param zoomin_size - the number of (coarse) pixels on either side of the edge to show
    :param pad - the number of (coarse) pixels ommitted on the edges of the fine region
    :param plot_size - the size of the panel to be added
    :param bottom - the location of the bottom of the panel to be added
    :param left - the location of the left of the panel to be added

    Typically you want to call this after having first called display_cov
    """

    C11 = cov[:G.nP, :G.nP]
    C22 = cov[G.nP:, G.nP:]
    C12 = cov[:G.nP, G.nP:]

    offset_end = G.offset + G.nP//G.window_size_ratio - pad
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


def plot_power_spec(G: ZoomConstrained, cov:np.ndarray, pad=0, errors=False):
    C11 = cov[:G.nP, :G.nP]
    C12 = cov[:G.nP, G.nP:]
    C22 = cov[G.nP:, G.nP:]

    if pad > 0:
        pad_highres = pad*G.pixel_size_ratio
        C22 = C22[pad_highres:-pad_highres, pad_highres:-pad_highres]
        C12 = C12[:, pad_highres:-pad_highres]

    k_nyq_1 = G.nP/2  # N pi / L, but L=1
    k_nyq_2 = G.nW * G.window_size_ratio/2

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

    actual_P1 = C11.diagonal()/k1[1]
    actual_P2 = C22.diagonal()/k2[1]
    target_P1 = 2.*G._get_variance_k(k1)*G.nP
    target_P2 = 2.*G._get_variance_k(k2)*G.nW/G.window_size_ratio

    if errors:
        p.plot(k1, actual_P1, "r")
        p.plot(k2, actual_P2, "g")
        p.legend(["Coarse", "Fine"], loc='lower left')
        p.semilogx()
    else:
        p.plot(k1, actual_P1,"r")
        p.plot(k2, actual_P2,"g")

        p.plot(G.k_low, G.C_low/G.k_low[1], "m:")
        p.plot(G.k_high, G.C_high/G.k_high[1], "y:")

        p.legend(["Coarse Actual", "Fine Actual", "Coarse Target", "Fine Target"], loc='lower left')
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
    C11 = cov[:G.nP, :G.nP]
    C12 = cov[:G.nP, G.nP:]
    C22 = cov[G.nP:, G.nP:]

    if pad>0:
        pad2 = pad*G.pixel_size_ratio
        C22  = C22[pad2:-pad2,pad2:-pad2]
        C12 = C12[:,pad2:-pad2]

    k_nyq_1 = G.nP # N pi / L, but L=1
    k_nyq_2 = G.nW*G.window_size_ratio

    U1 = fft_wrapper.unitary_fft_matrix(len(C11))
    U2 = fft_wrapper.unitary_fft_matrix(len(C22))

    # keep only positive frequency components
    U1 = U1[:len(U1)//2]
    U2 = U2[:len(U2)//2]

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
    p.imshow(C11.real, extent=[0,k_nyq_1+k_pixel_size_1,k_nyq_1+k_pixel_size_1,0], interpolation='none', vmin=vmin,vmax=vmax)

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

    C11 = cov[:G.nP,:G.nP]
    C22 = cov[G.nP:,G.nP:]
    C12 = cov[:G.nP,G.nP:]

    zoom_width = G.nP//G.window_size_ratio
    offset = G.offset
    nP = G.nP

    if pad>0:
        pixel_scale = (G.nW*G.window_size_ratio)//G.nP
        pad_fine = pad*pixel_scale
        C22 = C22[pad_fine:-pad_fine,pad_fine:-pad_fine]
        C12 = C12[:,pad_fine:-pad_fine]
        zoom_width-=pad*2
        offset+=pad

    p.imshow(C11,extent=[0,1.0,1.0,0],vmin=vmin,vmax=vmax,interpolation='none')
    if downgrade:
        zoom_fac = G.window_size_ratio*(G.nW//G.nP)
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

    zoom_window_start_coordinate = offset/nP
    zoom_window_end_coordinate = (offset+zoom_width)/nP

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
    k_nyq_1 = G.nP  # N pi / L, but L=1
    k_nyq_2 = G.nW * G.window_size_ratio

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

def overlay_grid(nP=256, nW=256, hires_window_scale=4, x_min=0.15, x_max=0.20):
    X = FastIdealizedZoomConstrained(lambda x:1, nP, nW, hires_window_scale)
    xP, xW = X.boundary_xs()
    xW = xW[(xW>x_min) & (xW<x_max)]
    xP = xP[(xP>x_min) & (xP<x_max)]
    yrange = p.ylim()
    p.vlines(xW,yrange[0],yrange[1],color='#eeeeee')
    p.vlines(xP,yrange[0],yrange[1],color='#cccccc')

def zoom_demo(nP=256, nW=256, hires_window_scale=4, hires_window_offset=8,plaw=-1.5,method=FilteredZoomConstrained,
              constraint_val=None, constraint_covec=None,
              no_random=False,pad=None,
              show_covec=False,errors=False,linewidth=None,
              constrain_potential=False,seed=None,
              splice_seed=None,diff_seed=None,
              splice_range=(100,200),plot_kwargs={}):

    if errors and not no_random:
        raise ValueError("To display errors in the convolved constraint, you must disable the random component")

    cov_this = functools.partial(powerlaw_covariance, plaw=plaw)

    X = method(cov_this, nP=nP, nW=nW, hires_window_scale=hires_window_scale, offset=hires_window_offset)
    if constraint_val is not None:
        X.add_constraint(constraint_val, constraint_covec, constrain_potential)

    if errors:
        Y = FastIdealizedZoomConstrained(cov_this, nP=nP, nW=nW, hires_window_scale=hires_window_scale, offset=hires_window_offset)
        if constraint_val is not None:
            Y.add_constraint(constraint_val, constraint_covec, constrain_potential)

    if pad is None:
        pad = X.get_default_plot_padding()

    if seed:
        np.random.seed(seed)

    if splice_seed:
        delta_P, delta_W = X.get_spliced_realization(seed, splice_seed, splice_range)
    else:
        delta_P, delta_W = X.realization(no_random=no_random,seed=seed)

    if diff_seed:
        delta_P_diff, delta_W_diff = X.realization(seed=diff_seed)
        delta_P-=delta_P_diff
        delta_W-=delta_W_diff


    if errors:
        delta_Ps, delta_Ws = Y.realization(no_random=no_random)
        delta_P-=delta_Ps
        delta_W-=delta_Ws
        if not (isinstance(errors, str) and "abs" in errors):
            delta_P/=abs(delta_Ws).max()
            delta_W/=abs(delta_Ws).max()


    if constraint_val is not None and (not constrain_potential):
        cov = FastIdealizedZoomConstrained(cov_this, nP, nW, hires_window_scale=hires_window_scale, offset=hires_window_offset).get_cov()[nP:,nP:]
        if constraint_val==0.0:
            cv_var_uncon = np.dot(constraint_covec,np.dot(cov,constraint_covec))
            cv_var_con = np.dot(constraint_covec, np.dot(X.get_cov()[nP:,nP:] , constraint_covec))+1e-10 #1e-10 noise level prevent NaN
            print("%s constraint value %.2f (target %.2f; rms %.2f; rms without constraint %.2f; suppression %.1f%%)"%(X.description,
                                                        np.dot(constraint_covec, delta_W),
                                                        constraint_val,
                                                        np.sqrt(cv_var_con),
                                                        np.sqrt(cv_var_uncon),
                                                        100.*(1.-np.sqrt(cv_var_con/cv_var_uncon))))
        else:
            print("%s constraint value %.2f (target %.2f)"%(X.description,
                                                        np.dot(constraint_covec, delta_W),
                                                        constraint_val))

    xP, xW = X.xs()

    if pad>0:
        hr_pad = pad*X.pixel_size_ratio
        xW = xW[hr_pad:-hr_pad]
        delta_W = delta_W[hr_pad:-hr_pad]

    line, = p.plot(xW, delta_W, label=X.description,linewidth=linewidth,**plot_kwargs)

    left_of_window = slice(0,hires_window_offset+pad+1)
    right_of_window = slice(hires_window_offset+nP//hires_window_scale-pad-1, None)


    for sl in left_of_window, right_of_window:
        plot_kwargs['color'] = line.get_color()
        p.plot(xP[sl], delta_P[sl], ":",linewidth=linewidth,**plot_kwargs)


    if show_covec:
        p.legend()
        p.twinx()
        p.plot(X.xs()[1], constraint_covec, color='#aaaaaa')

    if splice_seed:
        xP, xW = X.xs()
        y_range = max(abs(delta_W).max(),abs(delta_P).max())*2
        yl = p.ylim()
        p.fill_betweenx([-y_range, y_range], xW[splice_range[0]], xW[splice_range[1]-1], color='#eeeeee')
        p.ylim(*yl)

    return X

def compare_constraints(plaw=-1.0, velocity=False, errors=False, covector_width=50,
                        nP=256, nW=256, hires_window_offset=8, pad=None,
                        using_methods=[methods.traditional.TraditionalZoomConstrained, methods.filtered.FilteredZoomConstrained]):
    from . import constraint_vector, deriv_constraint_vector
    if velocity:
        covec = deriv_constraint_vector(covector_width, 256, 127)
        potential = True
    else:
        covec = constraint_vector(covector_width,256,127)
        potential = False

    if p.gca() is p.gcf().axes[0]:
        if velocity:
            p.title("Constrain velocity; $n=%.1f$"%plaw)
        else:
            p.title("Constrain density; $n=%.1f$"%plaw)

    val = 1.0
    random = False
    np.random.seed(5)

    for i,m in enumerate(using_methods):
        zoom_demo(no_random=not random, constraint_val=val, constraint_covec=covec,
                     method=m,
                     plaw=plaw, errors=errors, linewidth=1+len(using_methods)-i, constrain_potential=potential,
                     nP=nP, nW=nW, hires_window_offset=hires_window_offset,pad=pad)


    if not errors:
        zoom_demo(no_random=not random, constraint_val=val, constraint_covec=covec,
                     method=methods.idealized.IdealizedZoomConstrained,
                     plaw=plaw, errors=False, linewidth=1,
                     constrain_potential=potential,
                     nP=nP, nW=nW, hires_window_offset=hires_window_offset,pad=pad)
    p.legend()
    ax = p.gca()
    p.xlim(0.08, 0.38)
    if errors:
        p.ylim(-0.05, 0.05)
        p.ylabel("Fractional error")
        ax.yaxis.set_ticks([-0.04,-0.02,0.0,0.02,0.04])
        ax.yaxis.set_ticklabels(["$-4$%","$-2$%","0","2%","4%"])

    else:
        p.ylabel("Solution")
        ax.yaxis.set_ticks([0])
        ax.yaxis.set_ticklabels(["0"])

    p.xlabel("Position")

def compare_constraints_with_errors(*args, **kwargs):
    p.clf()
    f, (sub1, sub2) = p.subplots(2,1,sharex=True,squeeze=True,gridspec_kw={'hspace': 0},num=p.gcf().number)
    p.sca(sub1)
    compare_constraints(*args, errors=False, **kwargs)
    p.sca(sub2)
    compare_constraints(*args, errors=True, **kwargs)


def cov_zoom_demo(nP=256, nW=256,
                  hires_window_scale=4,
                  hires_window_offset=10,
                  estimate=False, Ntrials=2000,
                  plaw=-1.5, method=FilteredZoomConstrained,
                  pad=None,vmin=None,vmax=None,errors=False,
                  show_hh=True,one_element=None,subplot=False,
                  plot_type='real',
                  with_constraint=False,
                  inlay_plaws=[],
                  initialization_kwargs={}):
    """End-to-end construction and plotting function for understanding a zoom algorithm
    
    :param nP: The number of low-res pixels
    :param nW: The number of high-res pixels
    :param hires_window_scale: The size of the high-res window relative to the low-res box
    :param estimate: if True, use a covariance estimator instead of the exact algorithm (mainly useful for testing the exact algorithm is working)
    :param Ntrials: if estimate is True, the number of trials to use
    :param plaw: the power-law to be used for the underlying covariance matrix P(k) propto k^plaw
    :param method: the class to be used for constructing the zoom gaussian realisation
    :param pad: the real-space padding in low-resolution pixels to be used in the analysis, or None to use the method-recommended value
    :param vmin: passed to imshow routines
    :param vmax: passed to imshow routines
    :param errors: if True, show the errors as a fraction of the real-space variance instead of the actual output.
    :param show_hh: if True (default) show the high x high covariance (otherwise shows only low x low and high x low)
    :param one_element: if not None, show a test of the basic convolution algorithm by using the specified pixel only
                        in the white-noise construction. See method ZoomConstrained._iter_one_cov_element
    :param subplot: if True, place the plot into existing axes
    :param plot_type: one of 'real', 'combined' or 'pspec'. "Real" plots the real-space covariance matrix;
     "combined" plots both real-space and fourier-space covariance matrices. "pspec" plots the power spectrum.
    :param with_constraint: if True, apply a linear constraint/modification to the centre of the zoom window
    :param inlay_plaws: list of plaw values to produce an inlayed mini-diagram for
    :param initialization_kwargs: kwargs to pass to the target class initialiser
    :return: the class instance used to construct the plot
    """
    if not subplot:
        p.clf()

    cov_this = functools.partial(powerlaw_covariance, plaw=plaw)
    X = method(cov_this, nP=nP, nW=nW, hires_window_scale=hires_window_scale, offset=hires_window_offset, **initialization_kwargs)

    if pad is None:
        pad = X.get_default_plot_padding()

    if with_constraint:
        X.add_constraint(0.0)

    if estimate:
        cv_est = X.estimate_cov(Ntrials)
    else:
        cv_est = X.get_cov(one_element=one_element)

    if errors:
        Y = FastIdealizedZoomConstrained(cov_this, nP=nP, nW=nW, hires_window_scale=hires_window_scale, offset=hires_window_offset)
        if with_constraint:
            Y.add_constraint(0.0)

        true_cov = Y.get_cov(one_element=one_element)
        cv_est-=true_cov
        cv_est/=true_cov.max() # variance in hi-res region
        if vmin is None:
            vmin = -max(abs(cv_est[1, 1]), abs(cv_est[-nW // 2, -nW // 2]))
        if vmax is None:
            vmax = max(abs(cv_est[1, 1]), abs(cv_est[-nW // 2, -nW // 2]))

    else:
        cv_est/=cv_est.max()
        if vmin is None:
            vmin = -1.0
        if vmax is None:
            vmax = 1.0



    if plot_type=='combined':
        combined_plots(X, cv_est, pad, vmin=vmin, vmax=vmax)
    elif plot_type=='pspec':
        plot_power_spec(X, cv_est, pad, errors=errors)
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
            p.title(X.description+"; $n=%.1f$"%plaw)
        else:
            plot_covariance_slice(X, cv_est, one_element)

    return X


def compare_covariances(using_methods = [methods.ml.MLZoomConstrained,
                                         methods.traditional.TraditionalZoomConstrained,
                                         methods.filtered.FilteredZoomConstrained],
                        plaws=[-1.5],
                        **kwargs):
    """Compare covariances for MLZoomConstrained, TraditionalZoomConstrained and FilteredZoomConstrained"""

    p.clf()
    fig, axes = p.subplots(1, len(using_methods) + 1, num=p.gcf().number, squeeze=True,
                           gridspec_kw={'width_ratios': ([1.0] * len(using_methods)) + [0.1]})

    vmin = kwargs.setdefault('vmin',-0.05)
    vmax = kwargs.setdefault('vmax',0.05)
    errors = kwargs.setdefault('errors', True)

    kwargs['subplot'] = True
    firstplot = True
    for cl, ax in zip(using_methods, axes):
        p.sca(ax)
        kwargs['method']=cl
        instance = cov_zoom_demo(**kwargs)
        ax.yaxis.set_ticks([0.0,1.0])
        if firstplot:
            ax.yaxis.set_ticklabels(["0.0","1.0"])
            ax.yaxis.set_label_text("Pixel 2 location")
            firstplot = False
        else:
            ax.yaxis.set_label_text("")
            ax.yaxis.set_ticklabels([])
        ax.xaxis.set_label_text("Pixel 1 location")
        ax.xaxis.set_ticks([0,1])
        ax.xaxis.set_ticklabels(["0.0","1.0"])

    cbar = p.colorbar(cax=axes[-1])

    p.subplots_adjust(0.04, 0.14, 0.93, 0.95)

    if errors:
        cbar.set_label("Fractional error")
        cbar.set_ticks([vmin, 0, vmax])
        num_digits = str(int(-np.log10(vmax))-1)
        format = "$%."+num_digits+"f$%%"
        cbar.set_ticklabels([format % (vmin * 100),
                             "0",
                             format % (vmax * 100)])
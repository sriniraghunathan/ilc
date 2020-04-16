import numpy as np, sys, os, scipy as sc

h, k_B, c=6.626e-34,1.38e-23, 3e8
################################################################################################################

def get_foreground_power_spt(component, freq1=150, freq2=None, units='uk', lmax = None):
    """
    Foreground powers from George et al. 2015 results.

    Uses .sav file generated by Christain Reichardt.

    Parameters
    ----------
    component : str
        The foreground component to use. Must be one of
        'all', 'tSZ', 'kSZ', 'DG-Cl', 'DG-Po', 'RG', 'tSZ-CIB', 'Total', 'CMB'
    freq1 : int
        Frequency band. If `freq2` is specified, the cross-spectrum between
        the two frequencies will be returned. Otherwise autospectrum of freq1.
    freq2 : int, optional
        Frequency band for cross-spectrum with `freq1`
    units : str
        'k' or 'uk'. Note: default savfile is Dls in uK

    Returns
    -------
    fgnd_cls : array
        Power spectrum of `component` at specified frequency band.
    """
    components = [
        'all',
        'tSZ',
        'kSZ',
        'DG-Cl',
        'DG-Po',
        'RG',
        'tSZ-CIB',
        'Total',
        'CMB',
    ]
    if component not in components:
        raise ValueError(
            '{} not in list of possible foregrounds, must be one of {}'.format(
                component, components
            )
        )

    #filename = os.path.join(
    #    os.path.dirname(__file__), 'data/foregrounds/george_plot_bestfit_line.sav'
    #)
    filename = 'george_plot_bestfit_line.sav'
    from scipy.io import readsav

    data = readsav(filename)

    #from IPython import embed; embed()
    if freq2 is None:
        freq2 = freq1
    if freq1 == 90:
        freq1 = 95
    if freq2 == 90:
        freq2 = 95

    freqs = np.asarray(
        [(95, 95), (95, 150), (95, 220), (150, 150), (150, 220), (220, 220)]
    )
    dls_all = data['ml_dls'][(freqs[:, 0] == freq1) & (freqs[:, 1] == freq2)][0]
    labels = data['ml_dl_labels'].astype('str')
    el = np.asarray(data['ml_l'], dtype=int)

    if component == 'all':
        spec = el * 0.0
        for fg in components:
            if fg in ['all', 'tSZ-CIB', 'Total', 'CMB']:
                continue
            spec += dls_all[labels == fg][0]
    else:
        spec = dls_all[labels == component][0]

    # Changing Dls to Cls
    spec /= el * (el + 1.0) / 2.0 / np.pi
    if units.lower() == 'k':
        spec /= 1e12

    # Pad to l=0
    spec = np.concatenate((np.zeros(min(el)), spec))
    el = np.concatenate((np.arange(min(el)), el))

    if lmax is not None:
        el = el[:lmax]
        spec = spec[:lmax]

    return el, spec

################################################################################################################

def fn_dB_dT(nu, nu0 = None, temp = 2.725):
    nu *= 1e9
    x=h*nu/(k_B*temp)
    dBdT = x**4. * np.exp(x) / (np.exp(x)-1)**2.

    if nu0 is not None:
        nu0 *= 1e9
        x0=h*nu0/(k_B*temp)
        dBdT0 = x0**4 * np.exp(x0) / (np.exp(x0)-1)**2.
        return  dBdT / dbdT0
    else:
        return dBdT

def fn_BnuT(nu, temp = 2.725):
    nu *= 1e9
    x=h*nu/(k_B*temp)

    t1 = 2 * h * nu**3./ c**2.
    t2 = 1./ (np.exp(x)-1.)
    return t1 * t2

def coth(x):
    return (np.exp(x) + np.exp(-x)) / (np.exp(x) - np.exp(-x))

def compton_y_to_delta_Tcmb(freq1, freq2 = None, Tcmb = 2.73):

    """ad
    c.f:  table 1, sec. 3 of arXiv: 1303.5081; 
    table 8 of http://arxiv.org/pdf/1303.5070.pdf
    no relativistic corrections included.
    freq1, freq2 = frequencies in GHz to cover the bandpass
    freq2 = None will force freq1 to be the centre frequency
    """

    if freq1<1e3: freq1 = freq1 * 1e9

    if not freq2 is None:
        if freq2<1e3: freq2 = freq2 * 1e9
        freq = np.arange(freq1,freq2,delta_nu)
    else:
        freq = np.asarray([freq1])

    x = (h * freq) / (k_B * Tcmb)
    g_nu = x * coth(x/2.) - 4.

    return Tcmb * np.mean(g_nu)

def get_cl_dust(freq1, freq2, fg_model = 'george15', freq0 = 150, spec_index_dg_po = 1.505 - 0.077, spec_index_dg_clus = 2.51-0.2, Tcib = 20.):

    if fg_model:
        el, cl_dg_po_freq0 = get_foreground_power_spt('DG-Po', freq1 = freq0, freq2 = freq0)
        el, cl_dg_clus_freq0 = get_foreground_power_spt('DG-Cl', freq1 = freq0, freq2 = freq0)
        el_norm = 3000

    #conert to Dls
    Dls_fac = el * (el+1)/2/np.pi
    Dls_dg_po = Dls_fac * cl_dg_po_freq0
    Dls_dg_clus = Dls_fac * cl_dg_clus_freq0

    nr = ( fn_dB_dT(freq0) )**2.
    dr = fn_dB_dT(freq1) * fn_dB_dT(freq2)

    epsilon_nu1_nu2 = nr/dr

    bnu1 = fn_BnuT(freq1, temp = Tcib)
    bnu2 = fn_BnuT(freq2, temp = Tcib)
    bnu0 = fn_BnuT(freq0, temp = Tcib)

    etanu1_dg_po = ((1.*freq1*1e9)**spec_index_dg_po) * bnu1
    etanu2_dg_po = ((1.*freq2*1e9)**spec_index_dg_po) * bnu2
    etanu0_dg_po = ((1.*freq0*1e9)**spec_index_dg_po) * bnu0

    etanu1_dg_clus = ((1.*freq1*1e9)**spec_index_dg_clus) * bnu1
    etanu2_dg_clus = ((1.*freq2*1e9)**spec_index_dg_clus) * bnu2
    etanu0_dg_clus = ((1.*freq0*1e9)**spec_index_dg_clus) * bnu0

    Dls_dg_po = Dls_dg_po[el == el_norm][0] * epsilon_nu1_nu2 * (1.*etanu1_dg_po * etanu2_dg_po/etanu0_dg_po/etanu0_dg_po) * (el*1./el_norm)**2
    Dls_dg_clus = Dls_dg_clus[el == el_norm][0] * epsilon_nu1_nu2 * (1.*etanu1_dg_clus * etanu2_dg_clus/etanu0_dg_clus/etanu0_dg_clus) * (el*1./el_norm)**0.8

    cl_dg_po = Dls_dg_po / Dls_fac
    cl_dg_clus = Dls_dg_clus / Dls_fac

    cl_dg_po[np.isnan(cl_dg_po)] = 0.
    cl_dg_clus[np.isinf(cl_dg_clus)] = 0.

    return el, cl_dg_po, cl_dg_clus

def get_cl_tsz(freq1, freq2, freq0 = 150, fg_model = 'george15'):

    if fg_model:
        el, cl_tsz_freq0 = get_foreground_power_spt('DG-Po', freq1 = freq0, freq2 = freq0)

    tsz_fac_freq0 = compton_y_to_delta_Tcmb(freq0*1e9)
    tsz_fac_freq1 = compton_y_to_delta_Tcmb(freq1*1e9)
    tsz_fac_freq2 = compton_y_to_delta_Tcmb(freq2*1e9)

    scalefac = tsz_fac_freq1 * tsz_fac_freq2/ (tsz_fac_freq0**2.)

    cl_tsz = cl_tsz_freq0 * scalefac
    cl_tsz[np.isnan(cl_tsz)] = 0.
    cl_tsz[np.isinf(cl_tsz)] = 0.

    return el, cl_tsz

def get_cl_radio(freq1, freq2, freq0 = 150, fg_model = 'george15', spec_index_rg = -0.9):

    if fg_model:
        el, cl_rg_freq0 = get_foreground_power_spt('RG', freq1 = freq0, freq2 = freq0)
        el_norm = 3000

    #conert to Dls
    Dls_fac = el * (el+1)/2/np.pi
    Dls_rg = Dls_fac * cl_rg_freq0

    nr = ( fn_dB_dT(freq0) )**2.
    dr = fn_dB_dT(freq1) * fn_dB_dT(freq2)

    epsilon_nu1_nu2 = nr/dr

    Dls_rg = Dls_rg[el == el_norm][0] * epsilon_nu1_nu2 * (1.*freq1 * freq2/freq0/freq0)**spec_index_rg * (el*1./el_norm)**2

    cl_rg = Dls_rg / Dls_fac

    cl_rg[np.isnan(cl_rg)] = 0.

    return el, cl_rg

def get_cl_galactic(param_dict, component, freq1, freq2, bl1 = None, bl2 = None, pol = 1, which_gal_mask = 0):

    #fix: pol not working yet

    assert component in ['dust', 'sync']

    try:
        which_gal_mask = param_dict['which_gal_mask']
    except:
        pass

    if component == 'dust':
        cl_gal_dic_fname = param_dict['cl_gal_dic_dust_fname']
    elif component == 'sync':
        cl_gal_dic_fname = param_dict['cl_gal_dic_sync_fname']

    cl_gal_dic = np.load(cl_gal_dic_fname, allow_pickle = 1, encoding = 'latin1').item()['cl_dic'][which_gal_mask]

    try:
        cl_gal = cl_gal_dic[ (freq1, freq2) ]
    except:
        cl_gal = cl_gal_dic[ (freq2, freq1) ]

    if bl1 is not None:
        cl_gal /= bl1
    if bl2 is not None:
        cl_gal /= bl2

    el = np.arange( len(cl_gal) )

    return el, cl_gal

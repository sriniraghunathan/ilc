import numpy as np, sys, os, scipy as sc, healpy as H, foregrounds as fg
################################################################################################################

def get_covariance_dic(param_dict, freqarr, nl_dic = None, ignore_fg = []):

    #get the Cls for the baseline band
    el, cl_cmb = fg.get_foreground_power_george_2015('CMB', freq1 = param_dict['freq0'], freq2 = param_dict['freq0'])
    el, cl_ksz = fg.get_foreground_power_george_2015('kSZ', freq1 = param_dict['freq0'], freq2 = param_dict['freq0'])

    cl_dic = {}

    for freq1 in freqarr:
        for freq2 in freqarr:

            if (freq2, freq1) in cl_dic:
                cl_dic[(freq1, freq2)] = cl_dic[(freq2, freq1)]
                continue

            #get dust
            el,  cl_dg_po, cl_dg_clus = fg.get_cl_dust(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_dg_po = param_dict['spec_index_dg_po'], spec_index_dg_clus = param_dict['spec_index_dg_clus'], Tcib = param_dict['Tcib'])

            #get tsz
            el, cl_tsz = fg.get_cl_tsz(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'])

            #get radio
            el, cl_radio = fg.get_cl_radio(freq1, freq2, freq0 = param_dict['freq0'], fg_model = param_dict['fg_model'], spec_index_rg = param_dict['spec_index_rg'])

            if 'cmb' not in ignore_fg:
                cl = np.copy(cl_cmb)
            if 'ksz' not in ignore_fg:
                cl = cl + cl_ksz             
            if 'tsz' not in ignore_fg:
                cl = cl + cl_tsz
            if 'radio' not in ignore_fg:
                cl = cl + cl_radio
            if 'dust' not in ignore_fg:
                cl = cl + cl_dg_po + cl_dg_clus

            #make sure cl start from el=0 rather than el=10 which is the default for SPT G15 results
            lmin = min(el)
            cl = np.concatenate( (np.zeros(lmin), cl) )

            #noise auto power spectrum
            if nl_dic is not None:

                nl = nl_dic[freq1]                    
                if freq1 != freq2: 
                    nl = np.copy(nl) * 0.

                if len(cl) > len(nl):
                    cl = cl[:len(nl)]
                elif len(cl) < len(nl):
                    nl = nl[:len(cl)]

                el = np.arange(len(cl))

            else:
                nl = np.zeros(len(cl))
                #20191116 - fix this: there must be noise correlation in case of atmospheric noise
                print('\n\n\t\tfix me: there must be noise correlation in case of atmospheric noise')

            cl = cl + nl

            cl[np.isnan(cl)] = 0.

            ##loglog(cl); title('%s - %s' %(freq1, freq2)); show()

            cl_dic[(freq1, freq2)] = cl 


    return el, cl_dic  

################################################################################################################
def create_clmat(freqarr, elcnt, cl_dic):
    """
    freqarr  = array of frequency channel
    elcnt = \el index
    cl_dic = cl_cmb + cl_FG auto and cross spectra fr thefrequency channel
    """
    nc = len(freqarr)
    clmat = np.zeros( (nc, nc) )
    for ncnt1, freq1 in enumerate(freqarr):
        for ncnt2, freq2 in enumerate(freqarr):
            clmat[ncnt2, ncnt1] = cl_dic[(freq1, freq2)][elcnt]
    return clmat

################################################################################################################
def residual_power(param_dict, freqarr, el, cl_dic, final_comp = 'CMB', freqcalib_fac = None):

    acap = get_acap(freqarr, final_comp = final_comp, freqcalib_fac = freqcalib_fac)

    cl_cleaned = np.zeros( (len(el)) )
    for elcnt, el in enumerate(el):
        clmat = np.mat( create_clmat(freqarr, elcnt, cl_dic) )
        clinv = sc.linalg.pinv2(clmat)
        
        nr = 1.
        dr = np.dot( acap.T, np.dot(clinv, acap) )
        cl_cleaned[elcnt] = np.asarray(nr/dr).squeeze()

    if final_comp == 'CMB':
        el, cl_ini = fg.get_foreground_power_george_2015('CMB', freq1 = param_dict['freq0'])
    elif final_comp == 'kSZ':
        el, cl_ini = fg.get_foreground_power_george_2015('kSZ', freq1 = param_dict['freq0'])
    elif final_comp == 'tSZ' or final_comp == 'comptony':
        el, cl_ini = fg.get_foreground_power_george_2015('tSZ', freq1 = param_dict['freq0'])
        if final_comp == 'comptony':
            freqscale_fac = compton_y_to_delta_Tcmb(param_dict['freq0'] * 1e9)
            cl_ini = cl_ini/freqscale_fac

    cl_ini = cl_ini[:len(cl_cleaned)]

    cl_residual = cl_cleaned - cl_ini
    
    return cl_residual

################################################################################################################

def get_acap(freqarr, final_comp = 'CMB', freqcalib_fac = None):

    nc = len(freqarr)

    if freqcalib_fac is None: freqcalib_fac = np.ones(nc)

    if final_comp.lower() == 'cmb':
        freqscale_fac = np.ones(nc)

    elif final_comp.lower() == 'tsz' or final_comp.lower() == 'y':

        freqscale_fac = []
        for freq in sorted( freqarr ):
            freqscale_fac.append( compton_y_to_delta_Tcmb(freq * 1e9) )

        freqscale_fac = np.asarray( freqscale_fac )

        if final_comp.lower() == 'tsz': #tsz at 150 GHz
            freqscale_fac = freqscale_fac / freqscale_fac[1]

    acap = np.zeros(nc) + (freqscale_fac * freqcalib_fac) #assuming CMB is the same and calibrations factors are same for all channel
    acap = np.mat(acap).T #should be nc x 1

    return acap

################################################################################################################

def get_ilc_map(final_comp, map_dic, bl_dic, nside, lmax, cl_dic = None, nl_dic = None, lmin = 10, freqcalib_fac = None, ignore_fg = [], full_sky = 0, mapparams = None, apod_mask = None):

    """
    inputs:
    final_comp: 'CMB' or 'tsz' or 'y'

    map_dic: dicionary containing
    map_dic[95] = map_95 where map_95 is either a healpix map or a SPT3G map object or a simple flat-sky 2D numpy array
    map_dic[150] = map_150 
    map_dic[220] = map_220

    bl_dic = {} #dictionary containing Bl of different frequencies - either a 1d array for healpix maps or 2D array for SPT3G/flat-sky maps
    bl_dic[95], bl_dic[150], bl_dic[220]

    nside = healpix nside set to None for flat-sky maps

    lmax = lmax used for ILC (covariance)

    cl_dic = same as beam dic and contains Cl of different frequencies - either a 1d array for healpix maps or 2D array for SPT3G/flat-sky maps
    nl_dic = same as beam dic and contains Nl of different frequencies - either a 1d array for healpix maps or 2D array for SPT3G/flat-sky maps

    lmin = minimum \el to use

    freqcalib_fac = array containing calibration factors for different frequencies. for example: freqcalib_fac = [1., 1., 1.]

    """
    spt3g_maps = 0
    freqarr, maparr = [], []
    for keyname in sorted(map_dic):
        freqarr.append( keyname )
        curr_map = map_dic[keyname]
        if isinstance(curr_map, core.G3Frame):
            spt3g_maps = 1
            #get the apodisation mask before converting the map inot an array
            if str(apod_mask) == 'from_weight':
                apod_mask = apodmask.make_border_apodization(curr_map, radius_arcmin=90.0)
            curr_map = np.array( curr_map['T'] ) ## / core.G3Units.uK
            curr_map[np.isnan(curr_map)] = 0.
            curr_map[np.isinf(curr_map)] = 0.


        if apod_mask is not None:
            curr_map = curr_map * apod_mask

        maparr.append( curr_map )

    #get covariance
    if cl_dic is None:
        if spt3g_maps:
            el, cl_dic = get_spt3g_covariance_dic(map_dic, lmin, lmax)
        else:
            el, cl_dic = get_covariance_dic(freqarr, nl_dic = nl_dic, ignore_fg = ignore_fg)

    #get weights
    weightsarr = get_multipole_weightsarr(final_comp, freqarr, el, cl_dic, lmin, freqcalib_fac, ignore_fg)
    weightsarr_1d = np.copy(weightsarr)

    #convert weights to 2D if flat-sky
    if not full_sky:
        assert mapparams is not None
        weightsarr_2D = []
        for currW in weightsarr:
            el = np.arange(len(currW))
            currW_2D = flatsky.cl_to_cl2d(el, currW, mapparams) 
            weightsarr_2D.append(currW_2D)
        weightsarr = np.asarray( weightsarr_2D )


    #rebeaming
    rebeamarr = misc.rebeam( bl_dic )

    #modify weights to include rebeam
    weightsarr = weightsarr * rebeamarr
    '''
    plot(weightsarr[0], 'k-'); plot(weightsarr[1], 'r-'); plot(weightsarr[2], 'g-')
    weightsarr = weightsarr * rebeamarr
    plot(weightsarr[0], 'k--'); plot(weightsarr[1], 'r-'); plot(weightsarr[2], 'g-')
    show(); sys.exit()
    sys.exit()
    '''

    #get ilc map now
    ilc_map = apply_ilc_weightsarr(maparr, weightsarr, nside, lmax, full_sky = full_sky)

    return ilc_map, weightsarr_1d

################################################################################################################

def apply_ilc_weightsarr(maparr, weightsarr, nside, lmax, full_sky = 0, verbose = 0):

    '''
    clf()
    freqs = [95, 150]#, 220]
    colordic = {95: 'darkblue', 150: 'green', 220: 'darkred'}
    for frqcntr, freq in enumerate( freqs ):
        plot(weightsarr[frqcntr], color = colordic[freq], label = r'%s' %(freq))
    plot(np.sum(weightsarr, axis = 0), 'k', label = r'Sum')
    legend(loc = 1);show(); sys.exit()
    '''

    #get the ilc combined map now
    weighted_maparr = []
    for mm in range(len(maparr)):
        if full_sky:
            map_alm = H.map2alm(maparr[mm], lmax = lmax)
            curr_weight = weightsarr[mm][:lmax]
            
            map_alm_weighted = H.almxfl(map_alm, curr_weight)

            '''
            map_alm_weighted = np.zeros_like(map_alm)
            for el in range(len(curr_weight)):
                alm_inds = H.Alm.getidx(lmax, el, np.arange(el + 1))
                map_alm_weighted[alm_inds] =  curr_weight[el] * map_alm[alm_inds]
            '''
            #plot(map_alm_weighted.real, 'r');show();sys.exit()
            map_weighted = H.alm2map(map_alm_weighted, nside = nside, verbose = verbose, lmax = lmax)
            #clf();plot(curr_weight); show()
            #H.mollview(maparr[mm], sub = (1,2,1)); H.mollview(map_weighted, sub = (1,2,2)); show()
        else:
            curr_map = maparr[mm]
            weightsarr[mm][np.isnan(weightsarr[mm])]=0.
            weightsarr[mm][np.isinf(weightsarr[mm])]=0.
            map_weighted = np.fft.ifft2( np.fft.fft2(curr_map) * weightsarr[mm] ).real
        weighted_maparr.append(map_weighted)

    ilc_map = np.sum(weighted_maparr, axis = 0)

    return ilc_map

################################################################################################################

def get_multipole_weightsarr(final_comp, freqarr, el, cl_dic, lmin, freqcalib_fac, ignore_fg):

    acap = get_acap(freqarr, final_comp = final_comp, freqcalib_fac = freqcalib_fac)

    nc = len(freqarr)

    #get weightsarr
    if np.ndim(el) == 1:
        weightsarr = np.zeros( (nc, el.shape) )
        for elcnt, curr_el in enumerate( el ):
            if curr_el <= lmin: continue ## or el>=lmax: continue
            clmat = np.mat( create_clmat(freqarr, elcnt, cl_dic) )
            clinv = sc.linalg.pinv2(clmat)
            
            nr = np.dot(clinv, acap)
            dr = np.dot( acap.T, np.dot(clinv, acap) )

            weight_mat = np.asarray(nr/dr).squeeze()

            weightsarr[:, elcnt] = weight_mat
    else:
        #2d- PSDs
        from IPython import embed; embed()
        sys.exit()
        lx, ly = el
        for elcnt1, curr_lx in enumerate( lx ):
            for elcnt2, curr_ly in enumerate( ly ):
                curr_el = np.sqrt(curr_lx**2. + curr_ly**2.)
                if curr_el <= lmin: continue ## or el>=lmax: continue
            
                clmat = np.mat( create_clmat(freqarr, elcnt, cl_dic) )
                clinv = sc.linalg.pinv2(clmat)
            
            nr = np.dot(clinv, acap)
            dr = np.dot( acap.T, np.dot(clinv, acap) )

            weight_mat = np.asarray(nr/dr).squeeze()

            weightsarr[:, elcnt] = weight_mat

    return weightsarr

################################################################################################################

def get_spt3g_covariance_dic(map_dic, lmin, lmax, apod_mask = 'from_weight', return_2d = 1):

    freqarr = sorted(map_dic.keys())
    cl_dic = {}
    for cntr1, freq1 in enumerate( freqarr ):
        for cntr2, freq2 in enumerate( freqarr ):
            if (freq2, freq1) in cl_dic:
                cl_dic[(freq1, freq2)] = cl_dic[(freq2, freq1)]
                continue
            map1, map2 = map_dic[freq1], map_dic[freq2]

            if return_2d:
                mapps = map_analysis.calculate_powerspectra(map1, map2, apod_mask = apod_mask, return_2d = 1)['TT']

                mapps[np.isnan(mapps)] = 0.
                mapps[np.isinf(mapps)] = 0.
                cl_dic[(freq1, freq2)] = mapps 

                from IPython import embed; embed()

                el = mapps.get_lxly()

            else:
                cl = map_analysis.calculate_powerspectra(map1, map2, lmin = lmin, lmax = lmax, apod_mask = apod_mask)['TT']
                cl = np.concatenate( (np.zeros(lmin), cl) )
                cl[np.isnan(cl)] = 0.
                cl[np.isinf(cl)] = 0.
                cl_dic[(freq1, freq2)] = cl 
            
                el = np.arange(lmin, lmax)
                el = np.concatenate( (np.zeros(lmin), el) )

    return el, cl_dic

################################################################################################################

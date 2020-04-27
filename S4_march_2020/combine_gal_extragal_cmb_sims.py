import numpy as np, healpy as H, sys, glob, scipy as sc, warnings, os
sys.path.append('/data48/sri/git/tools/')
sys.path.append('/data48/sri/git/ilc/')

import tools, misc
import foregrounds as fg

reqd_freq = int( sys.argv[1] )

##############################################################################
opfolder = 'output'
if not os.path.exists(opfolder): os.system('mkdir %s' %(opfolder))
nside_out = 2048 ##4096
verbose = 0
lmax = 6000

produce_radio_skies = 0
if produce_radio_skies: #specs for radio sky sims - Gaussian reals from SPT-SZ G15 measurements
    lmin = 10
    freq0 = 150
    fg_model = 'george15'
    spec_index_rg = -0.9
    pol_frac_per_cent_radio = 0.03
##############################################################################
#specifications: 
#Noise values come from: https://cmb-s4.org/wiki/index.php/Expected_Survey_Performance_for_Science_Forecasting#Instrument_Definition 
#Deep and Wide Field (Hi-res) field
freq_dic = {
#freq: [beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P] 
#20: [11.0, None, None, None, None, None, None],
30: [7.3, 21.8, 471, 3.5, 30.8, 700, 1.4],
40: [5.5, 12.4, 428, 3.5, 17.6, 700, 1.4], 
95: [2.3, 2.0, 2154, 3.5, 2.9, 700, 1.4],
145: [1.5, 2.0, 4364, 3.5, 2.8, 700, 1.4],
220: [1.0, 6.9, 7334, 3.5, 9.8, 700, 1.4],
270: [0.8, 16.7, 7308, 3.5, 23.6, 700, 1.4],
}

##############################################################################
#modified beams come from https://docs.google.com/spreadsheets/d/1Q76Il6g33PeNhQEv1xkHe3PWpomxBPRcC7wgp7TQfus/edit#gid=26429016
mod_freq_beam_dic = {
#freq: modified_freq, modifed_beam_arcmins
20: [20, 10.0], 
30: [27, 7.4],
40: [39, 5.1],
95: [93, 2.2], 
145: [145, 1.4], 
220: [225, 1.0], 
270: [278, 0.9],
}
##############################################################################
'''
#corresponding frequencies in Caterina Umilta  (Galactic) / Yuuki Omori (Extragalactic) simulations
cumilta_yomori_freq_dic = {
#freq: freq_in_CUmilta_sims, freq_in_YOmori_sims
20: [20, 20], 
30: [27, 27], 
40: [39, 39], 
95: [93, 93],
145: [145, 145],
220: [225, 225],
270: [278, 278],    
}
'''
searchstr_gal_dust = 'sims_from_others/CUmilta/ampmod_maps/Ampmod_map_dust_lmax_5200_freq_xxx_rel0000.fits'
searchstr_gal_sync = 'sims_from_others/CUmilta/ampmod_maps/Ampmod_map_sync_lmax_5200_freq_xxx_rel0000.fits'
searchstr_extragal_cmb = 'sims_from_others/YOmori/cmbs4/mdpl2_sky_cmbs4_lcmbNG_tszNGmasked_kszNG_cibNG_xxxghz.fits'
##############################################################################
if produce_radio_skies:#foreground only to create Gaussian realisations of radio sky
    cl_radio_dic = {}
    for freq1 in sorted(freq_dic):
        el, cl_radio = fg.get_cl_radio(freq1, freq1, freq0 = freq0, fg_model = fg_model, spec_index_rg = spec_index_rg)
        cl_radio = np.concatenate( (np.zeros(lmin), cl_radio) )
        el = np.arange(len(cl_radio))

        cl_radio = cl_radio[:lmax]
        el = el[:lmax]

        cl_radio_dic[freq1] = cl_radio
##############################################################################
npix = H.nside2npix(nside_out)
if produce_radio_skies: random_seed = 100

for freq in sorted(freq_dic):

    if freq != reqd_freq: continue

    print('\n\t%s' %(freq))

    #opfname = '%s/%03d_nside%s.fits' %(opfolder, freq, nside_out)
    #if os.path.exists(opfname): continue

    #get specs
    beam_arcmins, white_noise_T, elknee_T, alphaknee_T, whitenoise_P, elknee_P, alphaknee_P = freq_dic[freq]

    #get sim filenames
    #freq_mod_cumilta, freq_mod_yomori =  cumilta_yomori_freq_dic[freq]
    freq_mod, beam_arcmins_mod = mod_freq_beam_dic[freq]
    gal_dust_fname = searchstr_gal_dust.replace('xxx', '%03d' %(freq_mod))
    gal_sync_fname = searchstr_gal_sync.replace('xxx', '%03d' %(freq_mod))
    extragal_cmb_fname = searchstr_extragal_cmb.replace('xxx', '%s' %(freq_mod)) ##'%03d' %(freq_mod))

    '''
    #read files and combine them into one
    fname_arr = [gal_dust_fname, gal_sync_fname, extragal_cmb_fname]
    hmap = np.zeros( (3, npix) )
    for fcnt, fname in enumerate( fname_arr ):
        print('\t\t\t%s' %fname)
        curr_hmap = H.read_map( fname, verbose = verbose, field = (0, 1, 2) ) #TQU
        curr_nside = H.get_nside(curr_hmap)
        if curr_nside != nside_out:
            print('\t\t\t\tdownsampling from %s to %s' %(curr_nside, nside_out))
            curr_hmap = H.ud_grade(curr_hmap, nside_out)
        if fname.find('YOmori')>-1: #then add beam
            print('\t\t\t\tadding beam: %s arcmins' %(beam_arcmins_mod))
            curr_hmap = H.smoothing(np.copy(curr_hmap), fwhm = np.radians(beam_arcmins_mod/60.), lmax = lmax, verbose = verbose)
        hmap += curr_hmap
    '''


    #read files and combine them into one
    for iter in range(2):
        if iter == 0:
            opfname = '%s/gal_%03d_nside%s.fits' %(opfolder, freq, nside_out)
            fname_arr = [gal_dust_fname, gal_sync_fname]#, extragal_cmb_fname]
        else:
            opfname = '%s/cmb_extragal_%03d_nside%s.fits' %(opfolder, freq, nside_out)
            fname_arr = [extragal_cmb_fname]


        hmap = np.zeros( (3, npix) )
        for fcnt, fname in enumerate( fname_arr ):
            print('\t\t\t%s' %fname)
            curr_hmap = H.read_map( fname, verbose = verbose, field = (0, 1, 2) ) #TQU
            curr_nside = H.get_nside(curr_hmap)
            if curr_nside != nside_out:
                print('\t\t\t\tdownsampling from %s to %s' %(curr_nside, nside_out))
                curr_hmap = H.ud_grade(curr_hmap, nside_out)
            if fname.find('YOmori')>-1: #then add beam
                print('\t\t\t\tadding beam: %s arcmins' %(beam_arcmins_mod))
                curr_hmap = H.smoothing(np.copy(curr_hmap), fwhm = np.radians(beam_arcmins_mod/60.), lmax = lmax, verbose = verbose)
            hmap += curr_hmap

        H.fitsfunc.write_map(opfname, hmap)



    '''
    opfname = '%s/radio_%03d_nside%s.fits' %(opfolder, freq, nside_out)
    hmap = np.zeros( (3, npix) )
    if produce_radio_skies: #create radio skies now from SPT-SZ measurements
        print('\t\t\tcreate radio Gaussian skies now from %s SPT-SZ measurements' %(fg_model))
        np.random.seed( random_seed ) #correalted foregrounds
        cl_list = [cl_radio_dic[freq], 0. , cl_radio_dic[freq]*(pol_frac_per_cent_radio**2.)/2., cl_radio_dic[freq]*(pol_frac_per_cent_radio**2.)/2.] #Old ordering for spt.uchicago.edu healpy: TT, TE, EE, BB
        radio_map = H.synfast(cl_list, nside = nside_out, lmax = lmax, verbose = verbose)
        radio_map = H.smoothing(np.copy(radio_map), fwhm = np.radians(beam_arcmins_mod/60.), lmax = lmax, verbose = verbose)
        hmap += radio_map
 
    H.fitsfunc.write_map(opfname, hmap)
    '''



sys.exit()

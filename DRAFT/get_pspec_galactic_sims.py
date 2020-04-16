################################################################################################################

def healpix_rotate_coords(hmap, coord):
    """
    coord = ['C', 'G'] to convert a map in RADEC to Gal.    
    """

    #get map pixel
    pixel = np.arange(len(hmap))

    #get angles in this map first
    nside = H.get_nside(hmap)
    angles = H.pix2ang(nside, pixel)

    #roate the angles to the desired new coordinate
    rotated_angles = H.Rotator(coord=coord)(*angles)

    #get the rotated pixel values
    rotated_pixel = H.ang2pix(nside, *rotated_angles)

    #initialise new map
    rot_hmap = np.zeros(len(pixel))

    #push the original map pixel to the new map (in the rotated pixel positions)
    rot_hmap[rotated_pixel] = hmap[pixel]

    return rot_hmap

############################################################
############################################################

import healpy as H, numpy as np, glob, sys, os

local = 1
if str(os.getcwd()).find('sri')>-1: local = 0


dust_or_sync = sys.argv[1] ##'sync' ##'dust'
lmax = 5200
nuarr = [20, 27, 39, 93, 145, 225, 278]
verbose = 0
nside = 2048

#lmax = 2000
#nside = 1024

log_file = 'tmp/pspec_%s.txt' %(dust_or_sync)
lf = open(log_file, 'w'); lf.close()

if local:
    sim_folder = 'galactic/CUmilta/ampmod_maps/'
else:    
    sim_folder = '/data48/sri/git/ilc/S4_march_2020/sims_from_others/CUmilta/ampmod_maps/'

opfname = '%s/cls_gal_%s.npy' %(sim_folder, dust_or_sync)

if not local:

    totthreads = 4
    os.putenv('OMP_NUM_THREADS',str(totthreads))

    #get filename prefix
    fname_pref = 'Ampmod_map_%s_lmax_5200_freq_xxx_rel0000.fits' %(dust_or_sync)

    map_dic = {}
    for nucntr, nu in enumerate( nuarr ):
        fname = '%s/%s' %(sim_folder, fname_pref)
        fname = fname.replace('xxx', '%03d' %(nu))

        logline = '\n\t%s\n' %fname
        lf = open(log_file,'a'); lf.writelines('%s\n' %(logline));lf.close()
        print(logline)

        currmap = H.read_map(fname, verbose = 0, field = (0,1,2))
        if H.get_nside(currmap) != nside:
            currmap = H.ud_grade(currmap, nside_out = nside)

        ##from IPython import embed; embed()

        map_dic[nu] = currmap

    #now get masks
    if (1):
        planck_mask_fname = '%s/HFI_Mask_GalPlane-apo0_2048_R2.00.fits' %(sim_folder)
        planck_mask = H.read_map(planck_mask_fname, verbose = verbose, field = (1,2,3))
        if H.get_nside(planck_mask) != nside:
            planck_mask = H.ud_grade(planck_mask, nside_out = nside)

        '''
        planck_mask = H.smoothing(np.copy(planck_mask), fwhm = np.radians(5.), lmax = lmax, verbose = verbose)
        thresh = 0.4
        for mask_iter in range(len(planck_mask)):
            planck_mask[mask_iter][planck_mask[mask_iter]<thresh] = 0.
            planck_mask[mask_iter][planck_mask[mask_iter]!=0] = 1.
        '''

        cmbs4_hit_map_fname = '%s/high_cadence_hits_el30_cosecant_modulation.fits' %(sim_folder)
        cmbs4_hit_map = H.read_map(cmbs4_hit_map_fname, verbose = verbose)
        cmbs4_hit_map[cmbs4_hit_map!=0] = 1.
        if H.get_nside(cmbs4_hit_map) != nside:
            cmbs4_hit_map = H.ud_grade(cmbs4_hit_map, nside_out = nside)

        
    tot_masks = len(planck_mask)

    logline = '\n\tget masks now\n'
    lf = open(log_file,'a'); lf.writelines('%s\n' %(logline));lf.close()
    print(logline)

    mask_arr = []
    for mask_iter in range(tot_masks):
        mask = np.copy(planck_mask[mask_iter])

        #simple rotation from gal to celestial
        mask = healpix_rotate_coords(mask, coord = ['G', 'C'])
        mask = H.smoothing(np.copy(mask), fwhm = np.radians(10.), lmax = lmax, verbose = verbose)
        thresh = 0.4
        mask[mask<thresh] = 0.
        mask[mask!=0] = 1.

        mask_arr.append( mask )   

    mask_arr = np.asarray(mask_arr)
    mask_arr = mask_arr * cmbs4_hit_map
    fsky_arr = np.mean(mask_arr, axis = 1)

    if (0):
        #from IPython import embed; embed()
        H.mollview(mask_arr[0], sub = (1,3,1)); H.mollview(mask_arr[1], sub = (1,3,2)); H.mollview(mask_arr[2], sub = (1,3,3)); show()
    
    resdic = {}
    resdic['fsky_arr'] = fsky_arr
    resdic['lmax'] = lmax
    resdic['cl_dic'] = {}
    for mask_iter in range(tot_masks):
        resdic['cl_dic'][mask_iter] = {}
        for nu1 in nuarr:
            for nu2 in nuarr:

                print(nu1, nu2)

                logline = '\n\tMask = %s: (%s,%s)\n' %(mask_iter, nu1, nu2)
                lf = open(log_file,'a'); lf.writelines('%s\n' %(logline));lf.close()
                print(logline)

                map1, map2 = map_dic[nu1], map_dic[nu2]

                curr_mask = mask_arr[mask_iter]
                fsky = np.mean(curr_mask)
                map1 = map1 * curr_mask
                map2 = map2 * curr_mask

                if (nu2, nu1) in resdic['cl_dic'][mask_iter]: continue
                curr_cl = H.anafast(map1, map2, lmax = lmax)
                curr_cl /= fsky
                resdic['cl_dic'][mask_iter][(nu1, nu2)] = curr_cl

    np.save(opfname, resdic)
    sys.exit()



#plot
#opfname = 'cl_websky.npy'
from pylab import *
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
#mpl.rcParams['hatch.linewidth'] = 0.5
rcParams['font.family'] = 'serif'
rc('text.latex',preamble=r'\usepackage{/Users/sraghunathan/.configs/apjfonts}')

resdic = np.load(opfname, allow_pickle = 1).item()
cl_arr = resdic['cl_arr']
fname_arr = resdic['fname_arr']
lmax = resdic['lmax']

spt_el = 3000.
spt_points_arr = [0, (9.6, 3.46), (66.1, 50.6), 2.9]
spt_points_err_arr = [(0., 0.), (0.36, 0.54), (3.1, 4.4), (1.3, 1.3)]

cl_mod_arr = np.zeros( (len(fname_arr)-1, lmax + 1) )
cl_mod_arr[:3] = np.copy(cl_arr)[:3]
cl_mod_arr[3] = np.sum(cl_arr[4:], axis = 0)


el = np.arange( lmax + 1 )
elnorm = 3000
elind = np.where(el == elnorm)
dl = cl_mod_arr[1] * el * ( el + 1)/ 2 / np.pi    
fac = ( np.sum(spt_points_arr[1]) / dl[elind] )**0.5
print(dl[elind], fac, dl[elind] * fac)

clf()
colorarr = ['navy', 'orangered', 'darkred', 'black']
ax = subplot(111, yscale = 'log')#, xscale = 'log')
for cntr, cl in enumerate( cl_mod_arr[:-1] ):
    if cntr == 0: continue
    lab = fname_arr[cntr].split('/')[-1].replace('_','-').replace('.fits', '').replace('0','').replace('nu','').upper()
    colorval = colorarr[cntr]
    if lab.lower().find('ksz')==-1:
        #lab = r'Unmasked: %s' %(lab) 
        lab = r'Masked + scaled: %s' %(lab) 

    ##cl = cl  * div_arr[cntr]**2.
    cl = ( cl ) * fac**2.
    el = np.arange( len(cl) )
    dl = cl * el * ( el + 1)/ 2 / np.pi    
    plot(el, dl, lw = 1., color = colorval, label = lab)
    spt_data_point = np.sum(spt_points_arr[cntr])
    spt_data_point_err = np.sqrt( spt_points_err_arr[cntr][0]**2. + spt_points_err_arr[cntr][1]**2. )
    print(spt_data_point, spt_data_point_err)
    errorbar(spt_el, spt_data_point, yerr = spt_data_point_err, color = colorval, marker = 'o', capsize = 2.)

    if cntr == 1 or cntr == 2:
        if cntr == 1:
            freq1 = 150
        elif cntr == 2:
            freq1 = 220

        el, cl_dg_po, cl_dg_clus = fg.get_cl_dust(freq1, freq1, freq0 = 150)
        dl_dg_po = cl_dg_po * el * ( el + 1)/ 2 / np.pi    
        dl_dg_clus = cl_dg_clus * el * ( el + 1)/ 2 / np.pi    
        dl_spt = dl_dg_po + dl_dg_clus
        spt_lwval = 2.
        plot(el, dl_spt, lw = spt_lwval, color = colorval,alpha = 0.5)#, label = r'SPT')
        plot(el, dl_dg_po, lw = spt_lwval, ls = '--', color = colorval, alpha = 0.5)#, label = r'SPT: Po')
        plot(el, dl_dg_clus, lw = spt_lwval, ls = ':' , color = colorval, alpha = 0.5)#, label = r'SPT: Clus')

if (0): #add sehgal
    sehgal_dic = np.load('cross_power_CIB_sehgal_SO_sims.fits.npy', allow_pickle = 1).item()
    sehgal_cl_dic = sehgal_dic['cl_dic']
    sehgal_lmax = sehgal_dic['lmax']
    sehgal_cl = sehgal_cl_dic[(217, 217)]

    sehgal_el = np.arange( len(sehgal_cl) )
    sehgal_dl = sehgal_cl * sehgal_el * ( sehgal_el + 1)/ 2 / np.pi    
    plot(sehgal_el, sehgal_dl, lw = 1., color = 'darkred', label = r'Sehgal 217 $\times$ 217', alpha = 0.5)

if (1):
    plot([], [], lw = spt_lwval, color = 'k',alpha = 0.5, label = r'SPT')
    plot([], [], lw = spt_lwval, ls = '--', color = 'k',alpha = 0.5, label = r'SPT: Po')
    plot([], [], lw = spt_lwval, ls = ':' , color = 'k',alpha = 0.5, label = r'SPT: Clus')
    #plot([], [], marker = '.', label = r'Masked: SPT (G15)', ls = 'None')

ylim(0.5, 2e3)
#plot([], [], marker = 'None', label = r'Unmasked Websky')
legend(loc = 2, fancybox = 1, fontsize = 12)
xlabel(r'Multipole $\ell$'); ylabel(r'D$_{\ell}$')
xlim(100., 4000)
show();


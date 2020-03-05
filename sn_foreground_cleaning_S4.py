
def coth(x):
        return (np.exp(x) + np.exp(-x)) / (np.exp(x) - np.exp(-x))

def compton_y_to_delta_Tcmb(nu):

        """ad
        c.f:  table 1, sec. 3 of arXiv: 1303.5081; 
        table 8 of http://arxiv.org/pdf/1303.5070.pdf
        no relativistic corrections included.
        nu_1, nu_2 = frequencies in GHz to cover the bandpass
        nu_2 = None will force nu_1 to be the centre frequency
        """

        x = h * nu / (k_B * Tcmb)
        g_nu = x * coth(x/2.) - 4.

        return Tcmb * np.mean(g_nu)

def fn_get_A_syn(nu, nu_ref = 150e9, beta = -3.1):
	return (nu/nu_ref)**beta

def fn_get_A_dust(nu, nu_ref = 150e9, beta = 1.59, Td = 19.6):
	t1 =  (nu/nu_ref)**(beta + 1)
	t2_nr = np.exp(h * nu_ref / (k_B * Td)) - 1.
	t2_dr = np.exp(h * nu / (k_B * Td)) - 1.

	t2 = t2_nr / t2_dr

	return t1 * t2

def fn_get_A_tSZ(nu, nu_ref = 150e9):

	ysz_Tsz_conv_fac_ref = compton_y_to_delta_Tcmb(nu_ref)
	ysz_Tsz_conv_fac = compton_y_to_delta_Tcmb(nu)

	return ysz_Tsz_conv_fac/ysz_Tsz_conv_fac_ref


import numpy as np
import modules.scl_cmb as scl_cmb, argparse
import scipy as sc
from pylab import *

sims = scl_cmb.simulations()

h=6.62607004e-34 #Planck constant in m2 kg / s
k_B=1.38064852e-23 #Boltzmann constant in m2 kg s-2 / K-1
Tcmb = 2.73 #Kelvin

'''
#pol = 1
which_exp = sys.argv[1]#'act'
pol = int(sys.argv[2])
'''

parser = argparse.ArgumentParser(description='')
parser.add_argument('-which_exp', dest='which_exp', action='store', help='which_exp', type=str, required = True) #some random number
parser.add_argument('-pol', dest='pol', action='store', help='pol', type=int, required = True)
parser.add_argument('-year', dest='year', action='store', help='year', type=float, default = -1.)
parser.add_argument('-totyears', dest='totyears', action='store', help='totyears', type=float, default = 5.) #only matters for 3G

#20190528
parser.add_argument('-s4_noise_scale_factor', dest='s4_noise_scale_factor', action='store', help='s4_noise_scale_factor', type=float, default = 1.)
parser.add_argument('-s4_beam_scale_factor', dest='s4_beam_scale_factor', action='store', help='s4_beam_scale_factor', type=float, default = 1.)


args = parser.parse_args()
which_exp = args.which_exp
pol = args.pol
year = args.year
totyears = args.totyears
s4_noise_scale_factor = args.s4_noise_scale_factor
s4_beam_scale_factor = args.s4_beam_scale_factor

if year == -1:
	yearscaling = -1
else:
	yearscaling = np.sqrt(totyears * 1./year)

print '\n\tExp = %s; Pol = %s\n' %(which_exp, pol)

if which_exp == 'spt':
	nuarr = [90e9,150e9, 220e9]
	#noisearr = [3,2,4]
	noisearr = [4.5,2.5,4.5]
	beamarr = [1.7,1.2,1.]
elif which_exp == 'CMB-HD':
	nuarr = [90e9,150e9,220e9]
	noisearr = [0.1, 0.1, 0.31]
	beamarr = [0.27, 0.17, 0.12]
elif which_exp == 'spt3g_july_2018':
	nuarr = [90e9,150e9, 220e9]
	noisearr = [3.,2.2,8.8]
	beamarr = [1.7,1.2,1.]
elif which_exp == 'spt3g_summerfield_2019_1500':
	nuarr = [90e9,150e9, 220e9]
	noisearr = [11.7, 10.0, 41.0]
	beamarr = [1.7,1.2,1.]
elif which_exp == 'spt3g_summerfield_2019_2500':
	nuarr = [90e9,150e9, 220e9]
	noisearr = [15.1, 12.9, 53.0]
	beamarr = [1.7,1.2,1.]
elif which_exp == 'spt3g_summerfield_2019_3000':
	nuarr = [90e9,150e9, 220e9]
	noisearr = [16.5, 14.2, 58.0]
	beamarr = [1.7,1.2,1.]
elif which_exp == 'spt3g_2019f2f_projectedB':
	nuarr = [90e9,150e9, 220e9]
	noisearr = [3.56,2.83,10.46]
	beamarr = [1.7,1.2,1.]
elif which_exp == 'spt3g_2017_2018':
	nuarr = [90e9,150e9, 220e9]
	noisearr_pol = np.asarray( [27., 20., 56.] )
	noisearr = noisearr_pol/1.414
	beamarr = [1.6,1.2,1.1] ##[1.7,1.2,1.]
elif which_exp == 'spt3g_2019':
	nuarr = [90e9,150e9, 220e9]
	noisearr_pol = np.asarray( [9.1, 7.6, 30.4] )
	noisearr = noisearr_pol/1.414
	beamarr = [1.6,1.2,1.1]
elif which_exp == 'spt3g_may_2018':
	nuarr = [90e9,150e9, 220e9]
	noisearr = [3.3,2.8,9.8]
	beamarr = [1.7,1.2,1.]
elif which_exp == 'spt3gdeep':
	nuarr = [90e9,150e9, 220e9]
	noisearr = [2.8,2.6,6.6]
	beamarr = [1.7,1.2,1.]
elif which_exp == 'spt3gwide':
	nuarr = [90e9,150e9, 220e9]
	noisearr = [3.6,3.3,8.5]
	beamarr = [1.7,1.2,1.]
elif which_exp == 'act':
	nuarr = [90e9,150e9, 220e9]
	noisearr = [8,7,25]
	beamarr = [2.2,1.4,1.]
	if pol: #remove 220
		nuarr = [90e9,150e9]
		noisearr = [8,7]
		beamarr = [2.2,1.4]
elif which_exp == 'planck':
	nuarr = [90e9,150e9, 220e9]
	noisearr = [40,26,100]
	beamarr = [10.,7.1,5.5]
elif which_exp == 'sptsz':
	nuarr = [90e9,150e9]
	noisearr = [40.,18.]
	beamarr = [1.7,1.2]
elif which_exp == 'SAdeep':
	nuarr = [90e9,150e9, 220e9]
	noisearr = [2.0, 2.0, 4.7]
	beamarr = [5.2,3.5,2.]
elif which_exp == 'SAwide':
	nuarr = [90e9,150e9, 220e9]
	noisearr = [5.5, 5.5, 20.0]
	beamarr = [5.2,3.5,2.]
	if pol: #remove 220
		nuarr = [90e9,150e9]
		noisearr = [5.5, 5.5]
		beamarr = [5.2,3.5]
elif which_exp == 'CMB-S4':
	nuarr = [90e9,150e9,220e9]
	#noisearr = [1.5, 1.5,2.5]
	noisearr = [1., 1., 1.]
	beamarr = [3.5,2.0,1.4]
elif which_exp == 'CMB-S4_2019':
	nuarr = [90e9,150e9,220e9]
	noisearr = [1.73, 1.73, 4.8]
	beamarr = [2.2,1.4,1.0]
elif which_exp == 'CMB-S4_sptlike_beam':
	nuarr = [90e9,150e9,220e9]
	#noisearr = [1.5, 1.5,2.5]
	noisearr = [1., 1., 1.]
	beamarr = [1.7,1.2,1.]
elif which_exp == 'CMB-S4_SOlike_beam':
	nuarr = [90e9,150e9,220e9]
	#noisearr = [1.5, 1.5,2.5]
	noisearr = [1., 1., 1.]
	beamarr = [5.2,3.5,2.]
elif which_exp == 'SO_SOlike_beam':
	nuarr = [90e9,150e9,220e9]
	#noisearr = [1.5, 1.5,2.5]
	noisearr = [4.5, 2.5, 4.5]
	beamarr = [5.2,3.5,2.]
elif which_exp == 'SO_sptlike_beam':
	nuarr = [90e9,150e9,220e9]
	noisearr = [4.5, 2.5, 4.5]
	beamarr = [1.7,1.2,1.]
elif which_exp == 'SO':
	nuarr = [90e9,150e9,220e9]
	noisearr = [4.5, 2.5, 4.5]
	beamarr = [3.5,2.0,1.4]
elif which_exp == 'CCAT_prime':
	nuarr = [90e9,150e9,220e9]
	noisearr = [6.5, 6.5, 6.5]
	beamarr = [2.2,1.4,1.]
elif which_exp == 'CCAT_prime_v2':
	nuarr = [90e9,150e9,220e9]
	noisearr = [1.,1.,2.]
	beamarr = [2.2,1.4,1.]
elif which_exp == 'PB_Anh':
	nuarr = [90e9,150e9]
	noisearr = [10., 10.]
	beamarr = [5.0,3.5]
elif which_exp == 'dummy':
	nuarr = [90e9, 150e9]
	noisearr = [55., 55.]
	beamarr = [5.0, 5.0]
elif which_exp == 'SO2019_goal':
	nuarr = [90e9,150e9,220e9]
	noisearr = [5.8, 6.3, 15.0]
	beamarr = [2.2,1.4,1.0]
elif which_exp == 'advact_for_so_delensing': 
	#scaled based on https://arxiv.org/pdf/1510.02809.pdf
	nuarr = [90e9,150e9,220e9]
	act_fsky, so_fsky = 0.4, 0.1
	scale_fsky_from_act  = ( act_fsky / so_fsky )**0.5
	noisearr = np.asarray( [8., 7., 25.] ) / scale_fsky_from_act
	beamarr = [2.2,1.4,0.9]
elif which_exp == 'S4_may_2019':
	nuarr = [90e9,150e9,220e9]
	noisearr = np.asarray( [1.5, 1.5, 4.8] ) * s4_noise_scale_factor
	beamarr = np.asarray( [2.2,1.4,1.0] ) * s4_beam_scale_factor
elif which_exp == 'S4':
    nuarr = [90e9,150e9,220e9]
    noisearr = [2.0, 2.0, 6.9]
    beamarr = [2.3, 1.5, 1.0]    

if pol: noisearr = np.asarray(noisearr) * 1.414
pol_frac_per_cent = 0.16/100. #4%
#print noisearr;sys.exit()

if yearscaling <> -1.:
	noisearr = np.array(noisearr) * yearscaling

nuarr = np.asarray( nuarr )

"""
Async = map(lambda x: fn_get_A_syn(x), nuarr)
Adust = map(lambda x: fn_get_A_dust(x), nuarr)
AtSZ = map(lambda x: fn_get_A_tSZ(x), nuarr)
Acmb = np.zeros(len(nuarr)) + 1.

A = np.column_stack( (Acmb, Async, AtSZ, Adust) )

print Async
print Adust
print AtSZ

nx = ny = 50
dx = dy = 0.5
simmapparams = [nx, ny, dx, dy]
NOISE =  sims.fn_get_noise(simmapparams, noisearr = noisearr)
BEAMS = map(lambda x: sims.fn_beam_stuff(simmapparams, use_beam =1, return_beam = 1, exp_beam = x)[0], beamarr)
for n in range(len(NOISE)):
	subplot(2,3,n+1);imshow(NOISE[n]);colorbar()
	NOISE[n] = np.fft.ifft2( np.fft.fft2(NOISE[n]) * BEAMS[n] ). real
	subplot(2,3,n+1+len(NOISE));imshow(NOISE[n]);colorbar()
show()


for i in range(nx):
	for j in range(ny):
		N = np.mat(NOISE[:,j,i])
		Ninv = np.linalg.inv(N)
		
		print N;sys.exit()

"""
############################################################################################################
############################################################################################################
############################################################################################################
#step 1: Create CMB maps for 3 different freqeuncies and add foregrounds, noise

boxsize = 600
dx = dy = 0.5
nx = int(boxsize/dx); ny = int(boxsize/dy)

simmapparams = [nx, ny, dx, dy]
log_file = 'tmp/logs_fg_removal.txt'
sims._ini_log_file(log_file)


# read and store param dict
paramfile = 'params/params_planck_r_0.0_2015_cosmo_lensed_LSS.txt'
params = np.recfromtxt(paramfile,usecols=[0],delimiter = '=')
paramvals = np.recfromtxt(paramfile,usecols=[1],delimiter = '=')
param_dict = {}
for p,pval in zip(params,paramvals):
        tmp = pval.strip()
        try:
                float(tmp)
                if tmp.find('.')>-1:
                        param_dict[p.strip()] = float(tmp)
                else:
                        param_dict[p.strip()] = int(tmp)
        except:
                if tmp == 'None':
                        param_dict[p.strip()] = None
                elif tmp[0] == '[':
                        #param_dict[p.strip()] = np.asarray( map(tmp[1:-1].split(',')) )
                        dummy = tmp[1:-1].split(',')[0]
                        try:
                                param_dict[p.strip()] = np.asarray( map(float, tmp[1:-1].split(',')) )
                        except:                         
                                param_dict[p.strip()] = tmp[1:-1].split(',')
                                for vcnt in range(len(param_dict[p.strip()])):
                                        param_dict[p.strip()][vcnt] = param_dict[p.strip()][vcnt].strip()
                else:
                        param_dict[p.strip()] = tmp.replace('\'','')

add_fg = 1
#if pol: add_fg  = 0
add_noise = 1
run_s2_dict = {}
if add_fg:
	run_s2_dict['fg_T'] = 'all_Gaussian'
sims.run_s2_dict = run_s2_dict
############################################
############################################
# read CAMB Dls
Dlfile_len = param_dict['Dlfile_len']
Dls_len = np.loadtxt(Dlfile_len,usecols=[0,1])
tqulen = Dls_len.shape[1] - 1
sims.tqulen = tqulen

########################################################################################################################
########################################################################################################################
########################################################################################################################

#initialize CMB sim stuffs
import time
show_plots = 0

"""
start = time.time()
CMB_MAPS = np.zeros( (len(nuarr), nx, ny) )
random_seed = 12345
Cls_dic = {}
for ncnt,(nuval, beamval, noiseval) in enumerate(zip(nuarr,beamarr,noisearr)):

	np.random.seed(random_seed) #all frequencies should have the same seed
	nuval = nuval/1e9
	
	param_dict['exp_beam'] = beamval
	param_dict['use_beam'] = 1
	sims._ini_params(param_dict) #just move this dic also to modules/scl_cmb.py for efficient usage
	sims.exp_beam = param_dict['exp_beam'] #analytic beam

	sims.Dls2map(Dls_len, simmapparams, CMB_outputscale = param_dict['CMB_outputscale'])
	sims.fn_beam_stuff(simmapparams)

	#loglog(Dls_len[:,0],Dls_len[:,1]);show()

	#create CMB map and add the respective frequency foreground using George et al. 2015 measurements
	SIMMAP = sims.fn_perform_lensing_nfw(simmapparams, 1, param_dict, use_beam = param_dict['use_beam'], no_lensing = 1, nu = nuval)

	#add noise and apodize map
	np.random.seed(int(nuval))
	SIMMAP = sims.fn_tf_noise(SIMMAP, simmapparams, 1, add_TF = param_dict['add_TF'], add_noise = add_noise, expnoiselevel = [noiseval], perform_apod = 0)
	CMB_MAPS[ncnt] = SIMMAP[0,0]

	#get auto power spectra
	Cls = sims.fn_plot_pow_spec(simmapparams,np.asarray( [CMB_MAPS[ncnt]] ))[0][0]
	Cls_dic[nuval] = Cls

	if show_plots:
		subplot(1,3,ncnt+1);imshow(CMB_MAPS[ncnt]);colorbar()

if show_plots: show()#;quit()

if show_plots: 
	clf();ax = subplot(111, xscale = 'log', yscale = 'log')
	for nuval in Cls_dic:
		Cls = Cls_dic[nuval]
		els = Cls[:,0]
		Dls_fac = els * (els+1) / 2 / np.pi
		plot(els, Cls[:,1] * Dls_fac, label = nuval)
	plot(Dls_len[:,0], Dls_len[:,1] * 1e12 * sims.Tcmb**2., 'k')
	legend(loc=1, fancybox=1);xlim(100,2e4);ylim(1e-2,1e4);show()

end = time.time()
print 'Time for simulated skies = %.1f' %((end-start))

#NOISE =  sims.fn_get_noise(simmapparams, noisearr = noisearr)
#BEAMS = map(lambda x: sims.fn_beam_stuff(simmapparams, use_beam =1, return_beam = 1, exp_beam = x)[0], beamarr)

############################################################################################################
############################################################################################################
############################################################################################################
#step2: calcualte auto and cross spectra between maps. Smooth all maps by the desired beam
start = time.time()
Cls_dic = {}
nuarr = nuarr/1e9
if show_plots:clf()


op_beam = beamarr[1]
maxinvBl_val = 1000.
op_Bl = sims.fn_beam_stuff(simmapparams, use_beam = 1, exp_beam = op_beam, return_beam = 1)

def fn_smooth_with_new_beam(MAP, expbeamval1, expbeamval2 = op_beam):

	Bl1 = sims.fn_beam_stuff(simmapparams, use_beam = 1, exp_beam = expbeamval1, return_beam = 1)
	if expbeamval1 > expbeamval2:
		Bl1[1./Bl1>maxinvBl_val] = 1./maxinvBl_val
	
	if expbeamval2 == op_beam:
		Bl2 = op_Bl
	else:
		Bl2 = sims.fn_beam_stuff(simmapparams, use_beam = 1, exp_beam = expbeamval2, return_beam = 1)
	
	mod_Bl = Bl2/Bl1
	mod_Bl[ mod_Bl == np.inf] = 0.
	mod_Bl[np.where(np.isnan(mod_Bl))] = 0.
	
	#clf();imshow(mod_Bl);colorbar();title('%s, %s' %(expbeamval1, expbeamval2));show()

	return np.fft.ifft2( np.fft.fft2(MAP) * mod_Bl ).real



for mcnt1, (MAP1, nuval1, beamval1) in enumerate(zip(CMB_MAPS, nuarr, beamarr)):
	if show_plots:ax = subplot(1,3,mcnt1+1, xscale = 'log', yscale = 'log')

	#subplot(121);imshow(MAP1);colorbar()
	MAP1 = fn_smooth_with_new_beam(MAP1, beamval1)
	#subplot(122);imshow(MAP1);colorbar();title(beamval1);show()#;sys.exit()
	
	for mcnt2, (MAP2, nuval2, beamval2) in enumerate(zip(CMB_MAPS, nuarr, beamarr)):

		MAP2 = fn_smooth_with_new_beam(MAP2, beamval2)

		cls = sims.fn_get_pow_spec(simmapparams, [MAP1,MAP2], binsize = None, K_to_uK = 1)
		Cls_dic[(nuval1, nuval2)] = cls
		els = cls[:,0]
		Dls_fac = els * (els+1) / 2 / np.pi
		if show_plots:plot(els, cls[:,1] * Dls_fac, label = '%s x %s' %(nuval1, nuval2))

	if show_plots:
		plot(Dls_len[:,0], Dls_len[:,1] * 1e12 * sims.Tcmb**2., 'k')
		legend(loc=1, fancybox=1);xlim(100,2e4)#;ylim(1,1e4)
if show_plots: show()
end = time.time()
print 'Time for power spectrum calculation = %.1f' %((end-start))
"""
############################################################################################################
############################################################################################################
############################################################################################################
#step2: get auto and cross spectra from George et al. 2015 results. Smooth all maps by the desired beam
start = time.time()
Cls_dic = {}
nuarr = nuarr/1e9
if show_plots:clf()

op_beam = beamarr[1]

def fn_get_Bl(beamval, els):

	fwhm_radians = np.radians(beamval/60.)
	sigma = fwhm_radians / np.sqrt(8. * np.log(2.))
	sigma2 = sigma ** 2
	Bl = np.exp(els * (els+1) * sigma2)

	return Bl

def fn_get_Nl(noiseval, els, beamval = op_beam, use_beam_window = 1, uk_to_K = 0): #here we will smooth noise by the output required beam

	if uk_to_K: noiseval = noiseval/1e6

	if use_beam_window:
		fwhm_radians = np.radians(beamval/60.)
		sigma = fwhm_radians / np.sqrt(8. * np.log(2.))
		sigma2 = sigma ** 2
		Bl = np.exp(els * (els+1) * sigma2)

	delta_T_radians = noiseval * np.radians(1./60.)
	Nl = np.tile(delta_T_radians**2., int(max(els)) + 1 )

	Nl = np.asarray( [Nl[int(el)] for el in els] )

	if use_beam_window: Nl *= Bl

	return Nl

if which_exp.find('spt')>-1:
	use_lknee = 1
else:
	use_lknee = 0
if use_lknee:	
	if not pol:
		elkneedic = {}
		elkneedic[150] = 1200.
		elkneedic[90] = 2200.
		elkneedic[220] = 2300.
	else:
		elkneedic = {}
		elkneedic[150] = 300.
		elkneedic[90] = 300.
		elkneedic[220] = 300.

show_plots = 0
for mcnt1, (nuval1, beamval1, noiseval1) in enumerate(zip(nuarr, beamarr, noisearr)):
	if show_plots:ax = subplot(1,3,mcnt1+1, yscale = 'log')#, xscale = 'log')

	for mcnt2, (nuval2, beamval2, noiseval2) in enumerate(zip(nuarr, beamarr, noisearr)):

		#print nuval1, nuval2
		if (nuval2, nuval1) in Cls_dic:
			Cls_dic[(nuval1, nuval2)] = Cls_dic[(nuval2, nuval1)]
			if show_plots:plot(els, Cls_dic[(nuval2, nuval1)] * Dls_fac, label = '%s x %s' %(nuval1, nuval2))
			continue

		if not pol:
			els, Cls = sims.fn_get_foreground_power('Total', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			'''
			#loglog(els, Cls_full, 'k')
			els, Cls = sims.fn_get_foreground_power('CMB', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			els, Cls_DGPo = sims.fn_get_foreground_power('DG-Po', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			els, Cls_DGclus = sims.fn_get_foreground_power('DG-Cl', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			els, Cls_RG = sims.fn_get_foreground_power('RG', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			els, Cls_tSZ = sims.fn_get_foreground_power('tSZ', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			els, Cls_kSZ = sims.fn_get_foreground_power('kSZ', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			els, Cls_tSZxcib = sims.fn_get_foreground_power('tSZ-CIB', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			els, Cls_others = sims.fn_get_foreground_power('cirrus, tSZ-RG, and other', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			Cls_FG_extragal = Cls_DGPo + Cls_DGclus + Cls_RG + Cls_tSZ + Cls_kSZ + Cls_others# + Cls_tSZxcib
			Cls = Cls + Cls_FG_extragal
			'''

			'''
			## els, Cls_kSZ = sims.fn_get_foreground_power('kSZ', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			## els, Cls_tSZxcib = sims.fn_get_foreground_power('tSZ-CIB', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			Cls_FG_extragal = Cls_DGPo + Cls_DGclus + Cls_RG## + Cls_tSZ + Cls_kSZ# + Cls_tSZxcib
			els, Cls_others = sims.fn_get_foreground_power('cirrus, tSZ-RG, and other', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			Cls = Cls + Cls_FG_extragal##  + Cls_others
			loglog(els, Cls, 'r--')
			show();sys.exit()
			'''

		else:
			els, Cls_DGPo = sims.fn_get_foreground_power('DG-Po', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			els, Cls_DGclus = sims.fn_get_foreground_power('DG-Cl', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			els, Cls_RG = sims.fn_get_foreground_power('RG', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)

			'''
			Cls_FG_extragal = Cls_DGPo + Cls_DGclus + Cls_RG
			Cls_FG_extragal = Cls_FG_extragal * pol_frac_per_cent

			print Cls_FG_extragal
			'''

			#3G sim numbers
			dg_pol_fraction = 0.035  # Planck: Bonavera et al. arXiv: 1705.10603
			rg_pol_fraction = 0.028  # ACTPol: arXiv: 1811.01854; SPTpol: Gupta et al. arXiv: 1907.02156

			Cls_FG_extragal = (Cls_DGPo + Cls_DGclus) * (dg_pol_fraction**2.) + (Cls_RG * rg_pol_fraction**2.)

			#print Cls_FG_extragal; sys.exit()

			els, Cls = sims.fn_get_foreground_power('CMB', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			Cls += Cls_FG_extragal

		els, Cls_CMB = sims.fn_get_foreground_power('CMB', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
		Cls_CMB *= 1e12

		Cls *= 1e12

		'''
		clf();
		nuarr = [90, 150, 220]
		for nucnt, nu in enumerate(nuarr):
			subplot(1,3,nucnt+1)
			arr = ['Total', 'CMB', 'DG-Po', 'DG-Cl', 'RG']
			for acnt, a in enumerate(arr):
				els, Cls = sims.fn_get_foreground_power(a, simmapparams, 1, just_Cls = 1, nu1 = nu, nu2 = nu)
				loglog(els, Cls, label  =a)

			els, Cls_DGPo = sims.fn_get_foreground_power('DG-Po', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			els, Cls_DGclus = sims.fn_get_foreground_power('DG-Cl', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)
			els, Cls_RG = sims.fn_get_foreground_power('RG', simmapparams, 1, just_Cls = 1, nu1 = nuval1, nu2 = nuval2)

			Cls_FG_extragal = Cls_DGPo + Cls_DGclus + Cls_RG
			loglog(els, Cls_FG_extragal, label  = 'EG')
			legend(loc = 1)
			title('george\_plot\_bestfit\_line.sav')
		show();### sys.exit()
		'''

		try:
			op_Bl
		except:
			op_Bl = fn_get_Bl(op_beam, els)

		Bl = fn_get_Bl(beamval2, els)

		#multiply by the required output beam
		#Cls = Cls * op_Bl**2.

		'''
		mod_Bl = op_Bl/Bl
		mod_Bl[mod_Bl == np.inf] = 0.
		mod_Bl[np.where(np.isnan(mod_Bl))] = 0.
		Cls = Cls * mod_Bl**2.
		clf();loglog(els, mod_Bl);show()#;quit()
		'''

		#if autospectrum: add noise convolved with the final output beam --> becuase we will filter all maps with the final beam before ILC
		if nuval1 == nuval2: 
			Nl = fn_get_Nl(noiseval2, els)
		else:
			Nl = np.zeros(len(Cls))

		Nl_futuristic = fn_get_Nl(1.7 * 1.414, els)

		#loglog(els, Nl)
		if use_lknee:
			if nuval1 == nuval2: 
				elknee = elkneedic[nuval1]
				Nl = Nl * (1. + (elknee/els) )
				Nl_futuristic = Nl_futuristic * (1. + (elknee/els) )
		#loglog(els, Nl);show();sys.exit()
		Cls = Cls + Nl

		if nuval1 == 150 and nuval2 == 150: 
			Nl150 = Nl


		Cls_dic[(nuval1, nuval2)] = Cls

		Dls_fac = els * (els+1) / 2 / np.pi
		if show_plots:plot(els, Cls * Dls_fac, label = '%s x %s' %(nuval1, nuval2));

	if show_plots:
		plot(Dls_len[:,0], Dls_len[:,1] * 1e12 * sims.Tcmb**2., 'k')
		legend(loc=1, fancybox=1);xlim(1,12e3);ylim(1.,1e4);#xlim(100,2e4)#;ylim(1,1e4)
if show_plots: show()#;quit()
end = time.time()
print 'Time for power spectrum calculation = %.1f' %((end-start))

############################################################################################################
############################################################################################################
############################################################################################################
#step3: get output power spectrum
nc = len(nuarr)
acap = np.zeros(nc) + 1. #assuming CMB is the same and calibrations factors are same for all channels
acap = np.mat(acap).T #should be nc x 1

def fn_create_Clmat(elcnt):
	Clmat = np.zeros( (nc, nc) )
	for ncnt1, nuval1 in enumerate(nuarr):
		for ncnt2, nuval2 in enumerate(nuarr):
			#Clmat[ncnt2, ncnt1] = Cls_dic[(nuval1, nuval2)][elcnt,1]
			Clmat[ncnt2, ncnt1] = Cls_dic[(nuval1, nuval2)][elcnt]
	return Clmat

CLEANED_POW_SPEC = np.zeros( (len(els)) )
for elcnt, el in enumerate(els):
	Clmat = np.mat( fn_create_Clmat(elcnt) )
	Clinv = sc.linalg.pinv2(Clmat)
	
	nr = 1.
	dr = np.dot( acap.T, np.dot(Clinv, acap) )
	final_pow = np.asarray(nr/dr).squeeze()

	CLEANED_POW_SPEC[elcnt] = final_pow

#RES_POW_SPEC = (Dls_len[:,1]*1e12) - CLEANED_POW_SPEC
RES_POW_SPEC = CLEANED_POW_SPEC - Cls_CMB

'''
#1/f component
elknee, alpha = 100., 1.
#RES_POW_SPEC_with_oneoverf = RES_POW_SPEC_with_oneoverf * (1 + (elknee/els)**2.)

l1, l2 = 400, 500
indsforfitting = np.where((els>=l1) & (els<=l2))
poly_fit_eq=np.polyfit(els[indsforfitting],RES_POW_SPEC[indsforfitting],3)
poly_fit=np.poly1d(poly_fit_eq)
indsforfitting = np.where((els<=l1))
RES_POW_SPEC[indsforfitting]=poly_fit(els[indsforfitting])

RES_POW_SPEC_with_oneoverf = np.copy(RES_POW_SPEC)
RES_POW_SPEC_with_oneoverf[:] = np.median(RES_POW_SPEC_with_oneoverf)
RES_POW_SPEC_with_oneoverf = RES_POW_SPEC_with_oneoverf * (1 + (elknee/els)**alpha)
'''

Dls_fac[:] = 1.

print noisearr, beamarr

if which_exp == 'CMB-HD' and pol:
	poly_deg = 1
	poly_fit_eq=np.polyfit(els, RES_POW_SPEC, poly_deg)
	poly_fit=np.poly1d(poly_fit_eq)
	RES_POW_SPEC = poly_fit(els)

if 1==1:#show_plots:
	clf(); ax = subplot(111, xscale = 'log', yscale = 'log')
	#plot(Dls_len[:,0], Dls_len[:,1] * 1e12 * sims.Tcmb**2., 'k', lw = 3.)
	plot(els, Cls_CMB  * Dls_fac, 'k', lw = 3.)
	plot(els, CLEANED_POW_SPEC * Dls_fac, 'gold', label = 'FG cleaned')
	plot(els, RES_POW_SPEC * Dls_fac, 'lime', label = 'Residual')
	if year <> -1.:
		plot(els, RES_POW_SPEC * Dls_fac / yearscaling, 'm', label = 'scaled from full')
	### plot(els, RES_POW_SPEC_with_oneoverf * Dls_fac, 'lightskyblue', label = 'Residual with 1/f')
	plot(els, Nl150 * Dls_fac, 'r', label = '150 GHz - White')
	plot(els, Nl_futuristic * Dls_fac, 'm', label = 'Futuristic - White')
	legend(loc=1, fancybox=1);xlim(1,12e3);#ylim(1e-3,1e4);#xlim(100,2e4)#;ylim(1,1e4)
	title('EXperiment = %s' %(which_exp.replace('_','\_')))
	show()
end = time.time()
print 'Time for get final power spectrum = %.1f' %((end-start))

import pickle, gzip
Nl_dic = {}
Nl_dic['Cl_CMB'] = [els, Cls_CMB]
Nl_dic['Cl_CMB_cleaned'] = [els, CLEANED_POW_SPEC]
Nl_dic['Cl_residual'] = [els, RES_POW_SPEC]
Nl_dic['channels'] = nuarr
Nl_dic['beams'] = beamarr
Nl_dic['noiselevels'] = noisearr
Nl_dic['beam_effective_arcmins'] = op_beam
Nl_dic['Bl_effective'] = op_Bl



if not pol:
	opf = 'data/foreground/cleaning/Cls_ILC_%s.pkl.gz' %(which_exp)
else:
	opf = 'data/foreground/cleaning/Cls_ILC_%s_pol.pkl.gz' %(which_exp)


if use_lknee:
	Nl_dic['lknee'] = elkneedic
	print elkneedic

	opf = opf.replace('.pkl.gz', 'with_1_f_noise.pkl.gz')	

if year<>-1.:
	opf = opf.replace('.pkl.gz', '_year_%s.pkl.gz' %(year))

if which_exp == 'S4_may_2019':
	opf = 'data/foreground/cleaning/%s/Cls_ILC.pkl.gz' %(which_exp)
	opf = opf.replace('.pkl.gz', '_beamx%s_noisex%s.pkl.gz' %(s4_beam_scale_factor, s4_noise_scale_factor))

	print opf


pickle.dump(Nl_dic, gzip.open(opf, 'wb'), protocol = 2)

quit()

############################################################################################################
############################################################################################################
############################################################################################################
#step3: get SMICA weights
"""
arXiv: 1303.5072 - appendix D; Eq. D2
"""
start = time.time()
nc = len(nuarr)
acap = np.zeros(nc) + 1. #assuming CMB is the same and calibrations factors are same for all channels
acap = np.mat(acap).T #should be nc x 1

def fn_create_Clmat(elcnt):
	Clmat = np.zeros( (nc, nc) )
	for ncnt1, nuval1 in enumerate(nuarr):
		for ncnt2, nuval2 in enumerate(nuarr):
			#Clmat[ncnt2, ncnt1] = Cls_dic[(nuval1, nuval2)][elcnt,1]
			Clmat[ncnt2, ncnt1] = Cls_dic[(nuval1, nuval2)][elcnt]
	return Clmat

W = np.zeros( (nc, len(els)) )
for elcnt, el in enumerate(els):
	Clmat = np.mat( fn_create_Clmat(elcnt) )
	Clinv = sc.linalg.pinv2(Clmat)
	
	nr = np.dot(Clinv, acap)
	dr = np.dot( acap.T, np.dot(Clinv, acap) )
	weight_mat = np.asarray(nr/dr).squeeze()

	W[:, elcnt] = weight_mat

end = time.time()
print 'Time for get weights = %.1f' %((end-start))

ax = subplot(111)
for ncnt, nuval in enumerate(nuarr):
	plot(els, W[ncnt], label = nuval)		

#ylim(-3.,3.)
legend(loc = 1, fancybox = 1)
show();quit()
############################################################################################################
############################################################################################################
############################################################################################################
#step4: get Nl for each channel and the final combination
maxel = int( max(els) + 1 )
ax = subplot(111, xscale = 'log', yscale = 'log')
use_beam_window = 1
colorarr =['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728','#9467bd', '#8c564b', '#e377c2', '#7f7f7f','#bcbd22', '#17becf']
Nl_tot = np.zeros( len(els) )
for ncnt, (nuval, noiseval, beamval) in enumerate(zip(nuarr, noisearr, beamarr)):

	if use_beam_window:
		fwhm_radians = np.radians(beamval/60.)
		sigma = fwhm_radians / np.sqrt(8. * np.log(2.))
		sigma2 = sigma ** 2
		Bl = np.exp(els * (els+1) * sigma2)

	delta_T_radians = noiseval * np.radians(1./60.)
	Nl = np.tile(delta_T_radians**2., maxel)

	Nl = np.asarray( [Nl[int(el)] for el in els] )


	if use_beam_window: Nl *= Bl
	Nl_weighted = W[ncnt]**2. * Nl

	#Nl_tot += (1./Nl_weighted)
	Nl_tot += Nl_weighted

	Dls_fac = els * (els+1) / 2 / np.pi
	plot(els, Nl * Dls_fac, label = nuval, color = colorarr[ncnt])
	plot(els, Nl_weighted * Dls_fac, color = colorarr[ncnt], marker = ',', ls = 'None')

#Nl_tot = 1./Nl_tot
import pickle, gzip
Nl_dic = {}
Nl_dic['Nl'] = [els, Nl_tot]
Nl_dic['channels'] = nuarr
Nl_dic['beams'] = beamarr
Nl_dic['noiselevels'] = noisearr
Nl_dic['beam_effective_arcmins'] = op_beam
Nl_dic['Bl_effective'] = op_Bl

"""
if not pol:
	pickle.dump(Nl_dic, gzip.open('data/foreground/cleaning/Nl_%s.pkl.gz' %(which_exp), 'wb'), protocol = 2)
else:
	pickle.dump(Nl_dic, gzip.open('data/foreground/cleaning/Nl_%s_pol_with_extragal_FG.pkl.gz' %(which_exp), 'wb'), protocol = 2)
"""
plot(Dls_len[:,0], Dls_len[:,1] * 1e12 * sims.Tcmb**2., 'k')
plot(els, Nl_tot * Dls_fac, color='darkred', lw = 3.)
xlim(500,1e4)
ylim(0.01,1e5)
title('Experiment = %s' %(which_exp))
legend(loc = 1, fancybox = 1)
show()

	





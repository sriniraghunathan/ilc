

import argparse, sys, numpy as np, scipy as sc
sys.path.append('/Users/sraghunathan/Research/SPTPol/analysis/git/tools/')
import flatsky, tools, misc
import ilc, foregrounds as fg


h=6.62607004e-34 #Planck constant in m2 kg / s
k_B=1.38064852e-23 #Boltzmann constant in m2 kg s-2 / K-1
Tcmb = 2.73 #Kelvin

parser = argparse.ArgumentParser(description='')
parser.add_argument('-paramfile', dest='paramfile', action='store', help='paramfile', type=str, default = 'params.ini')
parser.add_argument('-freqarr', dest='freqarr', action='store', type=int, nargs='+', default= [90, 150, 220, 270], help='freqarr')
parser.add_argument('-beamarr', dest='beamarr', action='store', type=float, nargs='+', default= [2.3, 1.5, 1.0, 0.8], help='freqarr')
parser.add_argument('-noisearr', dest='noisearr', action='store', type=float, nargs='+', default= [2.0, 2.0, 6.9, 16.7], help='noisearr')
parser.add_argument('-elkneearr', dest='elkneearr', action='store', type=float, nargs='+', default= [2154., 4364., 7334., 7308.], help='elkneearr')
parser.add_argument('-alphakneearr', dest='alphakneearr', action='store', type=float, nargs='+', default= [3.5, 3.5, 3.5, 3.5], help='alphakneearr')

args = parser.parse_args()
args_keys = args.__dict__
for kargs in args_keys:
    param_value = args_keys[kargs]
    if isinstance(param_value, str):
        cmd = '%s = "%s"' %(kargs, param_value)
    else:
        cmd = '%s = %s' %(kargs, param_value)
    exec(cmd)

freqarr = np.asarray( freqarr )
############################################################################################################
# read and store param dict
param_dict = misc.fn_get_param_dict(paramfile)
el = np.arange(param_dict['lmax'])

############################################################################################################
#get beam deconvolved noise nls

beam_noise_dic = {}
for (freq, beam, noise) in zip(freqarr, beamarr, noisearr):
    beam_noise_dic[freq] = [beam, noise]

elknee_dic = {}
for (freq, elknee, alphaknee) in zip(freqarr, elkneearr, alphakneearr):
    elknee_dic[freq] = [elknee, alphaknee]

nl_dic = {}
for freq in freqarr:
    beamval, noiseval = beam_noise_dic[freq]
    nl = misc.get_nl(noiseval, el, beamval, elknee_t=elknee_dic[freq][0], alpha_knee=elknee_dic[freq][1])
    nl[el<=param_dict['lmin']] = 0.
    nl[nl == 0.] = np.min(nl[nl!=0.])/1e3
    nl_dic[freq] = nl
    
print(nl_dic.keys())

############################################################################################################

#get the CMB, noise, and foreground covriance
try:
    ignore_fg = param_dict['ignore_fg']
except:
    ignore_fg = []
print(ignore_fg)
el, cl_dic = ilc.get_covariance_dic(param_dict, freqarr, nl_dic = nl_dic, ignore_fg = ignore_fg)
print(el)
sys.exit()
############################################################################################################

nc = len(freqarr)
acap = np.zeros(nc) + 1. #assuming CMB is the same and calibrations factors are same for all channels
acap = np.mat(acap).T #should be nc x 1

def fn_create_Clmat(elcnt):
	Clmat = np.zeros( (nc, nc) )
	for ncnt1, nuval1 in enumerate(freqarr):
		for ncnt2, nuval2 in enumerate(freqarr):
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
Nl_dic['channels'] = freqarr
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
nc = len(freqarr)
acap = np.zeros(nc) + 1. #assuming CMB is the same and calibrations factors are same for all channels
acap = np.mat(acap).T #should be nc x 1

def fn_create_Clmat(elcnt):
	Clmat = np.zeros( (nc, nc) )
	for ncnt1, nuval1 in enumerate(freqarr):
		for ncnt2, nuval2 in enumerate(freqarr):
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
for ncnt, nuval in enumerate(freqarr):
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
for ncnt, (nuval, noiseval, beamval) in enumerate(zip(freqarr, noisearr, beamarr)):

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
Nl_dic['channels'] = freqarr
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

	





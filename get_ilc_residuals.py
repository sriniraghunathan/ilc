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
#collect beam and noise into a dic
beam_noise_dic = {}
for (freq, beam, noise) in zip(freqarr, beamarr, noisearr):
    beam_noise_dic[freq] = [beam, noise]

############################################################################################################
#collect elknee and alpha into a dic
elknee_dic = {}
for (freq, elknee, alphaknee) in zip(freqarr, elkneearr, alphakneearr):
    elknee_dic[freq] = [elknee, alphaknee]
############################################################################################################
#get beam deconvolved noise nls
nl_dic = {}
for freq in freqarr:
    beamval, noiseval = beam_noise_dic[freq]
    nl = misc.get_nl(noiseval, el, beamval, elknee_t=elknee_dic[freq][0], alpha_knee=elknee_dic[freq][1])
    nl[el<=param_dict['lmin']] = 0.
    nl[nl == 0.] = np.min(nl[nl!=0.])/1e3
    nl_dic[freq] = nl
############################################################################################################
#get the CMB, noise, and foreground covriance
try:
    ignore_fg = param_dict['ignore_fg']
except:
    ignore_fg = []
print(ignore_fg)
el, cl_dic = ilc.get_covariance_dic(param_dict, freqarr, nl_dic = nl_dic, ignore_fg = ignore_fg)
############################################################################################################
#get the residual power now
cl_residual = ilc.residual_power(param_dict, freqarr, el, cl_dic, final_comp = 'CMB')#, freqcalib_fac = None)
print(cl_residual)
############################################################################################################

if 1==1:#show_plots:
    from pylab import *
    freq0, lmax = param_dict['freq0'], param_dict['lmax']
    el_, cl_cmb = fg.get_foreground_power_george_2015('CMB', freq1 = freq0, lmax = lmax)
    foregrounds_to_plot = ['kSZ', 'tSZ', 'DG-Po', 'DG-Cl', 'RG']

    clf(); 
    ax = subplot(111, xscale = 'log', yscale = 'log')
    plot(el, cl_cmb, 'gray', lw = 1., label = r'CMB')
    plot(el, cl_residual, 'black', lw = 2., label = r'Residual')
    for curr_fg in foregrounds_to_plot:
        el_, cl_curr_fg = fg.get_foreground_power_george_2015(curr_fg, freq1 = freq0, lmax = lmax)
        plot(el, cl_curr_fg, lw = 0.5, ls = '--', label = r'%s' %(curr_fg), alpha = 1.)
    for freq in freqarr:
        plot(el, nl_dic[freq], lw = 0.5, ls = '-', label = r'Noise: %s' %(freq))#, alpha = 0.5)
    legend(loc=1, fancybox=1, ncol = 4, fontsize = 8);
    xlim(10,1e4);#ylim(1e-3,1e4);
    xlabel(r'Multipole $\ell$')
    ylabel(r'$C_{\ell}$')
    show()
    





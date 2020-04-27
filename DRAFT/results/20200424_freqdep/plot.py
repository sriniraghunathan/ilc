import numpy as np, glob, sys, os
from pylab  import *
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
#mpl.rcParams['hatch.linewidth'] = 0.5
rcParams['font.family'] = 'serif'
rc('text.latex',preamble=r'\usepackage{/Users/sraghunathan/.configs/apjfonts}')

flist = glob.glob('*.npy')

colordic = {
    '93-145-225-278': 'black',
    '93-145-225': 'royalblue',
    '93-145-278': 'darkgreen',
    '93-145': 'darkred',
}
ax = subplot(111, yscale = 'log')#, xscale = 'log')
for fname in flist:
    dic = np.load(fname, allow_pickle = 1).item()
    el = dic['el']
    cl_res_TT = dic['cl_residual']['T']

    freqstr = str( fname.split('_')[1].strip() )

    colorval = colordic[freqstr]

    print(freqstr, colorval)
    
    plot(el, cl_res_TT, label = r'%s' %(freqstr), color = colorval)

    #print(fname)
legend(loc = 1, fancybox = 1, fontsize = 12)
title(r'S4-wide ILC curves: TT', fontsize = 14)
xlim(20,7000);ylim(1e-8,1e6);
xlabel(r'Multipole $\ell$', fontsize = 14)
ylabel(r'$C_{\ell}\ [\mu{\rm K}^{2}]$', fontsize = 14)
show()
sys.exit()
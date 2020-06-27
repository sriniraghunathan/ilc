import numpy as np
fname = 'S4_ilc_20203030.npy'

result_dic = np.load(fname, allow_pickle = 1).item()
print(result_dic.keys())
el = result_dic['el']
cl_residual = result_dic['cl_residual']
nl_TT = cl_residual['T']
nl_PP = cl_residual['P']
print(len(el), len(nl_TT), len(nl_PP))


#plot the curves
from pylab import *
from matplotlib import rc;rc('text', usetex=True);rc('font', weight='bold');matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']
import os
rc('text.latex',preamble=r'\usepackage{/Users/sraghunathan/.configs/apjfonts}')

clf(); 
fsval = 8
lwval = 0.75
subplots_adjust(wspace=0.1, hspace = 0.1)
xmin, xmax = 20, 7000
which_spec_arr = ['T', 'P']
for cntr, which_spec in enumerate( which_spec_arr ):
    ax = subplot(1,2,cntr+1, xscale = 'log', yscale = 'log')
    plot(el, cl_residual[which_spec], 'black', lw = 2., label = r'Residual')
    xlim(xmin, xmax);
    ylim(1e-8,1e6);
    xlabel(r'Multipole $\ell$')
    if cntr == 0: 
        ylabel(r'$C_{\ell}$')
        ##legend(loc = 1, fontsize = 6, ncol = 2, handlelength = 2., handletextpad = 0.1)
    else:
        setp(ax.get_yticklabels(which = 'both'), visible=False)
    for label in ax.get_xticklabels(): label.set_fontsize(fsval)
    for label in ax.get_yticklabels(): label.set_fontsize(fsval)
    title(r'%s' %(which_spec))
show()

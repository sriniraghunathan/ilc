import numpy as np
fname = 's4_lmin10_lmax5000_lmaxtt5000.npy'

result_dic = np.load(fname, allow_pickle = 1).item()
print(result_dic.keys())

#noise_T_uKarmins
noise_T_uKarmins = result_dic['noise_T_uKarmins'] ###will be 0 if an ILC file has been used to generate these curves.

#lmin, lmax, lmax_tt: are the maximum CMB multipole used. 
#If lmax_tt is different from lmax, then that means lmax for CMB-T has been restricted to lmax_tt
lmin, lmax, lmax_tt = result_dic['lmin'], result_dic['lmax'], result_dic['lmax_tt']
#print(lmin, lmax, lmax_tt)


#lensing noise curves for different estimators
#Nl_TT, Nl_EB are for TT, EB etc. 
#Nl_MV is the minimum variance from all estimators: Note that ET == TE and has been considered only once.
#Nl_MVpol is the minimum variance combination from pol-only estimators: EE and EB.

el, cl_kk = result_dic['els'], result_dic['cl_kk']
nl_tt, nl_eb, nl_mv, nl_mvpol = result_dic['Nl_TT'], result_dic['Nl_EB'], result_dic['Nl_MV'], result_dic['Nl_MVpol']
print(len(el), len(nl_tt), len(nl_eb))
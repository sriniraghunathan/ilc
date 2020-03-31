import numpy as np
fname = 'S4_ilc_20203030.npy'

result_dic = np.load(fname, allow_pickle = 1).item()
print(result_dic.keys())
el = result_dic['el']
cl_residual = result_dic['cl_residual']
nl_TT = cl_residual['T']
nl_PP = cl_residual['P']
print(len(el), len(nl_TT), len(nl_PP))
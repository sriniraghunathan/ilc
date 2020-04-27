import numpy as np, glob, sys

nside, lmax = 2048, 3500
t_only = 0

searchstr = 'cls_galactic_sims_xxxx_CUmilta_20200319_maskplanck_nside%s_lmax%s_mask?.npy' %(nside, lmax)

for which_comp in ['dust', 'sync']:
    curr_searchstr = searchstr.replace('xxxx', which_comp)
    flist = sorted( glob.glob(curr_searchstr) )

    opdic = {}
    opdic['cl_dic'] = {}
    fksy_arr = []
    for fcntr, f in enumerate( flist ):
        which_mask = int( f.split('_')[-1].replace('mask','').replace('.npy','') )
        curr_dic = np.load(f, allow_pickle = 1).item()

        opdic['cl_dic'][which_mask] = curr_dic['cl_dic'][0]

        lmax = curr_dic['lmax']
        fsky = curr_dic['fsky_arr'][0]
        fksy_arr.append( fsky )

        print(which_comp, which_mask, len( opdic['cl_dic'][which_mask].keys()), len( curr_dic['cl_dic'][0].keys()) )

    opdic['lmax'] = lmax
    opdic['fsky_arr'] = np.asarray( fksy_arr )

    opfname = f.replace('_mask%s.npy' %(which_mask), '.npy')
    print(opfname)
    np.save(opfname, opdic)

sys.exit()



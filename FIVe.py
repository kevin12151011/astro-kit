import time
import numpy as np 


def FIVe(data,v_grad,SNR_seed,r_grp,r_asso,num_grp,num_max):
    '''
    To find velocity-coherent structures (VCS) using FIVe algorithm 
        (Hacar et al, 2013).

    Inputs
        data: dict with:
                'pix2pc': conversion factor, from pixel to pc
                'points'  
                        'i-j-k': pixel coords (referring to img.fits) along 
                                    dec,ra,v

                                'coord': [ra,dec]
                                'I': intensity, in fits units
                                'sigma': noise level
                                'SNR'
                                'v': [km/s]

        v_grad: [km/s/pc], v gradient
        SNR_seed: SNR threshold of seeds
        r_grp: [pix], range for grouping
        r_asso: [pix], range for association
        num_group: minimal number of components in a group
        num_max: maximum number of components in each pixel

    Outputs
        Seed: list of seed points
        Group: list of groups of point names
        Iso: list of isolated points

    Notes:
        1. CAUSION: not advancing each group by one step in association!
    '''

    t0 = time.time()

    # constants
    pix2pc = data['pix2pc']

    # find seeds --------------------------------------------------------------

    print('Finding seeds...')

    Seed = []
    Iso = list(data['points'].keys())
    for name in data['points']:
        tmp = data['points'][name]

        i, j, k = [int(item) for item in name.split('-')]
        v = tmp['v'] # [km/s]

        flag_SNR = False
        if tmp['SNR']>=SNR_seed:
            flag_SNR = True 

        flag_cen = False
        num_valid = 0
        for ii in range(i-1,i+2):
            for jj in range(j-1,j+2):
                if ii==i and jj==j:
                    continue
                for kk in range(num_max):
                    name_neigh = '%d-%d-%d'%(ii,jj,kk)
                    if name_neigh in data['points']:
                        v_neigh = data['points'][name_neigh]['v']
                        v_grad_neigh = (abs(v-v_neigh)/
                                        (pix2pc*((i-ii)**2+(j-jj)**2)**.5))
                        if (v_grad_neigh<=v_grad and 
                            data['points'][name_neigh]['SNR']>=SNR_seed):
                            num_valid += 1
        if num_valid >=4:
            flag_cen = True

        if flag_SNR and flag_cen:
            Seed.append(name)
            Iso.pop(Iso.index(name))

    # FoF ---------------------------------------------------------------------

    Group = []
    Points = Seed.copy()
    L = len(Points)

    while Points:
        nxt = []
        nxt.append(Points.pop(0))
        p = 0    
        while p < len(nxt):
            name = nxt[p]
            i, j, k = [int(item) for item in name.split('-')]
            tmp = data['points'][name]
            v = tmp['v'] # [km/s]

            for ii in range(int(np.ceil(i-r_grp)),int(np.ceil(i+r_grp))):
                for jj in range(int(np.ceil(j-r_grp)),int(np.ceil(j+r_grp))):
                    if ii==i and jj==j:
                        continue
                    for kk in range(num_max):
                        name_neigh = '%d-%d-%d'%(ii,jj,kk)
                        if name_neigh not in Points:
                            continue
                        v_neigh = data['points'][name_neigh]['v']
                        v_grad_neigh = (abs(v-v_neigh)/
                                        (pix2pc*((i-ii)**2+(j-jj)**2)**.5))
                        if v_grad_neigh<=v_grad:
                            nxt.append(Points.pop(Points.index(name_neigh)))
            p += 1

        Group.append(nxt)
        print('FoF %.2f%%'%(100*(1-len(Points)/L)))

    # association -------------------------------------------------------------

    # CAUSION: not advancing each group by one step!
    for ind in range(len(Group)):
        print('Associating group %d/%d.'%(ind,len(Group)))

        group = Group[ind]
        p = 0
        while p < len(group):
            name = group[p]
            i, j, k = [int(item) for item in name.split('-')]
            tmp = data['points'][name]
            v = tmp['v'] # [km/s]

            for ii in range(int(np.ceil(i-r_asso)),int(np.ceil(i+r_asso))):
                for jj in range(int(np.ceil(j-r_asso)),int(np.ceil(j+r_asso))):
                    if ii==i and jj==j:
                        continue
                    for kk in range(num_max):
                        name_neigh = '%d-%d-%d'%(ii,jj,kk)
                        if name_neigh not in Iso:
                            continue
                        v_neigh = data['points'][name_neigh]['v']
                        v_grad_neigh = (abs(v-v_neigh)/
                                        (pix2pc*((i-ii)**2+(j-jj)**2)**.5))
                        if v_grad_neigh<=v_grad:
                            group.append(Iso.pop(Iso.index(name_neigh)))
            p += 1

        Group[ind] = group


    # drop small groups
    Group1 = []
    for group in Group:
        if len(group)>=num_grp:
            Group1.append(group)
        else:
            Iso.extend(group)

    t = time.time()-t0
    print('Time consumed for FIVe: %d s, %.1f m, %.4f h.'%(t,t/60,t/3600))

    return Seed, Group1, Iso




























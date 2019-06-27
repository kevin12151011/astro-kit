import time
import numpy as np 
from scipy.spatial import KDTree


def FIVe(Name,I,J,V,SNR,SNR_seed,v_grad,r_seed,r_grp,r_asso,frac_nei,num_grp):
    '''
    Purpose
    ---------
    To find velocity-coherent structures (VCS) using FIVe algorithm 
        (Hacar et al, 2013).


    Inputs
    --------
    Name: list of name_p in naturally sorted order
    I: 1darray, 'i' of the points
    J: 1darray
    V: 1darray
    SNR: 1darray, SNRs of points
    SNR_seed: SNR threshold of seeds
    v_grad: [km/s/pix], v gradient
    r_seed: [pix], radius for searching neighboring seeds
    r_grp: [pix], range for grouping
    r_asso: [pix], range for association
    frac_nei: minimal fraction of neighboring seeds of one seed
    num_grp: minimal number of components in a group


    Outputs
    ----------
    Seed: list of seed points
    Group: list of groups of point names
    Iso: list of isolated points
    '''


    t0 = time.time()

    # find seeds --------------------------------------------------------------

    print('Finding seeds...')

    ind_s0 = SNR >= SNR_seed
    Seed0 = Name[ind_s0] 

    # pop isolated seeds
    I_s0 = I[ind_s0] 
    J_s0 = J[ind_s0] 
    V_s0 = V[ind_s0] 
    tree_s0 = KDTree(np.transpose([I_s0,J_s0,V_s0]))
    nei_s0 = tree_s0.query_ball_point(np.transpose([I_s0,J_s0,V_s0]),r_seed)
    N_nei_s0 = np.array([len(i) for i in nei_s0])

    N_cr = frac_nei* 4/3*np.pi*r_seed**3
    ind_s = N_nei_s0 >= N_cr
    Seed = Seed0[ind_s]
    I_s = I_s0[ind_s]
    J_s = J_s0[ind_s]
    V_s = V_s0[ind_s]
    tree_s = KDTree(np.transpose([I_s,J_s,V_s]))

    # Iso
    Iso = []
    for name in Name:
        if name not in Seed:
            Iso.append(name)
    Iso = np.array(Iso)


    # FoF ---------------------------------------------------------------------

    print('FoFing...')

    Group = []
    Seed1 = list(Seed)
    L = len(Seed1)
    I_s1 = list(I_s)
    J_s1 = list(J_s)
    V_s1 = list(V_s)

    while Seed1:
        nxt = [Seed1.pop(0)]
        I_s1.pop(0)
        J_s1.pop(0)
        V_s1.pop(0)

        p = 0
        while p < len(nxt):
            print(p)
            name = nxt[p]
            ind = list(Seed).index(name) 
            i = I_s[ind]
            j = J_s[ind]
            v = V_s[ind]

            if not I_s1:
                break

            tree = KDTree(np.transpose([I_s1,J_s1,V_s1]))
            Index_n = sorted(tree.query_ball_point([i,j,v],r_grp),reverse=True)

            for index_n in Index_n:
                i_n = I_s1[index_n]
                j_n = J_s1[index_n]
                v_n = V_s1[index_n]
                if abs(v-v_n)/((i-i_n)**2+(j-j_n)**2)**.5 <= v_grad:
                    nxt.append(Seed1.pop(index_n))
                    I_s1.pop(index_n)
                    J_s1.pop(index_n)
                    V_s1.pop(index_n)
            p += 1

        Group.append(nxt)
        print('FoF %.2f%%'%(100*(1-len(Seed1)/L)))


    # association -------------------------------------------------------------

    # CAUSION: not advancing each group by one step!
    # for ind in range(len(Group)):
    #     print('Associating group %d/%d.'%(ind,len(Group)))

    #     group = Group[ind]
    #     p = 0
    #     while p < len(group):
    #         name = group[p]
    #         i, j, k = [int(item) for item in name.split('-')]
    #         tmp = data['points'][name]
    #         v = tmp['v_lsr'] # [km/s]

    #         for ii in range(int(np.ceil(i-r_asso)),int(np.ceil(i+r_asso))):
    #             for jj in range(int(np.ceil(j-r_asso)),int(np.ceil(j+r_asso))):
    #                 if ii==i and jj==j:
    #                     continue
    #                 for kk in range(num_max):
    #                     name_neigh = '%d-%d-%d'%(ii,jj,kk)
    #                     if name_neigh not in Iso:
    #                         continue
    #                     v_neigh = data['points'][name_neigh]['v_lsr']
    #                     v_grad_neigh = (abs(v-v_neigh)/
    #                                     (pix2pc*((i-ii)**2+(j-jj)**2)**.5))
    #                     if v_grad_neigh<=v_grad:
    #                         group.append(Iso.pop(Iso.index(name_neigh)))
    #         p += 1

    #     Group[ind] = group


    # # drop small groups
    # Group1 = []
    # for group in Group:
    #     if len(group)>=num_grp:
    #         Group1.append(group)
    #     else:
    #         Iso.extend(group)

    t = time.time()-t0
    print('Time consumed for FIVe: %d s, %.1f m, %.4f h.'%(t,t/60,t/3600))

    return Seed, Group




























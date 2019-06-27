import os
import shutil
import numpy as np
import montage_wrapper as mw


def re_project(orig_dir,Name_f0,proj_dir,Name_f1,para,exact=True,
    clear_files=True):
    '''
    To re-project fits files using montage_wrapper.

    Inputs
        orig_dir: directory of fits to be re-projected, should omit 
                  unnecessary '/'
        Name_f0: list of fits names before re-projection
        proj_dir: directory of re-projected fits, DO NOT LET proj_dir=orig_dir!
        Name_f1: list of fits names after re-projection, ordered like Name_f0
        para: list of parameters defining the framework of re-projected image,
              which includes:

                NAXIS1: integer
                NAXIS2: integer
                projection: 3-letter string, e.g. 'TAN'
                CRVAL1
                CRVAL2
                CRPIX1
                CRPIX2
                CD1_1
                CD1_2
                CD2_1
                CD2_2

        exact: whether the output shape exactly match the FITS header
        clear_files: whether to delete intermediate files
    
    Outputs
        re-projected fits in proj_dir

    Caveats
        1. This code is still in test.
        2. Only tested for 2D fits with one hdu.
        3. It's quite often that "exact=True" doesn't work, try to use mAdd. 
            Sometimes mAdd gives fits with zero file sizes. Maybe mAdd can't 
            handle single fits.
    '''


    # create a folder containing the original fits
    raw_dir = orig_dir+'/fits_orig'
    os.mkdir(raw_dir)
    for name_f in Name_f0:
        shutil.copy('%s/%s'%(orig_dir,name_f),raw_dir)

    # make images_table
    images_table = raw_dir+'/images_table.txt'
    mw.mImgtbl(raw_dir,images_table)

    # create header file
    (NAXIS1,NAXIS2,projection,CRVAL1,CRVAL2,CRPIX1,CRPIX2,CD1_1,CD1_2,CD2_1,
        CD2_2) = para
    f = open(raw_dir+'/header.txt','w')
    f.write('SIMPLE  = T\n')
    f.write('BITPIX  = -64\n')
    f.write('BUNIT  = none\n')
    f.write('NAXIS   = 2\n')
    f.write('NAXIS1  = %d\n'%NAXIS1)
    f.write('NAXIS2  = %d\n'%NAXIS2)
    f.write("CTYPE1  = 'RA---%s'\n"%projection)
    f.write("CTYPE2  = 'DEC--%s'\n"%projection)
    f.write('CRPIX1  = %d\n'%CRPIX1)
    f.write('CRPIX2  = %d\n'%CRPIX2)
    f.write('CRVAL1  = %f\n'%CRVAL1)
    f.write('CRVAL2  = %f\n'%CRVAL2)
    f.write('CD1_1   = %f\n'%CD1_1)
    f.write('CD1_2   = %f\n'%CD1_2)
    f.write('CD2_1   = %f\n'%CD2_1)
    f.write('CD2_2   = %f\n'%CD2_2)
    f.write('HISTORY =  By Yue Cao\n')
    f.write('END')
    f.close()

    # re-project
    stats_table = raw_dir+'/stats_table.txt'
    mw.mProjExec(images_table=images_table,template_header=raw_dir+'/header.txt',
        raw_dir=raw_dir,proj_dir=proj_dir,stats_table=stats_table,exact=exact)

    # delete intermediate files
    if clear_files:
        shutil.rmtree(raw_dir)
        for name_f in Name_f0:
            n_f = name_f.split('.')[0]
            os.remove('%s/hdu0_%s_area.fits'%(proj_dir,n_f))  

    # rename the re-projected fits
    for i in range(len(Name_f0)):
        os.rename('%s/hdu0_%s'%(proj_dir,Name_f0[i]),
            '%s/%s'%(proj_dir,Name_f1[i]))



# def re_project(orig_dir,proj_dir,header,exact=True,clear_files=True):
    '''
    To re-project fits files using montage_wrapper.

    Inputs
        orig_dir: directory of fits to be re-projected
        proj_dir: directory of re-projected fits 
        header: FITS header template to be used in generation of output FITS
        exact: whether the output shape exactly match the FITS header
        clear_files: whether to delete intermediate files
    
    Outputs
        re-projected fits in proj_dir

    Caveats
        1. This code is still in test, everything can go wrong.
        2. Only tested for 2D fits with one hdu.
        3. It's quite often that "exact=True" doesn't work, try to use mAdd. 
            Sometimes mAdd gives fits with zero file sizes. Maybe mAdd can't 
            handle single fits.
    '''

    # make images_table
    # images_table = orig_dir+'/images_table.txt'
    # mw.mImgtbl(orig_dir,images_table)

    # # re-project
    # stats_table = orig_dir+'/stats_table.txt'
    # mw.mProjExec(images_table=images_table,template_header=header,
    #     raw_dir=orig_dir,proj_dir=proj_dir,stats_table=stats_table,exact=exact)

    # # rename the re-projected fits
    # for name_f in os.listdir(proj_dir):
    #     if 'hdu0_' in name_f and '.fits' in name_f and '_area' not in name_f:
    #         os.rename('%s/%s'%(proj_dir,name_f),'%s/%s'%(proj_dir,name_f[5:]))

    # # delete intermediate files
    # if clear_files:
    #     os.remove(images_table)
    #     os.remove(stats_table)
    #     Name_f = os.listdir(proj_dir)
    #     for name_f in Name_f:
    #         if '_area.fits' in name_f:
    #             os.remove('%s/%s'%(proj_dir,name_f))  















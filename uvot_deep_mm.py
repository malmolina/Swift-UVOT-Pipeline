###########################################################################################
"""
uvot_deep_mm.py      Created By: Lea Hagen, Edited By: Mallory Molina    June 2018                   
This program has the basic framework of uvot_deep.py with the following additions:

1) Windowed frames are no longer included in the uvotimsum command (which caused fatal
errors in uvot_deep)
2) If LSS is misaligned, the program now aligns the image to the sky (counts) image to 
allow the program to continue to run smoothly
3) If ALL frames are windowed frames, the code logs the information in swift_uvot.log, 
which is in the directory that holds the original uvot_deep_mm.py copy 

This code includes the package reproject, which may need to be installed separately. 
Instructions on installation are on their website, which is linked to on the github page
"""
"""
VERSION 2.0      Created by: Mallory Molina    September 2021
Updates:
1) 2020 time-dependent throughput loss correction is applied.
2) 1x1 images are re-binned to 2x2 so they can be included
3) dead-time correction is applied
4) log now stores observations that are skipped for windowed frames (frame time != 0.0110322s)
5) log now stores observations that are skipped because they are not aspect corrected
6) observations without uat files (no star tracker) are skipped and stored in log
"""
###########################################################################################


#Import packages
import numpy as np
from astropy.io import fits
import glob
import math
import os
import subprocess
from reproject import reproject_exact
from config_uvot_mosaic import __ROOT__
from astropy.utils.data import get_pkg_data_filename
import pdb
import os.path
from astropy.time import Time
import astropy.units as u

def uvot_deep(main_dir,obs_dir, input_folders,output_prefix,filter_list=['w2','m2','w1','uu','bb','vv'],scattered_light=False):
    """
    For a set of UVOT images downloaded from HEASARC, do processing on each snapshot:
    * create a counts image
    * create exposure map
    * create LSS image
    * create mask image for bad pixels
    * create scattered light image (USE WITH CAUTION)

    Cases where an image will be skipped:
    * If imaging for a filter doesn't exist, it will be skipped, even if that filter name is in the input.
    * UVOT images are generally 2x2 binned.  If any images are unbinned, they will be skipped.
    * If a particular snapshot has no aspect correction, the astrometry is unreliable, so it will be skipped.
    
    The resulting counts images and exposure maps will be ready to use - corrected for
    LSS and bad pixels masked.  SSS has not yet been implemented, so keep that in mind,
    but it is unlikely to be an issue.


    Modeled off of Michael Siegel's code uvot_deep.pro


    Parameters
    ----------
    main_dir: string
        directory that holds the initial build of uvot_deep_mm.py

    obs_dir: string
        directory that holds the observations of interest

    input_folders : list of strings
        each item of the string is the 11-digit name of the folder downloaded from HEASARC

    output_prefix : string
        the prefix for output files (be sure to include an underscore or similar for readability)

    filter_list : list of strings
        some or all of ['w2','m2','w1','uu','bb','vv'] (default is all of them)

    scattered_light : boolean (default=False)
        choose whether to generate scattered light images - this is turned off until LMZH
        figures out which fits files should be used and writes understandable documentation


    Returns
    -------
    nothing

    """

    # full path to most recent teldef files
    caldb = os.environ['CALDB']
    teldef = {'uu':sorted(glob.glob(caldb+'/data/swift/uvota/bcf/teldef/*uu*'))[-1],
                  'bb':sorted(glob.glob(caldb+'/data/swift/uvota/bcf/teldef/*bb*'))[-1],
                  'vv':sorted(glob.glob(caldb+'/data/swift/uvota/bcf/teldef/*vv*'))[-1],
                  'w1':sorted(glob.glob(caldb+'/data/swift/uvota/bcf/teldef/*w1*'))[-1],
                  'm2':sorted(glob.glob(caldb+'/data/swift/uvota/bcf/teldef/*m2*'))[-1],
                  'w2':sorted(glob.glob(caldb+'/data/swift/uvota/bcf/teldef/*w2*'))[-1] }
                  

    # ------------------------
    # identify the filters in each snapshot
    # ------------------------

    # dictionary to hold filters that exist for each folder
    filter_exist = {key:[] for key in input_folders}
    
    for i in input_folders:

        # list all of the sky images
        sk_list = glob.glob(i + '/uvot/image/*_sk.img')

        # check that images exist
        if len(sk_list) == 0:
            print('No images found for input folder: ' + i)

        # grab the filter from the filename of each sky image
        for sk in sk_list:
            filter_name = sk[-9:-7]
            if filter_name in filter_list:
                filter_exist[i].append(filter_name)


    # ------------------------
    # go through each filter and build the images
    # ------------------------

    for filt in filter_list:

        # get the images that have observations in that filter
        obs_list = [im for im in filter_exist.keys() if filt in filter_exist[im]]

        # check that images exist
        if len(obs_list) == 0:
            print('No images found for filter: ' + filt)
            continue

        # dictionary to hold information about each image
        image_info = {'aspect_corr':[],'binning':[],'exposure':[],'frame_time':[],'extension':[],
                          'sk_image':[],'sk_image_corr':[],'exp_image':[],'exp_image_mask':[],
                          'lss_image':[],'mask_image':[],'sl_image':[] }

        # initialize HDUs to hold all of the extensions
        hdu_sk_all = fits.HDUList()
        hdu_ex_all = fits.HDUList()
        hdu_sl_all = fits.HDUList()
        
        for obs in obs_list:

            print('')
            print('*************************************************************')
            print('  observation ', obs, ', filter = ', filt)
            print('*************************************************************')
            print('')

            
            # --- 1. create the mask & LSS & scattered light images,
            #        and correct the counts and exposure images
            
            # counts image (labeled as sk)
            sk_image = obs+'/uvot/image/sw'+obs+'u'+filt+'_sk.img'
            # exposure image
            ex_image = obs+'/uvot/image/sw'+obs+'u'+filt+'_ex.img'
            # attitude file
            att_sat = obs+'/auxil/sw'+obs+'sat.fits'
            # make a more accurate attitude file
            att_uat = obs+'/auxil/sw'+obs+'uat.fits'
            corr_file = obs+'/uvot/hk/sw'+obs+'uac.hk'
            if not os.path.isfile(att_uat):
                cmd = 'uvotattcorr attfile=' + att_sat + ' corrfile=' + corr_file + ' outfile=' + att_uat + 'chatter=4'
                subprocess.run(cmd, shell=True)

            #If UAT file cannot be created, assume star tracker was not on, skip observation and store it in log file
            filename = 'swift_uvot.log'
            if not os.path.isfile(att_uat):
                if os.path.exists(main_dir+filename):
                    append_write = 'a' # append if already exists
                else:
                    append_write = 'w' # make a new file if not
                logfle = open(main_dir+filename,append_write)
                logfle.write(output_prefix+filt+', obsid '+str(obs)+' skipped, uat file not created'+'\n')
                logfle.close()
                continue


            # scattered light images
            if scattered_light:
                scattered_light(obs, filt, teldef[filt])

            # mask and bad pixel images (which also fixes the exposure map)
            mask_image(obs, filt, teldef[filt])
            
            # LSS images
            lss_image(obs, filt)

            # do corrections to sky images (LSS, mask) 
            corr_sk(obs, filt)


            # --- 2. assemble info about each extension in this observation

            with fits.open(sk_image) as hdu_sk:
                for i in range(1,len(hdu_sk)):
                    image_info['aspect_corr'].append(hdu_sk[i].header['ASPCORR'])
                    image_info['binning'].append(hdu_sk[i].header['BINX'])
                    image_info['exposure'].append(hdu_sk[i].header['EXPOSURE'])
                    image_info['frame_time'].append(hdu_sk[i].header['FRAMTIME'])
                    image_info['extension'].append(i)
                    image_info['sk_image'].append(sk_image)
                    image_info['sk_image_corr'].append(obs+'/uvot/image/sw'+obs+'u'+filt+'_sk_corr.img')
                    image_info['exp_image'].append(ex_image)
                    image_info['exp_image_mask'].append(obs+'/uvot/image/sw'+obs+'u'+filt+'_ex_mask.img')
                    image_info['lss_image'].append(obs+'/uvot/image/sw'+obs+'u'+filt+'.lss')
                    image_info['mask_image'].append(obs+'/uvot/image/sw'+obs+'u'+filt+'_mask.img')
                    image_info['sl_image'].append(obs+'/uvot/image/sw'+obs+'u'+filt+'.sl')
                    
            # --- 3. make one file with ALL OF THE EXTENSIONS

            hdu_sk_all= append_ext(hdu_sk_all, image_info['sk_image_corr'][-1], image_info)
            hdu_ex_all = append_ext(hdu_ex_all, image_info['exp_image_mask'][-1], image_info)
            hdu_er_all = append_ext(hdu_sk_all, image_info['sk_image_corr'][-1], image_info)
            if scattered_light:
                hdu_sl_all= append_ext(hdu_sl_all, image_info['sl_image'][-1], image_info)
                    

        #Write out skipped observations to log (frame time not 0.0110322s or not aspect-corrected)
        obsids,obsidxs=np.unique(np.array(image_info['sk_image']),return_index=True)
        aspcors=np.array(image_info['aspect_corr'])[obsidxs]
        ftimes=np.array(image_info['frame_time'])[obsidxs]
        filename = 'swift_uvot.log'
        if os.path.exists(main_dir+filename):
            append_write = 'a' # append if already exists
        else:
            append_write = 'w' # make a new file if not
        logfle = open(main_dir+filename,append_write)
        for i in range(0,len(aspcors)):
            if aspcors[i] != 'DIRECT':
                if aspcors[i] != 'UNICORR':
                        logfle.write(output_prefix+filt+', obsid '+str(obsids[i].split('/')[0])+'  skipped, not Aspect Corrected'+ '\n')
        for i in range(0,len(ftimes)):
            if ftimes[i] != 0.0110322:
                logfle.write(output_prefix+filt+', obsid '+str(obsids[i].split('/')[0])+'  skipped, frame time = '+str(ftimes[i])+'\n')
        logfle.close()

        # write out all of the combined extensions
        hdu_sk_all.writeto(output_prefix + filt + '_sk_all.fits', overwrite=True)
        hdu_ex_all.writeto(output_prefix + filt + '_ex_all.fits', overwrite=True)
        if scattered_light:
            hdu_sl_all.writeto(output_prefix + filt + '_sl_all.fits', overwrite=True)
            
                    
        # --- 4. stack all of the extensions together into one image

        print('')
        print('  ** stacking images')
        print('')
        

        # counts image
        cmd = 'uvotimsum ' + output_prefix + filt + '_sk_all.fits ' + \
              output_prefix + filt + '_sk.fits exclude=none clobber=yes'
        subprocess.run(cmd, shell=True)

        # exposure map
        cmd = 'uvotimsum ' + output_prefix + filt + '_ex_all.fits ' + \
              output_prefix + filt + '_ex.fits method=EXPMAP exclude=none clobber=yes'
        subprocess.run(cmd, shell=True)

        # make a count rate image too
        if os.path.isfile(output_prefix + filt + '_sk.fits'):
            #If all data was not windowed, proceed as planned
            with fits.open(output_prefix + filt + '_sk.fits') as hdu_sk, fits.open(output_prefix + filt + '_ex.fits') as hdu_ex:
                cr_hdu = fits.PrimaryHDU(data=hdu_sk[1].data/hdu_ex[1].data, header=hdu_sk[1].header)
                cr_hdu.writeto(output_prefix + filt + '_cr.fits', overwrite=True)

        #make an error image too
        if os.path.isfile(output_prefix + filt + '_sk.fits'):
            #If all data was not windowed, proceed as planned
            with fits.open(output_prefix + filt + '_sk.fits') as hdu_sk:
                err_hdu = fits.PrimaryHDU(header=hdu_sk[0].header)
                err_hdu.append(fits.ImageHDU(data=hdu_sk[1].data, header=hdu_sk[1].header))
                err_hdu.writeto(output_prefix + filt + '_sk_err.fits', overwrite=True)

        else:
            #Otherwise, delete blank fits files, and record note in log file
            filename = 'swift_uvot.log'
            if os.path.exists(main_dir+filename):
                append_write = 'a' # append if already exists
            else:
                append_write = 'w' # make a new file if not
            logfle = open(main_dir+filename,append_write)
            logfle.write(output_prefix+filt+'  -  Not processed'+ '\n')
            logfle.close()
            os.remove(obs_dir+output_prefix + filt + '_sk_all.fits')
            os.remove(obs_dir+output_prefix + filt + '_ex_all.fits')


def append_ext(hdu_all, new_fits_file, image_info):
    """
    append extenstions from new_fits_file to hdu_all (checking that the new
    extensions have an aspect correction and are 2x2 binned) after correcting for sensitivity loss.
    Parameters
    ----------
    hdu_all : HDU object
        the HDU that we're appending to
    new_fits_file : string
        name of the fits file that has extensions to append
    
    image_info : dict
        dictionary that has the extacted info about binning and aspect correction

    Returns
    -------
    hdu_all : HDU object
        the same as the input HDU, with the new extensions appended
    """
    
    with fits.open(new_fits_file) as hdu_new:

        # if this is the first image, copy over the primary header
        if len(hdu_all) == 0:
            hdu_all.append(fits.PrimaryHDU(header=hdu_new[0].header))

        wind,aspcor=[],[]
        # if frametime is 0.011032s (not windowed data) and the aspect corrections are ok, append the array
        for i in range(1,len(hdu_new)):
            dict_ind = i-len(hdu_new)
            filt=image_info['sk_image'][dict_ind].split('/')[-1][13:16]

            if ((image_info['binning'][dict_ind] == 1) | (image_info['binning'][dict_ind] == 2)) & \
                ((image_info['aspect_corr'][dict_ind] == 'DIRECT') | (image_info['aspect_corr'][dict_ind] == 'UNICORR'))& \
                (image_info['frame_time'][dict_ind] == 0.0110322):
                    #Rebin Data if necessary
                    if image_info['binning'][dict_ind] == 1:
                        hdu_new[i].data,hdu_new[i].header=rebin(hdu_new[i].data,hdu_new[i].header)

                    # Correct for Degradation of Detector
                    hdu_new[i].data = tdtl_corr(str(hdu_new[i].header['DATE-OBS']),hdu_new[i].data,filt)

                    #Correct for Read-Out (Dead) Time for full-frame observations
                    hdu_new[i].data = hdu_new[i].data*1.016

                    #Append Corrected File
                    hdu_all.append(fits.ImageHDU(data=hdu_new[i].data, header=hdu_new[i].header))
    
    return hdu_all

                        
def scattered_light(obs_folder, obs_filter, teldef_file):

    """
    Create scattered light images with the same orientation as the input snapshots

    Parameters
    ----------
    obs_folder : string
        the 11-digit name of the folder downloaded from HEASARC

    obs_filter : string
        one of the UVOT filters ['w2','m2','w1','uu','bb','vv']
    
    teldef_file : string
        full path+name for the teldef file


    Returns
    -------
    nothing

    """

    print('')
    print('  ** scattered light images')
    print('')

    # counts image (labeled as sk)
    sk_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_sk.img'
    # attitude files
    att_uat = obs_folder+'/auxil/sw'+obs_folder+'uat.fits'
    att_sat = obs_folder+'/auxil/sw'+obs_folder+'sat.fits'

    
    with fits.open(sk_image) as hdu_sk:
        
        # grab the ra/dec/roll from the sky file
        ra_pnt = str(hdu_sk[0].header['RA_PNT'])
        dec_pnt = str(hdu_sk[0].header['DEC_PNT'])
        roll_pnt = str(hdu_sk[0].header['PA_PNT'])

        # create HDU for the scattered light images
        hdu_sl = fits.HDUList()
        # copy over the primary header from the sky image
        hdu_sl.append(fits.PrimaryHDU(header=hdu_sk[0].header))

        # for each image extension, make the SL image
        for i in range(1,len(hdu_sk)):
 
            # create image
            skytime = '{:.7f}'.format( (hdu_sk[i].header['TSTART'] + hdu_sk[i].header['TSTOP'])/2 )
            cmd = 'swiftxform infile='+__ROOT__+'/scattered_light_images/scal_'+obs_filter+'_smooth_2x2.fits' + \
                  ' outfile=temp.sl attfile='+att_uat + ' teldeffile=' + teldef_file + ' method=AREA' + \
                  ' to=sky clobber=yes bitpix=-32 ra='+ra_pnt + ' dec='+dec_pnt + ' roll='+roll_pnt + \
                  ' skytime=MET:'+skytime
            subprocess.run(cmd, shell=True)
            
            # append it to the big fits file
            with fits.open('temp.sl') as hdu_sl_slice:
                hdu_sl.append(fits.ImageHDU(data=hdu_sl_slice[0].data, header=hdu_sl_slice[0].header))
                
            # delete the image
            os.remove('temp.sl')
            

        # write out all of the compiled SL images
        sl_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'.sl'
        hdu_sl.writeto(sl_image, overwrite=True)

    subprocess.run('rm .nfs*', shell=True)
    

def mask_image(obs_folder, obs_filter, teldef_file):

    """
    Create a bad pixel map, and use that to create mask images and masked exposure maps

    Parameters
    ----------
    obs_folder : string
        the 11-digit name of the folder downloaded from HEASARC

    obs_filter : string
        one of the UVOT filters ['w2','m2','w1','uu','bb','vv']
    
    teldef_file : string
        full path+name for the teldef file


    Returns
    -------
    nothing

    """
    
    print('')
    print('  ** mask images')
    print('')

    # counts image (labeled as sk)
    sk_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_sk.img'
    #exposure image
    ex_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_ex.img'
    # attitude files
    att_uat = obs_folder+'/auxil/sw'+obs_folder+'uat.fits'
    att_sat = obs_folder+'/auxil/sw'+obs_folder+'sat.fits'


    # make bad pixel map (detector coordinates)
    bad_pix = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'.badpix'
    subprocess.run('uvotbadpix infile='+sk_image + ' badpixlist=CALDB' + \
                   ' outfile='+bad_pix + ' clobber=yes', shell=True)

    # regenerate exposure maps
    # makes two images:
    # - mask image (wcs image)
    # - exposure map (wcs image) with bad pixels as NaN
    ex_image_new = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_ex_mask.img'
    mask_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_mask.img'
    cmd = 'uvotexpmap infile='+sk_image + ' outfile='+ex_image_new + ' maskfile='+mask_image + \
          ' badpixfile='+bad_pix + ' method=MEANFOV attfile='+att_uat + ' teldeffile='+teldef_file + \
          ' masktrim=25 clobber=yes'
    subprocess.run(cmd, shell=True)

    # create an exposure map with 0 at both the bad pixels and masked areas
    with fits.open(ex_image_new) as hdu_ex, fits.open(mask_image) as hdu_mask:

        # create HDU for the improved exposure map
        hdu_ex_new = fits.HDUList()
        # create HDU for the corresponding mask image
        hdu_mask_new = fits.HDUList()
        
        # copy over the primary headers
        hdu_ex_new.append(fits.PrimaryHDU(header=hdu_ex[0].header))
        hdu_mask_new.append(fits.PrimaryHDU(header=hdu_mask[0].header))

        # for each image extension, make the new images
        for i in range(1,len(hdu_ex)):

            new_mask_array = hdu_mask[i].data
            new_mask_array[np.isnan(hdu_ex[i].data)] = 0

            new_ex_array = hdu_ex[i].data
            new_ex_array[new_mask_array == 0] = 0
                       
            # append them to the big fits files
            hdu_ex_new.append(fits.ImageHDU(data=new_ex_array, header=hdu_ex[i].header))
            hdu_mask_new.append(fits.ImageHDU(data=new_mask_array, header=hdu_mask[i].header))

    # write out the new fits files
    hdu_ex_new.writeto(ex_image_new, overwrite=True)
    hdu_mask_new.writeto(mask_image, overwrite=True)
  

def lss_image(obs_folder, obs_filter):

    """
    Create LSS images for each snapshot

    Parameters
    ----------
    obs_folder : string
        the 11-digit name of the folder downloaded from HEASARC

    obs_filter : string
        one of the UVOT filters ['w2','m2','w1','uu','bb','vv']
    

    Returns
    -------
    nothing

    """

    print('')
    print('  ** LSS images')
    print('')

    # counts image (labeled as sk)
    sk_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_sk.img'
    #exposure image
    ex_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_ex.img'
    # attitude files
    att_uat = obs_folder+'/auxil/sw'+obs_folder+'uat.fits'
    att_sat = obs_folder+'/auxil/sw'+obs_folder+'sat.fits'


    # create LSS image
    lss_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'.lss'
    subprocess.run('uvotskylss infile='+sk_image + ' outfile='+lss_image + \
                   ' attfile='+att_uat +' clobber=yes', shell=True)
                   

def corr_sk(obs_folder, obs_filter):

    """
    Correct counts images for LSS and mask them

    counts_new = counts_old / lss * mask


    Parameters
    ----------
    obs_folder : string
        the 11-digit name of the folder downloaded from HEASARC

    obs_filter : string
        one of the UVOT filters ['w2','m2','w1','uu','bb','vv']
    

    Returns
    -------
    nothing

    """

    print('')
    print('  ** correcting sk images')
    print('')

    # counts image (labeled as sk)
    sk_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_sk.img'
    # LSS image
    lss_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'.lss'
    # mask image
    mask_image = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_mask.img'

    
    with fits.open(sk_image) as hdu_sk, fits.open(lss_image) as hdu_lss, fits.open(mask_image) as hdu_mask:

        # create HDU for the new counts image
        hdu_sk_new = fits.HDUList()
        # copy over the primary header
        hdu_sk_new.append(fits.PrimaryHDU(header=hdu_sk[0].header))

        # for each image extension, make the new image
        for i in range(1,len(hdu_sk)):
            #Test to make sure the LSS image is the same size and aligned with the sky image
            if len(hdu_lss[i].data) != len(hdu_sk[i].data) or len(hdu_lss[i].data[0]) != len(hdu_sk[i].data[0]):
                #If not, align images
                print('LSS Misaligned...')
                hd_test=fits.open((sk_image))[i]
                hdu1=fits.open((lss_image))[i]
                print('Aligning LSS to Sky Image')
                lss_test, footprint=reproject_exact(hdu1,hd_test.header)
                print('Aligned Images')
                # divide by lss and multiply by mask
                new_sk_array = hdu_sk[i].data / lss_test * hdu_mask[i].data
            else:
                # divide by lss and multiply by mask
                new_sk_array = hdu_sk[i].data / hdu_lss[i].data * hdu_mask[i].data



            # remove NaNs from dividing by 0
            new_sk_array[np.isnan(new_sk_array)] = 0
 
            # append to the big fits file
            hdu_sk_new.append(fits.ImageHDU(data=new_sk_array, header=hdu_sk[i].header))

    # write out the new fits file
    sk_image_corr = obs_folder+'/uvot/image/sw'+obs_folder+'u'+obs_filter+'_sk_corr.img'
    hdu_sk_new.writeto(sk_image_corr, overwrite=True)

def tdtl_corr(obs_date,imdata,obs_filter):

    """
    Calculate Time-Dependent Throughput Loss (TDTL) Correction

    Parameters
    ----------
    obs_date : float
        Observation Date
    
    imdata : array
        Uncorrected Data

    obs_filter : string
        one of the UVOT filters ['uw2','um2','uw1','uuu','ubb','uvv','uwh']


    Returns
    -------
    corrected data (image, exposure or scattered light)

    """
    print('')
    print('  ** TDTL Correction')
    print('')
    #Compile Information from Calibration File
    sens_file=fits.open('swusenscorr20041120v006.fits')
    if obs_filter == 'uvv':
        tdtl_sens=sens_file[1].data
    if obs_filter == 'ubb':
        tdtl_sens=sens_file[2].data
    if obs_filter == 'uuu':
        tdtl_sens=sens_file[3].data
    if obs_filter == 'uw1':
        tdtl_sens=sens_file[4].data
    if obs_filter == 'um2':
        tdtl_sens=sens_file[5].data
    if obs_filter == 'uw2':
        tdtl_sens=sens_file[6].data
    if obs_filter == 'uwh':
        tdtl_sens=sens_file[7].data
    sens_file.close()
    #Calculate Time between launch date and now 
    t0=Time('2005-01-01T00:00:00',format='fits',scale='tt')
    t1=Time(str(obs_date),format='fits',scale='tt')
    dti = t1 - t0
    dtsec = dti.sec
    idx = np.abs(dtsec -tdtl_sens['TIME']).argmin()
    if dtsec < tdtl_sens['TIME'][idx]:
        idx = idx-1 
    DT     = (dtsec - tdtl_sens['TIME'][idx])/31557600 #Years between start of interval and observation
    OFFSET = tdtl_sens['OFFSET'][idx]
    SLOPE  = tdtl_sens['SLOPE'][idx]

    #Apply TDTL correction
    corr_data = imdata * (1 + OFFSET) * (1 + SLOPE)**DT

    return corr_data

def rebin(obs_data,obs_header):
    """
    Rebin 1x1 binned data and update WCS

    Parameters
    ----------
    obs_data : 2D array
        1x1 binned data
    
    obs_header : dictionary
        observation header 

    Returns
    -------
    rebinned data, updated header and updated observing dictionary

    """
    
    print('')
    print('  ** Re-bin 1x1 Binned Data')
    print('')

    obs_bin = np.zeros((math.floor(obs_data.shape[0]/2),math.floor(obs_data.shape[1]/2)))

    
    #Confirm Arrays are even, if not discard last element in dimension
    if obs_data.shape[0] % 2 > 0:
        a0_max = len(obs_data)-1
    else:
        a0_max = len(obs_data)
    
    if obs_data.shape[1] % 2 > 0:
        a1_max = len(obs_data[0])-1
    else:
        a1_max = len(obs_data[0])

    #Rebin Data to 2x2 Bins
    for i in range(0,a0_max,2):
        for j in range(0,a1_max,2):
            binval = obs_data[i][j+1]+obs_data[i+1][j+1]+obs_data[i][j]+obs_data[i+1][j]
            obs_bin[int(i/2)][int(j/2)] = binval


    #Update Observation Header
    obs_header['BINX'] = 2
    obs_header['BINY'] = 2
    obs_header['CDELT1P'] = 2
    obs_header['CDELT2P'] = 2
    obs_header['NAXIS1'] = obs_bin.shape[1]
    obs_header['NAXIS2'] = obs_bin.shape[0]
    obs_header['CDELT1'] = obs_header['CDELT1']*2.
    obs_header['CDELT2'] = obs_header['CDELT2']*2.
    obs_header['CDELT1D'] = obs_header['CDELT1D']*2.
    obs_header['CDELT2D'] = obs_header['CDELT2D']*2.
    obs_header['CRPIX1'] = obs_header['CRPIX1']/2.
    obs_header['CRPIX2'] = obs_header['CRPIX2']/2.
    obs_header['CRPIX1P'] = obs_header['CRPIX1P']/2.
    obs_header['CRPIX2P'] = obs_header['CRPIX2P']/2.
    obs_header['CRPIX1D'] = obs_header['CRPIX1D']/2.
    obs_header['CRPIX2D'] = obs_header['CRPIX2D']/2.


    return(obs_bin,obs_header)




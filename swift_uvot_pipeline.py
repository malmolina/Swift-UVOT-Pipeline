###########################################################################################
"""
swift_uvot_pipeline.py      Created By: Mallory Molina    June 2018					  
This program uses uvot_deep.py created by Lea Hagen (and modeled off of Michael Siegel's 
code uvot_deep.pro) to create an automated pipeline creating mosaics of calibrated UVOT  
images for multiple observations. The data must NOT be windowed and must be 2x2 binned in
order to be included in the mosaic. Additional details can be found on the uvot_deep.py  
documentation. This code uses uvot_deep_mm.py which requires an additional python package 
reproject.py. Installation directions are given on their website.

This program requires an input file. The format of the input file is as follows
Line 1: directory that holds uvot_deep_mm.y and config_uvot_mosaic.py
Line 2: "All" or "None" - If "All", the individual frame filess will be saved. If "None", 
all individual frame files will be deleted
Lines 3 onwards: three column files with the directory to observations, prefix, and filters
Example:

/Users/userid/programs/
None
/Users/userid/data/obj1   obj1_    [w1,w2,m2]
/Users/userid/data/obj2   obj2_    [w1,w2,m2]
/Users/userid/data/obj3   obj3_    [w1,w2,m2]
/Users/userid/data/obj4   obj4_    [w1,w2,m2]

Use Python 3 to run this Code via the command line:
> python swift_uvot_pipeline.py

The code will then ask for the full path to the input file. After you input the filename, 
the code will complete the data reduction automatically and keep the requested files.
"""
"""
VERSION 2.0      Created by: Mallory Molina    September 2021
Updates:
1) Now requires 2020 TDTL calibration file
"""

"""
VERSION 3.0      Created by: Mallory Molina    February 2024
Updates:
1) Fixed Bugs in uvot_deep_mm.py
"""
###########################################################################################
import os
import glob
from shutil import copy2
import pathlib
from shutil import move
import importlib


def read_tbl(infile):
	"""This reads the input file for swift pipeline
	The first row is the directory holding uvot_deep.py and config_uvot_mosaic.py
	Second row denotes what files to keep. If "None", then only the reduced images will be saved. If "All", then everything will be kept
	The rest of the rows are the directory to the observations (tarred or not), the prefix, and filters of interest"""
	f=open(infile,'r')
	lines=f.readlines()
	f.close()
	dirs=[]
	prefix=[]
	filts=[]
	#store uvot_deep directory
	ud_dir=lines[0].replace('\n','')
	#Note whether to keep individual frames or not
	keep=lines[1].replace('\n','')
	lines=lines[2:]
	#create dictionary that holds uvot_deep parameters
	for row in lines:
		val=row.split()
		if val[0][-1] == '/':
			dirs.append(val[0])
		else:
			dirs.append(val[0]+'/')
		if val[1][-1] != '_':
			prefix.append(val[1]+'_')
		else:
			prefix.append(val[1])
		val[2]=val[2].replace('[','')
		val[2]=val[2].replace(']','')
		val2=val[2].split(',')
		filt=[]
		for i in range(0,len(val2)):
			filt.append(val2[i])
		filts.append(filt)
	obs_pars={'directories': dirs, 'prefixes':prefix, 'filters':filts}
	return ud_dir,keep,obs_pars
def filter_sort(f_array):
	"""This sorts the filters so it can be read correctly by uvot_deep"""
	filters=[]
	for flt in f_array:
		if 'w1' in flt:
			filters.append('w1')
		elif 'w2' in flt:
			filters.append('w2')
		elif 'm2' in flt:
			filters.append('m2')
		elif 'uu' in flt:
			filters.append('uu')
		elif 'bb' in flt:
			filters.append('bb')
		elif 'vv' in flt:
			filters.append('vv')
	return filters

#Start routine by entering full path to input file
infle=input('Enter Full path to input file: ')
uvot_dir,keep,obs_info=read_tbl(infle)
#Start routine
for i in range(0,len(obs_info['directories'])):
	#define main object directory, sort filters copy necessary files
	drct=obs_info['directories'][i]
	filts=filter_sort(obs_info['filters'][i])
	os.chdir(drct)
	copy2(uvot_dir+'uvot_deep_mm.py',drct+'uvot_deep_mm.py')
	copy2(uvot_dir+'config_uvot_mosaic.py',drct+'config_uvot_mosaic.py')
	copy2(uvot_dir+'swusenscorr20041120v006.fits',drct+'swusenscorr20041120v006.fits')
	#untar and/or unzip UVOT HEASARC files
	tfle=glob.glob('*.tar')
	if len(tfle) > 0:
		os.system('tar -zxvf '+tfle[0])
	dirs=glob.glob('*/')
	for j in range(0,len(dirs)):
		dirs[j]=dirs[j].replace('/','')
	for obs in dirs:
		agz_test=glob.glob(obs+'/auxil/*.gz')
		if len(agz_test) > 0:
			os.system('gunzip '+obs+'/auxil/*.gz')
			print(obs+' aux files unzipped')
		uvgz_test=glob.glob(obs+'/uvot/*/*.gz')
		if len(uvgz_test) > 0:
			os.system('gunzip '+obs+'/uvot/*/*.gz')
			print(obs+' uvot files unzipped')
	#import uvot_deep or reload
	if drct == obs_info['directories'][0]:
		import uvot_deep_mm
	else:
		importlib.reload(uvot_deep_mm)
	#run uvot_deep
	uvot_deep_mm.uvot_deep(uvot_dir,drct,dirs,obs_info['prefixes'][i],filts)
	#Store or delete frames as requested
	if keep == 'None':
		for j in range(0,len(dirs)):
			cmd='rm -rf '+dirs[j]+'/'
			os.system(cmd)
		if len(tfle) > 0:
			os.remove(dir1+manid+'/swift/'+tfle[0])
	#remove extranneous running files
	os.system('rm -rf __pycache__/')
	os.remove(drct+'config_uvot_mosaic.py')
	os.remove(drct+'uvot_deep_mm.py')
	os.remove(drct+'swusenscorr20041120v006.fits')




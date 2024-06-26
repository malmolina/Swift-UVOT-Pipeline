Swift UVOT Pipeline Patch
----------
This patch uses the `uvot_deep.py code` created by Lea Hagen (forked here). There are two major updates:

1) Fixes minor bugs in `uvot_deep.py` (see `uvot_deep_mm.py`) as described below, and
2) Creates a pipeline to automatically create mosaics of calibrated UVOT images for multiple objects using downloaded HEASARC data.

Files, software and python packages needed to run Swift UVOT Pipeline Patch
----------
Files:
1) `swift_uvot_pipeline.py`
2) `uvot_deep_mm.py`
3) `swusenscorr20041120v006.fits`
4) `config_uvot_mosaic.py`

Python packages:
1) astropy
2) reproject
Note that python 3 is required

Additional Software (see below for more detailed instructions):
1) HEASARC FTOOLS
2) HEASARC CALDB

`uvot_deep_mm.py` Update (February 2024)
----------
This patch fixes the following bugs:
1) The dead-time (i.e., read-out time) correction was originally applied to both count and exposure maps. It is now only applied to exposure maps
2) We fixed the indexing issue in the Time-Dependent Throughput Loss (TDTL) correction algorithm so that the correct offset and slope are now used in the TDTL calculation.
3) The TDTL correction was originally applied to both the count and exposure maps. It is now only applied to the counts maps.
4) The re-binning routine originally added all the measurements from each bin in the exposure maps, increasing the exposure time by a factor of 4 and artificially lowering the final count rate. We now take the average value for each new bin to preserve the correct exposure time.
5) In the original data, some edge pixels for the exposure maps were set to zero. When the data were rebinned, we did not account for this, resulting in the new, larger edge pixels to have an incorrect value. We have zeroed out all larger pixels affected by this issue.
6) The counts images were added twice during the <code>uvotimsum</code> procedure, causing an artificial inflation of the final count rate. Now all images are only added once.</li>


`uvot_deep_mm.py` Update (September 2021)
----------
This code is an updated version of the `uvot_deep.py` code created by Lea Hagen. The changes are listed below: 

1) Windowed frames are no longer included in the uvotimsum command (which caused fatal errors in uvot_deep)
2) If the large scale structure correction map (LSS image) is misaligned, the program now aligns the image to the sky (counts) image to allow the program to continue to run smoothly
3) If ALL frames in a given filter are windowed frames, the code logs the information in swift_uvot.log, which is in the directory that holds the original `uvot_deep_mm.py` copy. The final images for that filter are not made, and any empty files associated with that filter are deleted.
4) If a uat file cannot be created, the observation is skipped and the action is recorded in swift_uvot.log.
5) The time-dependent throughput loss is corrected for each frame before stacking.
6) The dead-time correction is applied to each frame before stacking.
7) 1x1 binned images are re-binned to 2x2 binning so they can be included in the final mosaic.
8) A counts error image is created

This allows for a smooth run of uvot_deep over multiple objects, regardless of data types or any issues with the LSS images.

In order to run this patch, uvot_deep must be appropriately set up (Including the HEASEARC FTOOLS and CALDB). Note that only the pre-compiled binary version is needed to run `swift_uvot_pipeline.py`. As a consequence, the original instructions are attached at the end. Additionally, the package `reproject` is required. This is a package that can be installed using `pip` in python3. 

Required packages: astropy, reproject

Python 3 is required.

Swift UVOT Pipeline Code
----------
This program uses an updated version of `uvot_deep.py` written by Lea Hagen (and modeled off of Michael Siegel's code `uvot_deep.pro`) to create an automated pipeline creating mosaics of calibrated UVOT images for multiple observations. Additional details can be found in the General UVOT-Mosaic documentation (see below). This code uses `uvot_deep_mm.py` which requires an additional python package`reproject`. Installation directions are given on their website, but it can be installed using `pip`.

This program requires an input file. The format of the input file is as follows:

Line 1: directory that holds `uvot_deep_mm.y`, `config_uvot_mosaic.py` and `swusenscorr20041120v006.fits'

Line 2: "All" or "None" - If "All", the individual frame files will be saved. If "None", 
all individual frame files will be deleted

Lines 3 onwards: three column files with the directory to observations, prefix, and filters
Example:
```
/Users/userid/programs/
None
/Users/userid/data/obj1   obj1_    [w1,w2,m2]
/Users/userid/data/obj2   obj2_    [w1,w2,m2]
/Users/userid/data/obj3   obj3_    [w1,w2,m2]
/Users/userid/data/obj4   obj4_    [w1,w2,m2]
```
Running the Swift UVOT Pipeline
----------
When running this code use `python` to initiate script in terminal using a python3 environment:
```
> python swift_uvot_pipeline.py
Enter Full path to input file: 
```
After you input the full path to the input file (such as `/Users/userid/input_files/input.dat`), the code will complete the data reduction automatically and keep the individual frames if requested. NOTE: Untarring or unzipping files is NO LONGER NECESSARY with this program. The only required steps are downloading the UVOT and Swift Auxiliary data from HEASARC (<https://heasarc.gsfc.nasa.gov/cgi-bin/W3Browse/swift.pl>) in whatever format is most convenient (tarred file or by using a script), creating the input file, and running `swift_uvot_pipeline.py`.


If the individual frames are kept, the tar file is also currently saved. If none of the individual frames are kept, the tar file holding the original download is deleted. 

Products of the Swift UVOT Pipeline
----------

The pipeline will create the following files: 

- `*_sk_all.fits`: each extension is a counts ("sky") image, in units of counts per pixel
- `*_sk.fits`: all extensions from `*_sk.fits` added together
- `*_sk_err.fits` : error image of `*_sk.fits` file, or the square root of the total counts image
- `*_ex_all.fits`: each extension is an exposure map, in units of seconds
- `*_ex.fits`: all extensions from `*_ex.fits` added together
- `*_cr.fits`: count rate image (`*_sk.fits` divided by `*_ex.fits`), in counts per second per pixel
- swift_uvot.log: A log file containing information on the objects not processed and observations excluded 

The code does not convert to magnitude or flux values. The AB magnitude system Zero points and flux conversions are found here: https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/swift/docs/uvot/uvot_caldb_AB_10wa.pdf

To convert from count rate (cr) to AB magnitude (m_AB), with the zero point ZP from the documentation, use the following equation:

`m_AB = ZP - 2.5*LOG10(cr)`

Similarly to convert from count rate (cr) to flux (f) using the conversion factor (C) from the documenation, use the following equation:

`f = C * cr`

Citing the Swift UVOT Pipeline
----------

If you use the Swift-UVOT-Pipeline, please use the reference for the most up-to-date version of the code which was presented in Molina et al. (2023) (https://ui.adsabs.harvard.edu/abs/2023ApJS..268...63M/abstract).

The code was first introduced in Molina et al. (2020b) (https://ui.adsabs.harvard.edu/abs/2020ApJS..251...11M/abstract).

-------------------------
-------------------------

UVOT MOSAIC General Instructions
----------
This code creates a calibrated mosaic from the individual UVOT snapshots downloaded from the archive.

The code and documentation are in active development.  Please consult @lea-hagen before using this.


How to use
----------

Required packages: astropy

Python 3 is required.

This is built around the UVOT processing tools that are part of HEASOFT/FTOOLS.  Instructions to download and install HEASOFT are here:
<https://heasarc.gsfc.nasa.gov/lheasoft/install.html>

You will also need the latest CALDB files.  Download/installation information is here:
<https://heasarc.gsfc.nasa.gov/docs/heasarc/caldb/install.html>

Installation of `uvot-mosaic`: Either download or clone the repository.  You can keep the code wherever you like as long as it's in your python path.


Running `uvot_deep.py` to combine and stack images
-------

Download the desired images from HEASARC (<https://heasarc.gsfc.nasa.gov/cgi-bin/W3Browse/swift.pl>) and ensure that you've chosen to download both the UVOT and Swift Auxiliary data.  The downloads will be organized in folders named with the Observation ID (e.g., 00037723002), which is a combination of the target ID (00037723) and segment (002).  All of the fits files will need to be unzipped (`gunzip */*.gz`, `gunzip */*/*.gz`, etc.). For all of the observations you wish to stack, put their folders in the same directory, and run `uvot_deep.py` from that directory.

Example: Download two observations of the edge of the M31 disk, with Obs IDs 00037723001 and 00037723002.  You will have a directory structure something like
```
~/example/00037723001
~/example/00037723001/auxil
~/example/00037723001/uvot
~/example/00037723001/uvot/hk
~/example/00037723001/uvot/image
~/example/00037723001/uvot/products
~/example/00037723002
~/example/00037723002/auxil
~/example/00037723002/uvot
~/example/00037723002/uvot/hk
~/example/00037723002/uvot/image
~/example/00037723002/uvot/products
```
From `~/example`, run
```
> import uvot_deep
> uvot_deep.uvot_deep(['00037723001','00037723002'], 'test_', ['w2','m2','w1'])
```
After a lot of verbosity from the UVOT tools, this will create several files in `~/example` for each filter (replace `ff` with the filter name).
- `test_ff_sk_all.fits`: each extension is a counts ("sky") image, in units of counts per pixel
- `test_ff_sk.fits`: all extensions from `test_ff_sk.fits` added together
- `test_ff_ex_all.fits`: each extension is an exposure map, in units of seconds
- `test_ff_ex.fits`: all extensions from `test_ff_ex.fits` added together
- `test_ff_cr.fits`: count rate image (`test_ff_sk.fits` divided by `test_ff_ex.fits`), in counts per second per pixel

It will also create files in `~/example/[obsid]/uvot/image` for each exposure, which you most likely won't need to look at.  But if you're curious, this is the list, where `[obsid]` is the Observation ID and `[ff]` is the filter.
- `sw[obsid]u[ff].badpix`: bad pixel map
- `sw[obsid]u[ff]_ex_mask.img`: masked exposure map
- `sw[obsid]u[ff].lss`: large scale sensitivity (LSS) map
- `sw[obsid]u[ff]_mask.img`: mask image
- `sw[obsid]u[ff]_sk_corr.img`: sky (counts) image, corrected for LSS and masked
- `sw[obsid]u[ff].sl`: scattered light image (assuming this option is enabled)


Running `offset_mosaic.py` to adjust background for individual snapshots
-------

The background values in UVOT images are known to change, likely due to scattered light from the Earth/sun/moon, but sometimes also from UV-bright sources in or near the field of view.  `offset_mosaic.py` is being written to do offsets between snapshots to better account for this.

Swift UVOT Pipeline Patch
----------
This patch uses the `uvot_deep.py code` created by Lea Hagen (forked here). There are two major updates:

1) Fixes minor bugs in `uvot_deep.py` (see `uvot_deep_mm.py`) as described below, and
2) Creates a pipeline to automatically create mosaics of calibrated UVOT images for multiple objects using downloaded HEASARC data.

`uvot_deep_mm.py` Update
----------
This code is basically the same as the `uvot_deep.py` code created by Lea Hagen, with three minor changes: 

1) Windowed frames are no longer included in the uvotimsum command (which caused fatal errors in uvot_deep)
2) If the large scale structure correction map (LSS image) is misaligned, the program now aligns the image to the sky (counts) image to allow the program to continue to run smoothly
3) If ALL frames in a given filter are windowed frames, the code logs the information in swift_uvot.log, which is in the directory that holds the original `uvot_deep_mm.py` copy. The final images for that filter are not made, and any empty files associated with that filter are deleted.

This allows for a smooth run of uvot_deep over multiple objects, regardless of data types or any issues with the LSS images.

In order to run this patch, uvot_deep must be appropriately set up (Including the HEASEARC FTOOLS and CALDB). As a consequence, the original instructions are attached at the end. Additionally, the package `reproject` is required. This is a program that can be installed using `pip` in python3. 

Required packages: astropy, reproject

Python 3 is required.

Swift UVOT Pipeline Code
----------
This program uses an updated version of `uvot_deep.py` written by Lea Hagen (and modeled off of Michael Siegel's code `uvot_deep.pro`) to create an automated pipeline creating mosaics of calibrated UVOT images for multiple observations. The data must NOT be windowed and must be 2x2 binned in order to be included in the mosaic. Additional details can be found in the General UVOT-Mosaic documentation (see below). This code uses `uvot_deep_mm.py` which requires an additional python package`reproject`. Installation directions are given on their website, but it can be installed using `pip`.

This program requires an input file. The format of the input file is as follows:

Line 1: directory that holds `uvot_deep_mm.y` and `config_uvot_mosaic.py`

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
When running this code use `ipython` to initiate script in terminal using a python3 environment:
```
> ipython swift_uvot_pipeline.py
Enter Full path to input file: 
```
After you input the full path to the input file (such as `/Users/userid/input_files/input.dat`), the code will complete the data reduction automatically and keep the individual frames if requested. NOTE: Untarring or unzipping files is NO LONGER NECESSARY with this program. The only required steps are downloading the UVOT and Swift Auxiliary data from HEASARC (<https://heasarc.gsfc.nasa.gov/cgi-bin/W3Browse/swift.pl>) in whatever format is most convenient (tarred file or by using a script), creating the input file, and running `swift_uvot_pipeline.py`.


If the individual frames are kept, the tar file is also currently saved. If none of the individual frames are kept, the tar file holding the original download is deleted. 

Current Updates (Under Construction)
----------
Current plans for updating this pipeline include creating more options for saving and deleting image frames, storing output files in a requested architecture by the user, and automatically downloading data from scripts stored from HEASARC.

If you would like other features not mentioned here, please contact @malmolina.

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

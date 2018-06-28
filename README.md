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

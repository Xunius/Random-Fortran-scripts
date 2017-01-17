# Random-Fortran-scripts

Some random fortran scripts for various purposes, with a f2py signature file for compiling to a Python module.

## pcmin

Compute the maximum wind speed and minimum central pressure of Tropical Cyclones, adopted from Prof. Emanuel's Fortran 77 code.

To build a python module:
```
f2py -c pcmin.pyf pcmin.f90
```
To use in python:

```
import pcmin
print pcmin.pcmin_mod.pcmin.__doc__
```

See Emanuel's original Fortran (77) and Matlab versions: ftp://texmex.mit.edu/pub/emanuel/TCMAX/


## regiongrow

Region grow algorithm by recursively searching a map from a given seed.

Two subroutines are included in the module:

- `regiongrow_real(data,seedy,seedx,val,eps,res,ny,nx,diagonal)`:

search a map `data` from a given seed denoted by `seedy` and `seedx`, and return 1s in `data` where `data` is in the range [`val` - `eps`, `val` + `eps`], and the rest of the returned map is filled up with 0s.

- `regiongrow_bin(data,seedy,seedx,val,res,ny,nx,diagonal)`:

binary region grow search. Return 1s where `data==val` and 0s otherwise.

For both, set `diagnoal` to 1 to treat diagnal cells (within the 3x3 element) as neighbours, set to 0 to exclude diagnoal cells.


To build a python module:
```
f2py -c regiongrow.pyf regiongrow.f90
```
To use in python:

```
import regiongrow
print regiongrow.regiongrow.regiongrow_real.__doc__
print regiongrow.regiongrow.regiongrow_bin.__doc__
```
**NOTE** beware of the differences in the indexing conventions of Python v.s. Fortran: Python indexing starts from 0, while Fortran starts from 1 (by default). So don't forget to increment the seed indices by 1 before calling the functions.


## conv2d

Perform 2D convolution and moving averages and assign missings when missing data amount exceeds a given threshold.

Two subroutines are included in the module:

- `convolve2d(slab,slabmask,kernel,kernelflag,max_missing,resslab,resmask,hs,ws)`:

convolve `slab` with kernel `kernel`. Missing values in `slab` are denoted by 1s in `slabmask`, where 0s denote valid points. On the contrary, 1s in `kernelflag` denote valid points in the kernel and 0s otherwise. This is to facilitate valid data count (e.g. a Gaussian kernel has near-zero values around the kernel edges where ambiguity may exist whether they should be treated as 0s). `max_missing` specifies the maximum tolerable ratio of missing values for each integration: E.g. `max_missing = 0.5`, then among the valid kernel points (`kernelflag==1`), at least 50 % of values in `data` has to be valid, otherwise the integration in is assigned a value 0 to `reslab`, and correspondingly `resmask` is assigned a value 1 to denote insufficient data. No filling or wrapping is done at the edges, the kernel is cropped when it extends outside `data`.

This is to complement similar function in `scipy.ndimage.convolve()`, which assigns a `nan` for all cells when any missing overlaps with the kernel (even when overlaps with 0s in kernel). The `astropy.convolution.convolve()` interpolates missings in the data slab, but doesn't allow control over the level of missing tolerance.

- `runmean2d(slab,slabmask,kernel,kernelflag,max_missing,resslab,resmask,hs,ws)`:

perform a 2d moving average of `data` using weights in a given `kernel`, and treat missings as in `convolve2d()`.

To build a python module:
```
f2py -c conv2d.pyf conv2d.f90
```
To use in python:

```
import conv2d
print conv2d.conv2d.convolve2d.__doc__
print conv2d.conv2d.runmean2d.__doc__
```

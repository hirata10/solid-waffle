"""
Back end routines for solid-waffle.

Functions
---------
get_version
    Version number of script
get_nside
    Function to get array size from format codes.
get_num_slices
    Gets the number of time slices in an image cube.
load_segment
    Function to load an image segment.
pyIRC_percentile
    Routine to get percentile cuts with a mask removed.
pyIRC_mean
    Routine to get mean with a mask removed.
ref_corr
    Get reference corrections from left & right pixel sets.
ref_array
    Get reference corrections from left & right pixel sets for a full list of files.
ref_array_onerow
    Similar to ref_array, but for one row being valid. Saves time if you only need that one row.
ref_array_block
    Extracts reference pixel data in a specified range of y values.
pixel_data
    Generate a 4D date cube containing information on a region of the detector.
gen_nl_cube
    Routine to get nonlinearity curve.
compute_gain_corr
    Gets the correction to the gain from using the full model versus the beta model.
compute_gain_corr_many
    Gets the correction to the gain from using the full model versus the beta model.
compute_xc_corr
    Gets the correction to the adjacent-pixel correlation from using the full model versus the beta model.
compute_xc_corr_many
    Gets the correction to the adjacent-pixel correlation from using the full model versus the beta model.
gain_alphacorr
    Routine to get IPC-corrected gain from properties of a difference image.
gain_alphabetacorr
    Get IPC+NL-corrected gain.
basic
    Basic characterization of a data cube.
corr_5x5
    Extracts 5x5 correlation matrix from light and dark data.
corrstats
    Routine to obtain statistical properties of a region of the detector across many time slices.
polychar
    Routine to characterize of a region of the detector across many time slices.
bfe
    Routines to compute the BFE coefficients.
hotpix
    Selects hot pixels.
hotpix_ipc
    Return IPC data from a list of hot pixels.
slidemed_percentile
    Sliding median function.
get_vmin_vmax
    Generates min and max range for a color bar based on inter-quartile range.

Classes
-------
IndexDictionary
    Table of indices.

"""

import copy
import sys
import warnings

import asdf
import fitsio
import numpy as np
import scipy
import scipy.ndimage
import scipy.stats
from astropy.io import fits
from scipy.signal import fftconvolve

import solid_waffle  # noqa: F401

from .ftsolve import (
    decenter,
    pad_to_N,
    solve_corr,
    solve_corr_many,
    solve_corr_vis,
)

# <== TESTING PARAMETERS ONLY ==>
#
# [these are false and should only be set to true for debugging purposes]
Test_SubBeta = False

# <== THESE FUNCTIONS DEPEND ON THE FORMAT OF THE INPUT FILES ==>

# Recall the naming convention for formatpars:
#    1 ...  999 = lab test data
# 1001 ... 1999 = dedicated test configurations
# 2001 ... 2999 = WFI flight-like data (ASDF)


def get_version():
    """Version number of script"""
    return globals()["solid_waffle"].__version__


def get_nside(formatpars):
    """
    Function to get array size from format codes.

    Parameters
    ----------
    formatpars : int
        Format code.

    Returns
    -------
    int
        The array size (including reference pixels).

    Notes
    -----
    This is 4096 for Roman, but we have codes that enable us to run the script on H2RG data.

    """

    # Lab tests
    if formatpars == 1:
        return 4096
    if formatpars == 2:
        return 2048
    if formatpars == 3:
        return 4096
    if formatpars == 4:
        return 4096
    if formatpars == 5:
        return 4096
    if formatpars == 6:
        return 4096
    if formatpars == 7:
        return 2048

    # Test configurations
    if formatpars == 1001:
        return 512
    if formatpars == 2002:
        return 512

    # asdf files
    if formatpars == 2001:
        return 4096


# Get number of time slices
def get_num_slices(formatpars, filename):
    """
    Gets the number of time slices in an image cube.

    Parameters
    ----------
    formatpars : int
        The format code.
    filename : str
        The file name.

    Returns
    -------
    int
        The number of time slices in that file.

    """

    # Switch based on input format
    if formatpars == 1 or formatpars == 2 or formatpars == 5 or formatpars == 1001:
        hdus = fits.open(filename)
        ntslice = int(hdus[0].header["NAXIS3"])
        hdus.close()
    elif formatpars == 3 or formatpars == 7:
        hdus = fits.open(filename)
        ntslice = len(hdus) - 1
        hdus.close()
    elif formatpars == 4 or formatpars == 6:
        hdus = fits.open(filename)
        ntslice = int(hdus[1].header["NAXIS3"])
        hdus.close()
    elif formatpars == 2001 or formatpars == 2002:
        af = asdf.open(filename)
        ntslice = np.shape(af["roman"]["data"])[0]
        af.close()
    else:
        raise ValueError("Invalid formatpars = {formatpars}")
    return ntslice


def load_segment(filename, formatpars, xyrange, tslices, verbose, use_fitsio=True):
    """
    Function to load an image segment.

    Parameters
    ----------
    filename : str
        Name of the source FITS file.
    formatpars : int
        Format code.
    xyrange : list of int
        List in the form ``[xmin,xmax,ymin,ymax]``.
        Numpy convention so the first row & column are 0, and excludes xmax and ymax.
    tslices : list of int
        Time slices to use (beginning slice is 1).
    verbose : bool
        Whether to print a lot of outputs.
    use_fitsio : bool, optional
        Whether to use the ``fitsio`` package for FITS files (usually no reason to
        turn this off --- False would force use of astropy instead).

    Returns
    -------
    np.array of float
        A 3D array of the data, shape (len(`tslices`), ymax-ymin, xmax-xmin).
        All formats are in "ramp slope negative" format (digital saturation is 0).

    Notes
    -----
    Floating-point return type was chosen instead of the native uint16 so that differences
    don't underflow to 65535 minus a small number. But integers are exactly represented.

    """

    if verbose:
        print("Reading:", filename)

    # Get dimensions of output cube
    nxuse = xyrange[1] - xyrange[0]
    nyuse = xyrange[3] - xyrange[2]
    ntslice_use = len(tslices)
    output_cube = np.zeros((ntslice_use, nyuse, nxuse))

    # Switch based on input format
    if formatpars == 1 or formatpars == 2 or formatpars == 5 or formatpars == 1001:
        if use_fitsio:
            fileh = fitsio.FITS(filename)
            N = get_nside(formatpars)
            for ts in range(ntslice_use):
                t = tslices[ts]
                if ts > 0 and tslices[ts] == tslices[ts - 1]:
                    output_cube[ts, :, :] = output_cube[ts - 1, :, :]  # don't read this slice again
                else:
                    output_cube[ts, :, :] = 65535 - np.array(
                        fileh[0][t - 1, xyrange[2] : xyrange[3], xyrange[0] : xyrange[1]]
                    )
            fileh.close()
        else:
            hdus = fits.open(filename)
            in_hdu = hdus[0]
            ntslice = in_hdu.data.shape[0]
            if verbose:
                print("input shape -> ", in_hdu.data.shape)
                print("number of slices =", ntslice, ", used =", ntslice_use)
            for ts in range(ntslice_use):
                t = tslices[ts]
                output_cube[ts, :, :] = (
                    65535 - in_hdu.data[t - 1, xyrange[2] : xyrange[3], xyrange[0] : xyrange[1]]
                )
            hdus.close()
    elif formatpars == 3 or formatpars == 7:
        if use_fitsio:
            fileh = fitsio.FITS(filename)
            N = get_nside(formatpars)
            for ts in range(ntslice_use):
                t = tslices[ts]
                output_cube[ts, :, :] = np.array(fileh[t][xyrange[2] : xyrange[3], xyrange[0] : xyrange[1]])
            fileh.close()
        else:
            raise ValueError("Error: non-fitsio methods not yet supported for formatpars=3 or 7")
    elif formatpars == 4:
        if use_fitsio:
            fileh = fitsio.FITS(filename)
            N = get_nside(formatpars)
            for ts in range(ntslice_use):
                t = tslices[ts]
                if ts >= 1 and t == tslices[ts - 1]:
                    output_cube[ts, :, :] = output_cube[ts - 1, :, :]  # asked for same slice again
                else:
                    output_cube[ts, :, :] = np.array(
                        fileh[1][0, t - 1, xyrange[2] : xyrange[3], xyrange[0] : xyrange[1]]
                    )
            fileh.close()
        else:
            raise ValueError("Error: non-fitsio methods not yet supported for formatpars=4")
    elif formatpars == 6:
        if use_fitsio:
            fileh = fitsio.FITS(filename)
            N = get_nside(formatpars)  # might need this in the future # noqa: F841
            for ts in range(ntslice_use):
                t = tslices[ts]
                if ts >= 1 and t == tslices[ts - 1]:
                    output_cube[ts, :, :] = output_cube[ts - 1, :, :]  # asked for same slice again
                else:
                    output_cube[ts, :, :] = 65535 - np.array(
                        fileh[1][0, t - 1, xyrange[2] : xyrange[3], xyrange[0] : xyrange[1]]
                    )
            fileh.close()
        else:
            raise ValueError("Error: non-fitsio methods not yet supported for formatpars=6")
    elif formatpars == 2001 or formatpars == 2002:
        fileh = asdf.open(filename)
        _ = get_nside(formatpars)  # might need in the future
        for ts in range(ntslice_use):
            t = tslices[ts]
            if ts >= 1 and t == tslices[ts - 1]:
                output_cube[ts, :, :] = output_cube[ts - 1, :, :]  # asked for the same slice again
            else:
                output_cube[ts, :, :] = 65535 - np.array(
                    fileh["roman"]["data"][t - 1, xyrange[2] : xyrange[3], xyrange[0] : xyrange[1]]
                )
        fileh.close()
    else:
        raise ValueError(f"Error! Invalid formatpars = {formatpars}")

    return output_cube


# <== FUNCTIONS BELOW HERE ARE INDEPENDENT OF THE INPUT FORMAT ==>


class IndexDictionary:
    """
    Table of indices.

    These are designed for consistency with the outputs lists of certain functions.


    Parameters
    ----------
    itype : int
        Index type (currently only accepts 0).

    Attributes
    ----------
    N : int
        Number of parameters.
    Nb : int
        Number of "basic" parameters.
    Nbb : int
        Number of "basic" + BFE parameters.
    p : int, optional
        Number of non-linearity coefficients (polynomial order is `p` + 1).
    ngood : int
        Column number for how many pixels in a sub-pixel are good.
    graw : int
        Column number for raw gain.
    gacorr : int
        Column number for IPC-corrected gain.
    g : int
        Column number for best estimate of the gain.
    alphaH : int
        Column number for horizontal IPC.
    alphaV : int
        Column number for vertical IPC.
    beta : int
        Column number for non-linearity parameter (in inverse electrons).
    I : int
        Column number for flat intensity in e/p/s.
    alphaD : int
        Column number for diagonal IPC.
    tCH : int
        Column number for horizontal pixel-pixel covariance.
    tCV : int
        Column number for vertical pixel-pixel covariance.
    ker0 : int, optional
        Column number for the (0,0) position in the BFE kernel.
    s : int, optional
        Size of BFE kernel (range on both axes is from -`s` to +`s`, inclusive).

    Methods
    -------
    addbfe
        Adds columns for BFE kernel.
    addhnl
        Adds columns for higher-order non-linearity kernel.

    Notes
    -----
    Additional columns:
    - If BFE is enabled:
      the column numbers for the BFE kernel (dx,dy) are `ker0 + dy*(2s+1) + dx`.
    - If non-linearity is enabled:
      The column numbers
      for the higher-order linearity terms are that the degree ``d`` term is in column
      ``Nbb + d-2``.

    It is not allowed to add non-linearity and then BFE.

    """

    def __init__(self, itype):
        # basic characterization parameter index list -- outputs for "basic" function
        if itype == 0:
            self.Nb = 11  # number of basic parameters
            #
            self.ngood = 0
            self.graw = 1
            self.gacorr = 2
            self.g = 3
            self.alphaH = 4
            self.alphaV = 5
            self.beta = 6
            self.I = 7
            self.alphaD = 8
            self.tCH = 9
            self.tCV = 10
            #
            # polychar indices: output is list [1, ***, residual]
            # intermediate terms are basic output [ind1:ind2]
            self.ind1 = 3
            self.ind2 = 9
            self.indp1 = 1
            self.indp2 = self.indp1 + self.ind2 - self.ind1
            #
            self.N = self.Nbb = self.Nb

    def addbfe(self, s):
        """
        Adds columns for the BFE kernel (the size is (2s+1, 2s+1), but flattened).

        Parameters
        ----------
        s : int
            Size of kernel (dx and dy range from -`s` to +`s`).

        Returns
        -------
        None

        """

        self.s = s
        self.ker0 = self.Nb + 2 * s * (s + 1)  # BFE (0,0) kernel index
        self.Nbb += (2 * s + 1) * (2 * s + 1)
        self.N += (2 * s + 1) * (2 * s + 1)
        if hasattr(self, "p"):
            print(self.p)
            raise ValueError("Should add BFE columns first, then non-linearity.")

    def addhnl(self, p):
        """
        Adds columns for higher-order non-linearity coefficients (p coefs, for total degree 1+p).

        Parameters
        ----------
        p : int
            Number of non-linearity coefficients (starts at degree 2 since 0 and 1 are linear,
            so the total degree is 1+p).

        Returns
        -------
        None

        """

        self.p = p
        self.N += p


def pyIRC_percentile(this_array, mask, perc, disc=True):
    """
    Routine to get percentile cuts with a mask removed.

    Parameters
    ----------
    this_array : np.array
        The array from which we want percentiles.
    mask : np.array of bool or int
        The mask (True or >=1 for good pixels, False or <=0 for bad).
        Must be the same size as `this_array`. (Usually they are the same
        shape, but this is not strictly necessary since they are flattened.)
    perc : float
        The desired percentile (between 0 and 100).
    disc : bool, optional
        Interpolate a percentile based on the input data being discrete?

    Returns
    -------
    float
        Masked percentile of `this_array`.

    """

    val = this_array.flatten()
    ma = mask.flatten()
    w = np.array([val[x] for x in np.where(ma > 0.5)])
    n = np.size(w)
    if disc:
        ctr = np.percentile(w, perc)
        n1 = np.count_nonzero(w < ctr - 0.499999)
        n2 = np.count_nonzero(w <= ctr + 0.499999)
        assert n1 <= n2
        if n1 == n2:
            return ctr
        dctr = (perc / 100.0 * n - (n1 + n2) / 2.0) / float(n2 - n1)
        return ctr + dctr
    # w -= np.modf(np.linspace(0,(1.+np.sqrt(5.))/2*(n-1), num=n))[0] - .5
    return np.percentile(w, perc)


def pyIRC_mean(this_array, mask):
    """
    Routine to get mean with a mask removed.

    Parameters
    ----------
    this_array : np.array
        The array from which we want percentiles.
    mask : np.array of bool or int
        The mask (True or >=1 for good pixels, False or <=0 for bad).
        Must be the same size as `this_array`. (Usually they are the same
        shape, but this is not strictly necessary since they are flattened.)

    Returns
    -------
    float
        The mean of the masked array.

    """

    val = this_array.flatten()
    ma = mask.flatten()
    w = np.array([val[x] for x in np.where(ma > 0.5)])
    return np.mean(w)


def ref_corr(filename, formatpars, yrange, tslices, verbose):
    """
    Get reference corrections from left & right pixel sets.

    Parameters
    ----------
    filename : str
        The input file name.
    formatpars : int
        The input format.
    yrange : list of int
        The minimum and maximum row number, ``yrange[0]<=y<yrange[1]``.
    tslices : list of int
        A list of the time slices to use in ascending order (first in the file is 1).
    verbose : bool
        Whether to talk a lot.

    Returns
    -------
    output_ref : list of np.array
        A list of reference pixel corrections in these rows. There are 2*``ntslices``+1
        entries, each a 1D numpy array (where ``ntslices`` is the length of `tslices`).

    Notes
    -----
    The contents of `output_ref` are as follows.
    In pseudocode where ``S[j]`` refers to the signal in time
    slice ``ntslices[j]``:

    - In all cases, the first ``ntslices`` elements are the median of the corresponding
      time slice, i.e., ``output_ref[j]`` is the median of ``S[j]``.

    - The next ``ntslices`` elements are the medians of the difference frames, i.e.,
      ``output_ref[ntslices+j]`` is the median of the difference of slices
      ``S[0]-S[j]``.

    - If there are at least 2 slices, then the last element is the median of the double difference
      ``(D[-2]-D[-1]) - (D[0]-D[1])``. Otherwise it is zero.
      (This is mainly used to measure curvature of the reference pixel ramp.)

    Both the left and right side reference pixels are used and medianed together.

    """

    # Side length of the array (needed to find reference pixel indexing)
    N = get_nside(formatpars)
    # Number of time slices
    ntslices = len(tslices)
    # Clear list
    output_ref = []

    # Build arrays of reference pixels
    my_array_band = load_segment(filename, formatpars, [0, N] + yrange, tslices, False)
    my_array_L = my_array_band[:, :, :4]
    my_array_R = my_array_band[:, :, -4:]
    # my_array_L = load_segment(filename, formatpars, [0,4]+yrange, tslices, False)
    # my_array_R = load_segment(filename, formatpars, [N-4,N]+yrange, tslices, False)
    my_array_LR = np.concatenate((my_array_L, my_array_R), axis=2)
    if verbose:
        print(N, my_array_LR.shape)

    for ts in range(ntslices):
        output_ref.append(np.median(my_array_LR[ts, :, :]))
    for ts in range(ntslices):
        diff_array = my_array_LR[0, :, :] - my_array_LR[ts, :, :]
        output_ref.append(np.median(diff_array))
    if ntslices > 1:
        diff_array = (
            my_array_LR[ntslices - 2, :, :]
            - my_array_LR[ntslices - 1, :, :]
            - (my_array_LR[0, :, :] - my_array_LR[1, :, :])
            * (tslices[-1] - tslices[-2])
            / float(tslices[1] - tslices[0])
        )
        output_ref.append(np.median(diff_array))
    else:
        output_ref.append(0)
    return output_ref


def ref_array(filelist, formatpars, ny, tslices, verbose):
    """
    Get reference corrections from left & right pixel sets for a full list of files.

    Parameters
    ----------
    filelist : list of str
        The list of input files to use.
    formatpars : int
        The format parameter.
    ny : int
        Number of groups to bin the rows for making a list of reference pixels
       (should be a power of 2).
    tslices : list of int
        A list of the time slices to use in ascending order (first in the file is 1).
    verbose : bool
        Whether to talk a lot.

    Returns
    -------
    output_array : np.array
        An array of the reference pixel information, including medians of differences.
        The shape is (``num_files``, `ny`, 2*``ntslice_use``+1), where ``num_files`` is the
        number of files in `filelist`; and ``ntslice_use`` is the number of time slices.

        The indexing is that ``output_array[i,iy,j]`` is the median reference pixel in file
        ``i``, in row group ``iy``, and in combination of time slices ``j`` (see notes since there
        are more than ``ntslice_use`` options for ``j``).

    Notes
    -----
    The combination of time slices (last axis of `output_array`) is as follows.
    In pseudocode where ``S[j]`` refers to the signal in time
    slice ``ntslices[j]``:

    - In all cases, the first ``ntslices`` elements are the median of the corresponding
      time slice, i.e., ``output_array[:,:,j]`` is the median of ``S[j]``.

    - The next ``ntslices`` elements are the medians of the difference frames, i.e.,
      ``output_array[:,:,ntslices+j]`` is the median of the difference of slices
      ``S[0]-S[j]``.

    - If there are at least 2 slices, then the last element ``output_array[:,:,-1]``
      is the median of the double difference
      ``(D[-2]-D[-1]) - (D[0]-D[1])``. Otherwise it is zero.
      (This is mainly used to measure curvature of the reference pixel ramp.)

    """

    num_files = len(filelist)
    ntslices = len(tslices)
    output_array = np.zeros((num_files, ny, 2 * ntslices + 1))

    dy = get_nside(formatpars) // ny
    for ifile in range(num_files):
        for iy in range(ny):
            ymin = dy * iy
            ymax = ymin + dy
            output_array[ifile, iy, :] = np.asarray(
                ref_corr(filelist[ifile], formatpars, [ymin, ymax], tslices, False)
            )
            if verbose:
                print(ifile, iy)
                print(output_array[ifile, iy, :])

    return output_array


def ref_array_onerow(filelist, formatpars, iy, ny, tslices, verbose):
    """
    Similar to ref_array, but for one row being valid. Saves time if you only need that one row.

    **Only use this function if you are absolutely sure you don't need the other rows!**

    Parameters
    ----------
    filelist : list of str
        The list of input files to use.
    formatpars : int
        The format parameter.
    iy : int
        Which row group (between 0 and `ny`-1) needs to be good.
    ny : int
        Number of groups to bin the rows for making a list of reference pixels
       (should be a power of 2).
    tslices : list of int
        A list of the time slices to use in ascending order (first in the file is 1).
    verbose : bool
        Whether to talk a lot.

    Returns
    -------
    output_array : np.array
        An array of the reference pixel information, including medians of differences.

    See Also
    --------
    ref_array
        See for description of the `output_array` (with the warning that only that one
        row group is going to be good).

    """

    num_files = len(filelist)
    ntslices = len(tslices)
    output_array = np.zeros((num_files, ny, 2 * ntslices + 1))
    dy = get_nside(formatpars) // ny
    for ifile in range(num_files):
        ymin = dy * iy
        ymax = ymin + dy
        output_array[ifile, iy, :] = np.asarray(
            ref_corr(filelist[ifile], formatpars, [ymin, ymax], tslices, False)
        )
        if verbose:
            print(ifile, iy)
            print(output_array[ifile, iy, :])
    return output_array


def ref_array_block(filelist, formatpars, yrange, tslices, verbose):
    """
    Extracts reference pixel data in a specified range of y values.

    Parameters
    ----------
    filelist : list of str
        The list of input files to use.
    formatpars : int
        The format parameter.
    yrange : list or tuple of int
        Which row range to extract; length 2 ( = ymin, ymax ), with usual Python
        convention (extracts ymin<=y<ymax).
    tslices : list of int
        A list of the time slices to use in ascending order (first in the file is 1).
    verbose : bool
        Whether to talk a lot.

    Returns
    -------
    output_array : np.array
        A 2D array, shape (``num_files``, 2*``ntslices``+1).

        The indexing is that ``output_array[i,j]`` is the median reference pixel in file
        ``i``, and in combination of time slices ``j`` (see notes since there
        are more than ``ntslice_use`` options for ``j``).

    Notes
    -----
    The combination of time slices (last axis of `output_array`) is as follows.
    In pseudocode where ``S[j]`` refers to the signal in time
    slice ``ntslices[j]``:

    - In all cases, the first ``ntslices`` elements are the median of the corresponding
      time slice, i.e., ``output_array[:,j]`` is the median of ``S[j]``.

    - The next ``ntslices`` elements are the medians of the difference frames, i.e.,
      ``output_array[:,ntslices+j]`` is the median of the difference of slices
      ``S[0]-S[j]``.

    - If there are at least 2 slices, then the last element ``output_array[:,-1]``
      is the median of the double difference
      ``(D[-2]-D[-1]) - (D[0]-D[1])``. Otherwise it is zero.
      (This is mainly used to measure curvature of the reference pixel ramp.)

    """

    num_files = len(filelist)
    ntslices = len(tslices)
    output_array = np.zeros((num_files, 2 * ntslices + 1))

    if len(yrange) < 2:
        raise ValueError(f"Error in ref_array_block: yrange = {yrange}")
    for ifile in range(num_files):
        ymin = yrange[0]
        ymax = yrange[1]
        output_array[ifile, :] = np.asarray(
            ref_corr(filelist[ifile], formatpars, [ymin, ymax], tslices, False)
        )
        if verbose:
            print(ifile)
            print(output_array[ifile, :])

    return output_array


def pixel_data(filelist, formatpars, xyrange, tslices, maskinfo, verbose):
    """
    Generate a 4D date cube containing information on a region of the detector.

    Parameters
    ----------
    filelist : list of str
        A list of the input file names.
    formatpars : int
        Which input format type to use.
    xyrange : list of int
        The rectangular region to pull, in [xmin,xmax,ymin,ymax] format.
        Python format, i.e., the first row and column are zero, and xmax and ymax are not included.
    tslices : list of int
        The time slices to use (first time slice is 1).
    maskinfo : list or tuple
        A list-like object with at least 2 elements, [``cut_offset``, ``do_mask``]. Here ``cut_offset``
        is the range around median to accept (default: 0.1, must be within 10% of median); and
        ``do_mask`` is a boolean on whether to do the masking.
        If we don't have at least 2 elements, defaults to [0.1, True].
    verbose : bool
        Whether to print lots of information.

    Returns
    -------
    output_array : np.array
        A 4D array. The shape is (``num_files``+1, ``ntslices``, dy, dx),
        where dy=ymax-ymin and dx=xmax-xmin are the sizes of the regions on the detector;
        ``ntslices`` is the number of time slices requested (length of `tslice`); and ``num_files``
        is the number of files (length of `filelist`). The last slice, ``output_array[-1,:,:,:]``
        is the mask (good=True).

    """

    # Masking parameters
    cut_offset = 0.1
    if len(maskinfo) >= 1:
        cut_offset = maskinfo[0]
    do_mask = True
    if len(maskinfo) >= 2:
        do_mask = maskinfo[1]

    num_files = len(filelist)
    ntslices = len(tslices)
    output_array = np.zeros((num_files + 1, ntslices, xyrange[3] - xyrange[2], xyrange[1] - xyrange[0]))

    for ifile in range(num_files):
        output_array[ifile, :, :, :] = load_segment(filelist[ifile], formatpars, xyrange, tslices, verbose)

    # Generate mean CDS image and consider the median
    mCDS = np.mean(output_array[0:num_files, 0, :, :], axis=0) - np.mean(
        output_array[0:num_files, -1, :, :], axis=0
    )
    mCDS_med = np.median(mCDS)
    if do_mask:
        a = (1.0 / mCDS_med) * mCDS
        goodmap = np.where(np.logical_and(a > 1 - cut_offset, a < 1 + cut_offset), 1, 0)
    else:
        goodmap = np.ones_like(mCDS)
    for f in range(num_files):
        for t in range(ntslices):
            goodmap *= np.where(output_array[f, t, :, :] > 0, 1, 0)
    if verbose:
        print("Median =", mCDS_med, "cut_offset =", cut_offset)
        print(goodmap)
        print(goodmap.shape)
    # Copy map of good pixels into the output
    for t in range(ntslices):
        output_array[num_files, t, :, :] = goodmap

    return output_array


def gen_nl_cube(filelist, formatpars, timeslice, ngrid, Ib, usemode, swi, verbose):
    """
    Routine to get nonlinearity curve.

    Parameters
    ----------
    filelist : list of str
        A list of the input file names.
    formatpars : int
        Which input format type to use.
    timeslice : int or list of int
        Which samples to use. If a list, uses ``timeslice[1]`` through ``timeslice[2]``,
        assuming reset at time
        ``timeslice[0]``. If an integer, does time slices 1 ... `timeslice`,
        assuming reset at slice 0.
    ngrid : list or tuple of int
        Number of cells on each axis; length 2, y first:``[ny,nx]`` or ``(ny,nx)``.
    Ib : variable
        Deprecated; can pass anything.
    usemode : str
        Either ``'dev'`` (deviation from beta fit) or ``'abs'`` (absolute -- zero of time is absolute).
    swi : class
        Column table.
    verbose : bool
        Whether to talk a lot.

    Returns
    -------
    output_array : np.array
        Reference corrected signal in DN; shape = (``nt``, ``ny``, ``nx``), where ``nt``
        is the number of time slices requested.
        This is the median within a file, and then we take the mean across files.
    fit_array : np.array
        The polynomial fit in DN; shape = (``nt``, ``ny``, ``nx``).
        Would be equal to `output_array` if the fit is perfect.
    deriv_array : np.array
        The derivative of the polynomial fit in DN/frame; shape = (``nt``, ``ny``, ``nx``).
    coefs_array : np.array, optional
        The polynomial coefficients for the ramps; shape (``order``+1, ``ny``, ``nx``).
        Order is ascending, i.e., constant term is ``coefs_array[0,:,:]``, then the linear
        term is ``coefs_array[1,:,:]``, etc. Only returned if `usemode` is ``'abs'``.

    """

    # Extract basic information
    nfiles = len(filelist)
    nx = ngrid[1]
    ny = ngrid[0]
    N = get_nside(formatpars)
    dx = N // nx
    dy = N // ny

    # Check whether we have a list or single slice
    if isinstance(timeslice, list):
        tref = timeslice[0]
        tmin = timeslice[1]
        tmax = timeslice[2]
        nt = tmax - tmin + 1
    else:
        tref = 0
        tmin = 1
        nt = tmax = timeslice

    output_array = np.zeros((nt, ny, nx))
    temp_array = np.zeros((nfiles, ny, nx))

    # order of polynomial fit per pixel
    my_order = min(5, tmax - tmin - 1)
    if usemode == "abs":
        my_order = swi.p

    if verbose:
        print("Nonlinear cube:")
        sys.stdout.write("  reference pixel extraction ...")
        sys.stdout.flush()

    # Extract reference information
    # Now ref_signal[ifile, iy, it] contains it_th time slice of group iy of ref pixels in file ifile
    ref_signal = ref_array(filelist, formatpars, ny, range(tmin, tmax + 1), False)

    if verbose:
        print("  done.")
        sys.stdout.write("Time slices:")
        sys.stdout.flush()

    # Now loop over times
    for t in range(tmin, tmax + 1):
        temp_array[:, :, :] = 0.0
        if verbose:
            sys.stdout.write(f" {t:2d}")
            sys.stdout.flush()
        for ifile in range(nfiles):
            val = load_segment(filelist[ifile], formatpars, [0, N, 0, N], [1, t], False)  # make 2D array
            valc = val[1, :, :] - val[0, :, :]
            for iy in range(ny):
                for ix in range(nx):
                    temp_array[ifile, iy, ix] = np.median(
                        valc[dy * iy : dy * (iy + 1), dx * ix : dx * (ix + 1)]
                    ) - (ref_signal[ifile, iy, t - tmin] - ref_signal[ifile, iy, 0])
        output_array[t - tmin, :, :] = -np.mean(temp_array, axis=0)
        # <-- note: we flipped the sign so that the signal is positive

    if verbose:
        print("")

    # Make fit and derivatives
    coefs_array = np.zeros((my_order + 1, ny, nx))
    fit_array = np.zeros_like(output_array)
    deriv_array = np.zeros_like(output_array)
    for iy in range(ny):
        for ix in range(nx):
            p = np.poly1d(
                np.polyfit(np.asarray(range(tmin, tmax + 1)) - tref, output_array[:, iy, ix], my_order)
            )
            q = np.poly1d.deriv(p)
            fit_array[:, iy, ix] = p(np.array(range(tmin, tmax + 1)) - tref)
            deriv_array[:, iy, ix] = q(np.array(range(tmin, tmax + 1)) - tref)
            coefs_array[: p.order + 1, iy, ix] = p.c[::-1]
    if usemode == "dev":
        return output_array, fit_array, deriv_array
    else:
        return output_array, fit_array, deriv_array, coefs_array


def compute_gain_corr(fit_array, deriv_array, Ib, tslices, reset_frame):
    """
    Gets the correction to the gain from using the full model versus the beta model.

    Parameters
    ----------
    fit_array : np.array of float
        Signal in DN for true curve (length ``tmax`` array, starting with frame 1).
    deriv_array : np.array of float
        Signal rate in DN/frame (length ``tmax`` array, starting with frame 1).
    Ib : float
        The charge per frame times beta (unitless).
    tslices : list of int
        Length 3 list: the time slices used for the quadratic fit in determining beta.
    reset_frame : int or float
        Reset frame (ideally int, but float would be OK if you want to represent an
        imperfect reset as we observe in Roman detectors).

    Returns
    -------
    float
        The ractional gain error, log( gain[full NL] / gain[est. quad.]),
        caused by using a beta-model for the nonlinearity curve instead of the full curve.

    See Also
    --------
    compute_gain_corr_many : Similar but for multiple superpixels at once.

    """

    # unpack time information
    ta = tslices[0] - reset_frame
    tb = tslices[1] - reset_frame
    td = tslices[2] - reset_frame
    # indices
    ina = tslices[0] - 1
    inb = tslices[1] - 1
    ind = tslices[2] - 1

    # We want the nonlinearity corrections (mean[td]-mean[tb])/(var[tad]-var[tab])
    # which, for a non-linearity curve f(t) with f(0)=0, f'(0)=1, is
    # [f(td)-f(tb)] / { [ta*(f'(td)-f'(ta))^2 + tad*f'(td)^2] - [ta*(f'(tb)-f'(ta))^2 + tab*f'(tb)^2] }
    # = [f(td)-f(tb)] / [ ta*(f'(td)^2-f'(tb)^2 - 2f'(ta)*(f'(td)-f'(tb)) )
    #                    + td*f'(td)^2 - ta*f'(td)^2 - tb*f'(tb)^2 + ta*f'(tb)^2 ]
    # = [f(td)-f(tb)] / [ -2 ta f'(ta) (f'(td)-f'(tb)) + td f'(td)^2 - tb f'(tb)^2 ]
    #
    # Let us call that factor e^{epsilon}
    # Now for the beta-model, f(t) = t - Ib * t^2
    # so e^{epsilon} = tbd [1 - Ib (tb+td) ] / [
    #                       4 Ib ta tbd (1 - 2 Ib ta) + td (1 - 2 Ib td)^2 - tb (1 - 2 Ib tb)^2 ]
    #                = [1 - Ib (tb+td) ] / [
    #                       4 Ib ta (1 - 2 Ib ta) + 1 - 4 Ib (tb+td) + Ib^2 (td^2 + tb td + tb^2) ]
    # to lowest order in Ib (what is in the notes):
    # epsilon ~ Ib (-4ta+3tb+3td)

    true_expepsilon = (fit_array[ind] - fit_array[inb]) / (
        -2 * ta * deriv_array[ina] * (deriv_array[ind] - deriv_array[inb])
        + td * deriv_array[ind] ** 2
        - tb * deriv_array[inb] ** 2
    )
    return np.log(true_expepsilon * deriv_array[0]) - Ib * (-4.0 * ta + 3.0 * tb + 3.0 * td)


def compute_gain_corr_many(fit_array, deriv_array, Ib, tslices, reset_frame, is_good):
    """
    Gets the correction to the gain from using the full model versus the beta model.

    This version works for many superpixels.

    Parameters
    ----------
    fit_array : np.array of float
        Signal in DN for true curve. Shape (``tmax``, ``ny``, ``nx``),
        where the first axis corresponds to the time stamp and the others are superpixel
        indices. Time stamps start at frame 1.
    deriv_array : np.array of float
        Signal rate in DN/frame, same shape as `fit_array`.
    Ib : np.array of float
        The charge per frame times beta (unitless). This is a 2D array, shape (``ny``, ``nx``).
    tslices : list of int
        Length 3 list: the time slices used for the quadratic fit in determining beta.
    reset_frame : int or float
        Reset frame (ideally int, but float would be OK if you want to represent an
        imperfect reset as we observe in Roman detectors).
    is_good : np.array of int or bool
        Mask array, shape (``ny``, ``nx``), True or 1 indicates good pixels,
        False or 0 indicates bad.

    Returns
    -------
    np.array of float
        The ractional gain error, log( gain[full NL] / gain[est. quad.]),
        caused by using a beta-model for the nonlinearity curve instead of the full curve.
        This is 2D with all the superpixels, shape (``ny``, ``nx``).

    See Also
    --------
    compute_gain_corr : Similar but for only one superpixel.

    """

    out_array = np.zeros_like(fit_array[0, :, :])
    ny = np.shape(fit_array)[1]
    nx = np.shape(fit_array)[2]
    for iy in range(ny):
        for ix in range(nx):
            if is_good[iy, ix] > 0.5:
                out_array[iy, ix] = compute_gain_corr(
                    fit_array[:, iy, ix], deriv_array[:, iy, ix], Ib[iy, ix], tslices, reset_frame
                )
    return out_array


def compute_xc_corr(fit_array, deriv_array, Ib, tslices, reset_frame):
    """
    Gets the correction to the adjacent-pixel correlation from using the full model versus the beta model.

    Parameters
    ----------
    fit_array : np.array of float
        Signal in DN for true curve (length ``tmax`` array, starting with frame 1).
    deriv_array : np.array of float
        Signal rate in DN/frame (length ``tmax`` array, starting with frame 1).
    Ib : float
        The charge per frame times beta (unitless).
    tslices : list of int
        Length 3 list: the time slices used for the quadratic fit in determining beta.
    reset_frame : int or float
        Reset frame (ideally int, but float would be OK if you want to represent an
        imperfect reset as we observe in Roman detectors).

    Returns
    -------
    float
        The correction to the adjacent-pixel correlation,
        (full correlation) / (beta model correlation)
        caused by using a beta-model for the nonlinearity curve instead of the full curve.

    See Also
    --------
    compute_xc_corr_many : Similar but for multiple superpixels at once.

    """

    # unpack time information
    ta = tslices[0] - reset_frame
    tb = tslices[1] - reset_frame
    # indices
    ina = tslices[0] - 1
    inb = tslices[1] - 1

    # We want the correction
    # f'(tb)^2 + ta/tab * (f'(tb)-f'(ta)) - (1 - 4 Ib tb)
    return (
        deriv_array[inb] ** 2 - ta / (tb - ta) * (deriv_array[inb] - deriv_array[ina]) ** 2
    ) / deriv_array[0] ** 2 - (1.0 - 4 * Ib * tb)


def compute_xc_corr_many(fit_array, deriv_array, Ib, tslices, reset_frame, is_good):
    """
    Gets the correction to the adjacent-pixel correlation from using the full model versus the beta model.

    This version works for many superpixels.

    Parameters
    ----------
    fit_array : np.array of float
        Signal in DN for true curve. Shape (``tmax``, ``ny``, ``nx``),
        where the first axis corresponds to the time stamp and the others are superpixel
        indices. Time stamps start at frame 1.
    deriv_array : np.array of float
        Signal rate in DN/frame, same shape as `fit_array`.
    Ib : np.array of float
        The charge per frame times beta (unitless). This is a 2D array, shape (``ny``, ``nx``).
    tslices : list of int
        Length 3 list: the time slices used for the quadratic fit in determining beta.
    reset_frame : int or float
        Reset frame (ideally int, but float would be OK if you want to represent an
        imperfect reset as we observe in Roman detectors).
    is_good : np.array of int or bool
        Mask array, shape (``ny``, ``nx``), True or 1 indicates good pixels,
        False or 0 indicates bad.

    Returns
    -------
    np.array of float
        The ractional gain error, log( gain[full NL] / gain[est. quad.]),
        caused by using a beta-model for the nonlinearity curve instead of the full curve.
        This is 2D with all the superpixels, shape (``ny``, ``nx``).

    See Also
    --------
    compute_xc_corr_many : Similar but for only one superpixel.

    """

    out_array = np.zeros_like(fit_array[0, :, :])
    ny = np.shape(fit_array)[1]
    nx = np.shape(fit_array)[2]
    for iy in range(ny):
        for ix in range(nx):
            if is_good[iy, ix] > 0.5:
                out_array[iy, ix] = compute_xc_corr(
                    fit_array[:, iy, ix], deriv_array[:, iy, ix], Ib[iy, ix], tslices, reset_frame
                )
    return out_array


def gain_alphacorr(graw, CH, CV, signal):
    """
    Routine to get IPC-corrected gain from properties of a difference image.

    Parameters
    ----------
    graw : float
        Uncorrected gain (e/DN)
    CH : float
        Horizontal correlation (DN^2)
    CV : float
        Vertical correlation (DN^2)
    signal : float
        Signal in this ramp (DN)

    Returns
    -------
    list of float
        If successful, returns a list of [gain, alphaH, alphaV]
        (with gain alpha-corrected and in e/DN).
        Returns an empty list if failed.

    """

    g = graw
    for _ in range(100):
        alphaH = CH * g / (2 * signal)
        alphaV = CV * g / (2 * signal)
        if alphaH + alphaV > 0.25:
            return []  # FAIL!
        g = graw * ((1 - 2 * (alphaH + alphaV)) ** 2 + 2 * (alphaH**2 + alphaV**2))
    return [g, alphaH, alphaV]


# Routine to get IPC+NL-corrected gain
#
# Inputs:
#   graw        = uncorrected gain (e/DN)
#   CH          = horizontal correlation (DN^2)
#   CV          = vertical correlation (DN^2)
#   signal      = signal in this ramp (DN)
#   frac_dslope = mean signal rate in (cd) / mean signal rate in (ab) - 1
#   times       = list of times [a,b,c,d] used, normalized to reference slice
#
# Output list:
#   gain g (alpha corr), e/DN
#   alphaH
#   alphaV
#   beta
#   current I (electrons per time slice)
#
# returns [] if failed
def gain_alphabetacorr(graw, CH, CV, signal, frac_dslope, times):
    """
    Get IPC+NL-corrected gain.

    Parameters
    ----------
    graw : float
        Uncorrected gain (e/DN)
    CH : float
        Horizontal correlation (DN^2)
    CV : float
        Vertical correlation (DN^2)
    signal : float
        Signal in this ramp (DN)
    frac_dslope : float
        Signal ratio for non-linearity, S_{cd}/S_{ab}-1 (unitless).
    times : list of int
        The time slices to use (length 4: [ta, tb, tc, td]).

    Returns
    -------
    list of float
        If successful, returns a list of [gain, alphaH, alphaV, beta, current]
        (with gain corrected and in e/DN; beta in 1/e; and current in e/frame).
        Returns an empty list if failed.

    Notes
    -----

    This is solving the following set of equations
    (see Hirata's brighter-fatter effect paper)::

      # in pseudocode
      graw = g * [ 1 + beta I (3tb+3td-4ta) ] / [ (1-4alpha)^2 + 2alphaH^2 + 2alphaV^2 ]
      CH = (2 I tad alphaH / g^2) [ 1 - 4alpha - 4 beta I td ]
      CV = (2 I tad alphaV / g^2) [ 1 - 4alpha - 4 beta I td ]
      signal = I tad [ 1 - beta I (ta+td) ] / g
      frac_dslope = - beta I (tc+td-ta-tb)

    """

    # Initial guess
    g = graw
    alpha = alphaH = alphaV = beta = 0
    I_ = signal * g / (times[3] - times[0])

    # Iterate
    # (100 iterations is overkill for this problem if alpha and beta are small)
    for _ in range(100):
        g = (
            graw
            * ((1 - 4 * alpha) ** 2 + 2 * (alphaH**2 + alphaV**2))
            / (1 + beta * I_ * (3 * (times[1] + times[3]) - 4 * times[0]))
        )
        if g < 1e-3:
            print("Gain did not converge")
            print("IN:", graw, CH, CV, signal, frac_dslope, times)
            print("STATUS:", g, alphaH, alphaV, alpha, I_, beta)
            raise RuntimeError("Gain did not converge")
        temp = (1 - 4 * alpha - 4 * beta * I_ * times[3]) * 2 * I_ * (times[3] - times[0]) / g**2
        alphaH = CH / temp
        alphaV = CV / temp
        if alphaH + alphaV > 0.25:
            return []  # FAIL!
        alpha = (alphaH + alphaV) / 2.0
        I_ = signal * g / (times[3] - times[0]) / (1 - beta * I_ * (times[3] + times[0]))
        beta = -frac_dslope / I_ / (times[2] + times[3] - times[0] - times[1])
        if np.fabs(beta) * I_ * (times[3] + times[0]) > 0.5:
            return []  # FAIL!

    return [g, alphaH, alphaV, beta, I_]


def basic(region_cube, dark_cube, tslices, lightref, darkref, ctrl_pars, verbose):
    """
    Basic characterization of a data cube.

    Parameters
    ----------
    region_cube : np.array
        4D array of the region of interest. shape: (num_files+1, nt, dy, dx),
        where num_files is the number of files, nt is the number of time slices, and
        (dy, dx) is the shape of the region on the SCA used. The last slice, ``region_cube[-1,:,:,:]``,
        is the mask (1 for good, 0 for bad).
    dark_cube : np.array
        Like `region_cube`, but for the darks. It is optional whether there is a separate mask
        (it isn't used, so it is OK if the first axis has length num_files or num_files+1).
    tslices : list of int
        List of the time slice numbers; length ``nt``.
    lightref : np.array
        Reference pixel table for correcting light exposures. shape = (num_files, 2*nt+1);
        the way the time axis is managed is described in :func:`ref_corr`.
    darkref : np.array
        Similar to `lightref`, but for the dark exposures.
    ctrl_pars : class
        Contains the control parameters as attributes (see Notes).
    verbose : bool
        Whether to print lots of information.

    Returns
    -------
    list
        The basic calibration parameters. The return information depends on whether ``ctrl_pars.full_corr``
        is True or False::

          # True (default):
          [number of good pixels, gain_raw, gain_acorr, gain_abcorr, aH, aV, beta, I, 0., tCH, tCV]
          # False:
          [number of good pixels, median, variance, tCH, tCV, tCD]

        Returns the null list [] if failed.

    Notes
    -----
    The `ctrl_pars` class contains the following attributes:

    - ``epsilon`` : float
      Fraction of data points to cut for computing correlations (default 0.01)
    - ``subtr_corr`` : bool
      Do mean subtraction for the IPC correlation? (default to True)
    - ``noise_corr`` : bool
      Do noise subtraction for the IPC correlation? (default to True)
    - ``reset_frame`` : int
      Reset frame (default to 0)
    - ``subtr_href`` : bool
      Horizontal reference pixel subtraction? (default to True)
    - ``full_corr`` : bool
      Which parameters to report? (default to True = standard basic pars; False = correlation data instead)
    - ``leadtrailSub`` : bool
      Perform lead-trail subtraction? (default to False)
    - ``g_ptile`` : float
      Percentile for inter-quantile range (default to 75)

    This includes a test so this won't crash if tslices[1]>=tslices[-1] but returns meaningful
    cross-correlation C_{abab} (everything else is nonsense in this case).

    """

    # Settings:
    newMeanSubMethod = True  # use False only for test/debug
    leadtrailSub = True  # subtract leading & trailing (by +/-4 pix) from horiz & vert correlations

    g_ptile = 75.0  # percentile use for inter-quantile range for variance (default: 75, giving standard IQR)

    # Extract basic parameters
    num_files = region_cube.shape[0] - 1
    nt = region_cube.shape[1]
    dy = region_cube.shape[2]
    dx = region_cube.shape[3]
    # npix = dx * dy
    if nt != len(tslices):
        raise RuntimeError("Error in pyirc.basic: incompatible number of time slices")
    if verbose:
        print("nfiles = ", num_files, ", ntimes = ", nt, ", dx,dy=", dx, dy)
    treset = 0
    if hasattr(ctrl_pars, "reset_frame"):
        treset = ctrl_pars.reset_frame

    # First get correlation parameters
    epsilon = 0.01
    if hasattr(ctrl_pars, "epsilon"):
        epsilon = ctrl_pars.epsilon
    subtr_corr = True
    if hasattr(ctrl_pars, "subtr_corr"):
        subtr_corr = ctrl_pars.subtr_corr
    noise_corr = True
    if hasattr(ctrl_pars, "noise_corr"):
        noise_corr = ctrl_pars.noise_corr
    if verbose:
        print("corr pars =", epsilon, subtr_corr, noise_corr)
    #

    # Reference pixel subtraction?
    subtr_href = True
    if hasattr(ctrl_pars, "subtr_href"):
        subtr_href = ctrl_pars.subtr_href

    # return full correlation information?
    full_corr = True
    if hasattr(ctrl_pars, "full_corr"):
        full_corr = ctrl_pars.full_corr

    # lead-trail subtraction for IPC correlations?
    if hasattr(ctrl_pars, "leadtrailSub"):
        leadtrailSub = ctrl_pars.leadtrailSub

    # quantile for variance?
    if hasattr(ctrl_pars, "g_ptile"):
        g_ptile = ctrl_pars.g_ptile

    # Get means and variances at the early and last slices
    # (i.e. 1-point information)
    gauss_iqr_in_sigmas = scipy.stats.norm.ppf(g_ptile / 100.0) * 2  # about 1.349 for g_ptile=75.
    box1 = region_cube[0:num_files, 0, :, :] - region_cube[0:num_files, 1, :, :]
    box2 = region_cube[0:num_files, 0, :, :] - region_cube[0:num_files, -1, :, :]
    box2Noise = dark_cube[0:num_files, 0, :, :] - dark_cube[0:num_files, -1, :, :]
    #
    if subtr_href:
        for f in range(num_files):
            if verbose:
                print(
                    "lightref.shape=",
                    lightref.shape,
                    "subtr ->",
                    lightref[f, nt + 1],
                    lightref[f, 2 * nt - 1],
                    darkref[f, 2 * nt - 1],
                )
            box1[f, :, :] -= lightref[f, nt + 1]
            box2[f, :, :] -= lightref[f, 2 * nt - 1]
            box2Noise[f, :, :] -= darkref[f, 2 * nt - 1]
    mean1 = np.mean(box1, axis=0)
    mean2 = np.mean(box2, axis=0)
    med1 = np.median(mean1)
    med2 = np.median(mean2)
    var1 = 0
    var2 = 0
    corr_mask = region_cube[-1, 0, :, :]
    for if1 in range(1, num_files):
        for if2 in range(if1):
            temp_box = box1[if1, :, :] - box1[if2, :, :]
            iqr1 = pyIRC_percentile(temp_box, corr_mask, g_ptile) - pyIRC_percentile(
                temp_box, corr_mask, 100 - g_ptile
            )
            temp_box = box2[if1, :, :] - box2[if2, :, :]
            iqr2 = pyIRC_percentile(temp_box, corr_mask, g_ptile) - pyIRC_percentile(
                temp_box, corr_mask, 100 - g_ptile
            )
            var1 += (iqr1 / gauss_iqr_in_sigmas) ** 2 / 2.0
            var2 += (iqr2 / gauss_iqr_in_sigmas) ** 2 / 2.0
            if verbose:
                print("Inner loop,", if1, if2, temp_box.shape)
    var1 /= num_files * (num_files - 1) / 2.0
    var2 /= num_files * (num_files - 1) / 2.0
    if var2 <= var1 and tslices[1] < tslices[-1]:
        return []  # FAIL!
    gain_raw = (med2 - med1) / (var2 - var1 + 1e-100)  # in e/DN
    # 1e-100 does nothing except to prevent an error when var1 and var2 are exactly the same

    # Correlations of neighboring pixels, in DN^2
    #
    tCH = tCV = tCD = 0
    for if1 in range(1, num_files):
        for if2 in range(if1):
            temp_box = box2[if1, :, :] - box2[if2, :, :]

            # Run through twice if we have noise, otherwise once
            nrun = 2 if noise_corr else 1
            for icorr in range(nrun):
                # clipping
                cmin = pyIRC_percentile(temp_box, corr_mask, 100 * epsilon)
                cmax = pyIRC_percentile(temp_box, corr_mask, 100 * (1 - epsilon))
                this_mask = np.where(np.logical_and(temp_box > cmin, temp_box < cmax), 1, 0) * corr_mask
                if np.sum(this_mask) < 1:
                    return []  # FAIL!
                # mean subtraction
                mean_of_temp_box = np.sum(temp_box * this_mask) / np.sum(this_mask)
                if subtr_corr and newMeanSubMethod:
                    temp_box -= mean_of_temp_box

                # Correlations in horizontal and vertical directions
                maskCV = np.sum(this_mask[:-1, :] * this_mask[1:, :])
                maskCH = np.sum(this_mask[:, :-1] * this_mask[:, 1:])
                CV = np.sum(this_mask[:-1, :] * this_mask[1:, :] * temp_box[:-1, :] * temp_box[1:, :])
                CH = np.sum(this_mask[:, :-1] * this_mask[:, 1:] * temp_box[:, :-1] * temp_box[:, 1:])
                if maskCH < 1 or maskCV < 1:
                    return []
                CH /= maskCH
                CV /= maskCV

                # diagonal directions
                if not full_corr:
                    maskCD1 = np.sum(this_mask[:-1, :-1] * this_mask[1:, 1:])
                    maskCD2 = np.sum(this_mask[:-1, 1:] * this_mask[1:, :-1])
                    CD1 = np.sum(
                        this_mask[:-1, :-1] * this_mask[1:, 1:] * temp_box[:-1, :-1] * temp_box[1:, 1:]
                    )
                    CD2 = np.sum(
                        this_mask[:-1, 1:] * this_mask[1:, :-1] * temp_box[:-1, 1:] * temp_box[1:, :-1]
                    )
                    if maskCD1 < 1 or maskCD2 < 1:
                        return []
                    CD1 /= maskCD1
                    CD2 /= maskCD2
                    CD = (CD1 + CD2) / 2.0

                if leadtrailSub:
                    maskCVx1 = np.sum(this_mask[:-1, :-4] * this_mask[1:, 4:])
                    maskCHx1 = np.sum(this_mask[:, :-5] * this_mask[:, 5:])
                    CVx1 = np.sum(
                        this_mask[:-1, :-4] * this_mask[1:, 4:] * temp_box[:-1, :-4] * temp_box[1:, 4:]
                    )
                    CHx1 = np.sum(this_mask[:, :-5] * this_mask[:, 5:] * temp_box[:, :-5] * temp_box[:, 5:])
                    if maskCHx1 < 1 or maskCVx1 < 1:
                        return []
                    CHx1 /= maskCHx1
                    CVx1 /= maskCVx1
                    maskCVx2 = np.sum(this_mask[:-1, 4:] * this_mask[1:, :-4])
                    maskCHx2 = np.sum(this_mask[:, :-3] * this_mask[:, 3:])
                    CVx2 = np.sum(
                        this_mask[:-1, 4:] * this_mask[1:, :-4] * temp_box[:-1, 4:] * temp_box[1:, :-4]
                    )
                    CHx2 = np.sum(this_mask[:, :-3] * this_mask[:, 3:] * temp_box[:, :-3] * temp_box[:, 3:])
                    if maskCHx2 < 1 or maskCVx2 < 1:
                        return []
                    CHx2 /= maskCHx2
                    CVx2 /= maskCVx2
                    CH -= (CHx1 + CHx2) / 2.0
                    CV -= (CVx1 + CVx2) / 2.0
                    #
                    # correction of the diagonal directions
                    if not full_corr:
                        maskCDx1 = np.sum(this_mask[:-1, :-5] * this_mask[1:, 5:])
                        maskCDx2 = np.sum(this_mask[:-1, :-3] * this_mask[1:, 3:])
                        maskCDx3 = np.sum(this_mask[1:, :-5] * this_mask[:-1, 5:])
                        maskCDx4 = np.sum(this_mask[1:, :-3] * this_mask[:-1, 3:])
                        CDx1 = np.sum(
                            this_mask[:-1, :-5] * this_mask[1:, 5:] * temp_box[:-1, :-5] * temp_box[1:, 5:]
                        )
                        CDx2 = np.sum(
                            this_mask[:-1, :-3] * this_mask[1:, 3:] * temp_box[:-1, :-3] * temp_box[1:, 3:]
                        )
                        CDx3 = np.sum(
                            this_mask[1:, :-5] * this_mask[:-1, 5:] * temp_box[1:, :-5] * temp_box[1:, 5:]
                        )
                        CDx4 = np.sum(
                            this_mask[1:, :-3] * this_mask[:-1, 3:] * temp_box[1:, :-3] * temp_box[1:, 3:]
                        )
                        if maskCDx1 < 1 or maskCDx2 < 1 or maskCDx3 < 1 or maskCDx4 < 1:
                            return []
                        CDx1 /= maskCDx1
                        CDx2 /= maskCDx2
                        CDx3 /= maskCDx3
                        CDx4 /= maskCDx4
                        CD -= (CDx1 + CDx2 + CDx3 + CDx4) / 4.0

                if subtr_corr and not newMeanSubMethod and not leadtrailSub:
                    CH -= mean_of_temp_box**2
                    CV -= mean_of_temp_box**2
                tCH += CH * (1 if icorr == 0 else -1)
                tCV += CV * (1 if icorr == 0 else -1)
                if not full_corr:
                    if subtr_corr and not newMeanSubMethod and not leadtrailSub:
                        CD -= mean_of_temp_box**2
                    tCD += CD * (1 if icorr == 0 else -1)

                if verbose:
                    print("pos =", if1, if2, "iteration", icorr, "cmin,cmax =", cmin, cmax)
                    print("Mask size", np.sum(this_mask), "correlations =", maskCH, maskCV, "data:", CH, CV)

                temp_box = box2Noise[if1, :, :] - box2Noise[if2, :, :]
                # end nested for loop
    #
    # Normalize covariances. Note that taking the difference of 2 frames doubled the covariance
    # matrix, so we have introduced cov_clip_corr
    xi = scipy.stats.norm.ppf(1 - epsilon)
    cov_clip_corr = (1.0 - np.sqrt(2.0 / np.pi) * xi * np.exp(-xi * xi / 2.0) / (1.0 - 2.0 * epsilon)) ** 2
    tCH /= num_files * (num_files - 1) * cov_clip_corr
    tCV /= num_files * (num_files - 1) * cov_clip_corr
    if not full_corr:
        tCD /= num_files * (num_files - 1) * cov_clip_corr

    # if we don't need full correlations, exit now
    if not full_corr:
        return [np.sum(this_mask), med2, var2, tCH, tCV, tCD]

    # Curvature information (for 2nd order NL coefficient)
    if tslices[-1] != tslices[-2]:
        if subtr_href:
            for f in range(num_files):
                box1[f, :, :] += lightref[f, nt + 1]
        boxD = (
            region_cube[0:num_files, -2, :, :]
            - region_cube[0:num_files, -1, :, :]
            - (tslices[-1] - tslices[-2]) / float(tslices[1] - tslices[0]) * box1
        )
        # difference map
        if subtr_href:
            for f in range(num_files):
                box1[f, :, :] -= lightref[f, nt + 1]
                boxD[f, :, :] -= (
                    (tslices[-1] - tslices[-2]) / float(tslices[1] - tslices[0]) * lightref[f, 2 * nt]
                )
        fac0 = fac1 = 0
        for if1 in range(num_files):
            box1R = box1[if1, :, :]
            boxDR = boxD[if1, :, :]
            c1min = pyIRC_percentile(box1R, corr_mask, 100 * epsilon)
            if c1min <= 0.5:
                c1min = 0.5  # should have no effect if successful, but prevents division by 0 if failure
            c1max = pyIRC_percentile(box1R, corr_mask, 100 * (1 - epsilon))
            cDmin = pyIRC_percentile(boxDR, corr_mask, 100 * epsilon)
            cDmax = pyIRC_percentile(boxDR, corr_mask, 100 * (1 - epsilon))
            this_file_mask = np.where(
                np.logical_and(
                    box1R > c1min, np.logical_and(box1R < c1max, np.logical_and(boxDR > cDmin, boxDR < cDmax))
                ),
                corr_mask,
                0,
            )
            fac0 += np.sum(this_file_mask * boxDR)
            fac1 += np.sum(this_file_mask * box1R)
        if fac1 < 0.5:
            return []  # FAIL!
        frac_dslope = fac0 / fac1 / ((tslices[-1] - tslices[-2]) / float(tslices[1] - tslices[0]))
    else:
        frac_dslope = 0.0
    if verbose:
        print("frac_dslope =", frac_dslope)

    if verbose:
        print("Group 1 ->", med1, var1)
        print("Group 2 ->", med2, var2)
        print("correlations in Group 2:", tCH, tCV)
        print("factors used: xi =", xi, ", cov_clip_corr =", cov_clip_corr)

    # Get alpha-corrected gains
    out = gain_alphacorr(gain_raw, tCH, tCV, med2)
    if tslices[1] >= tslices[-1] and len(out) < 1:
        return [
            np.sum(this_mask),
            gain_raw,
            gain_raw,
            gain_raw,
            0.0,
            0.0,
            0.0,
            med2 / gain_raw / (tslices[1] - tslices[0]),
            0.0,
            tCH,
            tCV,
        ]
    if len(out) < 1:
        return []  # FAIL!
    gain_acorr = out[0]
    aH = out[1]
    aV = out[2]

    if tslices[1] >= tslices[-1]:
        return [
            np.sum(this_mask),
            gain_raw,
            gain_acorr,
            gain_acorr,
            aH,
            aV,
            0.0,
            med2 / gain_acorr / (tslices[1] - tslices[0]),
            0.0,
            tCH,
            tCV,
        ]

    out = gain_alphabetacorr(gain_raw, tCH, tCV, med2, frac_dslope, [t - treset for t in tslices])
    if len(out) < 1:
        return []  # FAIL!
    gain_abcorr = out[0]
    aH = out[1]
    aV = out[2]
    beta = out[3]
    I_ = out[4]

    return [np.sum(this_mask), gain_raw, gain_acorr, gain_abcorr, aH, aV, beta, I_, 0.0, tCH, tCV]


def corr_5x5(region_cube, dark_cube, tslices, lightref, darkref, ctrl_pars, verbose):
    """
    Extracts 5x5 correlation matrix from light and dark data.

    Parameters
    ----------
    region_cube : np.array
        4D array of the region of interest. shape: (num_files+1, nt, dy, dx),
        where num_files is the number of files, nt is the number of time slices, and
        (dy, dx) is the shape of the region on the SCA used. The last slice, ``region_cube[-1,:,:,:]``,
        is the mask (1 for good, 0 for bad).
    dark_cube : np.array
        Like `region_cube`, but for the darks. It is optional whether there is a separate mask
        (it isn't used, so it is OK if the first axis has length num_files or num_files+1).
    tslices : list of int
        List of the time slice numbers; length ``nt``.
    lightref : np.array
        Reference pixel table for correcting light exposures. shape = (num_files, 2*nt+1);
        the way the time axis is managed is described in :func:`ref_corr`.
    darkref : np.array
        Similar to `lightref`, but for the dark exposures.
    ctrl_pars : class
        A class containing the control parameters as attributes. These are optional (but
        recommended); if specified, they follow the same format as in :func:`basic`.
    verbose : bool
        Whether to print lots of information.

    Returns
    -------
    list
        The return list has 5 entries. In what follows, "a", "b", and "d" correspond
        to time slices ``tslices[0]``, ``tslices[1]``, and ``tslices[-1]``:

        - Number of good pixels

        - Accumulated signal S_{bd} (in DN; median of mean method).

        - Estimated robust variance (from inter-quantile range) of S_{ab} (DN^2).

        - Estimated robust variance (from inter-quantile range) of S_{ad} (DN^2).

        - Correlation function out to 2 pixels of S_{ad}, i.e., C_{adad}, shape = (5,5)
          cenetered on (0,0).

    """

    # Settings:
    newMeanSubMethod = True  # use False only for test/debug
    leadtrailSub = True  # subtract leading & trailing (by +/-4 pix) from horiz & vert correlations

    g_ptile = 75.0  # percentile use for inter-quantile range for variance (default: 75, giving standard IQR)

    # Extract basic parameters
    num_files = region_cube.shape[0] - 1
    nt = region_cube.shape[1]
    dy = region_cube.shape[2]
    dx = region_cube.shape[3]
    # npix = dx * dy
    if nt != len(tslices):
        raise RuntimeError("Error in pyirc.corr_5x5: incompatible number of time slices")
    if verbose:
        print("nfiles = ", num_files, ", ntimes = ", nt, ", dx,dy=", dx, dy)
    # treset = 0
    # if hasattr(ctrl_pars, "reset_frame"):
    #     treset = ctrl_pars.reset_frame

    # First get correlation parameters
    epsilon = 0.01
    if hasattr(ctrl_pars, "epsilon"):
        epsilon = ctrl_pars.epsilon
    subtr_corr = True
    if hasattr(ctrl_pars, "subtr_corr"):
        subtr_corr = ctrl_pars.subtr_corr
    noise_corr = True
    if hasattr(ctrl_pars, "noise_corr"):
        noise_corr = ctrl_pars.noise_corr
    if verbose:
        print("corr pars =", epsilon, subtr_corr, noise_corr)
    #

    # Reference pixel subtraction?
    subtr_href = True
    if hasattr(ctrl_pars, "subtr_href"):
        subtr_href = ctrl_pars.subtr_href

    # lead-trail subtraction for IPC correlations?
    if hasattr(ctrl_pars, "leadtrailSub"):
        leadtrailSub = ctrl_pars.leadtrailSub

    # quantile for variance?
    if hasattr(ctrl_pars, "g_ptile"):
        g_ptile = ctrl_pars.g_ptile

    # Get means and variances at the early and last slices
    # (i.e. 1-point information)
    gauss_iqr_in_sigmas = scipy.stats.norm.ppf(g_ptile / 100.0) * 2  # about 1.349 for g_ptile=75.
    box1 = region_cube[0:num_files, 0, :, :] - region_cube[0:num_files, 1, :, :]
    box2 = region_cube[0:num_files, 0, :, :] - region_cube[0:num_files, -1, :, :]
    box2Noise = dark_cube[0:num_files, 0, :, :] - dark_cube[0:num_files, -1, :, :]
    #
    if subtr_href:
        for f in range(num_files):
            if verbose:
                print(
                    "lightref.shape=",
                    lightref.shape,
                    "subtr ->",
                    lightref[f, nt + 1],
                    lightref[f, 2 * nt - 1],
                    darkref[f, 2 * nt - 1],
                )
            box1[f, :, :] -= lightref[f, nt + 1]
            box2[f, :, :] -= lightref[f, 2 * nt - 1]
            box2Noise[f, :, :] -= darkref[f, 2 * nt - 1]
    mean1 = np.mean(box1, axis=0)
    mean2 = np.mean(box2, axis=0)
    # med1 = np.median(mean1)
    # med2 = np.median(mean2)
    med21 = np.median(mean2 - mean1)
    var1 = 0
    var2 = 0
    corr_mask = region_cube[-1, 0, :, :]

    C_shift_mean = np.zeros((dy, dx))
    tC_all = np.zeros((dy, dx))

    for if1 in range(1, num_files):
        for if2 in range(if1):
            temp_box = box1[if1, :, :] - box1[if2, :, :]
            iqr1 = pyIRC_percentile(temp_box, corr_mask, g_ptile) - pyIRC_percentile(
                temp_box, corr_mask, 100 - g_ptile
            )
            temp_box = box2[if1, :, :] - box2[if2, :, :]
            iqr2 = pyIRC_percentile(temp_box, corr_mask, g_ptile) - pyIRC_percentile(
                temp_box, corr_mask, 100 - g_ptile
            )
            var1 += (iqr1 / gauss_iqr_in_sigmas) ** 2 / 2.0
            var2 += (iqr2 / gauss_iqr_in_sigmas) ** 2 / 2.0
            if verbose:
                print("Inner loop,", if1, if2, temp_box.shape)

    var1 /= num_files * (num_files - 1) / 2.0
    var2 /= num_files * (num_files - 1) / 2.0
    if var2 <= var1 and tslices[1] < tslices[-1]:
        return []  # FAIL!

    # Correlations of neighboring pixels, in DN^2
    #
    for if1 in range(1, num_files):
        for if2 in range(if1):
            temp_box = box2[if1, :, :] - box2[if2, :, :]

            # Run through twice if we have noise, otherwise once
            nrun = 2 if noise_corr else 1
            if verbose:
                print("if1,if2=", if1, if2, " nrun: ", nrun)
            for icorr in range(nrun):
                # clipping
                cmin = pyIRC_percentile(temp_box, corr_mask, 100 * epsilon)
                cmax = pyIRC_percentile(temp_box, corr_mask, 100 * (1 - epsilon))
                this_mask = np.where(np.logical_and(temp_box > cmin, temp_box < cmax), 1, 0) * corr_mask
                if np.sum(this_mask) < 1:
                    return []  # FAIL!
                # mean subtraction
                mean_of_temp_box = np.sum(temp_box * this_mask) / np.sum(this_mask)
                if subtr_corr and newMeanSubMethod:
                    temp_box -= mean_of_temp_box

                # Correlations in all directions
                # masktmp = correlate2d(this_mask, this_mask,mode='same')
                # C_all = correlate2d(this_mask*temp_box, this_mask*temp_box, mode='same')
                dy2 = dy // 2
                dx2 = dx // 2
                masktmp = fftconvolve(this_mask, np.flip(this_mask), mode="full")[
                    dy2 : -dy2 + 1, dx2 : -dx2 + 1
                ]
                C_all = fftconvolve(this_mask * temp_box, np.flip(this_mask * temp_box), mode="full")[
                    dy2 : -dy2 + 1, dx2 : -dx2 + 1
                ]

                if np.any(masktmp < 1):
                    return []

                C_all /= masktmp

                if leadtrailSub:
                    C_pos_shift = np.zeros_like(C_all)
                    C_neg_shift = np.zeros_like(C_all)

                    C_pos_shift[:, :-8] = C_all[
                        :, 8:
                    ]  # values of the correlation matrix 8 columns to the right
                    C_neg_shift[:, 8:] = C_all[
                        :, :-8
                    ]  # values of the correlation matrix 8 columns to the left

                    # The 8 columns at the right edge just take the negative shift values,
                    # the 8 columns at the left edge just take the positive shift values,
                    # and in the middle the mean of the two shifts is computed:

                    C_shift_mean[:, 8:-8] = np.mean([C_pos_shift[:, 8:-8], C_neg_shift[:, 8:-8]], axis=0)
                    C_shift_mean[:, :8] = C_pos_shift[:, :8]
                    C_shift_mean[:, -8:] = C_neg_shift[:, -8:]

                    C_all = C_all - C_shift_mean

                # need to update the lines below to use C_all
                if subtr_corr and not newMeanSubMethod and not leadtrailSub:
                    C_all -= mean_of_temp_box**2

                tC_all += C_all * (1 if icorr == 0 else -1)

                if verbose:
                    print("pos =", if1, if2, "iteration", icorr, "cmin,cmax =", cmin, cmax)
                    # Below needs to be adjusted
                    # print ('Mask size', np.sum(this_mask), 'correlations =', maskCH, maskCV,
                    #        'data:', CH, CV)

                temp_box = box2Noise[if1, :, :] - box2Noise[if2, :, :]
                # end nested for loop

    #
    # Normalize covariances. Note that taking the difference of 2 frames doubled the covariance
    # matrix, so we have introduced cov_clip_corr
    xi = scipy.stats.norm.ppf(1 - epsilon)
    cov_clip_corr = (1.0 - np.sqrt(2.0 / np.pi) * xi * np.exp(-xi * xi / 2.0) / (1.0 - 2.0 * epsilon)) ** 2
    tC_all /= num_files * (num_files - 1) * cov_clip_corr

    # extract 5x5 matrix in the center of tC_all here:
    # hard-coded to return only 5x5 arrays, we should add option to specify
    # Find the "center" of this array
    c_y = dy - dy // 2
    c_x = dx - dx // 2
    tC_all_5x5 = tC_all[c_y - 3 : c_y + 2, c_x - 3 : c_x + 2]
    decenter_tC_all = decenter(tC_all_5x5)  # Might come in handy
    if verbose:
        print("tCH, tCV: ", decenter_tC_all[0, 1], decenter_tC_all[1, 0])

    # Return the correlations
    return [np.sum(this_mask), med21, var1, var2, tC_all_5x5]


def corrstats(lightfiles, darkfiles, formatpars, box, tslices, sensitivity_spread_cut, ctrl_pars):
    """
    Routine to obtain statistical properties of a region of the detector across many time slices.

    Parameters
    ----------
    lightfiles : list
        The list of "light" (flat field) files.
    darkfiles : list
        The list of dark files.
    formatpars : int
        The file format.
    box : list of int
        The box boundaries in the form [xmin, xmax, ymin, ymax], with the usual
        convention that pixels are included if xmin<=x<xmax and ymin<=y<ymax.
    tslices : list if int
        The list of time slices to use, with length at least 2. The first two entries
        are ``tmin`` and ``tmax`` (with tmin<=t<tmax). If the length is exactly 2, then
        all time slice combinations are computed (see Returns). Otherwise, the additional
        entries are "deltas". For example, [4,10,1,3] means that computations are done for
        4<=ti<tj<10, but only with combinations that have tj-ti equal to 1 or 3.
    sensitivity_spread_cut : float
        What percentage response from median flat to cut for identifying good pixels (typically 0.1).
    ctrl_pars : class
        A class containing the control parameters as attributes. These are optional (but
        recommended); if specified, they follow the same format as in :func:`basic`.

    Returns
    -------
    data : np.array
        Array of return information, shape = (nt, nt, 6), where nt is the number of time slices
        (so data for time slice pair ti,tj is in ``data[ti-tmin,tj-tmin,:]``). The fields on the last
        axis are:

        - Number of good pixels

        - Accumulated signal S_{ab} (in DN; median of mean method).

        - Variance of S_{ab} (DN^2).

        - Correlation function for horizontal pixels, C_{abab}(1,0) (DN^2).

        - Correlation function for vertical pixels, C_{abab}(0,1) (DN^2).

        - Correlation function for diagonal pixels, C_{abab}(+/-1,+/-1), averaged
          over the two diagonal directions (DN^2).

    """

    # make copy of ctrl_pars, but force 5th element to be False
    ctrl_pars2 = copy.copy(ctrl_pars)
    ctrl_pars2.full_corr = False

    tmin = tslices[0]
    tmax = tslices[1]
    nt = tmax - tmin
    # build cube of good pixels, medians, variances, correlations
    data = np.zeros((nt, nt, 6))
    # and get mask (last 'time' slice) -- only thing we are extracting from region_cube_X
    region_cube_X = pixel_data(
        lightfiles,
        formatpars,
        box[:4],
        [tmin, tmax - 1, tmax - 1, tmax - 1],
        [sensitivity_spread_cut, True],
        False,
    )

    # Get list of (good pix, median, var, cov_H, cov_V)
    for ti in range(nt - 1):
        for tj in range(ti + 1, nt):
            if tslices[2:] == [] or tj - ti in tslices[2:] or tj - ti == nt - 1:
                t1 = tmin + ti
                t2 = tmin + tj
                tarray = [t1, t2, t2, t2]
                lightref = ref_array_block(lightfiles, formatpars, box[2:4], tarray, False)
                darkref = ref_array_block(darkfiles, formatpars, box[2:4], tarray, False)
                if not ctrl_pars.subtr_href:
                    lightref[:, :] = 0.0
                    darkref[:, :] = 0.0
                region_cube = pixel_data(
                    lightfiles, formatpars, box[:4], tarray, [sensitivity_spread_cut, False], False
                )
                dark_cube = pixel_data(
                    darkfiles, formatpars, box[:4], tarray, [sensitivity_spread_cut, False], False
                )
                # switch to the mask from above
                region_cube[-1, :, :, :] = region_cube_X[-1, :, :, :]
                dark_cube[-1, :, :, :] = region_cube_X[-1, :, :, :]
                B = basic(region_cube, dark_cube, tarray, lightref, darkref, ctrl_pars2, False)
                if len(B) == 6:
                    data[ti, tj, :] = np.asarray(B)
                # print (t1, t2, data[ti,tj,:], len(B))

    return data


def polychar(
    lightfiles,
    darkfiles,
    formatpars,
    box,
    tslices,
    sensitivity_spread_cut,
    ctrl_pars,
    addInfo,
    swi,
    corrstats_data=None,
):
    """
    Routine to characterize of a region of the detector across many time slices.

    Parameters
    ----------
    lightfiles : list
        The list of "light" (flat field) files.
    darkfiles : list
        The list of dark files.
    formatpars : int
        The file format.
    box : list of int
        The box boundaries in the form [xmin, xmax, ymin, ymax], with the usual
        convention that pixels are included if xmin<=x<xmax and ymin<=y<ymax.
    tslices : list if int
        The list of time slices to use, with length at least 2. The first two entries
        are ``tmin`` and ``tmax`` (with tmin<=t<tmax). If the length is exactly 2, then
        all time slice combinations are computed (see Returns). Otherwise, the additional
        entries are "deltas". For example, [4,10,1,3] means that computations are done for
        4<=ti<tj<10, but only with combinations that have tj-ti equal to 1 or 3.
    sensitivity_spread_cut : float
        What percentage response from median flat to cut for identifying good pixels (typically 0.1).
    ctrl_pars : class
        A class containing the control parameters as attributes. These are optional (but
        recommended); see the Notes.
    addInfo : list
        Some additional information needed for IPNL corrections to the inferred gain and IPC data.
        The length may be 0 (null, no corrections), 2, or 3. If there is a correction, then the entries
        are:

        - ``addInfo[0]`` : str
          Either 'bfe' or 'nlipc' (which form of IPNL to assume dominant).

        - ``addInfo[1]`` : np.array
          BFE kernel, shape (2s+1, 2s+1), centered at 0, units in inverse electrons.

        - ``addInfo[2]`` : np.array, optional
          1D array of polynomial coefficients, needed if ``ctrl_pars.use_allorder`` is True.
          This is in DN-based units, starting with the quadratic coefficient (unit: DN^-1).

    swi : class
        Column table.
    corrstats_data : np.array, optional
        If given, saved data from :func:`corrstats` (saves time if alraedy computed).

    Returns
    -------
    list
       The list entries are [isgood (1 = good, 0 = bad), gain (e/DN), alpha_H (IPC), alpha_V (IPC),
       beta (1/e), Intensity (e/frame), alpha_D (IPC), change in alpha from previous iteration
       (residual)]. Returns the empty list [] in the event of a failure.

    Notes
    -----
    The `ctrl_pars` class contains the following attributes.
    They follow the same format as in :func:`basic`, except that ``use_allorder`` is added:

    - ``epsilon`` : float
      Fraction of data points to cut for computing correlations (default 0.01)
    - ``subtr_corr`` : bool
      Do mean subtraction for the IPC correlation? (default to True)
    - ``noise_corr`` : bool
      Do noise subtraction for the IPC correlation? (default to True)
    - ``reset_frame`` : int
      Reset frame (default to 0)
    - ``subtr_href`` : bool
      Horizontal reference pixel subtraction? (default to True)
    - ``full_corr`` : bool
      Which parameters to report? (default to True = standard basic pars; False = correlation data instead)
    - ``leadtrailSub`` : bool
      Perform lead-trail subtraction? (default to False)
    - ``g_ptile`` : float
      Percentile for inter-quantile range (default to 75)
    - ``use_allorder`` : bool
      Whether to use the full polynomial expansion for the non-linearity correction? (default to False)

    """

    # Check whether we have non-linearity information
    if hasattr(ctrl_pars, "use_allorder") and ctrl_pars.use_allorder and len(addInfo) < 3:
        print("Error: polychar: not enough fields in addInfo")
        return []

    # Check time range
    if len(tslices) < 4:
        print("Error: polychar: not enough data", tslices)
        return []
    if tslices[2] >= tslices[3] or tslices[3] >= tslices[1] - tslices[0] or tslices[1] - tslices[0] < 3:
        print("Error: polychar: invalid slices range", tslices)
        return []

    # Get correlation function data (including adjacent steps))
    if corrstats_data is None:
        data = corrstats(
            lightfiles, darkfiles, formatpars, box, tslices + [1], sensitivity_spread_cut, ctrl_pars
        )
    else:
        data = np.copy(corrstats_data)

    # check if this is good
    nt = tslices[1] - tslices[0]
    for ti in range(nt - 1):
        for tj in range(ti + 1, nt):
            if data[ti, tj, 0] == 0 and tj - ti in [1, tslices[2], tslices[3]]:
                return [0, 0, 0, 0, 0, 0]

    # Determine whether we are applying corrections
    applyCorr = False
    if len(addInfo) >= 2:
        applyCorr = True
        typeCorr = addInfo[0]
        ipnl = addInfo[1]
        sBFE = np.shape(ipnl)[0] // 2

    # Fit of differences as a function of slice number
    # slope = -2*beta*I^2/g
    # intercept = (I - beta I^2)/g
    npts = tslices[1] - tslices[0] - 1
    diff_frames = np.zeros(npts)
    for j in range(npts):
        diff_frames[j] = data[j, j + 1, 1]  # median from frame tslices[0]+j -> tslices[0]+j+1
    slopemed, icpt = np.linalg.lstsq(
        np.vstack([np.array(range(npts)) + tslices[0] - ctrl_pars.reset_frame, np.ones(npts)]).T,
        diff_frames,
        rcond=-1,
    )[0]
    # If using 'allorder', let's subtract out the higher-order terms:
    if hasattr(ctrl_pars, "use_allorder") and ctrl_pars.use_allorder:
        xr = np.array(range(npts)) + tslices[0] - ctrl_pars.reset_frame
        i = 100
        err = 10
        etarget = 1e-9 * np.abs(icpt)
        while i >= 0 and err > etarget:
            I__g = icpt - 0.5 * slopemed  # might use this in the future # noqa: F841
            diff_frames_reduced = diff_frames.copy()
            icpt_old = icpt
            slopemed_old = slopemed
            for j in range(3, swi.p + 1):
                diff_frames_reduced -= (
                    addInfo[2][j - 2] * ((xr + 1) ** j - xr**j) * (icpt - slopemed * 0.5) ** j
                )
            slopemed, icpt = np.linalg.lstsq(np.vstack([xr, np.ones(npts)]).T, diff_frames_reduced, rcond=-1)[
                0
            ]
            err = np.sqrt((icpt - icpt_old) ** 2 + (slopemed - slopemed_old) ** 2)
            if i == 0:
                print("higher order loop failed to converge {:12.5E} vs {:12.5E} (target)", err, etarget)
                return []

    # Difference of correlation functions
    #
    # Cdiff = I/g^2 * ((1-4a)^2 + 2aH^2 + 2aV^2) * t_{bd}
    # - 4(1-8a)beta I^2/g^2 * (t_{ad}t_d - t_{ab}t_b + (e-1)/2*t_{bd})
    # where e = npts2 is number of bins averaged together
    #
    # and horizontal and vertical cross-correlations
    # CH = 2 I t_{ab} / g^2 * ( 1-4a - 4 beta (I t_b + 1/2 + (e-1)/2*I) ) * aH
    # CV = 2 I t_{ab} / g^2 * ( 1-4a - 4 beta (I t_b + 1/2 + (e-1)/2*I) ) * aV
    #
    npts2 = tslices[1] - tslices[0] - tslices[3]
    Cdiff = CV = CH = CD = 0.0
    for j in range(npts2):
        Cdiff += data[j, j + tslices[3], 2] - data[j, j + tslices[2], 2]
        CH += data[j, j + tslices[3], 3]
        CV += data[j, j + tslices[3], 4]
        CD += data[j, j + tslices[3], 5]
    Cdiff /= npts2
    CH /= npts2
    CV /= npts2
    CD /= npts2

    # initialize with no IPC or NL
    alphaH = alphaV = alphaD = alpha = beta = 0.0
    da = 1.0
    # dummy initializations; these get over-written before they are used
    I_ = g = 1.0
    Cdiffcorr = 0.0
    iCycle = 0
    nCycle = 100
    while iCycle < nCycle:
        alphaH_old = alphaH
        alphaV_old = alphaV
        alphaD_old = alphaD
        g_old = g  # to track convergence

        # Get combination of I and gain from difference of correlation functions
        tbrack = (
            tslices[3] * (tslices[0] + tslices[3] - ctrl_pars.reset_frame)
            - tslices[2] * (tslices[0] + tslices[2] - ctrl_pars.reset_frame)
            + (npts2 - 1) / 2.0 * (tslices[3] - tslices[2])
        )
        I__g2 = (
            (Cdiff - Cdiffcorr + 4.0 * (1.0 - 8.0 * alpha) * beta * I_**2 / g**2 * tbrack)
            / (tslices[3] - tslices[2])
            / ((1.0 - 4 * alpha - 4 * alphaD) ** 2 + 2 * alphaH**2 + 2 * alphaV**2 + 4 * alphaD**2)
        )

        # Now use slopemed = -2 beta I^2/g, icpt = (I - beta I^2)/g, and I/g^2 to solve for I, beta, and g
        g = (icpt - slopemed / 2.0) / I__g2
        I_ = I__g2 * g**2
        beta = -g * slopemed / 2.0 / I_**2

        # Corrections to horiz. and vert. IPC
        #
        CHcorr = CVcorr = CDcorr = 0.0
        if applyCorr:
            if typeCorr.lower() == "bfe":
                CHcorr = (ipnl[sBFE, sBFE + 1] + ipnl[sBFE, sBFE - 1]) / 2.0 * (I_ / g * tslices[3]) ** 2
                CVcorr = (ipnl[sBFE + 1, sBFE] + ipnl[sBFE - 1, sBFE]) / 2.0 * (I_ / g * tslices[3]) ** 2
                CDcorr = (
                    (
                        ipnl[sBFE + 1, sBFE + 1]
                        + ipnl[sBFE + 1, sBFE - 1]
                        + ipnl[sBFE - 1, sBFE + 1]
                        + ipnl[sBFE - 1, sBFE - 1]
                    )
                    / 4.0
                    * (I_ / g * tslices[3]) ** 2
                )
                Cdiffcorr = ipnl[sBFE, sBFE] * (I_ / g) ** 2 * (tslices[3] ** 2 - tslices[2] ** 2)
            if typeCorr.lower() == "nlipc":
                CHcorr = (
                    (ipnl[sBFE, sBFE + 1] + ipnl[sBFE, sBFE - 1])
                    / 2.0
                    * (I_ / g) ** 2
                    * tslices[3]
                    * (tslices[0] + tslices[3] + (npts2 - 1) * 0.5)
                    * 2
                )
                CVcorr = (
                    (ipnl[sBFE + 1, sBFE] + ipnl[sBFE - 1, sBFE])
                    / 2.0
                    * (I_ / g) ** 2
                    * tslices[3]
                    * (tslices[0] + tslices[3] + (npts2 - 1) * 0.5)
                    * 2
                )
                CDcorr = (
                    (
                        ipnl[sBFE + 1, sBFE + 1]
                        + ipnl[sBFE + 1, sBFE - 1]
                        + ipnl[sBFE - 1, sBFE + 1]
                        + ipnl[sBFE - 1, sBFE - 1]
                    )
                    / 4.0
                    * (I_ / g) ** 2
                    * tslices[3]
                    * (tslices[0] + tslices[3] + (npts2 - 1) * 0.5)
                    * 2
                )
                Cdiffcorr = (
                    ipnl[sBFE, sBFE]
                    * (I_ / g) ** 2
                    * (
                        (tslices[0] + tslices[3]) * tslices[3]
                        - (tslices[0] + tslices[2]) * tslices[2]
                        + (tslices[3] - tslices[2]) * (npts2 - 1) * 0.5
                    )
                )

            # apply corrections from ftsolve
            if ctrl_pars.fullnl and typeCorr.lower() == "bfe":
                beta_cm = beta
                if ctrl_pars.use_allorder:
                    beta_cm = -addInfo[2] / g ** np.linspace(1, swi.p - 1, num=swi.p - 1)
                if Test_SubBeta:
                    beta_cm = beta
                t0 = tslices[0] - ctrl_pars.reset_frame
                CF_BigStep = solve_corr_many(
                    ipnl,
                    21,
                    I_,
                    g,
                    beta_cm,
                    0.0,
                    [t0, t0 + tslices[3], t0, t0 + tslices[3], npts2],
                    [alphaV, alphaH, alphaD],
                    [0.0, 0.0, 0.0],
                    sBFE,
                )
                CF_SmallStep = solve_corr_many(
                    ipnl,
                    21,
                    I_,
                    g,
                    beta_cm,
                    0.0,
                    [t0, t0 + tslices[2], t0, t0 + tslices[2], npts2],
                    [alphaV, alphaH, alphaD],
                    [0.0, 0.0, 0.0],
                    sBFE,
                )
                Cdiffcorr = (
                    CF_BigStep[sBFE, sBFE]
                    - CF_SmallStep[sBFE, sBFE]
                    - (
                        I_
                        / g**2
                        * ((1 - 4 * alpha - 4 * alphaD) ** 2 + 2 * alphaH**2 + 2 * alphaV**2 + 4 * alphaD**2)
                        * (tslices[3] - tslices[2])
                        - 4 * (1 - 8 * alpha) * beta * I_**2 / g**2 * tbrack
                    )
                )
                ad3 = tslices[0] + tslices[3] - ctrl_pars.reset_frame
                CHcorr = (CF_BigStep[sBFE, sBFE + 1] + CF_BigStep[sBFE, sBFE - 1]) / 2.0 - (
                    2.0
                    * I_
                    / g**2
                    * tslices[3]
                    * (1.0 - 4 * alpha - 4 * alphaD - 4 * beta * (I_ * ad3 + (npts2 - 1) * 0.5 * I_ + 0.5))
                    * alphaH
                    + 4.0 * I_ / g**2 * tslices[3] * alphaV * alphaD
                )
                CVcorr = (CF_BigStep[sBFE + 1, sBFE] + CF_BigStep[sBFE - 1, sBFE]) / 2.0 - (
                    2.0
                    * I_
                    / g**2
                    * tslices[3]
                    * (1.0 - 4 * alpha - 4 * alphaD - 4 * beta * (I_ * ad3 + (npts2 - 1) * 0.5 * I_ + 0.5))
                    * alphaV
                    + 4.0 * I_ / g**2 * tslices[3] * alphaH * alphaD
                )
                CDcorr = (
                    CF_BigStep[sBFE + 1, sBFE + 1]
                    + CF_BigStep[sBFE - 1, sBFE + 1]
                    + CF_BigStep[sBFE + 1, sBFE - 1]
                    + CF_BigStep[sBFE - 1, sBFE - 1]
                ) / 4.0 - (
                    2.0 * I_ / g**2 * tslices[3] * (1.0 - 4 * alpha - 4 * alphaD) * alphaD
                    + 2.0 * I_ / g**2 * tslices[3] * alphaH * alphaV
                )

        factor = (
            2.0
            * I__g2
            * tslices[3]
            * (
                1.0
                - 4.0 * alpha
                - 4.0 * alphaD
                - 4.0
                * beta
                * (I_ * (tslices[0] + tslices[3] - ctrl_pars.reset_frame + (npts2 - 1.0) / 2.0) + 0.5)
            )
        )
        factor_raw = 2.0 * I__g2 * tslices[3]
        alphaH = (CH - CHcorr - 2.0 * alphaV * alphaD * factor_raw) / factor
        alphaV = (CV - CVcorr - 2.0 * alphaH * alphaD * factor_raw) / factor
        alphaD = ((CD - CDcorr) / factor_raw - alphaH * alphaV) / (1.0 - 4.0 * alpha - 4.0 * alphaD)
        alpha = (alphaH + alphaV) / 2.0
        da = np.abs(alphaH_old - alphaH) + np.abs(alphaV_old - alphaV) + np.abs(alphaD_old - alphaD)
        dg = np.abs(g_old - g)
        iCycle += 1
        if iCycle < nCycle - 2 and da < 1e-8 and dg < 1e-8:
            iCycle = nCycle - 2  # fast exit from loop

    return [1, g, alphaH, alphaV, beta, I_, alphaD, da]


def bfe(region_cube, tslices, basicinfo, ctrl_pars_bfe, swi, verbose):
    """
    Routines to compute the BFE coefficients.

    Parameters
    ----------
    region_cube : np.array
        4D array of the region of interest. shape: (num_files+1, nt, dy, dx),
        where num_files is the number of files, nt is the number of time slices, and
        (dy, dx) is the shape of the region on the SCA used. The last slice, ``region_cube[-1,:,:,:]``,
        is the mask (1 for good, 0 for bad).
    tslices : list of int
        A list of at least 4 time slices. The first two and last two ("a", "b", "c", and "d") are
        used for the BFE determination.
    basicinfo : list
        Output from :func:`basic` (inclides gains, IPC, and non-linearity).
    ctrl_pars_bfe : class
        Parameters to control BFE determination; see Notes.
    swi : class
        Column table.
    verbose : bool
        Whether to talk a lot.

    Returns
    -------
    np.array
        The BFE kernel, shape (2s+1, 2s+1), Antilogus coefficients in inverse electrons.

    Notes
    -----
    The (optional) attributes in `ctrl_pars_bfe` are::

      - ``epsilon`` : float
        Fraction of data to cut in computing correlation coefficients (default to 0.01).

      - ``treset`` : int or float
        Reset frame (default to 0). Fractional values are possible.

      - ``BSub`` : bool
        Perform baseline subtraction? (default to True)

      - ``vis`` : bool
        Does this have visible light information (for quantum yield)? (default to False)

      - ``Phi`` : np.array
        2D quantum yield + charge diffusion kernel (only used if `ctrl_pars_bfe`.vis is True).
        This is omega/(1+omega)*p2, where 1+omega is the quantum yield (omega is the probability
        of getting 2 charges from 1 photon) and p2[s+dy,s+dx] is the pairwise probability that two charges
        generated at the same point land in the pixel (dx,dy). Here the kernel has shape (2s+1, 2s+1),
        centered at zero (so it is symmetric).

    """

    N = 21  # <-- size for ftsolve

    # Extract parameters from basicinfo
    gain = basicinfo[swi.g]
    aH = basicinfo[swi.alphaH]
    aV = basicinfo[swi.alphaV]
    beta = basicinfo[swi.beta]
    I_ = basicinfo[swi.I]

    # Extract basic parameters
    num_files = region_cube.shape[0] - 1
    nt = region_cube.shape[1]
    dy = region_cube.shape[2]
    dx = region_cube.shape[3]
    # npix = dx * dy
    if nt != len(tslices):
        raise RuntimeError("Error in basic: incompatible number of time slices")
    if verbose:
        print("nfiles = ", num_files, ", ntimes = ", nt, ", dx,dy=", dx, dy)
    treset = 0
    if hasattr(ctrl_pars_bfe, "treset"):
        treset = ctrl_pars_bfe.treset

    # for visible flats
    hasvis = False
    if hasattr(ctrl_pars_bfe, "vis") and ctrl_pars_bfe.vis:
        hasvis = True
        normPhi = np.sum(ctrl_pars_bfe.Phi)  # this is omega/(1+omega)
        omega = normPhi / (1 - normPhi)
        p2 = np.zeros_like(ctrl_pars_bfe.Phi)
        if np.abs(normPhi) > 1e-49:
            p2 = ctrl_pars_bfe.Phi / normPhi  # this prevents an exception if omega=0
        p2 = pad_to_N(p2, N)  # still centered

    # BFE kernel size:
    # sBFE = range; fsBFE = full size
    sBFE = swi.s
    fsBFE = 2 * sBFE + 1
    sBFE_out = sBFE
    fsBFE_out = fsBFE

    # replace beta with a scalar value if necessary
    # note beta[0] is now 2nd order coef (in DN^-1) is to be converted to beta (in e^-1) and has opposite sign
    if ctrl_pars_bfe.fullnl and ctrl_pars_bfe.use_allorder:
        beta = -beta / gain ** np.linspace(1, swi.p - 1, num=swi.p - 1)
    if ctrl_pars_bfe.fullnl and ctrl_pars_bfe.use_allorder and Test_SubBeta:
        beta = beta[0]

    # Baseline subtraction -- requires bigger box
    BSub = True
    if hasattr(ctrl_pars_bfe, "BSub"):
        BSub = ctrl_pars_bfe.BSub
    if BSub:
        sBFE = max(sBFE_out, 10)
        fsBFE = 2 * sBFE + 1
        pad = 5  # Number of pixels in corr. fcn. to take for the baseline on each side in each row

    # Cut fraction and correction
    epsilon = 0.01
    if hasattr(ctrl_pars_bfe, "epsilon"):
        epsilon = ctrl_pars_bfe.epsilon
    xi = scipy.stats.norm.ppf(1 - epsilon)
    cov_clip_corr = (1.0 - np.sqrt(2.0 / np.pi) * xi * np.exp(-xi * xi / 2.0) / (1.0 - 2.0 * epsilon)) ** 2

    # Build the two slices to correlate
    box1 = region_cube[0:num_files, 0, :, :] - region_cube[0:num_files, 1, :, :]
    box3 = region_cube[0:num_files, -2, :, :] - region_cube[0:num_files, -1, :, :]
    corr_mask = region_cube[-1, 0, :, :]

    # setup for BFE kernel
    numBFE = np.zeros((fsBFE, fsBFE))
    denBFE = np.zeros((fsBFE, fsBFE))

    # Loop over the flat pairs we are going to use
    for if1 in range(1, num_files):
        for if2 in range(if1):
            # Build slices and mask
            slice_ab = box1[if1, :, :] - box1[if2, :, :]
            slice_cd = box3[if1, :, :] - box3[if2, :, :]
            ab_min = pyIRC_percentile(slice_ab, corr_mask, 100 * epsilon)
            ab_max = pyIRC_percentile(slice_ab, corr_mask, 100 * (1 - epsilon))
            cd_min = pyIRC_percentile(slice_cd, corr_mask, 100 * epsilon)
            cd_max = pyIRC_percentile(slice_cd, corr_mask, 100 * (1 - epsilon))
            this_file_mask_ab = np.where(np.logical_and(slice_ab > ab_min, slice_ab < ab_max), corr_mask, 0)
            this_file_mask_cd = np.where(np.logical_and(slice_cd > cd_min, slice_cd < cd_max), corr_mask, 0)
            if verbose:
                print(
                    if1,
                    if2,
                    slice_ab.shape,
                    slice_cd.shape,
                    np.sum(this_file_mask_ab),
                    np.sum(this_file_mask_cd),
                )

            # Mean subtraction
            slice_ab -= pyIRC_mean(slice_ab, this_file_mask_ab)
            slice_cd -= pyIRC_mean(slice_cd, this_file_mask_cd)
            # Set masked values to zero
            slice_ab *= this_file_mask_ab
            slice_cd *= this_file_mask_cd

            # Now get the correlation function ...
            # format is: numerator and denominator of C_{abcd}(2*sBFE-i,2*sBFE-j)
            for j in range(fsBFE):
                for i in range(fsBFE):
                    abminX = 0
                    abmaxX = dx
                    abminY = 0
                    abmaxY = dy
                    if i >= sBFE:
                        abmaxX += sBFE - i
                    else:
                        abminX += sBFE - i
                    if j >= sBFE:
                        abmaxY += sBFE - j
                    else:
                        abminY += sBFE - j
                    cdminX = abminX + i - sBFE
                    cdmaxX = abmaxX + i - sBFE
                    cdminY = abminY + j - sBFE
                    cdmaxY = abmaxY + j - sBFE

                    # Add up contributions to the correlation function
                    denBFE[j, i] += np.sum(
                        this_file_mask_ab[abminY:abmaxY, abminX:abmaxX]
                        * this_file_mask_cd[cdminY:cdmaxY, cdminX:cdmaxX]
                    )
                    numBFE[j, i] += (
                        np.sum(
                            slice_ab[abminY:abmaxY, abminX:abmaxX] * slice_cd[cdminY:cdmaxY, cdminX:cdmaxX]
                        )
                        / 2.0
                    )
                    # division by 2 since differencing two images doubles the answer

    BFEK = numBFE / (denBFE + 1e-99)
    BFEK *= gain**2 / (I_**2 * (tslices[1] - tslices[0]) * (tslices[-1] - tslices[-2]) * cov_clip_corr)

    # Baseline subtraction
    if BSub:
        for j in range(fsBFE):
            rowBL = (np.mean(BFEK[j, 0:pad]) + np.mean(BFEK[j, -pad:])) / 2.0
            BFEK[j, :] -= rowBL

    # Implement cr_converge.
    if ctrl_pars_bfe.fullnl:
        avals = [basicinfo[swi.alphaV], basicinfo[swi.alphaH], basicinfo[swi.alphaD]]
        avals_nl = [0, 0, 0]
        sigma_a = 0
        tol = 1.0e-11  # Pick a tolerance below which the two Crs are considered equal
        fsBFE_out = 2 * sBFE_out + 1
        observed_Cr = BFEK[sBFE - sBFE_out : sBFE + sBFE_out + 1, sBFE - sBFE_out : sBFE + sBFE_out + 1]
        BFEK_model = np.zeros((fsBFE_out, fsBFE_out)) + 1e-15
        element_diff = 10
        iters = 0
        while element_diff > tol and iters <= 100:
            # Note: solve_corr takes centered things, decenters/calculates internally
            if hasvis:
                theory_Cr = solve_corr_vis(
                    BFEK_model,
                    N,
                    I_,
                    gain,
                    beta,
                    sigma_a,
                    [t - treset for t in tslices],
                    avals,
                    avals_nl,
                    sBFE_out,
                    omega,
                    p2,
                ) * ((gain**2) / (I_**2 * (tslices[1] - tslices[0]) * (tslices[-1] - tslices[-2])))
            else:
                theory_Cr = solve_corr(
                    BFEK_model, N, I_, gain, beta, sigma_a, [t - treset for t in tslices], avals, avals_nl
                ) * ((gain**2) / (I_**2 * (tslices[1] - tslices[0]) * (tslices[-1] - tslices[-2])))
            if np.isnan(theory_Cr).any():
                warnings.warn("BFE loop diverged, generated NaN")
                return np.zeros((fsBFE_out, fsBFE_out)) + np.nan
            difference = theory_Cr - observed_Cr
            element_diff = np.amax(abs(difference))
            BFEK_model -= difference[::-1, ::-1]
            iters += 1
            if verbose:
                print(iter, BFEK_model)
            if iters > 99:
                warnings.warn("WARNING: NL loop has iterated 100 times")
                return np.zeros((fsBFE_out, fsBFE_out)) + np.nan
        return BFEK_model

    else:
        # Corrections for classical non-linearity
        BFEK[sBFE, sBFE] += 2 * (1 - 4 * (aH + aV)) * beta
        if sBFE >= 1:
            BFEK[sBFE, sBFE + 1] += 4 * aH * beta
            BFEK[sBFE, sBFE - 1] += 4 * aH * beta
            BFEK[sBFE + 1, sBFE] += 4 * aV * beta
            BFEK[sBFE - 1, sBFE] += 4 * aV * beta
        return BFEK[sBFE - sBFE_out : sBFE + sBFE_out + 1, sBFE - sBFE_out : sBFE + sBFE_out + 1]


def hotpix(darkfiles, formatpars, tslices, pars, verbose):
    """
    Selects hot pixels.

    Parameters
    ----------
    darkfiles : list of str
        A list of the filenames of the dark exposures.
    formatpars : int
        The format code.
    tslices : list of int
        The time slices to read (the first slice is 1).
    pars : np.array or array-like
        Parameters controlling the hot pixel selection. These should be
        ``[Smin, Smax, stability, f_isolation]`` (see Notes for detailed meaning).
    verbose : bool
        Whether to print lots of information.

    Returns
    -------
    row, col : np.array of int
        The row and column values of the selected hot pixels. The two arrays have the same
        length; the ``i``th hot pixel is at ``[row[i],col[i]]``.

    Notes
    -----
    The hot pixels in the array must meet the following criteria:

    - The apparent brightness in time slices up through tslices[-1] is assessed.
      Pixels are required to have a dark signal in DN between ``Smin`` and ``Smax``.

    - Repeatable to within a top-to-bottom error of ``stability`` as a fraction of the
      maximum signal (e.g. 0.1 for 10% repeatability).

    - Isolation: if ``f_isolation``>0, rejects pixels with neighbors that are at least this many
      times as bright as this pixel itself (e.g. 0.1 for 10% isolation).

    """

    # Build array for the dark cube
    ndarks = len(darkfiles)
    N = get_nside(formatpars)
    cube = np.zeros((ndarks, N, N))
    for f in range(ndarks):
        CDS = load_segment(darkfiles[f], formatpars, [0, N, 0, N], [1, tslices[-1]], False)
        cube[f, :, :] = CDS[0, :, :] - CDS[1, :, :]

    # Extract information on the pixels
    this_hot = np.zeros((N, N))
    ave_cube = np.mean(cube, axis=0)
    d_cube = np.max(cube, axis=0) - np.min(cube, axis=0)
    if verbose:
        print("time slices for hot pixel analysis ->", tslices)
        print(ave_cube)
        print("->", ave_cube.shape)
        print(d_cube)
        print("->", d_cube.shape)

    this_hot = np.where(np.logical_and(ave_cube >= pars[0], ave_cube <= pars[1]), 1, 0)

    # Isolation cut
    if verbose:
        print("Start with", np.sum(this_hot), "pixels before isolation cut")
    if pars[3] > 0:
        C = 2
        M = np.ones((2 * C + 1, 2 * C + 1))
        M[C, C] = 0
        isolation_mask = scipy.ndimage.maximum_filter(ave_cube, footprint=M, mode="constant", cval=0)
        # Also avoid pixels that border on reference pixels
        this_hot[: 4 + C, :] = 0
        this_hot[-(4 + C) :, :] = 0
        this_hot[:, : 4 + C] = 0
        this_hot[:, -(4 + C) :] = 0
        this_hot *= np.where(isolation_mask <= pars[3] * ave_cube, 1, 0)

    if verbose:
        print("Start with", np.sum(this_hot), "pixels")
    for t in tslices[1:]:
        for f in range(ndarks):
            CDS = load_segment(darkfiles[f], formatpars, [0, N, 0, N], [1, t], False)
            cube[f, :, :] = CDS[0, :, :] - CDS[1, :, :]
        d_cube = np.max(cube, axis=0) - np.min(cube, axis=0)
        this_hot *= np.where(d_cube <= pars[2] * ave_cube, 1, 0)
    if verbose:
        print(np.sum(this_hot))

    return np.where(this_hot > 0)


def hotpix_ipc(y, x, darkfiles, formatpars, tslices, pars, verbose):
    """
    Return IPC data from a list of hot pixels.

    Parameters
    ----------
    y, x : np.array of int
        Tables of hot pixel coordinates to use (probably selected from :func:`hotpix`).
    darkfiles : list of str
        List of dark files to use.
    formatpars : int
        Format code for dark files.
    tslices : list of int
        List of time slices to report.
    pars : list, [np.array, bool]
        Parameters to control data selection. If not empty, `pars`[0] is a map of the
        non-linearity times gain (units: 1/DN) and `pars`[1] is a boolean for whether
        the non-linearity should be referenced to a median stack of initial image.
    verbose : bool
        Whether to talk a lot.

    Returns
    -------
    np.array
        Data cube with the signal in DN from the hot pixels. The shape is (npix, nt, 10), where npix is the
        number of hot pixels (length of `y` and `x`); and nt is the number of time stamps
        (length of `tslices`). The last index indicates the position: 0 = that pixel; 1 = right; 2 = up;
        3 = left; 4 = downl; 5 = upper-right; 6 = upper-left; 7 = lower-left; 8 = lower-right;
        9 = background (from the next 5x5-3x3=13 pixels out).

    """

    # Build array for the dark cube
    ndarks = len(darkfiles)
    N = get_nside(formatpars)
    cube = np.zeros((ndarks, N, N))

    nt = len(tslices)
    npix = len(x)
    data = np.zeros((npix, nt, 10))

    # offset table
    dx = [0, 1, 0, -1, 0, 1, -1, -1, 1]
    dy = [0, 0, 1, 0, -1, 1, 1, -1, -1]

    # Perform nonlinearity correction?
    do_nonlin = False
    if len(pars) >= 1 and (type(pars[0]) is np.ndarray):
        do_nonlin = True
        m = pars[0]
        beta_gain = np.zeros((N, N))
        (ny1, nx1) = np.shape(m)
        kx1 = N // nx1
        ky1 = N // ny1
        for i in range(nx1):
            for j in range(ny1):
                beta_gain[ky1 * j : ky1 * (j + 1), kx1 * i : kx1 * (i + 1)] = m[j, i]
        # now beta_gain is an NxN map of beta*gain
        if verbose:
            print("beta*gain =", beta_gain)
    # baseline for NL correction is median image?
    medbaseline_nonlin = False
    if len(pars) >= 2:
        medbaseline_nonlin = pars[1]

    # background mask
    bkmask = np.ones((5, 5))
    bkmask[1:4, 1:4] = 0.0
    fourmask = False
    if fourmask:
        bkmask[:, :] = 0.0
        bkmask[2, 0] = bkmask[2, 4] = bkmask[0, 2] = bkmask[4, 2] = 1.0
    if verbose:
        print("bkmask =", bkmask)
    # 16 ones and 9 zeros

    # now make data cube
    for jt in range(nt):
        for f in range(ndarks):
            CDS = load_segment(darkfiles[f], formatpars, [0, N, 0, N], [1, tslices[jt]], False)
            cube[f, :, :] = CDS[0, :, :] - CDS[1, :, :]
            if do_nonlin:
                # non-linearity correction, if turned on
                cube_corr = cube[f, :, :]
                if medbaseline_nonlin:
                    cube_corr = (
                        2 * scipy.ndimage.median_filter(CDS[0, :, :], size=3) - CDS[0, :, :] - CDS[1, :, :]
                    )
                cube[f, :, :] = cube[f, :, :] * (1.0 + beta_gain * cube_corr)
        medframe = np.median(cube, axis=0)
        if verbose:
            print("med", np.shape(medframe), jt)
        for jpix in range(npix):
            for jpos in range(9):
                x_ = x[jpix] + dx[jpos]
                y_ = y[jpix] + dy[jpos]
                data[jpix, jt, jpos] = medframe[y_, x_]
            data[jpix, jt, 9] = (
                25.0 / 16.0 * np.mean(bkmask * medframe[y[jpix] - 2 : y[jpix] + 3, x[jpix] - 2 : x[jpix] + 3])
            )
            if fourmask:
                data[jpix, jt, 9] *= 16.0 / 4.0

    return data


def slidemed_percentile(x, y, p, mrange=[-1, 1], niter=64, pivot="pos"):  # list OK # noqa: B006
    """
    Sliding median function.

    Takes in points specified by `x` and `y` arrays, and percentile `p`,
    and returns slope m such that p% of the data are below the line y = m*x.

    (If all values in `x` are positive, this is simply a percentile of `y`/`x`.
    But this version is not biased when the noise causes a few values to fluctuate
    negative.)

    Parameters
    ----------
    x, y : np.array of float
        Vectors of the same length.
    p : float
        Desired percentile.
    mrange : list, optional
        Slope range to search, length 2 ([mmin, mmax]).
    niter : int, optional
        Number of bisections to perform.
    pivot : str, optional
        The pivot is not used but is here for forward compatibility; right now it assumes
        that the pivot point of the distribution is at x>0 (hence ordering in the bisection).
        In a future release if we need to change this the functionality is there.

    Returns
    -------
    float
        The slope of the line meeting that percentile criterion.

    """

    m1 = mrange[0]
    m2 = mrange[1]

    for _ in range(niter):
        m = (m1 + m2) / 2.0
        if np.nanpercentile(y - m * x, p) > 0:
            m1 = m
        else:
            m2 = m
    return m


def get_vmin_vmax(mydata, qext):
    """
    Generates min and max range for a color bar based on inter-quartile range.

    Parameters
    ----------
    mydata : np.array
        The data to consider
    qext : float
        Number of interquartile ranges to extend (should be >=0, with 0
        corresponding to 25th through 75th percentile).

    Returns
    -------
    float, float
        The minimum and maximum of the scale.

    """

    Q1 = np.nanpercentile(mydata, 25)
    Q2 = np.nanpercentile(mydata, 75)
    return Q1 - (Q2 - Q1) * qext, Q2 + (Q2 - Q1) * qext

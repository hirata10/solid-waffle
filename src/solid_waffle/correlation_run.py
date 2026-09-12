"""
Routines to run the infrared flat correlations.

Classes
-------
EmptyClass
    Blank, can add attributes later.
Config
    Extracts configuration data from a multiline string.

Functions
---------
run_ir_all
  Runs the IR characterization.
run_vis_all
  Runs the visible characterization.

"""

import re
import shutil
import sys
import time

import matplotlib
import numpy as np

from . import ftsolve, pyirc

matplotlib.use("Agg")
import matplotlib.pyplot as plt

plt.switch_backend("agg")


class EmptyClass:
    """Blank, can add attributes later."""

    pass


class Config:
    """
    Extracts configuration data from a multiline string.

    Some attributes may not be present depending on the configuration.

    Parameters
    ----------
    cfg : str or list
        The configuration (as a single string or list of strings).
    visible_run : bool, optional
        Is this configured to run visible light correlations?
    verbose : bool, optional
        Whether to specify lots of information when reading the configuration.

    Attributes
    ----------
    visible_run : bool
        Is this configured to run visible light correlations?
    outstem : str
        The stem for output files.
    use_cmap : str
        The color bar range to use.
    mydet : str
        Name of the detector.
    lightfiles : list of str
        List of the flat field file names.
    darkfiles : list of str
        List of the dark exposure file names, same length as `lightfiles`.
    formatpars : int
        The format type for these images.
    nx, ny : int
        The number of bins to break the detector array into on the x and y axes.
    tslices : list of int
        Time slices to use for BFE correlations. Length 4, ascending.
    tslicesM2a, tslicesM2b, tslicesM3 : list of int
        Time slices to use for alternative BFE tests. Length 4, ascending.
    fullref : bool
        Apply full reference pixel corrrection?
    sensitivity_spread_cut : float
        The fractional deviation from the median flat at which to cut pixels.
    critfrac : float
        Fraction of the pixels that need to be good to keep a superpixel.
    mychar : str
        Type of characterization to use (options are "Basic" and "Advanced").
    hotpix : bool
        Whether to do the hot pixel analysis.
    hotpix_ADU_range : list of float
        Hot pixel selection range (length 4: ``[Smin, Smax, stability, f_isolation]``).
    ref_for_hotpix_is_autocorr : bool
        Plot hot pixel results relative to autocorrelation results?
    hotpix_logtspace : bool
        Space hot pixel samples logarithmically?
    hotpix_slidemed : bool
        Use sliding median method to get hot pixel IPC measuremenent?
    s_bfe : int
        Radius of BFE kernel.
    p_order : int
        If >0, builds non-linearity polynomial table through this order.
    basicpar : class
        Basic characterization parameters.
    bfepar : class
        BFE characterization parameters.
    narrowfig : bool
        Output figures in narrow format?
    maskX, maskY : list of int
        Super-pixel regions to mask out (can be empty).
    tchar1, tchar2 : int
        Time steps for advanced characterization.
    ncycle : int
        Number of iterations for advanced characterization.
    ipnltype : str
        Inter-pixel non-linearity model: 'bfe' or 'nlipc'.
    nlfit_ts, nlfit_te : int
        Range of time stamps for fitting the non-linearity curve.
    NTMAX : int
        Maximum time slice allowed.
    swi : class
        Column information.
    N : int
        Size of the SCA.
    dx, dy : int
        Size of super-pixels.
    npix : int
        Number of pixels in a super-pixel.
    full_info : np.array
        Calibration results, shape = (`ny`, `nx`, ...).
    is_good : np.array
        Good super-pixel map, shape = (`ny`, `nx`).
    lightref, darkref : np.array
        Reference pixel corrections for the flats and darks.
    mean_full_info, std_full_info : np.array
        Array mean and standard deviation of `full_info`.
    nlfit, nlder : np.array
        Array of the ramp fit and derivative, shape = (nt, `ny`, `nx`).
    used_2a, used_2b, used_3 : bool
        Whether the alternative methods were implemented.
    ntM2a, ntM2b, ntM3 : int
        Number of time stamps in each alternative method.
    Method2a_slopes, Method2b_slopes, Method3_slopes : np.array
        Alternative gain/IPC calculation slopes, shape = (`ny`, `nx`).
    Method2a_vals, Method2b_vals, Method3_vals : np.array
        Alternative gain/IPC calculation values, shape = (`ny`, `nx`, ...).
    tfmin, tfmax, tfminB, tfmaxB, tfmin3, tfmax3 : int
        Time stamps for Methods 2a, 2b, and 3.
    slope_2a_BFE, slope_2a_NLIPC : float
        Predicted slopes for the different alternative methods for different sources of IPNL (Method 2a).
    slope_2b_BFE, slope_2b_NLIPC : float
        Predicted slopes for the different alternative methods for different sources of IPNL (Method 2b).
    slope_3_beta, slope_3_BFE, slope_3_NLIPC : float
        Predicted slopes for the different alternative methods for different sources of IPNL (Method 3).
    PV2a, PV2b, PV3 : float
        Peak-to-valley errors of alternative method fits.
    hotX, hotY : np.array of int
        Coordinates of the hot pixels; same length.
    htsteps : list of int
        Time steps for the hot pixel analysis.
    hotpix_signal : np.array
        Signal levels of the hot pixels, shape = (len(`hotX`), len(`htsteps`))
    hotpix_alpha : np.array
        IPC estimated from the hot pixels, shape = (len(`hotX`), len(`htsteps`))
    ipcmed_x, ipcmed_y, ipcmed_yerr : np.array
        Signal, IPC, and IPC error for the hot pixel signal-dependent IPC method. 1D arrays.
    grid_alphaCorr, grid_alphaCorrErr, grid_alphaHot, grid_alphaHotErr : np.array
        The correlation- and hot pixel-based IPC in each of the 16 hexadecants of the detector array,
        and their uncertainties.
    ts_vis, te_vis, tchar1_vis, tchar2_vis : int
        The time stamps for visible characterization.
    has_visbfe : bool
        Visible BFE enabled?
    tslices_visbfe : list of int
        Time slices for visible BFE analysis.
    copy_ir_bfe : bool
        Assume IR BFE applies to visible? (May improve S/N.)
    vis_out_data : np.array
        The visible characterization data, shape = (`ny`, `nx`, 56)
    vis_col : dict
        Column mapping for the visible characterization.

    Methods
    -------
    __init__
        Constructor.
    fit_parameters
        Build general parameters (gain, IPC, NL) and BFE Method 1.
    generate_nonlinearity
        Generates non-linearity data.
    write_basic_figure
        Makes a PDF figure of the SCA parameter maps.
    alt_methods
        Methods 2a, 2b, and 3 trend computation.
    method_23_plot
      Method 2 and 3 characterization Multi-panel figure showing basic characterization.
    text_output
      Generate a text summary.
    hotpix_analysis
      Hot pixel analysis.
    hotpix_plots
      Makes the hot pixel plots.
    compute_vis_quantities
      Computations for the visible light characterization.
    vis_plots
      Make plots for the visible light characterization.

    """

    def __init__(self, cfg, visible_run=False, verbose=False):
        self.visible_run = visible_run

        self.outstem = "default_output"
        self.use_cmap = "gnuplot"

        self.mydet = ""
        self.lightfiles = []
        self.darkfiles = []
        if self.visible_run:
            self.vislightfiles = []
            self.visdarkfiles = []
        self.formatpars = 1
        self.nx = self.ny = 32
        self.tslices = [3, 11, 13, 21]
        self.tslicesM2a = []
        self.tslicesM2b = []
        self.tslicesM3 = []
        self.fullref = True
        self.sensitivity_spread_cut = 0.1
        self.critfrac = 0.75
        self.mychar = "Basic"
        self.hotpix = False
        self.ref_for_hotpix_is_autocorr = False
        self.hotpix_logtspace = False
        self.hotpix_slidemed = False

        self.swi = pyirc.IndexDictionary(0)  # initialize column table; will add more

        # order parameters
        self.s_bfe = 2  # order of BFE parameters
        self.p_order = (
            0  # non-linearity polynomial table coefficients (table at end goes through order p_order)
        )
        # set to zero to turn this off

        # Parameters for basic characterization
        basicpar = EmptyClass()
        basicpar.epsilon = 0.01
        basicpar.subtr_corr = True
        basicpar.noise_corr = True
        basicpar.reset_frame = 1
        basicpar.subtr_href = True
        basicpar.full_corr = True
        basicpar.leadtrailSub = False
        basicpar.g_ptile = 75.0
        basicpar.fullnl = False
        basicpar.use_allorder = False
        if self.visible_run:
            basicpar.vis_med_correct = False
        self.basicpar = basicpar

        # Parameters for BFE
        bfepar = EmptyClass()
        bfepar.epsilon = 0.01
        bfepar.treset = basicpar.reset_frame
        bfepar.blsub = True
        bfepar.fullnl = False
        if self.visible_run:
            bfepar.vis = True
            self.copy_ir_bfe = False
        self.bfepar = bfepar

        # Separate parameters for visible BFE?
        if self.visible_run:
            self.has_visbfe = False

        # Plotting parameters
        self.narrowfig = False

        # Read in information
        if isinstance(cfg, list):
            cfg = "".join(cfg)  # convert to string
        content = cfg.splitlines()
        is_in_light = is_in_dark = False
        if self.visible_run:
            is_in_vislight = is_in_visdark = False
        self.maskX = []  # list of regions to mask
        self.maskY = []
        for line in content:
            # Cancellations
            m = re.search(r"^[A-Z]+\:", line)
            if m:
                is_in_light = is_in_dark = False
                if self.visible_run:
                    is_in_vislight = is_in_visdark = False

            # Searches for files -- must be first given the structure of this script!
            if is_in_light:
                m = re.search(r"^\s*(\S.*)$", line)
                if m:
                    self.lightfiles += [m.group(1)]
            if is_in_dark:
                m = re.search(r"^\s*(\S.*)$", line)
                if m:
                    self.darkfiles += [m.group(1)]
            if self.visible_run:
                if is_in_vislight:
                    m = re.search(r"^\s*(\S.*)$", line)
                    if m:
                        self.vislightfiles += [m.group(1)]
                if is_in_visdark:
                    m = re.search(r"^\s*(\S.*)$", line)
                    if m:
                        self.visdarkfiles += [m.group(1)]

            # -- Keywords go below here --

            # Search for outputs
            m = re.search(r"^OUTPUT\:\s*(\S*)", line)
            if m:
                self.outstem = m.group(1)
            # Search for input files
            m = re.search(r"^LIGHT\:", line)
            if m:
                is_in_light = True
            m = re.search(r"^DARK\:", line)
            if m:
                is_in_dark = True
            if self.visible_run:
                m = re.search(r"^VISLIGHT\:", line)
                if m:
                    is_in_vislight = True
                m = re.search(r"^VISDARK\:", line)
                if m:
                    is_in_visdark = True

            # Format
            m = re.search(r"^FORMAT:\s*(\d+)", line)
            if m:
                self.formatpars = int(m.group(1))

            # Bin sizes
            m = re.search(r"^NBIN:\s*(\d+)\s+(\d+)", line)
            if m:
                self.nx = int(m.group(1))
                self.ny = int(m.group(2))

            # Characterization type (Basic or Advanced)
            m = re.search(r"^CHAR:\s*(\S+)", line)
            if m:
                self.mychar = m.group(1)
                if self.mychar.lower() == "advanced":
                    m = re.search(r"^CHAR:\s*(\S+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\S+)", line)
                    if m:
                        self.tchar1 = int(m.group(2))
                        self.tchar2 = int(m.group(3))
                        self.ncycle = int(m.group(4))
                        self.ipnltype = m.group(5)
                    else:
                        print("Error: insufficient arguments: " + line + "\n")
                        exit()

            if self.visible_run:
                # Visible time stamp range
                m = re.search(r"^VISTIME:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", line)
                if m:
                    self.ts_vis = int(m.group(1))
                    self.te_vis = int(m.group(2))
                    self.tchar1_vis = int(m.group(3))
                    self.tchar2_vis = int(m.group(4))
                #
                m = re.search(r"^VISMEDCORR", line)
                if m:
                    self.basicpar.vis_med_correct = True
                #
                # Visible BFE
                m = re.search(r"^VISBFETIME:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", line)
                if m:
                    self.tslices_visbfe = [int(m.group(x)) for x in range(1, 5)]
                    self.has_visbfe = True

                m = re.search(r"^COPYIRBFE", line)
                if m:
                    self.copy_ir_bfe = True

            # Time slices
            m = re.search(r"^TIME:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", line)
            if m:
                self.tslices = [int(m.group(x)) for x in range(1, 5)]
            m = re.search(r"^TIME2A:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", line)
            if m:
                self.tslicesM2a = [int(m.group(x)) for x in range(1, 5)]
            m = re.search(r"^TIME2B:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", line)
            if m:
                self.tslicesM2b = [int(m.group(x)) for x in range(1, 5)]
            m = re.search(r"^TIME3:\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)", line)
            if m:
                self.tslicesM3 = [int(m.group(x)) for x in range(1, 5)]
            #
            # reference time slice
            m = re.search(r"^TIMEREF:\s*(\d+)", line)
            if m:
                self.bfepar.treset = basicpar.reset_frame = float(m.group(1)) if "." in m.group(1) else int(m.group(1))

            # reference pixel subtraction
            m = re.search(r"^REF\s+OFF", line)
            if m:
                self.fullref = False

            # sensitivity spread cut
            m = re.search(r"^SPREAD:\s*(\S+)", line)
            if m:
                self.sensitivity_spread_cut = float(m.group(1))

            # variance parameters
            m = re.search(r"^QUANTILE:\s*(\S+)", line)
            if m:
                self.basicpar.g_ptile = float(m.group(1))
            # correlation parameters
            m = re.search(r"^EPSILON:\s*(\S+)", line)
            if m:
                self.bfepar.epsilon = basicpar.epsilon = float(m.group(1))
            m = re.search(r"^IPCSUB:\s*(\S+)", line)
            if m:
                self.basicpar.leadtrailSub = m.group(1).lower() in ["true", "yes"]

            # Other parameters
            m = re.search(r"^DETECTOR:\s*(\S+)", line)
            if m:
                self.mydet = m.group(1)
            m = re.search(r"^COLOR:\s*(\S+)", line)
            if m:
                self.use_cmap = m.group(1)

            # Classical non-linearity
            m = re.search(r"^NLPOLY:\s*(\S+)\s+(\S+)\s+(\S+)", line)
            if m:
                self.p_order = int(m.group(1))
                self.nlfit_ts = int(m.group(2))
                self.nlfit_te = int(m.group(3))

            m = re.search(r"^FULLNL:\s*(\S+)\s+(\S+)\s+(\S+)", line)
            if m:
                self.basicpar.fullnl = m.group(1).lower() in ["true", "yes"]
                self.bfepar.fullnl = m.group(2).lower() in ["true", "yes"]
                self.basicpar.use_allorder = m.group(3).lower() in ["true", "yes"]

            # Hot pixels (adu min, adu max, cut stability, cut isolation)
            m = re.search(r"^HOTPIX:\s*(\S+)\s+(\S+)\s+(\S+)\s+(\S+)", line)
            if m:
                self.hotpix = True
                self.hotpix_ADU_range = [float(m.group(x)) for x in range(1, 5)]
            #
            # change reference for hot pixels from last point to autocorr
            m = re.search(r"^HOTREF\s+AUTOCORR", line)
            if m:
                self.ref_for_hotpix_is_autocorr = True
            # log spacing for times?
            m = re.search(r"^HOTPIX\s+LOGTSPACE", line)
            if m:
                self.hotpix_logtspace = True
            # sliding median alpha method?
            m = re.search(r"^HOTPIX\s+SLIDEMED", line)
            if m:
                self.hotpix_slidemed = True

            # Mask regions by hand
            m = re.search(r"^MASK:\s*(\d+)\s+(\d+)", line)
            if m:
                self.maskX.append(int(m.group(1)))
                self.maskY.append(int(m.group(2)))

            # Control figures
            m = re.search(r"^NARROWFIG", line)
            if m:
                self.narrowfig = True

        # set up array size parameters
        self.swi.addbfe(self.s_bfe)
        self.swi.addhnl(self.p_order)
        print("Number of output field per superpixel =", self.swi.N)

        # Check number of slices available
        NTMAX = 16384
        for f in self.lightfiles + self.darkfiles:
            nt = pyirc.get_num_slices(self.formatpars, f)
            if nt < NTMAX:
                NTMAX = nt
        self.NTMAX = NTMAX

        # Copy basicpar parameters to bfebar
        self.bfepar.use_allorder = self.basicpar.use_allorder

        if verbose:
            print(f"Output will be directed to {self.outstem:s}*")
            print("Light files:", self.lightfiles)
            print("Dark files:", self.darkfiles)
            if self.visible_run:
                print("Visible light files:", self.vislightfiles)
                print('"Visible" dark files:', self.visdarkfiles)
            print("Time slices:", self.tslices, "max=", self.NTMAX)
            print("Mask regions:", self.maskX, self.maskY)

        # Checks
        if len(self.lightfiles) != len(self.darkfiles) or len(self.lightfiles) < 2:
            raise ValueError(
                f"Failed: {len(self.lightfiles):d} light files and {len(self.darkfiles):d} dark files"
            )

        # Additional parameters Size of a block
        self.N = pyirc.get_nside(self.formatpars)
        # Side lengths
        self.dx = self.N // self.nx
        self.dy = self.N // self.ny
        # Pixels in a block
        self.npix = self.dx * self.dy

        # Initializations
        self.full_info = np.zeros((self.ny, self.nx, self.swi.N))
        self.is_good = np.zeros((self.ny, self.nx))

    def fit_parameters(self, verbose=False):
        """
        Build general parameters (gain, IPC, NL) and BFE Method 1.

        Parameters
        ----------
        verbose : bool
            Whether to talk a lot.

        Returns
        -------
        None

        """

        # Make table of reference pixel corrections for Method 1
        if self.fullref:
            self.lightref = pyirc.ref_array(self.lightfiles, self.formatpars, self.ny, self.tslices, False)
            self.darkref = pyirc.ref_array(self.lightfiles, self.formatpars, self.ny, self.tslices, False)
        else:
            self.lightref = np.zeros((len(self.lightfiles), self.ny, 2 * len(self.tslices) + 1))
            self.darkref = np.zeros((len(self.darkfiles), self.ny, 2 * len(self.tslices) + 1))
        self.basicpar.subtr_href = self.fullref

        if self.p_order > 0:
            # now coefficients for the info table note that in 'abs' mode, the full_info[:,:,0] grid is not
            # actually used, it is just there for consistency of the format I moved this up here since we want
            # to have these coefficients before the main program runs
            nlcubeX, nlfitX, nlderX, pcoefX = pyirc.gen_nl_cube(
                self.lightfiles,
                self.formatpars,
                [self.basicpar.reset_frame, self.nlfit_ts, self.nlfit_te],
                [self.ny, self.nx],
                self.full_info[:, :, 0],
                "abs",
                self.swi,
                False,
            )
            # fill in
            for iy in range(self.ny):
                for ix in range(self.nx):
                    if pcoefX[1, iy, ix] != 0:
                        self.full_info[iy, ix, self.swi.Nbb] = -pcoefX[0, iy, ix] / pcoefX[1, iy, ix]
                        for ooo in range(2, self.swi.p + 1):
                            self.full_info[iy, ix, self.swi.Nbb + ooo - 1] = (
                                pcoefX[ooo, iy, ix] / pcoefX[1, iy, ix] ** ooo
                            )
                    else:
                        self.full_info[iy, ix, pyirc.swi.Nbb] = -1e49  # error code

        # Detector characterization data in a cube (basic characterization + BFE Method 1)
        # Stdout calls are a progress indicator
        #
        if verbose:
            print("Method 1, progress of calculation:")
            sys.stdout.write("|")
            for _ in range(self.ny):
                sys.stdout.write(" ")
            print("| <- 100%")
            sys.stdout.write("|")
        for iy in range(self.ny):
            if verbose:
                sys.stdout.write("*")
                sys.stdout.flush()
            for ix in range(self.nx):
                region_cube = pyirc.pixel_data(
                    self.lightfiles,
                    self.formatpars,
                    [self.dx * ix, self.dx * (ix + 1), self.dy * iy, self.dy * (iy + 1)],
                    self.tslices,
                    [self.sensitivity_spread_cut, True],
                    False,
                )
                dark_cube = pyirc.pixel_data(
                    self.darkfiles,
                    self.formatpars,
                    [self.dx * ix, self.dx * (ix + 1), self.dy * iy, self.dy * (iy + 1)],
                    self.tslices,
                    [self.sensitivity_spread_cut, False],
                    False,
                )
                info = pyirc.basic(
                    region_cube,
                    dark_cube,
                    self.tslices,
                    self.lightref[:, iy, :],
                    self.darkref[:, iy, :],
                    self.basicpar,
                    False,
                )
                if len(info) > 0:
                    self.is_good[iy, ix] = 1
                    thisinfo = info.copy()
                    if self.basicpar.use_allorder:
                        thisinfo[self.swi.beta] = self.full_info[
                            iy, ix, self.swi.Nbb + 1 : self.swi.Nbb + self.swi.p
                        ]
                    if self.mychar.lower() == "advanced":
                        boxpar = [self.dx * ix, self.dx * (ix + 1), self.dy * iy, self.dy * (iy + 1)]
                        tpar_polychar = [self.tslices[0], self.tslices[-1] + 1, self.tchar1, self.tchar2]
                        corrstats_data = pyirc.corrstats(
                            self.lightfiles,
                            self.darkfiles,
                            self.formatpars,
                            boxpar,
                            tpar_polychar + [1],
                            self.sensitivity_spread_cut,
                            self.basicpar,
                        )
                        for _ in range(self.ncycle):
                            bfeCoefs = pyirc.bfe(
                                region_cube, self.tslices, thisinfo, self.bfepar, self.swi, False
                            )
                            if np.isnan(bfeCoefs).any():
                                bfeCoefs = np.zeros((2 * self.swi.s + 1, 2 * self.swi.s + 1))
                                self.is_good[iy, ix] = 0
                            Cdata = pyirc.polychar(
                                self.lightfiles,
                                self.darkfiles,
                                self.formatpars,
                                boxpar,
                                tpar_polychar,
                                self.sensitivity_spread_cut,
                                self.basicpar,
                                [self.ipnltype, bfeCoefs, thisinfo[self.swi.beta]],
                                self.swi,
                                corrstats_data=corrstats_data,
                            )
                            info[self.swi.ind1 : self.swi.ind2] = np.asarray(
                                Cdata[self.swi.indp1 : self.swi.indp2]
                            )
                            thisinfo = info.copy()
                            if self.basicpar.use_allorder:
                                thisinfo[self.swi.beta] = self.full_info[
                                    iy, ix, self.swi.Nbb + 1 : self.swi.Nbb + self.swi.p
                                ]
                    bfeCoefs = pyirc.bfe(region_cube, self.tslices, thisinfo, self.bfepar, self.swi, False)
                    if np.isnan(bfeCoefs).any():
                        bfeCoefs = np.zeros((2 * self.swi.s + 1, 2 * self.swi.s + 1))
                        self.is_good[iy, ix] = 0
                    info += bfeCoefs[0 : 2 * self.swi.s + 1, 0 : 2 * self.swi.s + 1].flatten().tolist()
                else:
                    info = np.zeros(self.swi.Nbb).tolist()

                # save the information, putting in 0's if the super-pixel is not good.
                if len(info) == self.swi.Nbb:
                    self.full_info[iy, ix, : self.swi.Nbb] = np.array(info)
                if info[0] < self.npix * self.critfrac:
                    self.is_good[iy, ix] = 0
                    self.full_info[iy, ix, 1:] = 0  # wipe out this super-pixel

        if verbose:
            print("|")

        # Mask regions
        for mask_index in range(len(self.maskX)):
            ix = self.maskX[mask_index]
            iy = self.maskY[mask_index]
            self.is_good[iy, ix] = 0
            self.full_info[iy, ix, :] = 0  # wipe out this super-pixel

        # if a pixel was set to not good for any other reason
        for iy in range(self.ny):
            for ix in range(self.nx):
                if self.is_good[iy, ix] < 0.5:
                    self.full_info[iy, ix, :] = 0  # wipe out this super-pixel

        self.mean_full_info = np.mean(np.mean(self.full_info, axis=0), axis=0) / np.mean(self.is_good)
        self.std_full_info = np.sqrt(
            np.mean(np.mean(self.full_info**2, axis=0), axis=0) / np.mean(self.is_good)
            - self.mean_full_info**2
        )
        if verbose:
            print(self.full_info.shape)
            print("Number of good regions =", np.sum(self.is_good))
            print("Mean info from good regions =", self.mean_full_info)
            print("Stdv info from good regions =", self.std_full_info)
            print("")

    def generate_nonlinearity(self, write_to_file=True):
        """
        Generates non-linearity data.

        Parameters
        ----------
        write_to_file : bool, optional
            Write the non-linearity table to a file as well?

        Returns
        -------
        None

        """

        # Non-linearity cube
        ntSub = self.tslices[-1]
        nlcube, self.nlfit, self.nlder = pyirc.gen_nl_cube(
            self.lightfiles,
            self.formatpars,
            ntSub,
            [self.ny, self.nx],
            self.full_info[:, :, self.swi.beta] * self.full_info[:, :, self.swi.I],
            "dev",
            self.swi,
            False,
        )
        if write_to_file:
            with open(self.outstem + "_nl.txt", "w") as thisOut:
                for iy in range(self.ny):
                    for ix in range(self.nx):
                        thisOut.write(
                            "{:3d} {:3d} {:1d} {:9.6f} {:9.6f}".format(
                                iy,
                                ix,
                                int(self.is_good[iy, ix]),
                                self.full_info[iy, ix, self.swi.beta]
                                * self.full_info[iy, ix, self.swi.g]
                                * 1e6,
                                self.full_info[iy, ix, self.swi.g],
                            )
                        )
                        for it in range(ntSub):
                            thisOut.write(f" {nlcube[it,iy,ix]:7.1f}")
                        thisOut.write("\n")

    def write_basic_figure(self):
        """Makes a PDF figure of the SCA parameter maps."""

        # these show up so much let's save them
        dx = self.dx
        dy = self.dy

        # Multi-panel figure showing basic characterization
        ar = float(self.nx / self.ny)
        spr = 2.2
        matplotlib.rcParams.update({"font.size": 8})
        F = plt.figure(figsize=(7, 9))
        S = F.add_subplot(3, 2, 1)
        S.set_title(r"Good pixel map (%)")
        S.set_xlabel(f"Super pixel X/{dx:d}")
        S.set_ylabel(f"Super pixel Y/{dy:d}")
        im = S.imshow(
            self.full_info[:, :, 0] * 100 / self.npix,
            cmap=self.use_cmap,
            aspect=ar,
            interpolation="nearest",
            origin="lower",
            vmin=100 * self.critfrac,
            vmax=100,
        )
        F.colorbar(im, orientation="vertical")
        S = F.add_subplot(3, 2, 2)
        S.set_title(r"Gain map $g$ (e/DN)")
        S.set_xlabel(f"Super pixel X/{dx:d}")
        S.set_ylabel(f"Super pixel Y/{dy:d}")
        svmin, svmax = pyirc.get_vmin_vmax(self.full_info[:, :, self.swi.g], spr)
        im = S.imshow(
            self.full_info[:, :, self.swi.g],
            cmap=self.use_cmap,
            aspect=ar,
            interpolation="nearest",
            origin="lower",
            vmin=svmin,
            vmax=svmax,
        )
        F.colorbar(im, orientation="vertical")
        S = F.add_subplot(3, 2, 3)
        S.set_title(r"IPC map $\alpha$ (%)")
        S.set_xlabel(f"Super pixel X/{dx:d}")
        S.set_ylabel(f"Super pixel Y/{dy:d}")
        svmin, svmax = pyirc.get_vmin_vmax(
            (self.full_info[:, :, self.swi.alphaH] + self.full_info[:, :, self.swi.alphaV]) / 2.0 * 100.0, spr
        )
        im = S.imshow(
            (self.full_info[:, :, self.swi.alphaH] + self.full_info[:, :, self.swi.alphaV]) / 2.0 * 100.0,
            cmap=self.use_cmap,
            aspect=ar,
            interpolation="nearest",
            origin="lower",
            vmin=svmin,
            vmax=svmax,
        )
        F.colorbar(im, orientation="vertical")
        S = F.add_subplot(3, 2, 4)
        S.set_title(r"Non-linearity map $\beta$ (ppm/e)")
        S.set_xlabel(f"Super pixel X/{dx:d}")
        S.set_ylabel(f"Super pixel Y/{dy:d}")
        svmin, svmax = pyirc.get_vmin_vmax(self.full_info[:, :, self.swi.beta] * 1e6, spr)
        im = S.imshow(
            self.full_info[:, :, self.swi.beta] * 1e6,
            cmap=self.use_cmap,
            aspect=ar,
            interpolation="nearest",
            origin="lower",
            vmin=svmin,
            vmax=svmax,
        )
        F.colorbar(im, orientation="vertical")
        S = F.add_subplot(3, 2, 5)
        S.set_title(r"Charge $It_{n,n+1}$ (e):")
        S.set_xlabel(f"Super pixel X/{dx:d}")
        S.set_ylabel(f"Super pixel Y/{dy:d}")
        svmin, svmax = pyirc.get_vmin_vmax(self.full_info[:, :, self.swi.I], spr)
        im = S.imshow(
            self.full_info[:, :, self.swi.I],
            cmap=self.use_cmap,
            aspect=ar,
            interpolation="nearest",
            origin="lower",
            vmin=svmin,
            vmax=svmax,
        )
        F.colorbar(im, orientation="vertical")
        S = F.add_subplot(3, 2, 6)
        S.set_title(r"IPNL $[K^2a+KK^\prime]_{0,0}$ (ppm/e):")
        S.set_xlabel(f"Super pixel X/{dx:d}")
        S.set_ylabel(f"Super pixel Y/{dy:d}")
        svmin, svmax = pyirc.get_vmin_vmax(self.full_info[:, :, self.swi.ker0] * 1e6, spr)
        im = S.imshow(
            self.full_info[:, :, self.swi.ker0] * 1e6,
            cmap=self.use_cmap,
            aspect=ar,
            interpolation="nearest",
            origin="lower",
            vmin=svmin,
            vmax=svmax,
        )
        F.colorbar(im, orientation="vertical")
        F.set_tight_layout(True)
        F.savefig(self.outstem + "_multi.pdf")
        plt.close(F)

    def alt_methods(self, verbose=False):
        """
        Methods 2a, 2b, and 3 trend computation.

        Parameters
        ----------
        verbose : bool, optional
            Whether to print a lot to the output.

        Returns
        -------
        None

        """

        # Method 2a
        #
        self.used_2a = False
        if len(self.tslicesM2a) != 4 or self.tslicesM2a[-1] <= self.tslicesM2a[-2]:
            if verbose:
                print(
                    "Error: tslicesM2a =", self.tslicesM2a, "does not have length 4 or has insufficient span."
                )
                print("Skipping Method 2a ...")
        else:
            # Proceed to implement Method 2a
            self.used_2a = True
            if verbose:
                print("Alt. time slices (Method 2a): ", self.tslicesM2a)
            self.tfmin = self.tslicesM2a[2]
            self.tfmax = self.tslicesM2a[3]
            self.ntM2a = self.tfmax - self.tfmin + 1
            if verbose:
                print("Method 2a, progress of calculation:")
                sys.stdout.write("|")
                for _ in range(self.ny):
                    sys.stdout.write(" ")
                print("| <- 100%")
                sys.stdout.write("|")
            self.Method2a_slopes = np.zeros((self.ny, self.nx))
            self.Method2a_vals = np.zeros((self.ny, self.nx, self.ntM2a))
            lngraw = np.zeros(self.ntM2a)
            for iy in range(self.ny):
                if verbose:
                    sys.stdout.write("*")
                    sys.stdout.flush()
                for ix in range(self.nx):
                    if self.is_good[iy, ix] == 1:
                        for t in range(self.ntM2a):
                            temp_tslices = [
                                self.tslicesM2a[0],
                                self.tslicesM2a[1],
                                self.tslicesM2a[1],
                                self.tfmin + t,
                            ]
                            if self.fullref:
                                lightref = pyirc.ref_array_onerow(
                                    self.lightfiles, self.formatpars, iy, self.ny, temp_tslices, False
                                )
                                darkref = pyirc.ref_array_onerow(
                                    self.darkfiles, self.formatpars, iy, self.ny, temp_tslices, False
                                )
                            region_cube = pyirc.pixel_data(
                                self.lightfiles,
                                self.formatpars,
                                [self.dx * ix, self.dx * (ix + 1), self.dy * iy, self.dy * (iy + 1)],
                                temp_tslices,
                                [self.sensitivity_spread_cut, True],
                                False,
                            )
                            dark_cube = pyirc.pixel_data(
                                self.darkfiles,
                                self.formatpars,
                                [self.dx * ix, self.dx * (ix + 1), self.dy * iy, self.dy * (iy + 1)],
                                temp_tslices,
                                [self.sensitivity_spread_cut, False],
                                False,
                            )
                            info = pyirc.basic(
                                region_cube,
                                dark_cube,
                                temp_tslices,
                                lightref[:, iy, :],
                                darkref[:, iy, :],
                                self.basicpar,
                                False,
                            )
                            self.Method2a_vals[iy, ix, t] = lngraw[t] = np.log(info[1])
                        # Build least squares fit
                        mS, cS = np.linalg.lstsq(
                            np.vstack([np.array(range(self.ntM2a)), np.ones(self.ntM2a)]).T, lngraw, rcond=-1
                        )[0]
                        self.Method2a_slopes[iy, ix] = mS / self.full_info[iy, ix, self.swi.I]
            if verbose:
                print("|")
                print(
                    "Mean slope d[ln graw]/d[I td] at fixed ta,tb =",
                    np.mean(self.is_good * self.Method2a_slopes) / np.mean(self.is_good),
                )
                print("")
            # Predicted slopes
            self.slope_2a_BFE = (
                3 * self.mean_full_info[self.swi.beta]
                - (1 + 4 * self.mean_full_info[self.swi.alphaH] + 4 * self.mean_full_info[self.swi.alphaV])
                * self.mean_full_info[self.swi.ker0]
            )
            self.slope_2a_NLIPC = (
                3 * self.mean_full_info[self.swi.beta]
                - 2
                * (1 + 4 * self.mean_full_info[self.swi.alphaH] + 4 * self.mean_full_info[self.swi.alphaV])
                * self.mean_full_info[self.swi.ker0]
            )

        # Method 2b
        #
        self.used_2b = False
        if len(self.tslicesM2b) != 4 or self.tslicesM2b[-1] <= self.tslicesM2b[-2]:
            print("Error: tslicesM2b =", self.tslicesM2b, "does not have length 4 or has insufficient span.")
            print("Skipping Method 2b ...")
        else:
            # Proceed to implement Method 2b
            self.used_2b = True
            if verbose:
                print("Alt. time slices (Method 2b): ", self.tslicesM2b)
            self.tfminB = self.tslicesM2b[2]
            self.tfmaxB = self.tslicesM2b[3]
            self.ntM2b = self.tfmaxB - self.tfminB + 1
            if verbose:
                print("Method 2b, progress of calculation:")
                sys.stdout.write("|")
                for _ in range(self.ny):
                    sys.stdout.write(" ")
                print("| <- 100%")
                sys.stdout.write("|")
            self.Method2b_slopes = np.zeros((self.ny, self.nx))
            self.Method2b_vals = np.zeros((self.ny, self.nx, self.ntM2b))
            lngraw = np.zeros(self.ntM2b)
            for iy in range(self.ny):
                if verbose:
                    sys.stdout.write("*")
                    sys.stdout.flush()
                for ix in range(self.nx):
                    if self.is_good[iy, ix] == 1:
                        for t in range(self.ntM2b):
                            temp_tslices = [
                                self.tslicesM2b[0] + t,
                                self.tslicesM2b[1] + t,
                                self.tslicesM2b[1] + t,
                                self.tslicesM2b[2] + t,
                            ]
                            if self.fullref:
                                lightref = pyirc.ref_array_onerow(
                                    self.lightfiles, self.formatpars, iy, self.ny, temp_tslices, False
                                )
                                darkref = pyirc.ref_array_onerow(
                                    self.darkfiles, self.formatpars, iy, self.ny, temp_tslices, False
                                )
                            region_cube = pyirc.pixel_data(
                                self.lightfiles,
                                self.formatpars,
                                [self.dx * ix, self.dx * (ix + 1), self.dy * iy, self.dy * (iy + 1)],
                                temp_tslices,
                                [self.sensitivity_spread_cut, True],
                                False,
                            )
                            dark_cube = pyirc.pixel_data(
                                self.darkfiles,
                                self.formatpars,
                                [self.dx * ix, self.dx * (ix + 1), self.dy * iy, self.dy * (iy + 1)],
                                temp_tslices,
                                [self.sensitivity_spread_cut, False],
                                False,
                            )
                            info = pyirc.basic(
                                region_cube,
                                dark_cube,
                                temp_tslices,
                                lightref[:, iy, :],
                                darkref[:, iy, :],
                                self.basicpar,
                                False,
                            )
                            self.Method2b_vals[iy, ix, t] = lngraw[t] = np.log(info[1])
                        # Build least squares fit
                        mS, cS = np.linalg.lstsq(
                            np.vstack([np.array(range(self.ntM2b)), np.ones(self.ntM2b)]).T, lngraw, rcond=-1
                        )[0]
                        self.Method2b_slopes[iy, ix] = mS / self.full_info[iy, ix, self.swi.I]
            if verbose:
                print("|")
                print(
                    "Mean slope d[ln graw]/d[I tb] at fixed tab,tad =",
                    np.mean(self.is_good * self.Method2b_slopes) / np.mean(self.is_good),
                )
                print("")
            # Predicted slopes
            self.slope_2b_BFE = 2 * self.mean_full_info[self.swi.beta]
            self.slope_2b_NLIPC = (
                2 * self.mean_full_info[self.swi.beta]
                + 2
                * (1 + 4 * self.mean_full_info[self.swi.alphaH] + 4 * self.mean_full_info[self.swi.alphaV])
                * self.mean_full_info[self.swi.ker0]
            )

        # Method 3
        #
        self.used_3 = False
        if len(self.tslicesM3) != 4 or self.tslicesM3[-1] <= self.tslicesM3[-2]:
            if verbose:
                print(
                    "Error: tslicesM3 =", self.tslicesM3, "does not have length 4 or has insufficient span."
                )
                print("Skipping Method 3 ...")
        else:
            # Proceed to implement Method 3
            self.used_3 = True
            if verbose:
                print("Alt. time slices (Method 3): ", self.tslicesM3)
            self.tfmin3 = self.tslicesM3[2]
            self.tfmax3 = self.tslicesM3[3]
            self.ntM3 = self.tfmax3 - self.tfmin3 + 1
            if verbose:
                print("Method 3, progress of calculation:")
                sys.stdout.write("|")
                for _ in range(self.ny):
                    sys.stdout.write(" ")
                print("| <- 100%")
                sys.stdout.write("|")
            self.Method3_slopes = np.zeros((self.ny, self.nx))
            self.Method3_vals = np.zeros((self.ny, self.nx, self.ntM3))
            # Method3_alphas = np.zeros((self.ny,self.nx,self.ntM3))
            CCraw = np.zeros(self.ntM3)
            for iy in range(self.ny):
                if verbose:
                    sys.stdout.write("*")
                    sys.stdout.flush()
                for ix in range(self.nx):
                    if self.is_good[iy, ix] == 1:
                        for t in range(self.ntM3):
                            temp_tslices = [
                                self.tslicesM3[0],
                                self.tslicesM3[1],
                                self.tslicesM3[0],
                                self.tfmin3 + t,
                            ]
                            if self.fullref:
                                lightref = pyirc.ref_array_onerow(
                                    self.lightfiles, self.formatpars, iy, self.ny, temp_tslices, False
                                )
                                darkref = pyirc.ref_array_onerow(
                                    self.darkfiles, self.formatpars, iy, self.ny, temp_tslices, False
                                )
                            region_cube = pyirc.pixel_data(
                                self.lightfiles,
                                self.formatpars,
                                [self.dx * ix, self.dx * (ix + 1), self.dy * iy, self.dy * (iy + 1)],
                                temp_tslices,
                                [self.sensitivity_spread_cut, True],
                                False,
                            )
                            dark_cube = pyirc.pixel_data(
                                self.darkfiles,
                                self.formatpars,
                                [self.dx * ix, self.dx * (ix + 1), self.dy * iy, self.dy * (iy + 1)],
                                temp_tslices,
                                [self.sensitivity_spread_cut, False],
                                False,
                            )
                            info = pyirc.basic(
                                region_cube,
                                dark_cube,
                                temp_tslices,
                                lightref[:, iy, :],
                                darkref[:, iy, :],
                                self.basicpar,
                                False,
                            )
                            self.Method3_vals[iy, ix, t] = CCraw[t] = (
                                (info[self.swi.tCH] + info[self.swi.tCV])
                                / 2.0
                                * self.full_info[iy, ix, self.swi.g] ** 2
                                / (self.full_info[iy, ix, self.swi.I] * (temp_tslices[-1] - temp_tslices[0]))
                            )
                            # Method3_alphas[iy,ix,t] = (info[pyirc.swi.alphaH]+info[pyirc.swi.alphaV])/2.
                        # Build least squares fit
                        mS, cS = np.linalg.lstsq(
                            np.vstack([np.array(range(self.ntM3)), np.ones(self.ntM3)]).T, CCraw, rcond=-1
                        )[0]
                        self.Method3_slopes[iy, ix] = mS / self.full_info[iy, ix, self.swi.I]
            if verbose:
                print("|")
                print(
                    "Mean slope d[g^2/(Itad) Cadj,ad]/d[I td] at fixed ta,tb =",
                    np.mean(self.is_good * self.Method3_slopes) / np.mean(self.is_good),
                )
                print("")
            # Predicted slopes
            ave = (
                self.mean_full_info[self.swi.ker0 - 1]
                + self.mean_full_info[self.swi.ker0 + 1]
                + self.mean_full_info[self.swi.ker0 - (2 * self.swi.s + 1)]
                + self.mean_full_info[self.swi.ker0 + (2 * self.swi.s + 1)]
            ) / 4.0
            self.slope_3_beta = (
                -4
                * (self.mean_full_info[self.swi.alphaH] + self.mean_full_info[self.swi.alphaV])
                * self.mean_full_info[self.swi.beta]
            )
            self.slope_3_BFE = (
                -4
                * (self.mean_full_info[self.swi.alphaH] + self.mean_full_info[self.swi.alphaV])
                * self.mean_full_info[self.swi.beta]
                + ave
            )
            self.slope_3_NLIPC = (
                -4
                * (self.mean_full_info[self.swi.alphaH] + self.mean_full_info[self.swi.alphaV])
                * self.mean_full_info[self.swi.beta]
                + ave * 2.0
            )

        # Non-linearity corrections, Methods 2 and 3:
        if verbose:
            print("Non-linearity correction tables:")
        if self.used_2a:
            if verbose:
                print("2a:")
            vec = []
            for t in range(self.tslicesM2a[2], self.tslicesM2a[3] + 1):
                offsets = pyirc.compute_gain_corr_many(
                    self.nlfit,
                    self.nlder,
                    self.full_info[:, :, self.swi.I] * self.full_info[:, :, self.swi.beta],
                    [self.tslicesM2a[0], self.tslicesM2a[1], t],
                    self.basicpar.reset_frame,
                    self.is_good,
                )
                if verbose:
                    print(t, np.mean(offsets * self.is_good) / np.mean(self.is_good))
                vec += [np.mean(offsets * self.is_good) / np.mean(self.is_good)]
            self.PV2a = max(vec) - min(vec)
            if verbose:
                print("PV: ", self.PV2a)
        if self.used_2b:
            if verbose:
                print("2b:")
            vec = []
            dt1 = self.tslicesM2b[1] - self.tslicesM2b[0]
            dt2 = self.tslicesM2b[2] - self.tslicesM2b[0]
            for t in range(self.tslicesM2b[0], self.tslicesM2b[3] - self.tslicesM2b[2] + 1):
                offsets = pyirc.compute_gain_corr_many(
                    self.nlfit,
                    self.nlder,
                    self.full_info[:, :, self.swi.I] * self.full_info[:, :, self.swi.beta],
                    [t, t + dt1, t + dt2],
                    self.basicpar.reset_frame,
                    self.is_good,
                )
                if verbose:
                    print(t, np.mean(offsets * self.is_good) / np.mean(self.is_good))
                vec += [np.mean(offsets * self.is_good) / np.mean(self.is_good)]
            self.PV2b = max(vec) - min(vec)
            if verbose:
                print("PV: ", self.PV2b)
        if self.used_3:
            if verbose:
                print("3:")
            vec = []
            for t in range(self.tslicesM3[2], self.tslicesM3[3] + 1):
                offsets = pyirc.compute_xc_corr_many(
                    self.nlfit,
                    self.nlder,
                    self.full_info[:, :, self.swi.I] * self.full_info[:, :, self.swi.beta],
                    [self.tslicesM3[0], t],
                    self.basicpar.reset_frame,
                    self.is_good,
                )
                alpha3 = (self.full_info[:, :, self.swi.alphaH] + self.full_info[:, :, self.swi.alphaV]) / 2.0
                offsets *= 2.0 * alpha3 * (1.0 - 4 * alpha3)
                if verbose:
                    print(t, np.mean(offsets * self.is_good) / np.mean(self.is_good))
                vec += [np.mean(offsets * self.is_good) / np.mean(self.is_good)]
            self.PV3 = max(vec) - min(vec)
            if verbose:
                print("PV: ", self.PV3)
                print("Method 3 implied slopes =", self.slope_3_beta, self.slope_3_BFE, self.slope_3_NLIPC)
                print("")

    def method_23_plot(self):
        """
        Method 2 and 3 characterization Multi-panel figure showing basic characterization.
        """

        matplotlib.rcParams.update({"font.size": 8})
        F = plt.figure(figsize=(3.5, 9))
        if self.used_2a:
            S = F.add_subplot(3, 1, 1)
            S.set_title(r"Raw gain vs. interval duration")
            S.set_xlabel(r"Signal level $It_{" + f"{self.tslicesM2a[0]:d}" + r",d}$ [ke]")
            S.set_ylabel(r"$\ln g^{\rm raw}_{" + f"{self.tslicesM2a[0]:d},{self.tslicesM2a[1]:d}" + r",d}$")
            SX = [
                np.mean(self.is_good * self.full_info[:, :, self.swi.I] * myt) / np.mean(self.is_good) / 1.0e3
                for myt in range(self.tfmin - self.tslicesM2a[0], self.tfmax + 1 - self.tslicesM2a[0])
            ]
            SY = [
                np.mean(self.is_good * self.Method2a_vals[:, :, t]) / np.mean(self.is_good)
                for t in range(self.ntM2a)
            ]
            SS = []  # std. dev. on the mean
            for t in range(self.ntM2a):
                SS += [
                    np.sqrt(
                        (
                            np.mean(self.is_good * self.Method2a_vals[:, :, t] ** 2) / np.mean(self.is_good)
                            - SY[t] ** 2
                        )
                        / (np.sum(self.is_good) - 1)
                    )
                ]
            xc = np.mean(np.array(SX))
            yc = np.mean(np.array(SY))
            S.set_xlim(min(SX) - 0.05 * (max(SX) - min(SX)), max(SX) + 0.05 * (max(SX) - min(SX)))
            xr = np.arange(min(SX), max(SX), (max(SX) - min(SX)) / 256.0)
            S.errorbar([xc], [min(SY)], yerr=[self.PV2a / 2.0], marker=",", color="k", ls="None")
            S.text(xc + 0.05 * (max(SX) - min(SX)), min(SY), "sys nl", color="k")
            S.errorbar(SX, SY, yerr=SS, marker="x", color="r", ls="None")
            S.plot(xr, yc + (xr - xc) * self.slope_2a_BFE * 1e3, "g--", label="pure BFE")
            S.plot(xr, yc + (xr - xc) * self.slope_2a_NLIPC * 1e3, "b-", label="pure NL-IPC")
            S.legend(loc=2)
        if self.used_2b:
            S = F.add_subplot(3, 1, 2)
            S.set_title(r"Raw gain vs. interval center")
            S.set_xlabel(r"Signal level $It_{" + f"{self.tslicesM2b[0]:d}" + r",a}$ [ke]")
            S.set_ylabel(
                r"$\ln g^{\rm raw}_{"
                + f"a,a+{self.tslicesM2b[1]-self.tslicesM2b[0]:d},a+{self.tslicesM2b[2]-self.tslicesM2b[0]:d}"
                + r"}$"
            )
            SX = [
                np.mean(self.is_good * self.full_info[:, :, self.swi.I] * myt) / np.mean(self.is_good) / 1.0e3
                for myt in range(self.ntM2b)
            ]
            SY = [
                np.mean(self.is_good * self.Method2b_vals[:, :, t]) / np.mean(self.is_good)
                for t in range(self.ntM2b)
            ]
            SS = []  # std. dev. on the mean
            for t in range(self.ntM2b):
                SS += [
                    np.sqrt(
                        (
                            np.mean(self.is_good * self.Method2b_vals[:, :, t] ** 2) / np.mean(self.is_good)
                            - SY[t] ** 2
                        )
                        / (np.sum(self.is_good) - 1)
                    )
                ]
            xc = np.mean(np.array(SX))
            yc = np.mean(np.array(SY))
            S.set_xlim(min(SX) - 0.05 * (max(SX) - min(SX)), max(SX) + 0.05 * (max(SX) - min(SX)))
            xr = np.arange(min(SX), max(SX), (max(SX) - min(SX)) / 256.0)
            S.errorbar([xc], [min(SY)], yerr=[self.PV2b / 2.0], marker=",", color="k", ls="None")
            S.text(xc + 0.05 * (max(SX) - min(SX)), min(SY), "sys nl", color="k")
            S.errorbar(SX, SY, yerr=SS, marker="x", color="r", ls="None")
            S.plot(xr, yc + (xr - xc) * self.slope_2b_BFE * 1e3, "g--", label="pure BFE")
            S.plot(xr, yc + (xr - xc) * self.slope_2b_NLIPC * 1e3, "b-", label="pure NL-IPC")
            S.legend(loc=2)
        if self.used_3:
            S = F.add_subplot(3, 1, 3)
            S.set_title(r"CDS ACF vs. signal")
            S.set_xlabel(r"Signal level $It_{" + f"{self.tslicesM3[0]:d}" + r",d}$ [ke]")
            S.set_ylabel(
                r"$g^2C_{"
                + f"{self.tslicesM3[0]:d}"
                + r"d"
                + f"{self.tslicesM3[0]:d}"
                + r"d}(\langle1,0\rangle)/[It_{"
                + f"{self.tslicesM3[0]:d}"
                + r"d}]$"
            )
            SX = [
                np.mean(self.is_good * self.full_info[:, :, self.swi.I] * myt) / np.mean(self.is_good) / 1.0e3
                for myt in range(self.tfmin3 - self.tslicesM3[0], self.tfmax3 + 1 - self.tslicesM3[0])
            ]
            SY = [
                np.mean(self.is_good * self.Method3_vals[:, :, t]) / np.mean(self.is_good)
                for t in range(self.ntM3)
            ]
            SS = []  # std. dev. on the mean
            for t in range(self.ntM3):
                SS += [
                    np.sqrt(
                        (
                            np.mean(self.is_good * self.Method3_vals[:, :, t] ** 2) / np.mean(self.is_good)
                            - SY[t] ** 2
                        )
                        / (np.sum(self.is_good) - 1)
                    )
                ]
            xc = np.mean(np.array(SX))
            yc = np.mean(np.array(SY))
            S.set_xlim(min(SX) - 0.05 * (max(SX) - min(SX)), max(SX) + 0.05 * (max(SX) - min(SX)))
            xr = np.arange(min(SX), max(SX), (max(SX) - min(SX)) / 256.0)
            S.errorbar([xc], [min(SY)], yerr=[self.PV3], marker=",", color="k", ls="None")
            S.text(xc + 0.05 * (max(SX) - min(SX)), min(SY) + self.PV3, "sys nl", color="k")
            S.errorbar(SX, SY, yerr=SS, marker="x", color="r", ls="None")
            S.plot(xr, yc + (xr - xc) * self.slope_3_BFE * 1e3, "g--", label="pure BFE")
            S.plot(xr, yc + (xr - xc) * self.slope_3_NLIPC * 1e3, "b-", label="pure NL-IPC")
            S.plot(xr, yc + (xr - xc) * self.slope_3_beta * 1e3, "k:", label="beta only")
            S.legend(loc=2)
        F.set_tight_layout(True)
        F.savefig(self.outstem + "_m23.pdf")
        plt.close(F)

    def text_output(self):
        """Generate a text summary.

        Parameters
        ----------
        None

        Returns
        -------
        str
            The output text table.

        """

        # Text output
        thisOut = ""

        # Print information in the file header
        thisOut += f"# This summary created at {time.asctime(time.localtime(time.time())):s}\n"
        thisOut += f"# Uses pyirc v{pyirc.get_version()}\n"
        thisOut += "# Detector: " + self.mydet + "\n"
        thisOut += "#\n# Files used:\n"
        thisOut += "# Light:\n"
        for f in self.lightfiles:
            thisOut += f"#   {f:s}\n"
        thisOut += "# Dark:\n"
        for f in self.darkfiles:
            thisOut += f"#   {f:s}\n"
        thisOut += f"# Input format {self.formatpars:d}\n"
        thisOut += "# Time slices:"
        for t in self.tslices:
            thisOut += f" {t:3d}"
        thisOut += "\n"
        thisOut += "# Mask: " + str(self.maskX) + "," + str(self.maskY) + "\n"
        thisOut += "#\n"
        thisOut += f"# Cut on good pixels {100*self.sensitivity_spread_cut:7.4f}% deviation from median\n"
        thisOut += (
            f"# Dimensions: {self.nx:3d}(x) x {self.ny:3d}(y) super-pixels,"
            f"{int(np.sum(self.is_good)):4d} good\n"
        )
        thisOut += f"# Frame number corresponding to reset: {self.basicpar.reset_frame}\n"
        thisOut += f"# Reference pixel subtraction for linearity: {str(self.fullref):s}\n"
        thisOut += f"# Quantile for variance computation = {self.basicpar.g_ptile:9.6f}%\n"
        thisOut += f"# Clipping fraction epsilon = {self.basicpar.epsilon:9.7f}\n"
        thisOut += "# Lead-trail subtraction for IPC correlation = " + str(self.basicpar.leadtrailSub) + "\n"
        thisOut += "# Characterization type: " + self.mychar + "\n"
        if self.mychar.lower() == "advanced":
            thisOut += f"#   dt = {self.tchar1:d},{self.tchar2:d}, ncycle = {self.ncycle:d}\n"
        thisOut += "# Non-linearity settings: "
        thisOut += "basicpar.fullnl={:s} bfepar.fullnl={:s} basicpar.use_allorder={:s}\n".format(
            str(self.basicpar.fullnl), str(self.bfepar.fullnl), str(self.basicpar.use_allorder)
        )
        thisOut += f"# BFE Method 1\n#   Baseline subtraction = {str(self.bfepar.blsub):s}\n"
        thisOut += f"# BFE Method 2a\n#   Enabled = {str(self.used_2a):s}\n"
        thisOut += f"# BFE Method 2b\n#   Enabled = {str(self.used_2b):s}\n"
        thisOut += f"# BFE Method 3\n#   Enabled = {str(self.used_3):s}\n"
        thisOut += f"# Hot pixels\n#   Enabled = {str(self.hotpix):s}\n"
        if self.hotpix:
            thisOut += f"#   Parameters = {str(self.hotpix_ADU_range):s}\n"
            if self.ref_for_hotpix_is_autocorr:
                thisOut += "#   ref for delta alpha = autocorr\n"
            else:
                thisOut += "#   ref for delta alpha = last frame used\n"
            if self.hotpix_slidemed:
                thisOut += "#   median method = sliding\n"
            else:
                thisOut += "#   median method = normal\n"
        thisOut += "# Associated figures:\n"
        thisOut += "#   {:s}\n".format(self.outstem + "_multi.pdf")
        thisOut += "#   {:s}\n".format(self.outstem + "_m23.pdf")
        thisOut += "#   {:s}\n".format(self.outstem + "_hotipc.pdf")
        thisOut += "#\n"
        thisOut += "# Columns:\n"
        col = 1
        thisOut += f"# {col:3d}, X (super pixel grid)\n"
        col += 1
        thisOut += f"# {col:3d}, Y (super pixel grid)\n"
        col += 1
        thisOut += f"# {col:3d}, number of good pixels\n"
        col += 1
        thisOut += f"# {col:3d}, raw gain (e/DN)\n"
        col += 1
        thisOut += f"# {col:3d}, alpha-corrected gain (e/DN)\n"
        col += 1
        thisOut += f"# {col:3d}, alpha,beta-corrected gain (e/DN)\n"
        col += 1
        thisOut += f"# {col:3d}, IPC alpha horizontal\n"
        col += 1
        thisOut += f"# {col:3d}, IPC alpha vertical\n"
        col += 1
        thisOut += f"# {col:3d}, nonlinearity beta (e^-1)\n"
        col += 1
        thisOut += f"# {col:3d}, charge per time slice (e)\n"
        col += 1
        thisOut += f"# {col:3d}, IPC alpha diagonal (if computed; otherwise all 0s)\n"
        col += 1
        thisOut += f"# {col:3d}, C_H at slices {self.tslices[0]:d},{self.tslices[-1]:d} (DN^2)\n"
        col += 1
        thisOut += f"# {col:3d}, C_V at slices {self.tslices[0]:d},{self.tslices[-1]:d} (DN^2)\n"
        col += 1
        for jb in range(2 * self.swi.s + 1):
            for ib in range(2 * self.swi.s + 1):
                thisOut += (
                    f"# {col:3d}, BFE kernel K^2a (+NL-IPC)"
                    f" at ({ib-self.swi.s:2d},{jb-self.swi.s:2d}) (e^-1)\n"
                )
                col += 1
        if self.swi.p > 0:
            thisOut += f"# {col:3d}, time intercept\n"
            col += 1
            for ip in range(2, self.swi.p + 1):
                thisOut += f"# {col:3d}, additional non-linearity coefficient, order {ip:d} (DN^-{ip-1:d})\n"
                col += 1
        if self.used_2a:
            thisOut += f"# {col:3d}, Method 2a slope (e^-1)\n"
            col += 1
        if self.used_2b:
            thisOut += f"# {col:3d}, Method 2b slope (e^-1)\n"
            col += 1
        if self.used_3:
            thisOut += f"# {col:3d}, Method 3 slope (e^-1)\n"
            col += 1
        thisOut += "#\n"
        # Now make the data table
        for iy in range(self.ny):
            for ix in range(self.nx):
                # Print the column first, then row (normal human-read order,
                # note this is the reverse of internal Python)
                thisOut += f"{ix:3d} {iy:3d}"
                for col in range(np.shape(self.full_info)[-1]):
                    thisOut += f" {self.full_info[iy,ix,col]:14.7E}"
                if self.used_2a:
                    thisOut += f" {self.Method2a_slopes[iy,ix]:14.7E}"
                if self.used_2b:
                    thisOut += f" {self.Method2b_slopes[iy,ix]:14.7E}"
                if self.used_3:
                    thisOut += f" {self.Method3_slopes[iy,ix]:14.7E}"
                thisOut += "\n"

        return thisOut

    def hotpix_analysis(self, verbose=False):
        """
        Hot pixel analysis.

        Parameters
        ----------
        verbose : bool, optional
            Whether to talk a lot.

        Returns
        -------
        str
            Hot pixel report.

        """

        # exit if the hot pixel analysis is not enabled.
        if not self.hotpix:
            return ""

        if verbose:
            print("Start hot pixels ...")
        self.hotY, self.hotX = pyirc.hotpix(
            self.darkfiles, self.formatpars, range(1, self.NTMAX), self.hotpix_ADU_range, True
        )
        if verbose:
            print("Number of pixels selected:", len(self.hotX))  # only printed for de-bugging -> , len(hotY)
        dtstep = 5  # <-- right now this is hard coded
        self.htsteps = range(1, self.NTMAX, dtstep)
        if self.hotpix_logtspace:
            self.htsteps = [1]
            for k in range(1, 12):
                if 2**k < self.NTMAX - 1:
                    self.htsteps += [2**k]
                if k >= 2 and 5 * 2 ** (k - 2) < self.NTMAX - 1:
                    self.htsteps += [5 * 2 ** (k - 2)]
                if 3 * 2 ** (k - 1) < self.NTMAX - 1:
                    self.htsteps += [3 * 2 ** (k - 1)]
                if k >= 2 and 7 * 2 ** (k - 2) < self.NTMAX - 1:
                    self.htsteps += [7 * 2 ** (k - 2)]
            self.htsteps += [self.NTMAX - 1]
        beta_gain = self.full_info[:, :, self.swi.beta] * self.full_info[:, :, self.swi.g]
        if verbose:
            print(beta_gain)
        hotcube = pyirc.hotpix_ipc(
            self.hotY, self.hotX, self.darkfiles, self.formatpars, self.htsteps, [beta_gain, False], True
        )
        nhstep = len(self.htsteps)
        if verbose:
            print("number of time steps ->", nhstep)
        fromcorr_alpha = np.zeros(len(self.hotX))
        hotpix_alpha = np.zeros((len(self.hotX), nhstep))
        hotpix_alpha_num = np.zeros((len(self.hotX), nhstep))
        hotpix_alpha_den = np.zeros((len(self.hotX), nhstep))
        hotpix_alphaD = np.zeros((len(self.hotX), nhstep))
        hotpix_signal = np.zeros((len(self.hotX), nhstep))
        #
        # generate and write hot pixel information
        thisOut = ""
        for jpix in range(len(self.hotX)):
            iy = self.hotY[jpix] // self.dy
            ix = self.hotX[jpix] // self.dx
            fromcorr_alpha[jpix] = (
                self.full_info[iy, ix, self.swi.alphaH] / 2.0 + self.full_info[iy, ix, self.swi.alphaV] / 2.0
            )
            thisOut += f"{self.hotX[jpix]:4d} {self.hotY[jpix]:4d} {fromcorr_alpha[jpix]:8.6f}"
            for t in range(1, nhstep):
                R = (np.mean(hotcube[jpix, t, 1:5]) - hotcube[jpix, t, -1]) / (
                    hotcube[jpix, t, 0] - hotcube[jpix, t, -1]
                )
                S = (np.mean(hotcube[jpix, t, 5:9]) - hotcube[jpix, t, -1]) / (
                    hotcube[jpix, t, 0] - hotcube[jpix, t, -1]
                )
                hotpix_alpha[jpix, t] = R / (1.0 + 4 * R + 4 * S)
                hotpix_alpha_num[jpix, t] = R
                hotpix_alpha_den[jpix, t] = 1.0 + 4 * R + 4 * S
                hotpix_alphaD[jpix, t] = S / (1.0 + 4 * R + 4 * S)
                hotpix_signal[jpix, t] = hotcube[jpix, t, 0] - hotcube[jpix, t, -1]
                thisOut += " {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f}".format(
                    hotcube[jpix, t, 0],
                    np.mean(hotcube[jpix, t, 1:5]),
                    hotcube[jpix, t, -1],
                    (hotcube[jpix, t, 1] + hotcube[jpix, t, 3] - hotcube[jpix, t, 2] - hotcube[jpix, t, 4])
                    / 4.0,
                    np.mean(hotcube[jpix, t, 5:9]),
                )
            thisOut += "\n"

        # report median levels
        if verbose:
            print("IPC relative to nominal (signal, median, uncert):")
        ipcmed_x = np.zeros(nhstep)
        ipcmed_y = np.zeros(nhstep)
        ipcmed_yerr = np.zeros(nhstep)
        delta = 0.5 / np.sqrt(len(self.hotX))
        for t in range(1, nhstep):
            my_y = hotpix_alpha[:, t] - hotpix_alpha[:, -1]
            if self.ref_for_hotpix_is_autocorr:
                my_y = hotpix_alpha[:, t] - fromcorr_alpha
            ipcmed_x[t] = np.nanpercentile(hotpix_signal[:, t], 50.0)
            if self.hotpix_slidemed:
                X_ = hotpix_alpha_den[:, t]
                Y_ = hotpix_alpha_num[:, t] - fromcorr_alpha * hotpix_alpha_den[:, t]
                ipcmed_y[t] = pyirc.slidemed_percentile(X_, Y_, 50)
                ipcmed_yerr[t] = (
                    pyirc.slidemed_percentile(X_, Y_, 50.0 + 100 * delta)
                    - pyirc.slidemed_percentile(X_, Y_, 50.0 - 100 * delta)
                ) / 2.0
            else:
                ipcmed_y[t] = np.nanpercentile(my_y, 50.0)
                ipcmed_yerr[t] = (
                    np.nanpercentile(my_y, 50.0 + 100 * delta) - np.nanpercentile(my_y, 50.0 - 100 * delta)
                ) / 2.0
            print(f"{ipcmed_x[t]:10.2f} {ipcmed_y[t]:9.6f} {ipcmed_yerr[t]:9.6f}")
        print("")
        print(
            "median alphaD:",
            "{:9.6f} {:9.6f}".format(
                np.nanpercentile(hotpix_alphaD[:, -1], 50.0),
                (
                    np.nanpercentile(hotpix_alphaD[:, -1], 50.0 + 100 * delta)
                    - np.nanpercentile(hotpix_alphaD[:, -1], 50.0 - 100 * delta)
                )
                / 2.0,
            ),
        )

        # bigger grid for IPC comparisons
        NG = 4
        grid_alphaCorr = np.zeros((NG, NG))
        grid_alphaCorrErr = np.zeros((NG, NG))
        grid_alphaHot = np.zeros((NG, NG))
        grid_alphaHotErr = np.zeros((NG, NG))
        if self.ny % NG == 0 and self.nx % NG == 0:
            stepx = self.nx // NG
            stepy = self.ny // NG
            sp = self.N // NG
            for ix in range(NG):
                for iy in range(NG):
                    # bin the auto-correlation measurements
                    suba = self.full_info[iy * stepy : (iy + 1) * stepy, ix * stepx : (ix + 1) * stepx, :]
                    mya = (suba[:, :, self.swi.alphaH] + suba[:, :, self.swi.alphaV]) / 2.0
                    pmask = suba[:, :, 0] > 0
                    n = mya[pmask].size
                    if n > 1:
                        grid_alphaCorr[iy, ix] = mya[pmask].mean()
                        grid_alphaCorrErr[iy, ix] = mya[pmask].std() / np.sqrt(n - 1)
                    #
                    # bin the hot pixel measurements -- use final time slice!
                    u = hotpix_alpha[np.logical_and(self.hotY // sp == iy, self.hotX // sp == ix), -1]
                    if u.size > 1:
                        grid_alphaHot[iy, ix] = np.nanpercentile(u, 50)
                        delta = 0.5 / np.sqrt(u.size)
                        grid_alphaHotErr[iy, ix] = (
                            np.nanpercentile(u, 50.0 + 100 * delta) - np.nanpercentile(u, 50.0 - 100 * delta)
                        ) / 2.0
        if verbose:
            print(grid_alphaCorr, grid_alphaCorrErr, grid_alphaHot, grid_alphaHotErr)

        # save information
        self.hotpix_signal = hotpix_signal
        self.hotpix_alpha = hotpix_alpha
        self.ipcmed_x = ipcmed_x
        self.ipcmed_y = ipcmed_y
        self.ipcmed_yerr = ipcmed_yerr
        self.grid_alphaCorr = grid_alphaCorr
        self.grid_alphaCorrErr = grid_alphaCorrErr
        self.grid_alphaHot = grid_alphaHot
        self.grid_alphaHotErr = grid_alphaHotErr

        return thisOut

    def hotpix_plots(self):
        """Makes the hot pixel plots."""

        if not self.hotpix:
            return

        nhstep = len(self.htsteps)

        # hot pixel plots
        matplotlib.rcParams.update({"font.size": 8})
        F = plt.figure(figsize=(3.5, 6)) if self.narrowfig else plt.figure(figsize=(7, 6))
        #
        # hot pixel locations
        S = F.add_subplot(2, 1, 1) if self.narrowfig else F.add_subplot(2, 2, 1)
        S.set_title("hot pixel locations: " + self.mydet)
        S.set_xlim(0, self.N)
        S.set_ylim(0, self.N)
        S.set_aspect("equal")
        S.xaxis.set_ticks(np.linspace(0, self.N, num=5))
        S.yaxis.set_ticks(np.linspace(0, self.N, num=5))
        S.grid(True, color="g", linestyle="-")
        S.scatter(self.hotX + 0.5, self.hotY + 0.5, s=3, marker=".", color="r")
        #
        # these two panels only in the full version, not the narrow version
        if not self.narrowfig:
            # hot pixel level
            S = F.add_subplot(2, 2, 2)
            S.set_xlabel(r"Signal level $S_{1," + f"{self.htsteps[-1]:d}" + "}$ [DN]")
            S.set_ylabel(r"IPC $\alpha$ [%]")
            SX = self.hotpix_signal[:, -1]
            SY = self.hotpix_alpha[:, -1] / 0.01
            S.set_title(r"IPC $\alpha$ for hot pixels")
            S.set_xlim(
                0.95 * (self.htsteps[-1] - 1) / (self.NTMAX - 1.0) * self.hotpix_ADU_range[0],
                1.05 * self.hotpix_ADU_range[1],
            )
            S.set_ylim(0, 4.0)
            S.grid(True, color="g", linestyle="-")
            S.scatter(SX, SY, s=3, marker=".", color="r")
            #
            # dependence on signal level
            S = F.add_subplot(2, 2, 3)
            S.set_xlabel(r"Signal level $S_{1,b}$ [DN]")
            S.set_ylabel(r"IPC $\alpha(S_{1,b})-\alpha(S_{1," + f"{self.htsteps[-1]:d}" + "})$ [%]")
            for t in range(1, nhstep):
                SXa = self.hotpix_signal[:, t]
                SYa = (self.hotpix_alpha[:, t] - self.hotpix_alpha[:, -1]) / 0.01
                if t == 1:
                    SX = SXa
                    SY = SYa
                else:
                    print(t, np.shape(SX), np.shape(SXa))
                    SX = np.concatenate((SX, SXa))
                    SY = np.concatenate((SY, SYa))
            S.set_title(r"IPC signal dependence $\Delta\alpha$")
            S.set_xlim(0.0, self.hotpix_ADU_range[1])
            S.set_ylim(-1.5, 1.5)
            S.xaxis.set_ticks(np.linspace(0, self.hotpix_ADU_range[1], num=6))
            S.yaxis.set_ticks(np.linspace(-1.5, 1.5, num=11))
            S.grid(True, color="g", linestyle="-", linewidth=0.5)
            S.scatter(SX, SY, s=0.25, marker="+", color="r")
            S.errorbar(
                self.ipcmed_x[1:],
                self.ipcmed_y[1:] / 0.01,
                yerr=self.ipcmed_yerr[1:] / 0.01,
                ms=2,
                marker="o",
                color="k",
                ls="None",
            )
        #
        # comparison with auto-correlations
        S = F.add_subplot(2, 1, 2) if self.narrowfig else F.add_subplot(2, 2, 4)
        S.set_title(r"hot pixels vs. autocorr. IPC $\alpha$")
        scale_test = np.concatenate((self.grid_alphaCorr.flatten(), self.grid_alphaHot.flatten()))
        smin = 0.92 * np.min(scale_test[scale_test > 0]) / 0.01
        smax = 1.08 * np.max(scale_test) / 0.01
        S.set_xlim(smin, smax)
        S.set_ylim(smin, smax)
        S.set_aspect("equal")
        S.set_xlabel(r"autocorrelation $\alpha$ [%]")
        S.set_ylabel(r"hot pixel $\alpha$ [%]")
        S.grid(True, color="g", linestyle="-", linewidth=0.5)
        S.errorbar(
            self.grid_alphaCorr.flatten() / 0.01,
            self.grid_alphaHot.flatten() / 0.01,
            xerr=self.grid_alphaCorrErr.flatten() / 0.01,
            yerr=self.grid_alphaHotErr.flatten() / 0.01,
            ms=1,
            marker="o",
            color="k",
            capsize=1,
            ls="None",
        )
        xr = np.linspace(0, 4, num=65)
        S.plot(xr, xr, "r-")

        F.set_tight_layout(True)
        F.savefig(self.outstem + "_hotipc.pdf")
        plt.close(F)

    def compute_vis_quantities(self, ir_output=None, verbose=False):
        """
        Computations for the visible light characterization.

        Parameters
        ----------
        ir_output : str, optional
            The data from the IR output characterization as a string; if not specified,
            tries to load from the file.
        verbose : bool, optional
            Whether to talk a lot.

        Returns
        -------
        None

        """

        # reference pixel subtraction flag
        self.basicpar.subtr_href = self.fullref

        # more allocations
        my_dim = self.swi.N
        self.full_info = np.zeros((self.ny, self.nx, my_dim))
        self.is_good = np.zeros((self.ny, self.nx))

        info_from_ir = np.loadtxt(self.outstem + "_summary.txt") if ir_output is None else ir_output
        for j in range(my_dim):
            self.full_info[:, :, j] = info_from_ir[:, j + 2].reshape((self.ny, self.nx))
        self.is_good = np.where(self.full_info[:, :, self.swi.g] > 1e-49, 1, 0)

        if verbose:
            print("Number of good regions =", np.sum(self.is_good))
            print("Lower-left corner ->", self.full_info[0, 0, :])

        if self.p_order == 0:
            raise ValueError("Error: did not include polynomial order")

        # Get Ie
        Ie = np.zeros((self.ny, self.nx))
        Ie_alt = np.zeros((self.ny, self.nx))
        Ie_alt2 = np.zeros((self.ny, self.nx))

        if verbose:
            print("computing Ie using", self.ts_vis, self.te_vis)
        nlcubeX, nlfitX, nlderX, pcoefX = pyirc.gen_nl_cube(
            self.vislightfiles,
            self.formatpars,
            [self.basicpar.reset_frame, self.ts_vis, self.te_vis],
            [self.ny, self.nx],
            self.full_info[:, :, 0],
            "abs",
            self.swi,
            False,
        )
        for iy in range(self.ny):
            for ix in range(self.nx):
                if pcoefX[1, iy, ix] != 0:
                    t = np.linspace(
                        self.ts_vis - self.basicpar.reset_frame,
                        self.te_vis - self.basicpar.reset_frame,
                        self.te_vis - self.ts_vis + 1,
                    )
                    Signal = np.zeros(self.te_vis - self.ts_vis + 1)
                    for ae in range(self.swi.p + 1):
                        Signal += pcoefX[ae, iy, ix] * t**ae
                    # iterative NL correction
                    LinSignal = np.copy(Signal)
                    for _ in range(32):
                        LS2 = np.copy(LinSignal)
                        LinSignal = np.copy(Signal)
                        LS2 += (
                            (LinSignal[-1] - LinSignal[0])
                            / (self.te_vis - self.ts_vis)
                            * (self.ts_vis - self.basicpar.reset_frame)
                        )
                        for o in range(2, self.swi.p + 1):
                            LinSignal -= self.full_info[iy, ix, self.swi.Nbb + o - 1] * LS2**o
                    Ie[iy, ix] = pcoefX[1, iy, ix] * self.full_info[iy, ix, self.swi.g]
                    Ie_alt[iy, ix] = (
                        (LinSignal[-1] - LinSignal[0])
                        / (self.te_vis - self.ts_vis)
                        * self.full_info[iy, ix, self.swi.g]
                    )
                    Sab = Signal[-1] - Signal[0]
                    Ie_alt2[iy, ix] = self.full_info[iy, ix, self.swi.g] * Sab / (self.te_vis - self.ts_vis)
                    beta_in_e = -self.full_info[
                        iy, ix, self.swi.Nbb + 1 : self.swi.Nbb + self.swi.p
                    ] / self.full_info[iy, ix, self.swi.g] ** np.linspace(
                        1, self.swi.p - 1, num=self.swi.p - 1
                    )  # in e , -
                    for _ in range(32):
                        btcorr = 0
                        for j in range(2, self.swi.p + 1):
                            btcorr += beta_in_e[j - 2] * Ie_alt2[iy, ix] ** (j - 1) * (t[-1] ** j - t[0] ** j)
                        Ie_alt2[iy, ix] = (
                            self.full_info[iy, ix, self.swi.g] * Sab / (self.te_vis - self.ts_vis - btcorr)
                        )
                else:
                    self.is_good[iy, ix] = 0  # error

        # we use the alt2 method
        Ie[:, :] = Ie_alt2

        # get vis:IR Ie ratio information
        vis_ir_ratio = Ie / self.full_info[:, :, self.swi.I]
        vis_ir_ratio_good = vis_ir_ratio[self.is_good > 0.5]
        if verbose:
            print("VIS:IR ratio information: ", np.shape(vis_ir_ratio_good))
            print("min, max =", np.amin(vis_ir_ratio_good), np.amax(vis_ir_ratio_good))
            print(
                "percentiles (5th,50th,95th)",
                np.percentile(vis_ir_ratio_good, 5),
                np.percentile(vis_ir_ratio_good, 50),
                np.percentile(vis_ir_ratio_good, 95),
            )
            print("")

        # Allocate space for visible information
        vis_bfek = np.zeros((self.ny, self.nx, 5, 5))
        vis_Phi = np.zeros((self.ny, self.nx, 5, 5))
        # omega and charge diffusion covariance
        QYomega = np.zeros((self.ny, self.nx))
        cdCov = np.zeros((self.ny, self.nx, 3))
        cdNiter = np.zeros((self.ny, self.nx))

        # Get correlation functions in each block
        nvis = self.te_vis - self.ts_vis - self.tchar2_vis + 1
        if verbose:
            print("Visible flat correlation functions, progress of calculation:")
            sys.stdout.write("|")
            for _ in range(self.ny):
                sys.stdout.write(" ")
            print("| <- 100%")
            sys.stdout.write("|")
        for iy in range(self.ny):
            if verbose:
                sys.stdout.write("*")
                sys.stdout.flush()
            if self.fullref:
                tslices0 = np.asarray(
                    [self.ts_vis, self.ts_vis + self.tchar1_vis, self.ts_vis + self.tchar2_vis]
                )
                lightref_array = []
                darkref_array = []
                for k in range(nvis):
                    tslicesk = (tslices0 + k).tolist()
                    lightref_array.append(
                        pyirc.ref_array(self.vislightfiles, self.formatpars, self.ny, tslicesk, False)
                    )
                    darkref_array.append(
                        pyirc.ref_array(self.vislightfiles, self.formatpars, self.ny, tslicesk, False)
                    )
            for ix in range(self.nx):
                if self.is_good[iy, ix] > 0.5:
                    # pull out basic parameters
                    basicinfo = self.full_info[iy, ix, : self.swi.Nb].tolist()
                    basicinfo[self.swi.I] = Ie[iy, ix]
                    basicinfo[self.swi.beta] = self.full_info[
                        iy, ix, self.swi.Nbb + 1 : self.swi.Nbb + self.swi.p
                    ]  # in DN, +
                    beta_in_e = -basicinfo[self.swi.beta] / basicinfo[self.swi.g] ** np.linspace(
                        1, self.swi.p - 1, num=self.swi.p - 1
                    )  # in e , -

                    tslices0 = np.asarray(
                        [self.ts_vis, self.ts_vis + self.tchar1_vis, self.ts_vis + self.tchar2_vis]
                    )
                    # initialize vector to stack correlation matrices:
                    corr_stack = []
                    for k in range(nvis):
                        tslicesk = (tslices0 + k).tolist()
                        region_cube = pyirc.pixel_data(
                            self.vislightfiles,
                            self.formatpars,
                            [self.dx * ix, self.dx * (ix + 1), self.dy * iy, self.dy * (iy + 1)],
                            tslicesk,
                            [self.sensitivity_spread_cut, True],
                            False,
                        )
                        dark_cube = pyirc.pixel_data(
                            self.visdarkfiles,
                            self.formatpars,
                            [self.dx * ix, self.dx * (ix + 1), self.dy * iy, self.dy * (iy + 1)],
                            tslicesk,
                            [self.sensitivity_spread_cut, False],
                            False,
                        )
                        if self.fullref:
                            lightref = lightref_array[k]
                            darkref = darkref_array[k]
                        else:
                            lightref = np.zeros((len(self.vislightfiles), self.ny, 2 * len(tslicesk) + 1))
                            darkref = np.zeros((len(self.visdarkfiles), self.ny, 2 * len(tslicesk) + 1))
                        info = pyirc.corr_5x5(
                            region_cube,
                            dark_cube,
                            tslicesk,
                            lightref[:, iy, :],
                            darkref[:, iy, :],
                            self.basicpar,
                            False,
                        )

                        corr_matrix = info[4]
                        var1 = info[2]
                        var2 = info[3]
                        # center of corr_matrix is element (2, 2) of the numpy array
                        corr_matrix[2][2] = var2 - var1

                        # median corrections to the central array of the auto-correlation matrix
                        # (so we multiply the measured variance by the measured/predicted median,
                        # this would perfectly correct for errors in Ie if the detector were exactly linear)
                        med21 = info[1]
                        predictmed = (
                            tslicesk[2]
                            * Ie[iy, ix]
                            * (
                                1.0
                                - np.sum(
                                    beta_in_e
                                    * (tslicesk[2] * Ie[iy, ix])
                                    ** np.linspace(1, self.swi.p - 1, num=self.swi.p - 1)
                                )
                            )
                            - tslicesk[1]
                            * Ie[iy, ix]
                            * (
                                1.0
                                - np.sum(
                                    beta_in_e
                                    * (tslicesk[1] * Ie[iy, ix])
                                    ** np.linspace(1, self.swi.p - 1, num=self.swi.p - 1)
                                )
                            )
                        ) / basicinfo[self.swi.g]
                        if self.basicpar.vis_med_correct:
                            corr_matrix[2][2] /= med21 / predictmed

                        corr_stack.append(corr_matrix)
                        # end loop over k

                    corr_mean = np.mean(corr_stack, axis=0)
                    # corr_mean is the v vector of eq. 34

                    # now get the cube of data for BFE
                    region_cube = pyirc.pixel_data(
                        self.vislightfiles,
                        self.formatpars,
                        [self.dx * ix, self.dx * (ix + 1), self.dy * iy, self.dy * (iy + 1)],
                        self.tslices,
                        [self.sensitivity_spread_cut, True],
                        False,
                    )

                    # iterate to solve BFE, Phi

                    np2 = 2
                    self.bfepar.Phi = np.zeros((2 * np2 + 1, 2 * np2 + 1))
                    self.bfepar.Phi[np2, np2] = 1.0e-12  # initialize to essentially zero
                    if self.copy_ir_bfe:
                        bfek_ir = self.full_info[iy, ix, self.swi.Nb : self.swi.Nbb].reshape(
                            (2 * np2 + 1, 2 * np2 + 1)
                        )
                        bfek = np.copy(bfek_ir)
                    else:
                        bfek = pyirc.bfe(region_cube, self.tslices, basicinfo, self.bfepar, self.swi, False)
                    tol = 1e-9
                    diff = 1
                    count = 0
                    NN = 21

                    while np.max(np.abs(diff)) > tol:
                        ts_vis_ref = self.ts_vis - self.basicpar.reset_frame
                        tslices_vis = [
                            ts_vis_ref,
                            ts_vis_ref + self.tchar2_vis,
                            ts_vis_ref,
                            ts_vis_ref + self.tchar2_vis,
                            nvis,
                        ]
                        tslices_vis1 = [
                            ts_vis_ref,
                            ts_vis_ref + self.tchar1_vis,
                            ts_vis_ref,
                            ts_vis_ref + self.tchar1_vis,
                            nvis,
                        ]
                        normPhi = np.sum(self.bfepar.Phi)  # this is omega/(1+omega)
                        omega = normPhi / (1 - normPhi)
                        p2 = self.bfepar.Phi / normPhi
                        sigma_a = 0.0
                        avals = [
                            basicinfo[self.swi.alphaV],
                            basicinfo[self.swi.alphaH],
                            basicinfo[self.swi.alphaD],
                        ]  # (aV, aH, aD)
                        truecorr = ftsolve.solve_corr_vis_many(
                            bfek,
                            NN,
                            basicinfo[self.swi.I],
                            basicinfo[self.swi.g],
                            beta_in_e,
                            sigma_a,
                            tslices_vis,
                            avals,
                            omega=omega,
                            p2=p2,
                        )
                        # if count==0:
                        #  print(tslices_vis, p2, truecorr)
                        truecorr[2, 2] = (
                            truecorr
                            - ftsolve.solve_corr_vis_many(
                                bfek,
                                NN,
                                basicinfo[self.swi.I],
                                basicinfo[self.swi.g],
                                beta_in_e,
                                sigma_a,
                                tslices_vis1,
                                avals,
                                omega=omega,
                                p2=p2,
                            )
                        )[2][2]
                        diff = (
                            basicinfo[self.swi.g] ** 2
                            / (2 * basicinfo[self.swi.I] * self.tchar2_vis)
                            * (corr_mean - truecorr)
                        )
                        diff[2, 2] = (
                            basicinfo[self.swi.g] ** 2
                            / (2 * basicinfo[self.swi.I] * (self.tchar2_vis - self.tchar1_vis))
                            * (corr_mean[2, 2] - truecorr[2, 2])
                        )
                        self.bfepar.Phi += 0.5 * (
                            diff + np.flip(diff)
                        )  # force symmetrization here to avoid instability

                        # update BFE
                        if self.copy_ir_bfe:
                            bfek = np.copy(bfek_ir)
                        else:
                            bfek = pyirc.bfe(
                                region_cube, self.tslices, basicinfo, self.bfepar, self.swi, False
                            )
                        count += 1

                        if count > 100:
                            if verbose:
                                print(
                                    f"100 iterations of BFE/Phi solver reached,"
                                    f" diff={np.max(np.abs(diff)):0.6f}"
                                )
                            break

                    # save information
                    vis_bfek[iy, ix, :, :] = bfek
                    vis_Phi[iy, ix, :, :] = self.bfepar.Phi
                    op2 = ftsolve.op2_to_pars(self.bfepar.Phi)
                    QYomega[iy, ix] = op2[0]
                    cdCov[iy, ix, 0] = op2[1]
                    cdCov[iy, ix, 1] = op2[2]
                    cdCov[iy, ix, 2] = op2[3]
                    cdNiter[iy, ix] = op2[-1]

                    # end loop over super-pixels
        if verbose:
            print("|")
            print("")

            # Now get ready to write information
            print("Mean BFE kernel:")
            print(np.mean(vis_bfek, axis=(0, 1)))
            print("Mean Phi kernel:")
            print(np.mean(vis_Phi, axis=(0, 1)))
            print("sigma Phi kernel:")
            print(np.std(vis_Phi, axis=(0, 1)))
            print("Charge diffusion parameters:")
            print(ftsolve.op2_to_pars(np.mean(vis_Phi, axis=(0, 1))))

        # put all information into a gigantic array
        vis_out_data = np.zeros((self.ny, self.nx, 56))
        vis_out_data[:, :, :25] = vis_bfek.reshape(self.ny, self.nx, 25)
        vis_out_data[:, :, 25:50] = vis_Phi.reshape(self.ny, self.nx, 25)
        vis_out_data[:, :, 50] = QYomega
        vis_out_data[:, :, 51:54] = cdCov
        vis_out_data[:, :, 54] = Ie
        vis_out_data[:, :, 55] = cdNiter
        ncol = 56
        #
        # now we have in each super-pixel, 55 "columns" of data
        # columns  0 .. 24 are the visible BFE kernel in e^-1
        #     (order: dy=-2 dx=-2; dy=-2 dx=-1; dy=-2 dx=0; ...)
        # columns 25 .. 49 are the visible Phi kernel (order: dy=-2 dx=-2; dy=-2 dx=-1; dy=-2 dx=0; ...)
        # column 50 is the quantum yield omega parameter
        # column 51 is Cxx charge diffusion in pixels^2
        # column 52 is Cxy charge diffusion in pixels^2
        # column 53 is Cyy charge diffusion in pixels^2
        # column 54 is visible current Ie (e per frame)
        # column 55 is number of iterations in p2 kernel
        self.vis_col = {
            "BFE00": 12,
            "Phi00": 37,
            "QYomega": 50,
            "Cxx": 51,
            "Cxy": 52,
            "Cyy": 53,
            "Ie": 54,
            "cdNiter": 55,
        }

        mean_vis_out_data = np.mean(np.mean(vis_out_data, axis=0), axis=0) / np.mean(self.is_good)
        std_vis_out_data = np.sqrt(
            np.mean(np.mean(vis_out_data**2, axis=0), axis=0) / np.mean(self.is_good) - mean_vis_out_data**2
        )
        if verbose:
            print("")
            print(vis_out_data.shape)
            print("Number of good regions =", np.sum(self.is_good))
            print("column, mean, stdev, stdev on the mean:")
            for k in range(ncol):
                print(
                    f"{k:2d} {mean_vis_out_data[k]:12.5E} {std_vis_out_data[k]:12.5E}"
                    f" {std_vis_out_data[k]/np.sqrt(np.sum(self.is_good)-1):12.5E}"
                )
            print("")

        # save to file and class
        np.savetxt(self.outstem + "_visinfo.txt", vis_out_data.reshape(self.ny * self.nx, ncol))
        self.vis_out_data = vis_out_data

    def vis_plots(self):
        """Make plots for the visible light characterization."""

        # Saving some figures of these quantities:
        matplotlib.rcParams.update({"font.size": 12})
        num_bins = 30
        F = plt.figure(figsize=(8, 6))
        S = F.add_subplot(2, 2, 1)
        S.hist(
            self.vis_out_data[:, :, self.vis_col["QYomega"]].ravel(), bins=np.linspace(0, 0.1, num=num_bins)
        )
        S.set_xlabel(r"$\omega$")

        S = F.add_subplot(2, 2, 2)
        S.hist(self.vis_out_data[:, :, self.vis_col["Ie"]].ravel(), bins=num_bins)
        S.set_xlabel(r"$I_e$")

        S = F.add_subplot(2, 2, 3)
        S.hist(
            self.vis_out_data[:, :, self.vis_col["cdNiter"]].ravel(), bins=np.linspace(0, 100, num=num_bins)
        )
        S.set_xlabel(r"Number of iterations")

        S = F.add_subplot(2, 2, 4)
        S.hist(
            self.vis_out_data[:, :, self.vis_col["Cxx"]].ravel(),
            num_bins,
            histtype="step",
            label=r"$C_{xx}$",
            linewidth=1.5,
            linestyle="-",
        )
        S.hist(
            self.vis_out_data[:, :, self.vis_col["Cxy"]].ravel(),
            num_bins,
            histtype="step",
            label=r"$C_{xy}$",
            linewidth=1.5,
            linestyle="--",
        )
        S.hist(
            self.vis_out_data[:, :, self.vis_col["Cyy"]].ravel(),
            num_bins,
            histtype="step",
            label=r"$C_{yy}$",
            linewidth=1.5,
            linestyle="-.",
        )
        S.set_xlabel(r"Charge diffusion component in pixels$^2$")
        S.legend(loc="upper right", fontsize=12, frameon=False)
        F.set_tight_layout(True)
        F.savefig(self.outstem + "_vis_hist.pdf", bbox_inches="tight")
        plt.close(F)

        F = plt.figure(figsize=(15, 8))
        S = F.add_subplot(2, 3, 1)
        S.set_title(r"$\omega$")
        S.set_xlabel(f"Super pixel X/{self.dx:d}")
        S.set_ylabel(f"Super pixel Y/{self.dy:d}")
        im = S.imshow(self.vis_out_data[:, :, self.vis_col["QYomega"]], cmap=self.use_cmap, origin="lower")
        F.colorbar(im, orientation="vertical")

        S = F.add_subplot(2, 3, 2)
        S.set_title(r"$I_e$")
        S.set_xlabel(f"Super pixel X/{self.dx:d}")
        # S.set_ylabel('Super pixel Y/{:d}'.format(self.dy))
        im = S.imshow(self.vis_out_data[:, :, self.vis_col["Ie"]], cmap=self.use_cmap, origin="lower")
        F.colorbar(im, orientation="vertical")

        S = F.add_subplot(2, 3, 3)
        S.set_title(r"Number of iterations")
        S.set_xlabel(f"Super pixel X/{self.dx:d}")
        # S.set_ylabel('Super pixel Y/{:d}'.format(self.dy))
        im = S.imshow(self.vis_out_data[:, :, self.vis_col["cdNiter"]], cmap=self.use_cmap, origin="lower")
        F.colorbar(im, orientation="vertical")

        S = F.add_subplot(2, 3, 4)
        S.set_title(r"$C_{xx}$")
        S.set_xlabel(f"Super pixel X/{self.dx:d}")
        S.set_ylabel(f"Super pixel Y/{self.dy:d}")
        im = S.imshow(self.vis_out_data[:, :, self.vis_col["Cxx"]], cmap=self.use_cmap, origin="lower")
        F.colorbar(im, orientation="vertical")

        S = F.add_subplot(2, 3, 5)
        S.set_title(r"$C_{xy}$")
        S.set_xlabel(f"Super pixel X/{self.dx:d}")
        # S.set_ylabel('Super pixel Y/{:d}'.format(self.dy))
        im = S.imshow(self.vis_out_data[:, :, self.vis_col["Cxy"]], cmap=self.use_cmap, origin="lower")
        F.colorbar(im, orientation="vertical")

        S = F.add_subplot(2, 3, 6)
        S.set_title(r"$C_{yy}$")
        S.set_xlabel(f"Super pixel X/{self.dx:d}")
        # S.set_ylabel('Super pixel Y/{:d}'.format(self.dy))
        im = S.imshow(self.vis_out_data[:, :, self.vis_col["Cyy"]], cmap=self.use_cmap, origin="lower")
        F.colorbar(im, orientation="vertical")

        # F.set_tight_layout(True)
        F.savefig(self.outstem + "_vis_matrices.pdf", bbox_inches="tight")
        plt.close(F)


def run_ir_all(infile):
    """
    Runs the IR characterization.

    Parameters
    ----------
    infile : str
        The input file.

    Returns
    -------
    None

    """

    # Get configuration and copy to the output.
    with open(infile) as f:
        cf = Config(f.readlines(), verbose=True)
    shutil.copyfile(infile, cf.outstem + "_config.txt")
    cf.fit_parameters(verbose=True)
    cf.generate_nonlinearity(write_to_file=True)
    cf.write_basic_figure()
    cf.alt_methods(verbose=True)
    cf.method_23_plot()
    with open(cf.outstem + "_summary.txt", "w") as f:
        f.write(cf.text_output())
    s = cf.hotpix_analysis(verbose=True)
    with open(cf.outstem + "_hot.txt", "w") as f:
        f.write(s)
    cf.hotpix_plots()


def run_vis_all(infile, run_ir_first=True):
    """
    Runs the visible characterization.

    Parameters
    ----------
    infile : str
        The input file.
    run_ir_first : bool, optional
        Run the IR characterization first (turn off if you already ran it!).

    Returns
    -------
    None

    """

    if run_ir_first:
        run_ir_all(infile)

    # Get configuration. Note this is a new configuration instance!
    with open(infile) as f:
        cf = Config(f.readlines(), visible_run=True, verbose=True)

    cf.compute_vis_quantities(verbose=True)
    cf.vis_plots()


### <== MAIN PROGRAM BELOW HERE ==> ###

if __name__ == "__main__":
    """Main run."""

    run_ir_all(sys.argv[1])

|badge1| |badge2|

.. |badge1| image:: https://codecov.io/gh/Roman-HLIS-Cosmology-PIT/solid-waffle/graph/badge.svg

.. |badge2| image:: https://github.com/Roman-HLIS-Cosmology-PIT/solid-waffle/actions/workflows/smoke-test.yml/badge.svg

solid-waffle
============

Overview
--------

**solid-waffle** is a repository designed for analysis of correlations in HxRG flats, with an assortment of other features. The repository contains scripts, configuration files, and utilities to support the characterization of HxRG detectors.

*Nota bene*\ :  These scripts are not ready for "black box" use in the sense that unanticipated defects in different detectors can cause issues in the stability of the characterization.  Some work and testing is required before given results can be trusted.  Manual masking is required in some cases, and we are working towards getting this code more streamlined and user-friendly in that regard. Please approach with caution and reach out to the `contacts <#contacts>`_ listed below for help.

See the `readthedocs page <https://solid-waffle.readthedocs.io/en/latest/>`_.

Table of Contents
-----------------


#. `Project Structure <#project-structure>`_

#. `Documentation <#documentation>`_

#. `Files Description <#files-description>`_

#. `Directories Overview <#directories-overview>`_

#. `How to Get Started <#how-to-get-started>`_

#. `Contacts <#contacts>`_

#. `References <#references>`_

Project Structure
-----------------

The repository is organized as follows:

Documentation
-------------

Some available help on:

* `How to run a solid-waffle characterization script <docs/ScriptInformation.rst>`_.

* `Noise and linearity scripts <docs/noise_linearity.rst>`_.

* `Single pixel reset script <docs/SPR.rst>`_.

* `Converting ASDF files <docs/asdf.rst>`_.

**Note**: If you want the old (non-module) interface, you can use the scripts in ``legacy_scripts/``. But we anticipate that most new code won't use that.

Files Description
-----------------


src/solid_waffle/correlation_run.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


* **Purpose**\ : Main script to call that takes a config file as an argument and runs the analysis.
* **Key Features**\ :

  * Detector characterization at basic and advanced levels, including nonlinearity and Brighter-Fatter Effect (BFE)
  * Hot pixel identification and calculation of interpixel capacitance (IPC)
  * Summary statistics in superpixels
  * Visualization such as maps of gain, IPC, BFE over superpixels
  * *If visible light data is supplied as well*: Estimates for charge diffusion covariance and quantum yield, also summarized and visualized in outputs

src/solid_waffle/pyirc.py
^^^^^^^^^^^^^^^^^^^^^^^^^

* **Purpose**\ : Main utility script containing many functions to load in data, calculate statistics, and estimate corrections.
* **Details**\ : Includes functions that will read the data in, compute reference pixel corrections, compute gain, do basic characterization, calculate correlations of neighboring pixels, compute BFE coefficients, and more.

src/solid_waffle/ftsolve.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **Purpose**\ : Routines for modeling Fourier-domain correlations across time (Feudenburg et al. 2020).

src/solid_waffle/spr_reduce.py
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* **Purpose**\ : Analysis of single pixel reset data.
* **Key Features**\ : 

  * Extraction of spatially varying IPC maps from single pixel reset data.

Directories Overview
--------------------

sample_configs
^^^^^^^^^^^^^^

* **Description**\ : Sample configuration files.
* **Key Files**\ :

  * example_config_wfirst_h4rg_18237: A sample configuration file which is simplified to the basics. This config requires only flat files, dark files, a format code, time frames, and an output location to be specified.

  * config.vis-sim, config.vis1: Example configuration files for visible light characterization.

  * ex_sim_config: Example configuration for the flat simulator. Run with ``python -m solid_waffle.flat_simulator.simulate_flat ex_sim_config``

src/solid_waffle/flat_simulator
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


* **Description**\ : Scripts and utilities for simulating flat field images.
* **Key Files**\ :

  * simulate_flat.py: Main script to generate a simulated flat field.

  * detector_function.py: Some utility functions called by simulate_flat.py

notebooks
^^^^^^^^^


* **Description**\ : Jupyter notebooks for development, analysis, and visualization. None of these are particularly well described currently.
* **Key Notebooks**\ :

  * plotting_tests.ipynb: Shows some visualizations of outputs.

projects
^^^^^^^^

* **Description**\ : Other projects related to detector characterization that aren't in the main workflow right now (and aren't part of the installable module).
* **Key Sub-Folders**\ :

  * crnl_tools: Tools for the "direct method" count-rate non-linearity ground test data.

  * speckle: Tools for speckle fringe data.


Contacts
--------

To communicate about this repository please reach out to:


* Chris Hirata (hirata.10 at osu dot edu)
* Ami Choi (ami.choi at nasa dot gov)

References
----------

For more detailed background on the concepts and methods used in this project, please refer to:

* Main correlation analysis projects for flats and darks:

  * Hirata, C. & Choi, A. (2020). *Brighter-fatter Effect in Near-infrared Detectors. I. Theory of Flat Autocorrelations*. Publications of the Astronomical Society of the Pacific, Volume 132, Issue 1007, id. 014501 `Link to abstract <https://ui.adsabs.harvard.edu/abs/2020PASP..132a4501H/abstract>`_
  * Choi, A. & Hirata, C. (2020). *Brighter-fatter Effect in Near-infrared Detectors. II. Autocorrelation Analysis of H4RG-10 Flats*. Publications of the Astronomical Society of the Pacific, Volume 132, Issue 1007, id. 014502 `Link to abstract <https://ui.adsabs.harvard.edu/abs/2020PASP..132a4502C/abstract>`_
  * Freudenburg, J., Givans, J. et al. (2020). *Brighter-fatter Effect in Near-infrared Detectors—III. Fourier-domain Treatment of Flat Field Correlations and Application to WFIRST*. Publications of the Astronomical Society of the Pacific, Volume 132, Issue 1013, id.074504 `Link to abstract <https://ui.adsabs.harvard.edu/abs/2020PASP..132g4504F/abstract>`_
  * Givans, J. et al. (2022). *Quantum Yield and Charge Diffusion in the Nancy Grace Roman Space Telescope Infrared Detectors*. Publications of the Astronomical Society of the Pacific, Volume 134, Issue 1031, id.014001 `Link to abstract <https://ui.adsabs.harvard.edu/abs/2022PASP..134a4001G/abstract>`_

* Speckle field analysis projects:

  * Hirata, C. & Merchant, C. (2022). *Pixel Centroid Characterization with Laser Speckle and Application to the Nancy Grace Roman Space Telescope Detector Arrays*. Publications of the Astronomical Society of the Pacific, Volume 134, Issue 1041, id.115001 `Link to abstract <https://ui.adsabs.harvard.edu/abs/2022PASP..134k5001H/abstract>`_
  * Macbeth, E., Laliotis, K. et al. (2025), in prep.

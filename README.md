# solid-waffle

## Overview
**solid-waffle** is a repository designed for analysis of correlations in HxRG flats, with an assortment of other features. The repository contains scripts, configuration files, and utilities to support the characterization of HxRG detectors.

*Nota bene*:  These scripts are not ready for "black box" use in the sense that unanticipated defects in different detectors can cause issues in the stability of the characterization.  Some work and testing is required before given results can be trusted.  Manual masking is required in some cases, and we are working towards getting this code more streamlined and user-friendly in that regard. Please approach with caution and reach out to the [contacts](#contacts) listed below for help.

## Table of Contents
1. [Project Structure](#project-structure)
2. [Files Description](#files-description)
   - [ScriptInformation.txt](#scriptinformationtxt)
   - [test_run.py](#test_runpy)
   - [test_run_vis.py](#test_run_vispy)
   - [pyirc.py](#pyircpy)
   - [example_config_wfirst_h4rg_18237](#example_config_wfirst_h4rg_18237)
3. [Directories Overview](#directories-overview)
   - [flat_simulator](#flat_simulator)
   - [notebooks](#notebooks)
4. [How to Get Started](#how-to-get-started)
5. [Contacts](#contacts)
6. [References](#references)

## Project Structure
The repository is organized as follows:


## Files Description

### ScriptInformation.txt
- **Purpose**: Reference information on how to build and run a solid-waffle script. 
- **Details**: Includes sections on how to write a configuration file containing the file names to be analyzed, their format, and analysis options. Also describes different available options for formats (e.g., FITS cube with ascending ramps, and a given extension structure with NAXIS values). There are many different analysis options to choose from.

### test_run.py
- **Purpose**: Main script to call that takes a config file as an argument and runs the analysis.
- **Key Features**:
  - Detector characterization at basic and advanced levels, including nonlinearity and Brighter-Fatter Effect (BFE)
  - Hot pixel identification and calculation of interpixel capacitance (IPC)
  - Summary statistics in superpixels
  - Visualization such as maps of gain, IPC, BFE over superpixels

### test_run_vis.py
- **Purpose**: Similar to `test_run.py` but designed to handle flats and darks in the visible part of the spectrum.
- **Key Features**:
  - If a flag is turned on, the first step of this script is to run `test_run.py` to generate initial processing and results from the infrared data
  - Estimates for charge diffusion covariance and quantum yield, also summarized and visualized in outputs

### pyirc.py
- **Purpose**: Main utility script containing many functions to load in data, calculate statistics, and estimate corrections.
- **Details**: Includes functions that will read the data in, compute reference pixel corrections, compute gain, do basic characterization, calculate correlations of neighboring pixels, compute BFE coefficients, and more.

### example_config_wfirst_h4rg_18237
- **Purpose**: A sample configuration file which is simplified to the basics
- **Details**: This config requires only flat files, dark files, a format code, time frames, and an output location to be specified. `ScriptInformation.txt` contains the information necessary to build a more complex config.

## Directories Overview

### flat_simulator
- **Description**: Scripts and utilities for simulating flat field images.
- **Key Files**:
  - simulate_flat.py: Main script to generate a simulated flat field.
  - ex_sim_config: Example config. Run with `python simulate_flat.py ex_sim_config`
  - detector_function.py: Some utility functions called by simulate_flat.py
  - write_config_flatsim.py: Utility script to quickly generate config files that can be used on the commandline or in a bash script.

### notebooks
- **Description**: Jupyter notebooks for development, analysis, and visualization. None of these are particularly well described currently.
- **Key Notebooks**:
  - plotting_tests.ipynb: Shows some visualizations of outputs.

## How to Get Started
Running the code requires python modules such as numpy, scipy, astropy, matplotlib, and fitsio, although note that we have not extensively tested this across different python versions and cannot guarantee everything will run smoothly.

0. The code expects input flats and darks in FITS format. Specifics on the expected format are described in `ScriptInformation.txt`. Aside from the array dimensions, the script does not use information from the FITS headers, only the image data. If you find that none of the available options work for the format your data is in, you will need to make a new format and associated format code in `pyirc.py` (modifying the functions `get_nside`, `get_num_slices`, and `load_segment` all in the first part of `pyirc.py`).
1. Start with a simple version of the configuration, such as provided in example_config_wfirst_h4rg_18237. You will replace the two placeholder files in the LIGHT section with your flat files, and the ones in the DARK section with dark files. You will also change the FORMAT parameter to the format code that matches your data format as described in `ScriptInformation.txt`. You may also need to adjust the TIME inputs depending on how many frames your files contain.
2. `python test_run.py <yourconfigfile>`
3. Output files will appear in the directory specified in the OUTPUT line of the config.

## Contacts
To communicate about this repository please reach out to:
* Chris Hirata (hirata.10 at osu dot edu)
* Ami Choi (ami.choi at nasa dot gov)

## References
For more detailed background on the concepts and methods used in this project, please refer to:

- Hirata, C. & Choi, A. (2020). *Brighter-fatter Effect in Near-infrared Detectors. I. Theory of Flat Autocorrelations*. Publications of the Astronomical Society of the Pacific, Volume 132, Issue 1007, id. 014501 [Link to abstract](https://ui.adsabs.harvard.edu/abs/2020PASP..132a4501H/abstract)
- Choi, A. & Hirata, C. (2020). *Brighter-fatter Effect in Near-infrared Detectors. II. Autocorrelation Analysis of H4RG-10 Flats*. Publications of the Astronomical Society of the Pacific, Volume 132, Issue 1007, id. 014502 [Link to abstract](https://ui.adsabs.harvard.edu/abs/2020PASP..132a4502C/abstract)
- Freudenburg, J., Givans, J. et al. (2020). *Brighter-fatter Effect in Near-infrared Detectors—III. Fourier-domain Treatment of Flat Field Correlations and Application to WFIRST*. Publications of the Astronomical Society of the Pacific, Volume 132, Issue 1013, id.074504 [Link to abstract](https://ui.adsabs.harvard.edu/abs/2020PASP..132g4504F/abstract)
- Givans, J. et al. (2022). *Quantum Yield and Charge Diffusion in the Nancy Grace Roman Space Telescope Infrared Detectors*. Publications of the Astronomical Society of the Pacific, Volume 134, Issue 1031, id.014001 [Link to abstract](https://ui.adsabs.harvard.edu/abs/2022PASP..134a4001G/abstract)

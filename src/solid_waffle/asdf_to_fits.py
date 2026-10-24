import glob
import os
import sys

import asdf
import numpy as np
from astropy.io import fits

# fields to copy into the FITS header (for Level 1 files)
USEFIELDS_L1 = {
    "ephemeris": {
        "ephemeris_reference_frame": ("E_REF_FR", "Ephemeris Reference Frame"),
        "spatial_x": ("E_SP_X", "X Spatial Coordinate of Roman (km)"),
        "spatial_y": ("E_SP_Y", "Y Spatial Coordinate of Roman (km)"),
        "spatial_z": ("E_SP_Z", "Z Spatial Coordinate of Roman (km)"),
        "time": ("E_TIME", "UTC Time of Ephemeris Information (MJD)"),
        "type": ("E_TYPE", "Ephemeris Type"),
        "velocity_x": ("E_VEL_X", "X Component of Roman Velocity (km/s)"),
        "velocity_y": ("E_VEL_Y", "Y Component of Roman Velocity (km/s)"),
        "velocity_z": ("E_VEL_Z", "Z Component of Roman Velocity (km/s)"),
    },
    "exposure": {
        "effective_exposure_time": ("EFF_EXPT", "Effective Exposure Time (s)"),
        "end_time": ("T_END", "Exposure End Time (UTC)"),
        "engineering_quality": ("QUALITY", "Engineering Data Quality"),
        "exposure_time": ("T_EXP", "Exposure Time (s)"),
        "hga_move": ("HGA_MOVE", "High Gain Antenna (HGA) Movement"),
        "ma_table_id": ("MATB_ID", "Multi-Accumulation Table Identifier"),
        "ma_table_name": ("MATB_NAM", "Name of the Multi-Accumulation Table"),
        "ma_table_number": ("MATB_NMB", "Multi-Accumulation Table Identification Number"),
        "nresultants": ("NRESLTNT", "Number of Resultants"),
        "read_pattern": ("RDP", "Read Pattern"),
        "start_time": ("T_START", "Exposure Start Time (UTC)"),
        "truncated": ("TRUNC", "Truncated MA Table"),
        "type": ("WFI_FLAT", "Exposure Type"),
    },
    "guide_star": {
        "epoch": ("GS_EPOCH", "Guide Star Coordinates Epoch"),
        "guide_mode": ("GS_MODE", "Guide Window Mode"),
        "guide_star_id": ("GS_ID", "Guide Star Identifier"),
        "guide_window_id": ("GS_ID", "Guide Window Identifier"),
        "window_xstart": ("GS_XMIN", "Guide Window X Start Position (pixels)"),
        "window_xstop": ("GS_XMAX", "Guide Window X Stop Position (pixels)"),
        "window_ystart": ("GS_YMIN", "Guide Window Y Start Position (pixels)"),
        "window_ystop": ("GS_YMAX", "Guide Window Y Start Position (pixels)"),
    },
    "instrument": {
        "data_encoding_offset": ("DATAENC", "Data Encoding Offset"),
        "detector": ("DETECTOR", "Detector"),
        "name": ("INSTRUME", "Instrument Name"),
        "optical_element": ("OPTELMNT", "Element Name"),
    },
    "program": {
        "category": ("PRG_CAT", "Program Category"),
        "investigator_name": ("PRG_PI", "Principal Investigator Name"),
        "science_category": ("PRG_SCI", "Science Category"),
        "subcategory": ("PRG_SUBC", "Program Subcategory"),
        "title": ("PRG_TITL", "Proposal Title"),
    },
    "rcs": {
        "active": ("RCS_ON", "Status of the Relative Calibration System (RCS)"),
        "bank": ("RCS_BANK", "Light Emitting Diode (LED) Bank Selection"),
        "counts": ("RCS_CTS", "Light Emitting Diode (LED) Flux (DN)"),
        "electronics": ("RCS_ELEC", "Relative Calibration System (RCS) Electronics Side"),
        "led": ("RCS_LED", "Light Emitting Diode (LED) Passband"),
    },
}


def main(input_dir=None, output_dir=None, fmatch="*.asdf", format="wfi_tvac"):
    """
    Converts all asdf files in a directory to fits files in a new directory

    Parameters
    ----------
    input_dir : str, optional
        The input directory; defaults to the current directory.
    output_dir : str, optional
        The output directory; defaults to a copy of the current directory name with
        ``_fits_converted`` appended.
    fmatch : str, optional
        Regular expression string to match filenames (glob-style); default is to take
        all asdf files.
    format : str, optional
        The input data format. Options are:

        - **"wfi_tvac_rst"** : WFI TVAC data; include the reset/read in the file.
        - **"wfi_tvac"** : WFI TVAC data; do **not** include the reset/read in the file.
        - **"flight_eng"** : Flight data, diagnostic mode; include the reset/read in the file.


    Returns
    -------
    None

    """

    current_dir = os.getcwd()
    parent_dir = os.path.dirname(current_dir)
    current_dir_name = os.path.basename(current_dir)

    # where to get and put files
    if input_dir is None:
        input_dir = current_dir
    if output_dir is None:
        output_dir = os.path.join(parent_dir, current_dir_name + "_fits_converted")

    os.makedirs(output_dir, exist_ok=True)
    asdf_files = glob.glob(os.path.join(input_dir, fmatch))
    print("ASDF files found: ")
    for fn in asdf_files:
        print(fn)
    for fn in asdf_files:
        base = os.path.basename(fn)
        new_file_path = os.path.join(output_dir, base.replace(".asdf", "_asdf_to.fits"))
        if os.path.exists(new_file_path):
            print(f"Skipping (already exists): {new_file_path}")
            continue
        with asdf.open(fn) as f:
            n = np.shape(f["roman"]["data"])[0]
            k = 1 if format in ["wfi_tvac_rst", "flight_eng"] else 0  # number of additional frames
            alldata = np.zeros((n + k, 4096, 4224), dtype=np.uint16)
            data = alldata[k:, :, :]
            data[:, :, :4096] = f["roman"]["data"]
            data[:, :, -128:] = f["roman"]["amp33"]
            if format == "wfi_tvac_rst":
                alldata[0, :, :4096] = f["roman"]["reset_reads"]
                alldata[0, :, -128:] = f["roman"]["amp33_reset_reads"]
            if format == "flight_eng":
                alldata[0, :, :4096] = f["roman"]["reference_read"]
                alldata[0, :, -128:] = f["roman"]["reference_amp33"]
                _offset = (
                    alldata[0].astype(np.int32) - f["roman"]["meta"]["instrument"]["data_encoding_offset"]
                )
                for j in range(n):
                    data[j, :, :] = np.clip(_offset + data[j].astype(np.int32), 0, 65535)
            h = fits.PrimaryHDU(alldata)
            h.header["TFRAME"] = f["roman"]["meta"]["exposure"]["frame_time"]
            h.header["TGROUP"] = f["roman"]["meta"]["exposure"]["frame_time"]

            # extra headers
            for metakey in USEFIELDS_L1:
                if metakey in f["roman"]["meta"]:
                    x = f["roman"]["meta"][metakey]
                    for d in USEFIELDS_L1[metakey]:
                        if d in x:
                            val = x[d]
                            if isinstance(val, list):
                                for j in range(len(val)):
                                    h.header[USEFIELDS_L1[metakey][d][0] + f"{j+1:04d}"] = (
                                        str(val[j]),
                                        USEFIELDS_L1[metakey][d][1] + f" {j+1:04d}",
                                    )
                            else:
                                if not isinstance(val, (int, float, bool, str)):
                                    val = str(val)
                                h.header[USEFIELDS_L1[metakey][d][0]] = (val, USEFIELDS_L1[metakey][d][1])

            # metadata
            ca = list(str(asdf.util.load_yaml(str(fn))))
            h2 = fits.BinTableHDU.from_columns(
                [
                    fits.Column(
                        name="config",
                        format="1A",
                        array=ca,
                    )
                ],
            )
            h2.header["EXTNAME"] = "CONFIG"
            h2.header["ORIGFILE"] = str(fn).split("/")[-1]
            fits.HDUList([h, h2]).writeto(new_file_path, overwrite=True)


if __name__ == "__main__":
    main(format=sys.argv[1] if len(sys.argv) > 1 else "wfi_tvac")

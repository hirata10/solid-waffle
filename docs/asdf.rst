ASDF file support in solid-waffle
#################################

Most Roman data are being distributed in ASDF format. While solid-waffle is capable of reading these directly, it has some scripts that access large numbers of parts of data files, and currently the FITS readers are much faster for this. Therefore, we expect most users will want to convert the data to FITS format (the scripts below convert to `format=1` in solid-waffle convention).

To do so, we provide the ``asdf_to_fits`` utility. This can be used as follows:

.. code-block:: python

    from solid_waffle import asdf_to_fits

    asdf_to_fits.main(
        input_dir="myindata",  # will read from the myindata/ directory
        output_dir="myoutdata",  # will write to the myoutdata/ directory
        format="flight_eng"  # which format; here flight engineering mode
    )

The other commonly used formats would be ``wfi_tvac`` (for the Instrument-level TVAC data) or ``wfi_tvac_rst`` (if you want to include the reset-read frames as well).

The above function will make FITS versions of the data: e.g., a file ``myindata/file1.asdf`` gets mapped to ``myoutdata/file1_asdf_to.fits``.

You may also supply the "fmatch" argument in glob format if you want to only convert some files:

.. code-block:: python

    asdf_to_fits.main(
        input_dir="myindata",  # will read from the myindata/ directory
        output_dir="myoutdata",  # will write to the myoutdata/ directory
        fmatch=f"files0?.asdf",  # will convert files08.asdf, but not files10.asdf
        format="flight_eng"  # which format; here flight engineering mode
    )

The output FITS file has a 3D data cube in the Primary HDU with shape (nreads, 4096, 4224) (the "4224" is because the amp33 data is 
appended to the right side). The TFRAME and TGROUP keywords are also set.

Level 1 metadata is now saved in the FITS header. So an example from a flat field in diagnostic mode would be::

 SIMPLE  =                    T / conforms to FITS standard                      
 BITPIX  =                   16 / array data type                                
 NAXIS   =                    3 / number of array dimensions                     
 NAXIS1  =                 4224                                                  
 NAXIS2  =                 4096                                                  
 NAXIS3  =                   14                                                  
 EXTEND  =                    T                                                  
 BSCALE  =                    1                                                  
 BZERO   =                32768                                                  
 TFRAME  =              3.16247                                                  
 TGROUP  =              3.16247                                                  
 E_REF_FR= 'EME2000 '           / Ephemeris Reference Frame                      
 E_SP_X  =   132775219.96924339 / X Spatial Coordinate of Roman (km)             
 E_SP_Y  =    63420012.52327047 / Y Spatial Coordinate of Roman (km)             
 E_SP_Z  =    27394272.53116869 / Z Spatial Coordinate of Roman (km)             
 E_TIME  =    61334.78680555556 / UTC Time of Ephemeris Information (MJD)        
 E_TYPE  = 'PREDICTED'          / Ephemeris Type                                 
 E_VEL_X =  -14.469271499190212 / X Component of Roman Velocity (km/s)           
 E_VEL_Y =    24.22814950783202 / Y Component of Roman Velocity (km/s)           
 E_VEL_Z =   10.580585299178246 / Z Component of Roman Velocity (km/s)           
 EFF_EXPT=             41.11211 / Effective Exposure Time (s)                    
 T_END   = '2026-10-21T18:53:29.288' / Exposure End Time (UTC)                   
 QUALITY = 'OK      '           / Engineering Data Quality                       
 T_EXP   =             41.11211 / Exposure Time (s)                              
 HGA_MOVE=                    F / High Gain Antenna (HGA) Movement               
 MATB_ID = 'SCI9010 '           / Multi-Accumulation Table Identifier            
 MATB_NAM= 'IM_DIAGNOSTIC'      / Name of the Multi-Accumulation Table           
 MATB_NMB=                 9010 / Multi-Accumulation Table Identification Number 
 NRESLTNT=                   14 / Number of Resultants                           
 RDP0001 = '[1]     '           / Read Pattern 0001                              
 RDP0002 = '[2]     '           / Read Pattern 0002                              
 RDP0003 = '[3]     '           / Read Pattern 0003                              
 RDP0004 = '[4]     '           / Read Pattern 0004                              
 RDP0005 = '[5]     '           / Read Pattern 0005                              
 RDP0006 = '[6]     '           / Read Pattern 0006                              
 RDP0007 = '[7]     '           / Read Pattern 0007                              
 RDP0008 = '[8]     '           / Read Pattern 0008                              
 RDP0009 = '[9]     '           / Read Pattern 0009                              
 RDP0010 = '[10]    '           / Read Pattern 0010                              
 RDP0011 = '[11]    '           / Read Pattern 0011                              
 RDP0012 = '[12]    '           / Read Pattern 0012                              
 RDP0013 = '[13]    '           / Read Pattern 0013                              
 T_START = '2026-10-21T18:52:45.013' / Exposure Start Time (UTC)                 
 TRUNC   =                    F / Truncated MA Table                             
 WFI_FLAT= 'WFI_FLAT'           / Exposure Type                                  
 GS_EPOCH= '?       '           / Guide Star Coordinates Epoch                   
 GS_MODE = 'WIM-TRACK'          / Guide Window Mode                              
 GS_ID   = '01063060060010070011' / Guide Window Identifier                      
 GS_XMIN =                 2902 / Guide Window X Start Position (pixels)         
 GS_XMAX =                 2917 / Guide Window X Stop Position (pixels)          
 GS_YMIN =                 3622 / Guide Window Y Start Position (pixels)         
 GS_YMAX =                 3637 / Guide Window Y Start Position (pixels)         
 DATAENC =                 4000 / Data Encoding Offset                           
 DETECTOR= 'WFI18   '           / Detector                                       
 INSTRUME= 'WFI     '           / Instrument Name                                
 OPTELMNT= 'DARK    '           / Element Name                                   
 PRG_CAT = 'COM     '           / Program Category                               
 PRG_PI  = 'Switzer, Eric'      / Principal Investigator Name                    
 PRG_SCI = 'None    '           / Science Category                               
 PRG_SUBC= 'WFI     '           / Program Subcategory                            
 PRG_TITL= 'CAR-216: RCS post-decon cooldown monitor' / Proposal Title           
 RCS_ON  =                    T / Status of the Relative Calibration System (RCS)
 RCS_BANK= '2       '           / Light Emitting Diode (LED) Bank Selection      
 RCS_CTS =                 1500 / Light Emitting Diode (LED) Flux (DN)           
 RCS_ELEC= 'A       '           / Relative Calibration System (RCS) Electronics S
 RCS_LED = '6       '           / Light Emitting Diode (LED) Passband            
 END                                                                             

There is also a "CONFIG" HDU that saves metadata, so you can write::

  >>> f = fits.open("OUT/im01_asdf_to.fits")
  >>> f["CONFIG"].header

      XTENSION= 'BINTABLE'           / binary table extension                         
      BITPIX  =                    8 / array data type                                
      NAXIS   =                    2 / number of array dimensions                     
      NAXIS1  =                    1 / length of dimension 1                          
      NAXIS2  =                 1188 / length of dimension 2                          
      PCOUNT  =                    0 / number of group parameters                     
      GCOUNT  =                    1 / number of groups                               
      TFIELDS =                    1 / number of table fields                         
      TTYPE1  = 'config  '                                                            
      TFORM1  = '1A      '                                                            
      EXTNAME = 'CONFIG  '                                                            
      ORIGFILE= 'im01.asdf'                                                           
      END                                                                             

  >>> origdict = yaml.safe_load("".join(f["CONFIG"].data["config"]))
  >>> origdict

      {
          'asdf_library': {'author': 'TheASDFDevelopers', 'homepage': 'http://github.com/asdf-format/asdf',
              'name': 'asdf', 'version': '5.3.1'},
          'history': {
              'extensions': [{'extension_class': 'asdf.extension._manifest.ManifestExtension',
                  'extension_uri': 'asdf://asdf-format.org/core/extensions/core-1.6.0',
                  'manifest_software': {'name': 'asdf_standard', 'version': '1.5.0'},
                  'software': {'name': 'asdf', 'version': '5.3.1'}}]},
          'roman': {
              'amp33': {'source': 0, 'datatype': 'uint16', 'byteorder': 'little', 'shape': [4, 4096, 128],
                  'offset': 34611200, 'strides': [34603008, 8448, 2]},
              'amp33_reference_read': {'source': 3, 'datatype': 'uint16', 'byteorder': 'little', 'shape': [1, 4096, 128]},
              'amp33_reset_reads': {'source': 4, 'datatype': 'uint16', 'byteorder': 'little', 'shape': [1, 4096, 128]},
              'data': {'source': 0, 'datatype': 'uint16', 'byteorder': 'little', 'shape': [4, 4096, 4096],
                  'offset': 34603008, 'strides': [34603008, 8448, 2]},
              'meta': {'exposure': {'frame_time': 3.15}},
              'reference_read': {'source': 1, 'datatype': 'uint16', 'byteorder': 'little', 'shape': [1, 4096, 4096]},
              'reset_reads': {'source': 2, 'datatype': 'uint16', 'byteorder': 'little', 'shape': [1, 4096, 4096]}}
      }

(I added some newlines in the outputs for readability on the GitHub viewer.)

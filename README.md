# mwPyRttov
If you need to build your own version, below the steps to build RTTOV and the Python wrapper

1) RTTOV must be downloaded from https://nwp-saf.eumetsat.int/ (registration is needed)

2) To use the Python wrapper you must have f2py (part of NumPy) installed. The wrapper is compatible only with Python3. 

3) create a new directory for the package
   $ mkdir rttov13

4) Copy the compressed RTTOV package file (rttov.tar.xz) into this new directory and extract the contents:
   $ cp rttov.tar.xz rttov13
   $ cd rttov13
   $ tar xf rttov131.tar.xz

5) RTTOV can be compiled without any external dependencies. However to make use of all functionality you should edit the file 
   build/Makefile.local with the location of your HDF5 installation (and optionally other external libraries) before compiling.
   For instance, HDF5 library is needed to access to the emissivity atlas

6) use the interactive script (recommended) to compile:

   $ cd src
   $ ../build/rttov_compile.sh
   
	example of interactive steps:

	Specify required compiler flag file (leave blank for default: gfortran)
	> gfortran-openmp

	Specify installation directory relative to top-level RTTOV directory (leave blank for default: ./)
	>

	Checking ../build/Makefile.local for user-specified HDF5 library...
	...found FFLAGS_HDF5 and LDFLAGS_HDF5: compiling with HDF5 library

	Checking ../build/Makefile.local for user-specified NetCDF library...
	...did not find FFLAGS_NETCDF or LDFLAGS_NETCDF: compiling without NETCDF library
	If you want to compile with NETCDF you must specify this in ../build/Makefile.local first

	Checking ../build/Makefile.local for user-specified LAPACK library...
	...did not find FFLAGS_LAPACK or LDFLAGS_LAPACK: compiling with lapack.f included in RTTOV

	Testing for f2py using 'f2py --fcompiler=gnu95'...
	...f2py detected: do you want to compile the Python wrapper and RTTOV GUI? (y/n)
	> y

	RTTOV COMPILATION SUMMARY 
	
	Compiling with flags           : gfortran-openmp
	Compiling in directory         : ./

	RTTOV features available:
	HDF5 coefficient I/O           : y
	Emissivity/BRDF atlases        : y
	C/C++ wrapper                  : y
	Python wrapper                 : y
	RTTOV GUI                      : y
	HTFRTC netCDF files            : n

	Compiling with user LAPACK lib : n

	Regenerating Makefiles using:
	$ ../build/Makefile.PL RTTOV_HDF=1 RTTOV_F2PY=1 RTTOV_USER_LAPACK=0

	Compiling RTTOV using:
	$ make ARCH=gfortran-openmp INSTALLDIR=./ 

	OK to continue and compile RTTOV? (y/n)
	> y

7) test that the build is successful (run some tests)

   $ cd ../rttov_test
   $ ./test_rttov13.sh ARCH=gfortran-openmp

8) test that that the python wrapper works OK (run a python example)

   $ cd ../wrapper
   $ python pyrttov_visirscatt_example.py
   
# Additional requirements

Mandatory for RTTOV-SCATT:
in the directory 'rtcoef_rttov13/hydrotable' the hydrotable*.dat file must be available for the specific MW sensor (e.g. hydrotable_gcom-w_amsr2.dat). 
hydrotable*.dat files can be downloaded from https://nwp-saf.eumetsat.int/site/software/rttov/download

Optional: to use the MW surface emissivity atlas (RTTOV must be build including the HDF5 library).
TELSEM atlas data must be located in the directory 'emis_data' and they can be downloaded from https://nwp-saf.eumetsat.int/site/software/rttov/download

# How to run
As an example to use this framework you can have a look at run_mw_rttov13v0.py or example_run_rttov.ipynb. 
Generally what you need:

-  The Python wrapper build as described above
-  The RTTOV rttov_wrapper_f2py.so library must be in your $PYTHONPATH
-  NWP data (1D vector), interpolated at the satellite location
-  You can choose to run RTTOV or RTTOV-SCATT for AMSR2 or SSMIS
-  mw_rttov_cfg.py: some parameters to configure for RTTOV (here you may also add addiitonal satellite)  
-  What you get after runing the run_rttov function:
  For each selected sensor's channel:
    
	simTb:       simulated brightness temperature

	simTb_noAtm: (This was implemented for the SIC algorithm) simulated brightness temperature with no Atm (e.g.: u10=0, v10=0, profQ=0, profCLWC=0, profCIWC=0, profCRWC=0, profCSWC=0, profCC=0)

	surfEmis:    surface emissivty used in the radiative transfer calculations

	strChsList:  list of channels (string that identifies GHz + pol)

	emis_param_dict: internal RTTOV data structure

	RTTOV-SCATT

		{'EmisTermsDownCld','EmisTermsUpCld','EmisTermsTauCld','EmisTermsDownClr','EmisTermsUpClr','EmisTermsTauClr','EmisTermsBsfc','EmisTermsCfrac'}
   
  	RTTOV
   
    		{'Rad2Down','Rad2Up','Rad2DnClear','Rad2UpClear','Rad2Surf','Rad2ReflDnClear'}


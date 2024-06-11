###############################################################################
# Author: Fabrizio Baordo (DMI)

# Follow the steps described in the README to generate the rttov_wrapper_f2py.so library (this example was based using RTTOV v13).

# This module can be used to:
    
# - simulate microwave brightness temperatures (BT) for different sensors 
# - run sensitivity study to evaluate RTTOV and RTTOV-SCATT 
# - run sensitivity study to evaluate the impact of the surface emissvity on the simulated BT 
# - evaluate the magnitude of the atm correction which is used for the OSISAF sic algorithm

###############################################################################

import sys
# The pyrttov package and the RTTOV rttov_wrapper_f2py.so library must be in your $PYTHONPATH
# Add location of pyrttov and rttov_wrapper_f2py.so to the python path
# WARNING: set RTTOV version (rttovVer) in mw_rttov_cfg.py 
sys.path.append('/home/fab/rttov130/wrapper')

import os
import xarray as xr
from datetime import datetime
from rttov13v0_interface import run_rttov

import logging

LOG = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='[%(levelname)s: %(asctime)s: %(name)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

############
# Inputs
############

# WARNING: mw_rttov_cfg.py contains 
# - some variables which can be used to tune RTTOV  
# - the logic to handle the satellite sensor according to the variable strSatId (e.g. channels to be simulated, satZenAngle, rttov files to get, constant ice emis)

# String to identify the satellite to simulate.
# The following strings are implemented in mw_rttov_cfg.py: amsr2, ssmis_f16 ssmis_f17, ssmis_f18
# To add any other satellite add a string and the necessary logic in mw_rttov_cfg.py 
# WARNING: in this implemenatation we consider a constant satZenAngle which is defined in the mw_rttov_cfg.py according to the sensor
strSatId = 'amsr2'

topPath    = '/home/fab/Desktop/ieee_rtm_corr_sic/nwp_swath_era5/20190225'
swath_time = '201902250529'

# Collocated NWP surface variables
ncSurfNwp2swath = os.path.join(topPath,swath_time, 'surf.era5TOamsr2.'+swath_time+'.nc')
# Collocated NWP model level variables
ncProfNwp2swath = os.path.join(topPath,swath_time, 'ml.era5TOamsr2.'+swath_time+'.nc')

# Expected RTTOV path which are in your RTTOV package. Within this paths we expect necessary input files to run RTTTOV  
# Mandatory: where to find rttov coef files for the mw sensor which must be simulated (e.g. rtcoef_gcom-w_1_amsr2.dat)
rttovCoefDir = '/home/fab/rttov132/rtcoef_rttov13/rttov13pred54L'
# Mandatory when rttovScatt = True: where to find hydrotable file for the mw sensor which must be simulated (e.g. hydrotable_gcom-w_amsr2.dat)
rttovHydDir  = '/home/fab/rttov132/rtcoef_rttov13/hydrotable'
# Mandatory when useEmisAtlas = True: where to find expected MW emissvity atlas files (hdf5) 
rttovAtlasEmisDir = '/home/fab/rttov132/emis_data'

# False: simulate atmospheric absorption (T, Q, CLWC are only used as NWP profile) 
# True: use RTTOV-SCATT to active the calculation of scattering (CIWC, CRWC, CSWC and CC are also mandatory as NWP profiles)
rttovScatt = True

# Emis is computed as Emis = (1-sic)*EmisOcean + sic*EmisIce, where
# EmisOcean is extracted from RTTOV (e.g. FATSEM emis) 
# EmisIce is:  
#    useEmisAtlas = False: use a constant values (defined in mw_rttov_cfg.py) or (to be implemented) any other estimate/model     
#    useEmisAtlas = True: use TELSEM MW emissvity atlas
useEmisAtlas = False

# active calculation of BT with 'no Atmosphere' (necessary for SIC algorithm)
# if True simTb_noAtm is populated, otherwise return empty list
calbBTnoAtm = False

###############################################################################

# NWP surface & model level variables
# Model To observation sapce interpolation: we collocate NWP variables 
# on the satellite swath (pyresample is used). For every observation location (lat/lon), 
# we have the corresponding model value. NWP AN/FC date/time is simply the closest 
# to the central time of the satellite scan

# Surface
# lats (1D, dim=nobs): latitude of observation location
# lons (1D, dim=nobs): langitude of observation location
# sic (1D, dim=nobs): sea ice cocnetration [0-1]
# tsk (1D, dim=nobs): skin temperature [K]
# t2m (1D, dim=nobs): 2m temperature [K]
# u10 (1D, dim=nobs): 10m wind u-component [m s-1]
# v10 (1D, dim=nobs): 10m wind v-component [m s-1]
# surfPre (1D, dim=nobs): surface pressure [Pa]
# Model level ('Profile')
# profT (2D, dim=(nobs,nlevFull)): temperature [K] on model level
# profQ (2D, dim=(nobs,nlevFull)): specific humidity [kg kg-1] on model level
# profCLWC (2D, dim=(nobs,nlevFull)): cloud liquid water content [kg kg-1] on model level
# profCIWC (2D, dim=(nobs,nlevFull)): cloud ice water content [kg kg-1] on model level
# profCRWC (2D, dim=(nobs,nlevFull)): specific rain water content [kg kg-1] on model level
# profCSWC (2D, dim=(nobs,nlevFull)): specific snow water content [kg kg-1] on model level
# profCC (2D, dim=(nobs,nlevFull)): fraction of Cloud Cover [0-1] on model level
# hPaFullLevels (2D, dim=(nobs,nlevFull)): full model level [hPa]
# hPaHalfLevels(2D, dim=(nobs,nlevHalf)): half model level [hPa]

# When RTTOV-SCATT is used (rttovScatt=True) the following variables are mandatory: 
# hPaHalfLevels, profCIWC, profCRWC, profCSWC, profCC

##########
# Outputs 
##########

# Using RTTOV, we calculate the simulated BT and the reference simulated BT with 
# no atmosphere so that we can calculate the atmospheric correction to apply to the observed BT. 
# For assessement, we also return the surface emissvity which is used in the 
# radiative transfer calculations. 

# For each selected sensor's channel:
    
# simTb:       simulated brightness temperature
# simTb_noAtm: simulated brightness temperature - no Atm: u10=0, v10=0, profQ=0, profCLWC=0, profCIWC=0, profCRWC=0, profCSWC=0, profCC=0
# surfEmis:    surface emissivty used in the radiative transfer calculations 
# strChsList:  list of channels (string that identifies GHz + pol) 
# emis_param_dict: internal RTTOV data structure
#################
# Load Variables
#################
# Collocated NWP surface variables
surfNwp2swath = xr.open_dataset(ncSurfNwp2swath)
# Collocated NWP model level variables
profNwp2swath = xr.open_dataset(ncProfNwp2swath)

# Surface
lats = surfNwp2swath['lat'].data
lons = surfNwp2swath['lon'].data
sic  = surfNwp2swath['sic'].data
skt  = surfNwp2swath['skt'].data
t2m  = surfNwp2swath['t2m'].data
u10  = surfNwp2swath['u10'].data
v10  = surfNwp2swath['v10'].data
sp   = surfNwp2swath['sp'].data

# sic must be [0-1]
if sic.max() > 1:
    sic = sic/100.

# Model Level
profT = profNwp2swath['t'].data
profQ = profNwp2swath['q'].data
profCLWC = profNwp2swath['clwc'].data
profCIWC = profNwp2swath['ciwc'].data
profCRWC = profNwp2swath['crwc'].data
profCSWC = profNwp2swath['cswc'].data
profCC   = profNwp2swath['cc'].data
hPaFullLevels = profNwp2swath['pf'].data
hPaHalfLevels = profNwp2swath['ph'].data

strYYYYMMDDHHMMs = surfNwp2swath.sensor_scan_start_time
strYYYYMMDDHHMMe = surfNwp2swath.sensor_scan_end_time

scanStartTime = datetime.strptime(strYYYYMMDDHHMMs, '%Y%d%m%H%M')
scanEndTime   = datetime.strptime(strYYYYMMDDHHMMe, '%Y%d%m%H%M')
centralSwathTime = scanStartTime + (scanEndTime - scanStartTime)/2 

# if available, provide satZenAngle at every obs location (np.array, 1D, dim=nobs)
# otherwise, satZenAngle = [0] and a constant value will be used according to the satellite sensor
satZenAngle = [0]

LOG.info("*** Run RTTOV *** ")
LOG.info('Inputs:')
LOG.info(' - sensor: {}'.format(strSatId))
LOG.info(' - nc nwp surface variables: {}'.format(ncSurfNwp2swath))
LOG.info(' - nc nwp model level variables: {}'.format(ncProfNwp2swath))
LOG.info(' - rttovScatt: {}'.format(rttovScatt))
LOG.info(' - useEmisAtlas: {}'.format(useEmisAtlas))
LOG.info(' - calbBTnoAtm: {}'.format(calbBTnoAtm))
LOG.info(' - rttov coef dir: {}'.format(rttovCoefDir))
LOG.info(' - rttov hydrotable dir: {}'.format(rttovHydDir))
LOG.info(' - rttov mw emis atlas dir: {}'.format(rttovAtlasEmisDir))

simTb, simTb_noAtm, surfEmis, strChsList, emis_param_dict = run_rttov(lats,
                                                                      lons,
                                                                      sic,
                                                                      skt,
                                                                      t2m,
                                                                      u10,
                                                                      v10,
                                                                      sp,
                                                                      profT,
                                                                      profQ,
                                                                      profCLWC,
                                                                      profCIWC,
                                                                      profCRWC,
                                                                      profCSWC,
                                                                      profCC,
                                                                      hPaFullLevels,
                                                                      hPaHalfLevels,
                                                                      strSatId,
                                                                      satZenAngle,
                                                                      centralSwathTime,
                                                                      rttovCoefDir,
                                                                      rttovHydDir,
                                                                      rttovAtlasEmisDir,
                                                                      rttovScatt,
                                                                      useEmisAtlas,
                                                                      calbBTnoAtm) 

LOG.info('RTM calculations done!')




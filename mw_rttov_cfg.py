###############################################################################
# This configuration is for
# - some variables which can be used to tune RTTOV  
# - handling the logic for the satellite sensor according to the variable strSatId (e.g. channels to be simulated, satZenAngle, rttov files to get, constant ice emis)
###############################################################################

rttovVer = 'v13.0'

# To check the RTTOV options which are configured before calling the forward model (e.g. rttov runDirect()) 
logRTTOVOptions = False

# if (sicFgThr > 0) tune the surface emis to use according to the ocean/seaice surface:
# where fgSic > sicFgThr:
#    Emis = (1-sic)*EmisOcean + sic*EmisIce
# otherwise
#    Emis = EmisOcean
sicFgThr = 0

###############
# RTTOV options 
###############

# Units of gas in the input profile  
# 0=ppmv over dry air; 1=kg/kg & 2=ppmv over moist air 
# kg/kg*10^6 = ppmv
GasUnits = 2

# if IcePolarisation > 0, enables the approximation for polarised scattering in RTTOV-SCATT simulations
# Default 1.4
IcePolarisation = 1.40

# This is to include contribution of cloud liquid water in RT calculations (NWP profile T,Q and CLW).
# By default, RTTOV ignores any CLW content above 322hPa. This limit is
# specified in the 'CLWCloudTop' option allowing you to modify it if you wish
CLWdata = True

# Activate lambertian option (not active for RTTOV-SCATT)
# DoLambertian: If true, activate treatment of surface as Lambertian instead of
# specular reflector for downwelling emitted radiance (default = false).
DoLambertian = True
# LambertianFixedAngle: If true, the Lambertian downwelling radiance is computed using a
# fixed effective angle. If false the angle is derived from a parameterisation computed 
# from the total atmospheric optical depth (default = true).
LambertianFixedAngle = True

#RTTOV interpolation mode (NWP pressure level variables)
#1:     Default, original RTTOV interpolation method,
#		Jacobians may show oscillations. Reasonable choice
#		when user levels are sparse compared to coef levels
#2:		Log-linear on optical. May be beneficial 
#		in execution time for direct-model depths
#		calculations, but not suitable for TL/AD/K
#3:		Similar to mode 1, but with somewhat reduced
#		oscillations. Reasonable choice when user levels are
#		sparse compared to coef levels
#4:		No oscillations, but most computationally expensive
#		method. Reasonable choice when user levels are dense
#		compared to coef levels
#5: 	No oscillations, but Jacobians may show small
#		artefacts due to interpolation, only slightly more
#		expensive than mode 1. Reasonable choice when user
#		levels are dense compared to coef levels.
InterpMode = 3

# use parallel computation 
Nthreads = 4

# Store variables in rttov data structure
StoreRad = True
# Store variables in rttov data structure (not actitve for rttov-scatt)
StoreRad2 = True
# store trans terms (not actitve for rttov-scatt)
StoreTrans = False  # memory problem on the laptop
# store emis term (active for rttov-scatt)
StoreEmisTerms = True    

# on/off checks and logs 
Verbose        = False
VerboseWrapper = False
# This must be False to run RTTOV with no atm
DoCheckinput   = False

# Emis over ocean: 0 = TESSEM2, 6 = Latest FASTEM
# only from v13.2: FastemVersion = 7--> for SURFEM-Ocean
FastemVersion = 6

###############################################################
# To process a new sensor with RTTOV add info in fucntion below
###############################################################
def get_mw_sensor_cfg(strSatId):

    # AMSR2
    if (strSatId == 'amsr2'):

        # hydrotable (used only by RTTOV-SCATT)
        hyd_file = 'hydrotable_gcom-w_amsr2.dat'
        # coef file        
        rttov_coef_file = 'rtcoef_gcom-w_1_amsr2.dat'
        # 0/1 active channels for RTM calculations (channels of interest for SIC: CH7: 18.7V, CH8: 18.7H; CH11: 36.5V; CH12: 36.5H, CH13: 89V; CH14: 89H)
        #activeChs =  [0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1]
        activeChs =  [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 1, 1, 1, 1]
        # string to identify active channels 
        #strChannels =  ['', '', '', '', '', '', '18.7V', '18.7H', '', '', '36.5V', '36.5H', '89V', '89H']
        strChannels =  ['6.9V', '6.9H', '7.3V', '7.3H', '10.65V', '10.65H', '18.7V', '18.7H', ' ',' ','36.5V', '36.5H', '89V', '89H']
        # ice emissivity, const values as those used by RT Wentz & Meissner
        iceEmis = [0.92, 0.85, 0.92, 0.85, 0.92, 0.85, 0.92, 0.85, -1, -1, 0.78, 0.7, 0.75, 0.68]
        # constant value for sat zen angle
        satZenAngle = 55.0
        
    # SSMIS16, SSMIS17, SSMIS18
    elif (strSatId == 'ssmis_f16') or (strSatId == 'ssmis_f17') or (strSatId == 'ssmis_f18'):

        # hydrotable (used only by RTTOV-SCATT)
        hyd_file = 'hydrotable_dmsp_ssmis.dat'
        # coef file
        rttov_coef_file = "rtcoef_dmsp_{}_ssmis.dat".format(strSatId.split('_')[1][1:3])
        # 0/1 active channels for RTM calculations (channels of interest for SIC: CH12: 19.35H, CH13: 19.35V; CH15: 37H, CH16: 37V; CH17: 91.65V, CH18: 91.65H)
        activeChs = [0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  1,  1,  1,  1,  0,  0,  0,  0,  0,  0]
        # string to identify active channels 
        strChannels = ['',  '',  '',  '',  '',  '',  '',  '',  '',  '',  '',  '19.35H', '19.35V',  '',  '37H', '37V',  '91.65V', '91.65H',  '',  '',  '',  '',  '',  '']
        # ice emissivity, const values as those used by RT Wentz & Meissner
        iceEmis = [-1,  -1,  -1, -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  0.90, 0.95, -1,  0.88, 0.93, 0.80, 0.75,  -1,  -1,  -1,  -1,  -1, -1]
        # constant value for sat zen angle
        satZenAngle = 53.1

    else:
        raise Exception("Logic for sensor {} not implemented!".format(strSatId))

    return rttov_coef_file, hyd_file, activeChs, strChannels, iceEmis, satZenAngle




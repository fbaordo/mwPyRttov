# Python wrapper for RTTOV
import pyrttov
import numpy as np
import mw_rttov_cfg as cfg
import scipy.constants as const

import logging

LOG = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='[%(levelname)s: %(asctime)s: %(name)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

##########################
# setting for this module
##########################
speedl = const.speed_of_light*100.

# if (sicThr > 0) identify ocean/seaice surface such as:

# where fgSic > sicThr:
#    Emis = (1-sic)*EmisOcean + sic*EmisIce
# otherwise
#    Emis = EmisOcean
sicThr = cfg.sicFgThr

#rttov version
rttovVer = cfg.rttovVer

############
# FUNCTIONS
############

def check_rttov_options(mwSatRttov):

    allOptionsList = dir(mwSatRttov.Options)
    LOG.info("**** Check RTTOV Options ****")
    for nO in range(0,len(allOptionsList)):
        strOption = allOptionsList[nO]
        if ('_' not in strOption) and ('defineStrOptions' not in strOption):
            strCheckOption = str(getattr(mwSatRttov.Options, strOption))
            LOG.info("%-20s --> %-6s" % (strOption, strCheckOption))
    LOG.info("*****************************")

    return None
########################################################
# Interface to call RTTOV and compute the BT (T+Q+CLWC)
########################################################
def calc_bt(nwp_rttov_lats,
            nwp_rttov_lons,
            nwp_rttov_fg_sic,
            nwp_rttov_skTemp,
            nwp_rttov_2mTemp,
            nwp_rttov_u10,
            nwp_rttov_v10,
            nwp_rttov_surfPre,
            nwp_rttov_T,
            nwp_rttov_Q,
            nwp_rttov_CLWC,
            modelFullLevhPa,
            strSatId,
            swathDate,
            satZenAngle,
            rttovCoefDir,
            rttovAtlasEmisDir,
            emisToRttov,
            satAzimuthAngle=None):

    # get sensor dependent rtCoef&channels
    rtCoefFile, rtHydFile, chsToKeepList, strChsList, emisIceChsList, constSatZenAngle = cfg.get_mw_sensor_cfg(strSatId)

    #########################
    # Subset channels
    #########################
    chsToProcessList = []
    chsStringList    = []

    for nCh in range(0, len(chsToKeepList)):

        if (chsToKeepList[nCh]):

            chsToProcessList.append(nCh+1)
            chsStringList.append(strChsList[nCh])

    #To be consistent with the calculation of surface emissvity
    if (sicThr > 0):
        indexIce   = (nwp_rttov_fg_sic > sicThr)
        indexOcean = (nwp_rttov_fg_sic <= sicThr)
    else:
        indexIce   = (nwp_rttov_fg_sic > 0)
        indexOcean = (nwp_rttov_fg_sic == 0)
                

    # Check possible negative or too small humidity value in NWP profile
    # --> put a very small value so RTTOV does not complain for 'exceeded allowed minimun'
    check_qprofile = (nwp_rttov_Q < 0.1000E-010)
    if ((check_qprofile == True).sum() != 0):
        nwp_rttov_Q[check_qprofile] = 0.1000E-09

    # ------------------------------------------------------------------------
    # Set up the profile data NWP
    # ------------------------------------------------------------------------

    # Pressure levels expected in hPa ordered from top atm level (python index 0) to surface (python last index)

    # Declare an instance of Profiles
    nlevels   = modelFullLevhPa.shape[1]
    nprofiles = modelFullLevhPa.shape[0]
    myProfiles = pyrttov.Profiles(nprofiles, nlevels)

    allZeros = np.asarray([0]*nprofiles)
    allOnes  = np.asarray([1]*nprofiles)

    # Profiles: 2D array [nprofiles, nlevels]
    myProfiles.GasUnits = cfg.GasUnits
    myProfiles.P = modelFullLevhPa
    myProfiles.T = nwp_rttov_T
    myProfiles.Q = nwp_rttov_Q
    myProfiles.CLW = nwp_rttov_CLWC

    # Surface: 2D array [nprofiles, nVars]

    # Populate SurfGeom data structure: Real [nprofiles][3] (latitude, longitude, elevation)
    myelevs = allZeros  # elevation is in km
    surfgeom = np.array([nwp_rttov_lats, nwp_rttov_lons, myelevs], dtype=np.float64)
    myProfiles.SurfGeom = surfgeom.transpose()

    # Populate DateTimes data structure: Integer [nprofiles][6] (Not  really necessary for MW)
    # The full date will be used to calculate the TOA solar irradiance for solar-affected simulations.
    # dts = [timestamp.year, timestamp.month, timestamp.day, timestamp.hour, timestamp.minute, timestamp.second]

    mydateTimeYYYY = np.array([swathDate.year]*nprofiles, dtype=np.int32)
    mydateTimeMM   = np.array([swathDate.month]*nprofiles, dtype=np.int32)
    mydateTimeDD   = np.array([swathDate.day]*nprofiles, dtype=np.int32)
    mydateTimeHH   = np.array([swathDate.hour]*nprofiles, dtype=np.int32)
    mydateTimeMM   = np.array([swathDate.minute]*nprofiles, dtype=np.int32)
    mydateTimeSS   = np.array([swathDate.second]*nprofiles, dtype=np.int32)

    datetimes = np.array([mydateTimeYYYY,
                          mydateTimeMM,
                          mydateTimeDD,
                          mydateTimeHH,
                          mydateTimeMM,
                          mydateTimeSS], dtype=np.int32)

    myProfiles.DateTimes = datetimes.transpose()

    # Populate Angles data structure: Real [nprofiles][4]: satzen (mandatory), satazi, sunzen, sunazi
    # satzen: Local satellite zenith angle (degrees) (mandatory)
    # satazi: Local satellite azimuth angle (0-360; measured clockwise, east=90deegre           (not mandatory)
    # sunzen: Local solar zenith angle (degrees), solar radiation only included up to 85deegre  (not mandatory)
    # sunazi: Local solar azimuth angle (0-360; measured clockwise, east=90deegre               (not mandatory)
    # NOTE: To turn on solar radiation for visible and near-IR channel,
    #       solar zenith angle, and satellite and solar azimuth angles must also be specified

    if len(satZenAngle) == 1:
        allSatZenithAngle = np.array([satZenAngle[0]]*nprofiles, dtype=np.float64)
    else:
        allSatZenithAngle = satZenAngle
    
    if satAzimuthAngle is not None:
        allSatAzimuthAngle = satAzimuthAngle   
    else:
        allSatAzimuthAngle = allZeros   

    angles = np.array([allSatZenithAngle, allSatAzimuthAngle, allZeros, allZeros], dtype=np.float64)
    myProfiles.Angles = angles.transpose()

    # Populate SurfType data structure:  Integer [nprofiles][2] (surftype, watertype)
    # surftype: (0=land, 1=sea, 2=seaice) watertype 0=fresh water, 1=ocean water.

    # set all ocean points
    rttovSurfType = np.asarray([1]*nprofiles, dtype=np.int32)

    # repalce ocean with ice
    rttovSurfType[indexIce] = 2

    surftype = np.array([rttovSurfType,
                         allOnes], dtype=np.int32)

    myProfiles.SurfType = surftype.transpose()

    # Populate S2m data structure: Real [nprofiles][6] (s2m%p, s2m%t, s2m%q, s2m%u, s2m%v, s2m%wfetc)

    windfetch = np.asarray([100000]*nprofiles)  # Wind fetch, value 100000m for open ocean – only used by sea surface solar BRDF model

    # if we have 2m dewpoint. Calculate: 2m specific humidity kg/kg
    # where: e = vapor pressure in mb; Td = dew point in deg C; p = surface pressure in mb (=hPa); q = specific humidity in kg/kg.
    #e_mb  = [6.112*np.exp( (17.67*(Td-273.15))/((Td-273.15) + 243.5) ) for Td in my2mDewTemp]
    #my2mQ = [0.001*(0.622 * e_mb[i])/(p - (0.378 * e_mb[i])) for i,p in enumerate(mySurfPres)]

    # for now simply get the value of the model level closest to surface
    my2mQ = nwp_rttov_Q[:,modelFullLevhPa.shape[1]-1]

    s2m = np.array([nwp_rttov_surfPre/100.,    # Surface pressure (hPa)
                    nwp_rttov_2mTemp,         # 2m temperature
                    my2mQ,                    # 2m Q
                    nwp_rttov_u10,            # 10m wind u component
                    nwp_rttov_v10,            # 10m wind v component
                    windfetch], dtype=np.float64)

    myProfiles.S2m = s2m.transpose()

    # Populate skin data structure: Real [nprofiles][9] (skinTemp, salinity, snow_fraction, foam_fraction, fastem(1:5))
    # - salinity: use constant value of 35 psu (used by FASTEM 4-6, TESSEM2)
    # - snow_fraction: Surface snow coverage fraction (0-1) (used by IR emissivity atlas)
    # - foam_fraction: used by FASTEM if mwSatRttov.Options.SupplyFoamFraction = True
    # - fastem(1:5): FASTEM for land/sea-ice surface types --> not necessary

    allSalinity = np.asarray([35]*nprofiles)

    skin = np.array([nwp_rttov_skTemp,
                     allSalinity,
                     allZeros,
                     allZeros,
                     allZeros,
                     allZeros,
                     allZeros,
                     allZeros,
                     allZeros], dtype=np.float64)

    myProfiles.Skin = skin.transpose()

    # ------------------------------------------------------------------------
    # Set up Rttov instances for an instrument
    # ------------------------------------------------------------------------
    mwSatRttov = pyrttov.Rttov()

    # Config RTTOV options

    # To include contribution of cloud liquid water in RT calculations
    mwSatRttov.Options.CLWData = cfg.CLWdata
    # To change treatment of specularity (together with surfemisrefl[3,:,:])
    mwSatRttov.Options.DoLambertian = cfg.DoLambertian
    mwSatRttov.Options.LambertianFixedAngle = cfg.LambertianFixedAngle
    # emis over ocean
    mwSatRttov.Options.FastemVersion   = cfg.FastemVersion
    # Interpolation to RTTOV pressure level
    mwSatRttov.Options.InterpMode = cfg.InterpMode
    # parallel calculation
    mwSatRttov.Options.Nthreads = cfg.Nthreads
    # store radiances
    mwSatRttov.Options.StoreRad  = cfg.StoreRad
    mwSatRttov.Options.StoreRad2 = cfg.StoreRad2    
    # store trans terms
    mwSatRttov.Options.StoreTrans = cfg.StoreTrans        
    # verbose on(off)
    mwSatRttov.Options.Verbose        = cfg.Verbose
    mwSatRttov.Options.VerboseWrapper = cfg.VerboseWrapper

    # Load the instruments
    try:
        mwSatRttov.FileCoef = '{}/{}'.format(rttovCoefDir,rtCoefFile)
        mwSatRttov.loadInst(chsToProcessList)
    except pyrttov.RttovError as e:
        LOG.error("Error loading instrument(s): {!s}".format(e))
        raise Exception("Error loading instrument(s): {!s}".format(e))

    # sensor dependednt
    nchans    = mwSatRttov.Nchannels
    #satHeight = pyrttov.rtwrap.rttov_get_coef_val_r0(mwSatRttov.InstId,'FC_SAT_HEIGHT')[1]
    ff_cwn    = mwSatRttov.WaveNumbers
    freq_ghz  = ff_cwn * speedl * 1.0E-09

    LOG.info("MW sensor: {}". format(strSatId))
    LOG.info("Rttov coef file: {}".format(mwSatRttov.FileCoef))
    LOG.info("Number of channels selected for RTTOV: {}".format(nchans))
    LOG.info("- List of selected channels: {}".format(chsToProcessList))
    LOG.info("- Freq (GHz): {}".format(freq_ghz))
    LOG.info("Gas units: {}".format(cfg.GasUnits))

    # Associate the profiles with each Rttov instance
    mwSatRttov.Profiles = myProfiles

    # surfemisrefl: Set up the surface emissivity/reflectance arrays and associate with the Rttov objects
    # 0 = control claculation of surface emis (-1 --> calcemis == TRUE)
    # 1 = control calculation of reflectance (-1 --> calcrefl == TRUE )
    # 2 = to provide diffuse reflectance values to rttov below 3um where calcrefl is False
    # 3 = to specify the surface specularity which is used by RTTOV only if the option DoLambertian is True
    # 4 (only from version v13.2): the surfemisrefl(4,:,:) array is used to specify the per-channel effective skin temperatures which are used with the use_tskin_eff option

    # Default: surface is treated as fully specular reflector --> specularity = 1 (DoLambertian == false)
    # fully Lambertian --> specularity = 0

    # For profiles/channels where an internal sea surface emissivity model is being used (i.e. FASTEM or TESSEM2 in
    # the MW) the Lambertian option is not valid and so is not applied for these channels.
    # When activated the Lambertian option is applied for sea surfaces where calcemis(:) is false and for land and sea-ice
    # surfaces regardless of calcemis(:).

    #Surface emissivity/reflectance arrays must be initialised *before every call to RTTOV*
    if rttovVer == 'v13.2':
        surfemisrefl = np.zeros((5,nprofiles,nchans), dtype=np.float64)
    else:
        surfemisrefl = np.zeros((4,nprofiles,nchans), dtype=np.float64)

    mwSatRttov.SurfEmisRefl = surfemisrefl

    # Surface emissivity/reflectance arrays must be initialised *before every call to RTTOV*
    # Negative values will cause RTTOV to supply emissivity/BRDF values
    surfemisrefl[:,:,:] = -1.

    surfemisrefl[0,:,:] = emisToRttov

    # in case we activate the Lambertian option
    # Default: surface is treated as fully specular reflector --> specularity = 1 (DoLambertian == false)
    # fully Lambertian --> specularity = 0
    if (cfg.DoLambertian):
        surfemisrefl[3,indexIce,:]   = 0
        surfemisrefl[3,indexOcean,:] = 1

    # Check rttov options
    if (cfg.logRTTOVOptions):
        check_rttov_options(mwSatRttov)

    try:
        LOG.info("****************************************************************")        
        LOG.info("RTTOV runDirect started (compute BT for NWP profile: T+Q+CLWC)")
        LOG.info("CLWC data: {}".format(cfg.CLWdata))
        LOG.info("Tot number of profiles: {}".format(nprofiles))
        LOG.info("Tot number of threads: {}".format(cfg.Nthreads))
        LOG.info("Do Lambertian: {}".format(cfg.DoLambertian))
        mwSatRttov.runDirect()
        LOG.info("RTTOV runDirect finished")
        LOG.info("****************************************************************")
    except pyrttov.RttovError as e:
        LOG.error("Error running RTTOV direct model: {!s}".format(e))
        raise Exception("Error running RTTOV direct model: {!s}".format(e))

    allSimulatedBT = mwSatRttov.Bt[:, :]
    #all parameters from StoreRad2 = True 
    
    emis_param_dict = {'Rad2Down':        mwSatRttov.Rad2Down,
                       'Rad2Up':          mwSatRttov.Rad2Up,
                       'Rad2DnClear':     mwSatRttov.Rad2DnClear,
                       'Rad2UpClear':     mwSatRttov.Rad2UpClear,
                       'Rad2Surf':        mwSatRttov.Rad2Surf,
                       'Rad2ReflDnClear': mwSatRttov.Rad2ReflDnClear}

    return allSimulatedBT, chsStringList, emis_param_dict

########################################################
# Interface to call RTTOV-SCATT and compute the BT
########################################################
def calc_bt_rscatt(nwp_rttov_lats,
                   nwp_rttov_lons,
                   nwp_rttov_fg_sic,
                   nwp_rttov_skTemp,
                   nwp_rttov_2mTemp,
                   nwp_rttov_u10,
                   nwp_rttov_v10,
                   nwp_rttov_surfPre,
                   nwp_rttov_T,
                   nwp_rttov_Q,
                   nwp_rttov_CLWC,
                   nwp_rttov_CIWC,
                   nwp_rttov_CRWC,
                   nwp_rttov_CSWC,
                   nwp_rttov_CC,
                   modelFullLevhPa,
                   modelHalfLevhPa,
                   strSatId,
                   swathDate,
                   satZenAngle,
                   rttovCoefDir,
                   rttovHydDir,
                   rttovAtlasEmisDir,
                   emisToRttov,
                   satAzimuthAngle=None):

    # get sensor dependent rtCoef&channels
    rtCoefFile, rtHydFile, chsToKeepList, strChsList, emisIceChsList, constSatZenAngle = cfg.get_mw_sensor_cfg(strSatId)

    #########################
    # Subset channels
    #########################
    chsToProcessList = []
    chsStringList    = []

    for nCh in range(0, len(chsToKeepList)):

        if (chsToKeepList[nCh]):
            chsToProcessList.append(nCh+1)
            chsStringList.append(strChsList[nCh])

    #To be consistent with the calculation of surface emissvity
    if (sicThr > 0):
        indexIce   = (nwp_rttov_fg_sic > sicThr)
        indexOcean = (nwp_rttov_fg_sic <= sicThr)
    else:
        indexIce   = (nwp_rttov_fg_sic > 0)
        indexOcean = (nwp_rttov_fg_sic == 0)

    ################################################################
    # Now prepare and run RTTOV-SCATT
    ################################################################

    # ------------------------------------------------------------------------
    # Set up the profile data NWP
    # ------------------------------------------------------------------------
    
    # Pressure levels expected in hPa ordered from top atm level (python index 0) to surface (python last index)

    # Declare an instance of Profiles
    nlevels   = modelFullLevhPa.shape[1]
    nprofiles = modelFullLevhPa.shape[0]

    myProfiles = pyrttov.ProfilesScatt(nprofiles, nlevels)

    allZeros = np.asarray([0]*nprofiles)

    # bottom level, Force to be identical to surface pressure otherwise RTTOV-SCATT complains
    modelHalfLevhPa[:,nlevels] = nwp_rttov_surfPre[:]/100.

    # Profiles: 2D array [nprofiles, nlevels]
    myProfiles.GasUnits = cfg.GasUnits
    myProfiles.P  = modelFullLevhPa
    myProfiles.Ph = modelHalfLevhPa

    myProfiles.T = nwp_rttov_T
    myProfiles.Q = nwp_rttov_Q

    myProfiles.HydroFrac = nwp_rttov_CC
    myProfiles.Clw       = nwp_rttov_CLWC
    myProfiles.Ciw       = nwp_rttov_CIWC
    myProfiles.Snow      = nwp_rttov_CSWC
    myProfiles.Rain      = nwp_rttov_CRWC

    # Equivalent code using the flexible hydro interface
    #myProfiles.HydroFrac1 = nwp_rttov_CC
    #myProfiles.Hydro4 = nwp_rttov_CLWC
    #myProfiles.Hydro5 = nwp_rttov_CIWC
    #myProfiles.Hydro2 = nwp_rttov_CSWC
    #myProfiles.Hydro1 = nwp_rttov_CRWC

    #######################################
    # Surface: 2D array [nprofiles, nVars] - Same as before
    #######################################

    # Populate SurfGeom data structure: Real [nprofiles][3] (latitude, longitude, elevation)
    myelevs = allZeros  # elevation is in km
    surfgeom = np.array([nwp_rttov_lats, nwp_rttov_lons, myelevs], dtype=np.float64)
    myProfiles.SurfGeom = surfgeom.transpose()

    # Populate DateTimes data structure: Integer [nprofiles][6] (Not  really necessary for MW)
    # The full date will be used to calculate the TOA solar irradiance for solar-affected simulations.
    # dts = [timestamp.year, timestamp.month, timestamp.day, timestamp.hour, timestamp.minute, timestamp.second]

    mydateTimeYYYY = np.array([swathDate.year]*nprofiles, dtype=np.int32)
    mydateTimeMM   = np.array([swathDate.month]*nprofiles, dtype=np.int32)
    mydateTimeDD   = np.array([swathDate.day]*nprofiles, dtype=np.int32)
    mydateTimeHH   = np.array([swathDate.hour]*nprofiles, dtype=np.int32)
    mydateTimeMM   = np.array([swathDate.minute]*nprofiles, dtype=np.int32)
    mydateTimeSS   = np.array([swathDate.second]*nprofiles, dtype=np.int32)

    datetimes = np.array([mydateTimeYYYY,
                          mydateTimeMM,
                          mydateTimeDD,
                          mydateTimeHH,
                          mydateTimeMM,
                          mydateTimeSS], dtype=np.int32)

    myProfiles.DateTimes = datetimes.transpose()

    # RTTOV-SCATT:
    # Populate Angles data structure: Real [nprofiles][2]: satzen (mandatory), satazi
    # satzen: Local satellite zenith angle (degrees) (mandatory)
    # satazi: Local satellite azimuth angle (0-360; measured clockwise, east=90deegre           (not mandatory)

    if len(satZenAngle) == 1:
        allSatZenithAngle = np.array([satZenAngle[0]]*nprofiles, dtype=np.float64)
    else:
        allSatZenithAngle = satZenAngle

    if satAzimuthAngle is not None:
        allSatAzimuthAngle = satAzimuthAngle   
    else:
        allSatAzimuthAngle = allZeros   

    angles = np.array([allSatZenithAngle, allSatAzimuthAngle], dtype=np.float64)
    myProfiles.Angles = angles.transpose()        
    
    # RTTOV-SCATT:
    # Populate SurfType data structure:  Integer [nprofiles] (surftype)
    # surftype: (0=land, 1=sea, 2=seaice)

    # set all ocean points
    rttovSurfType = np.asarray([1]*nprofiles, dtype=np.int32)

    # repalce ocean with ice
    rttovSurfType[indexIce] = 2

    myProfiles.SurfType = rttovSurfType

    # RTTOV-SCATT:
    # Populate S2m data structure: Real [nprofiles][5] (s2m%p, s2m%t, s2m%q, s2m%u, s2m%v)

    # if we have 2m dewpoint. Calculate: 2m specific humidity kg/kg
    # where: e = vapor pressure in mb; Td = dew point in deg C; p = surface pressure in mb (=hPa); q = specific humidity in kg/kg.
    #e_mb  = [6.112*np.exp( (17.67*(Td-273.15))/((Td-273.15) + 243.5) ) for Td in my2mDewTemp]
    #my2mQ = [0.001*(0.622 * e_mb[i])/(p - (0.378 * e_mb[i])) for i,p in enumerate(mySurfPres)]

    # for now simply get the value of the model level closest to surface
    my2mQ = nwp_rttov_Q[:,modelFullLevhPa.shape[1]-1]

    s2m = np.array([nwp_rttov_surfPre/100.,    # Surface pressure (hPa)
                    nwp_rttov_2mTemp,         # 2m temperature
                    my2mQ,                    # 2m Q
                    nwp_rttov_u10,            # 10m wind u component
                    nwp_rttov_v10], dtype=np.float64)

    myProfiles.S2m = s2m.transpose()

    # RTTOV-SCATT:
    # Populate skin data structure: Real [nprofiles][8] (skinTemp, salinity, foam_fraction, fastem(1:5))
    # - salinity: use constant value of 35 psu (used by FASTEM 4-6, TESSEM2)
    # - foam_fraction: used by FASTEM if mwSatRttov.Options.SupplyFoamFraction = True
    # - fastem(1:5): FASTEM for land/sea-ice surface types --> not necessary

    allSalinity = np.asarray([35]*nprofiles)

    skin = np.array([nwp_rttov_skTemp,
                     allSalinity,
                     allZeros,
                     allZeros,
                     allZeros,
                     allZeros,
                     allZeros,
                     allZeros], dtype=np.float64)

    myProfiles.Skin = skin.transpose()

    # ------------------------------------------------------------------------
    # Set up RttovScatt instance
    # ------------------------------------------------------------------------

    mwSatRttov = pyrttov.RttovScatt()

    # Config RTTOV options

    # By default, this flag is false, and the effective cloud fraction is
    # calculated internally in RTTOV-SCATT, otherwise you must supply your own
    mwSatRttov.Options.LuserCfrac = False
    # To change treatment of specularity (together with surfemisrefl[3,:,:])
    mwSatRttov.Options.DoLambertian = cfg.DoLambertian
    mwSatRttov.Options.LambertianFixedAngle = cfg.LambertianFixedAngle
    # emis over ocean
    mwSatRttov.Options.FastemVersion   = cfg.FastemVersion
    # Interpolation to RTTOV pressure level
    mwSatRttov.Options.InterpMode = cfg.InterpMode
    # parallel calculation
    mwSatRttov.Options.Nthreads = cfg.Nthreads
    # store radiances
    mwSatRttov.Options.StoreRad = cfg.StoreRad
    # store emis term (for calculating emis ret)
    mwSatRttov.Options.StoreEmisTerms = cfg.StoreEmisTerms    
    # if IcePolarisation > 0, enables the approximation for polarised scattering in RTTOV-SCATT simulations
    mwSatRttov.Options.IcePolarisation = cfg.IcePolarisation
    # verbose on(off)
    mwSatRttov.Options.Verbose        = cfg.Verbose
    mwSatRttov.Options.VerboseWrapper = cfg.VerboseWrapper


    # Load the instruments: RTTOV-SCATT loads all the channels
    try:
        mwSatRttov.FileCoef       = '{}/{}'.format(rttovCoefDir, rtCoefFile)
        mwSatRttov.FileHydrotable = '{}/{}'.format(rttovHydDir, rtHydFile)
        mwSatRttov.loadInst()
    except pyrttov.RttovError as e:
        LOG.error("Error loading instrument(s): {!s}".format(e))
        raise Exception("Error loading instrument(s): {!s}".format(e))

    # sensor dependednt
    nchans    = mwSatRttov.Nchannels
    #satHeight = pyrttov.rtwrap.rttov_get_coef_val_r0(mwSatRttov.InstId,'FC_SAT_HEIGHT')[1]
    ff_cwn    = mwSatRttov.WaveNumbers
    freq_ghz  = ff_cwn * speedl * 1.0E-09

    LOG.info("MW sensor: {}". format(strSatId))
    LOG.info("Rttov coef file: {}".format(mwSatRttov.FileCoef))
    LOG.info("Num of channels (RTTOV-SCATT loads all sensor s channels): {}".format(nchans))
    LOG.info("- Freq (GHz): {}".format(freq_ghz))
    LOG.info("Gas units: {}".format(cfg.GasUnits))

    # Associate the profiles with each Rttov instance
    mwSatRttov.Profiles = myProfiles
    
    # Assign surface emissvity
    #surfemis = np.zeros((nprofiles,nchans), dtype=np.float64)
    surfemis = np.zeros((nprofiles,len(chsToProcessList)), dtype=np.float64)
    mwSatRttov.SurfEmis = surfemis

    surfemis[:,:] = emisToRttov[:,:]

    # Check rttov options
    if (cfg.logRTTOVOptions):
        check_rttov_options(mwSatRttov)

    try:
        LOG.info("****************************************************************")        
        LOG.info("RTTOV-SCATT runDirect started")
        LOG.info("RTTOV-SCATT hydrotable: {}".format(mwSatRttov.FileHydrotable))        
        LOG.info("Tot number of profiles: {}".format(nprofiles))
        LOG.info("Tot number of threads: {}".format(cfg.Nthreads))
        LOG.info("Do Lambertian (not active for RTTOV-SCATT): {}".format(cfg.DoLambertian))        
        mwSatRttov.runDirect(chsToProcessList)
        LOG.info("RTTOV-SCATT runDirect finished")
        LOG.info("****************************************************************")
    except pyrttov.RttovError as e:
        LOG.error("Error running RTTOV-SCATT: {!s}".format(e))
        raise Exception("Error running RTTOV-SCATT: {!s}".format(e))

    allSimulatedBT = mwSatRttov.Bt[:,:]
    # all-sky clear column (basically RTTOV T+Q)
    #allSimulatedBT_clr = mwSatRttov.BtClear[:, :]
    
    #all parameters from StoreEmisTerms = True 
    
    emis_param_dict = {'EmisTermsDownCld': mwSatRttov.EmisTermsDownCld,
                       'EmisTermsUpCld':   mwSatRttov.EmisTermsUpCld,
                       'EmisTermsTauCld':  mwSatRttov.EmisTermsTauCld,
                       'EmisTermsDownClr': mwSatRttov.EmisTermsDownClr,
                       'EmisTermsUpClr':   mwSatRttov.EmisTermsUpClr,
                       'EmisTermsTauClr':  mwSatRttov.EmisTermsTauClr,
                       'EmisTermsBsfc':    mwSatRttov.EmisTermsBsfc,
                       'EmisTermsCfrac':   mwSatRttov.EmisTermsCfrac}
    
    
    return allSimulatedBT, chsStringList, emis_param_dict

#####################################################################
# Interface to call RTTOV and calculate BT with 'no atm'
#####################################################################
def calc_bt_no_atm(nwp_rttov_lats,
                   nwp_rttov_lons,
                   nwp_rttov_fg_sic,
                   nwp_rttov_skTemp,
                   nwp_rttov_2mTemp,
                   nwp_rttov_u10,
                   nwp_rttov_v10,
                   nwp_rttov_surfPre,
                   nwp_rttov_T,
                   nwp_rttov_Q,
                   modelFullLevhPa,
                   strSatId,
                   swathDate,
                   satZenAngle,
                   rttovCoefDir,
                   rttovAtlasEmisDir,
                   emisToRttov,
                   satAzimuthAngle=None):

    # get sensor dependent rtCoef&channels
    rtCoefFile, rtHydFile, chsToKeepList, strChsList, emisIceChsList, constSatZenAngle = cfg.get_mw_sensor_cfg(strSatId)

    #########################
    # Subset channels
    #########################
    chsToProcessList = []
    chsStringList    = []

    for nCh in range(0, len(chsToKeepList)):

        if (chsToKeepList[nCh]):

            chsToProcessList.append(nCh+1)
            chsStringList.append(strChsList[nCh])


    #To be consistent with the calculation of surface emissvity
    if (sicThr > 0):
        indexIce   = (nwp_rttov_fg_sic > sicThr)
        indexOcean = (nwp_rttov_fg_sic <= sicThr)
    else:
        indexIce   = (nwp_rttov_fg_sic > 0)
        indexOcean = (nwp_rttov_fg_sic == 0)


    ##########################################
    # Do again RT for no atm (only T profile)
    ##########################################
    # ------------------------------------------------------------------------
    # Set up the profile data NWP
    # ------------------------------------------------------------------------

    # Pressure levels expected in hPa ordered from top atm level (python index 0) to surface (python last index)

    # Declare an instance of Profiles
    nlevels   = modelFullLevhPa.shape[1]
    nprofiles = modelFullLevhPa.shape[0]
    myProfiles = pyrttov.Profiles(nprofiles, nlevels)

    allZeros = np.asarray([0]*nprofiles)
    allOnes  = np.asarray([1]*nprofiles)

    nwp_rttov_zero_Q = np.zeros((nprofiles, nlevels))
    nwp_rttov_zero_Q[:,:] = 0.1000E-10

    # In case we want to try to force also temp profile to zero
    #nwp_rttov_zero_T = np.zeros((nprofiles, nlevels))
    #nwp_rttov_zero_T[:,:] = 0.1000E-10

    # Profiles: 2D array [nprofiles, nlevels]
    myProfiles.GasUnits = cfg.GasUnits
    myProfiles.P = modelFullLevhPa
    myProfiles.T = nwp_rttov_T
    
    ###################
    # Zero Q profile !
    ###################
    myProfiles.Q = nwp_rttov_zero_Q

    #######################################
    # Surface: 2D array [nprofiles, nVars] 
    #######################################

    # Populate SurfGeom data structure: Real [nprofiles][3] (latitude, longitude, elevation)
    myelevs = allZeros  # elevation is in km
    surfgeom = np.array([nwp_rttov_lats, nwp_rttov_lons, myelevs], dtype=np.float64)
    myProfiles.SurfGeom = surfgeom.transpose()

    # Populate DateTimes data structure: Integer [nprofiles][6] (Not  really necessary for MW)
    # The full date will be used to calculate the TOA solar irradiance for solar-affected simulations.
    # dts = [timestamp.year, timestamp.month, timestamp.day, timestamp.hour, timestamp.minute, timestamp.second]

    mydateTimeYYYY = np.array([swathDate.year]*nprofiles, dtype=np.int32)
    mydateTimeMM   = np.array([swathDate.month]*nprofiles, dtype=np.int32)
    mydateTimeDD   = np.array([swathDate.day]*nprofiles, dtype=np.int32)
    mydateTimeHH   = np.array([swathDate.hour]*nprofiles, dtype=np.int32)
    mydateTimeMM   = np.array([swathDate.minute]*nprofiles, dtype=np.int32)
    mydateTimeSS   = np.array([swathDate.second]*nprofiles, dtype=np.int32)

    datetimes = np.array([mydateTimeYYYY,
                          mydateTimeMM,
                          mydateTimeDD,
                          mydateTimeHH,
                          mydateTimeMM,
                          mydateTimeSS], dtype=np.int32)

    myProfiles.DateTimes = datetimes.transpose()

    # Populate Angles data structure: Real [nprofiles][4]: satzen (mandatory), satazi, sunzen, sunazi
    # satzen: Local satellite zenith angle (degrees) (mandatory)
    # satazi: Local satellite azimuth angle (0-360; measured clockwise, east=90deegre           (not mandatory)
    # sunzen: Local solar zenith angle (degrees), solar radiation only included up to 85deegre  (not mandatory)
    # sunazi: Local solar azimuth angle (0-360; measured clockwise, east=90deegre               (not mandatory)
    # NOTE: To turn on solar radiation for visible and near-IR channel,
    #       solar zenith angle, and satellite and solar azimuth angles must also be specified

    if len(satZenAngle) == 1:
        allSatZenithAngle = np.array([satZenAngle[0]]*nprofiles, dtype=np.float64)
    else:
        allSatZenithAngle = satZenAngle
        
    if satAzimuthAngle is not None:
        allSatAzimuthAngle = satAzimuthAngle   
    else:
        allSatAzimuthAngle = allZeros   

    angles = np.array([allSatZenithAngle, allSatAzimuthAngle, allZeros, allZeros], dtype=np.float64)
    myProfiles.Angles = angles.transpose()

    # Populate SurfType data structure:  Integer [nprofiles][2] (surftype, watertype)
    # surftype: (0=land, 1=sea, 2=seaice) watertype 0=fresh water, 1=ocean water.

    # set all ocean points
    rttovSurfType = np.asarray([1]*nprofiles)

    # Associciate again surface type to the profile
    rttovSurfType[indexIce] = 2

    surftype = np.array([rttovSurfType,
                         allOnes], dtype=np.int32)

    myProfiles.SurfType = surftype.transpose()

    # Populate S2m data structure: Real [nprofiles][6] (s2m%p, s2m%t, s2m%q, s2m%u, s2m%v, s2m%wfetc)

    windfetch = np.asarray([100000]*nprofiles)  # Wind fetch, value 100000m for open ocean – only used by sea surface solar BRDF model

    # if we have 2m dewpoint. Calculate: 2m specific humidity kg/kg
    # where: e = vapor pressure in mb; Td = dew point in deg C; p = surface pressure in mb (=hPa); q = specific humidity in kg/kg.
    #e_mb  = [6.112*np.exp( (17.67*(Td-273.15))/((Td-273.15) + 243.5) ) for Td in my2mDewTemp]
    #my2mQ = [0.001*(0.622 * e_mb[i])/(p - (0.378 * e_mb[i])) for i,p in enumerate(mySurfPres)]

    # for now simply get the value of the model level closest to surface
    my2mQ = nwp_rttov_zero_Q[:,modelFullLevhPa.shape[1]-1]
    
    ##############
    # Zero Winds
    ##############
    s2m = np.array([nwp_rttov_surfPre/100.,    # Surface pressure (hPa)
                    nwp_rttov_2mTemp,          # 2m temperature
                    my2mQ,                     # 2m Q
                    allZeros,                  # 10m wind u component
                    allZeros,                  # 10m wind v component
                    windfetch], dtype=np.float64)

    myProfiles.S2m = s2m.transpose()

    # Populate skin data structure: Real [nprofiles][9] (skinTemp, salinity, snow_fraction, foam_fraction, fastem(1:5))
    # - salinity: use constant value of 35 psu (used by FASTEM 4-6, TESSEM2)
    # - snow_fraction: Surface snow coverage fraction (0-1) (used by IR emissivity atlas)
    # - foam_fraction: used by FASTEM if mwSatRttov.Options.SupplyFoamFraction = True
    # - fastem(1:5): FASTEM for land/sea-ice surface types --> not necessary

    allSalinity = np.asarray([35]*nprofiles)

    skin = np.array([nwp_rttov_skTemp,
                     allSalinity,
                     allZeros,
                     allZeros,
                     allZeros,
                     allZeros,
                     allZeros,
                     allZeros,
                     allZeros], dtype=np.float64)

    myProfiles.Skin = skin.transpose()

    # ------------------------------------------------------------------------
    # Set up Rttov instances for an instrument
    # ------------------------------------------------------------------------
    mwSatRttov = pyrttov.Rttov()

    # Config RTTOV options

    # Make sure we do not include contribution of cloud liquid water in RT calculations
    mwSatRttov.Options.CLWData = False
    # To change treatment of specularity (together with surfemisrefl[3,:,:])
    mwSatRttov.Options.DoLambertian = cfg.DoLambertian
    mwSatRttov.Options.LambertianFixedAngle = cfg.LambertianFixedAngle
    # emis over ocean
    mwSatRttov.Options.FastemVersion   = cfg.FastemVersion
    # Interpolation to RTTOV pressure level
    mwSatRttov.Options.InterpMode = cfg.InterpMode
    # parallel calculation
    mwSatRttov.Options.Nthreads = cfg.Nthreads
    # store radiances
    mwSatRttov.Options.StoreRad = cfg.StoreRad
    # verbose on(off)
    mwSatRttov.Options.Verbose        = cfg.Verbose
    mwSatRttov.Options.VerboseWrapper = cfg.VerboseWrapper
    # this must be false otherwise RTTOV complains because we do not have Q profile
    mwSatRttov.Options.DoCheckinput   = cfg.DoCheckinput

    # Load the instruments
    try:
        mwSatRttov.FileCoef = '{}/{}'.format(rttovCoefDir,rtCoefFile)
        mwSatRttov.loadInst(chsToProcessList)
    except pyrttov.RttovError as e:
        LOG.error("Error loading instrument(s): {!s}".format(e))
        raise Exception("Error loading instrument(s): {!s}".format(e))

    nchans = mwSatRttov.Nchannels
    #satHeight = pyrttov.rtwrap.rttov_get_coef_val_r0(mwSatRttov.InstId,'FC_SAT_HEIGHT')[1]
    ff_cwn    = mwSatRttov.WaveNumbers
    freq_ghz  = ff_cwn * speedl * 1.0E-09

    LOG.info("MW sensor: {}". format(strSatId))
    LOG.info("Rttov coef file: {}".format(mwSatRttov.FileCoef))
    LOG.info("Number of channels selected for RTTOV: {}".format(nchans))
    LOG.info("- List of selected channels: {}".format(chsToProcessList))
    LOG.info("- Freq (GHz): {}".format(freq_ghz))
    
    # Associate the profiles with each Rttov instance
    mwSatRttov.Profiles = myProfiles

    if rttovVer == 'v13.2':
        surfemisrefl = np.zeros((5,nprofiles,nchans), dtype=np.float64)
    else:
        surfemisrefl = np.zeros((4,nprofiles,nchans), dtype=np.float64)
        
    mwSatRttov.SurfEmisRefl = surfemisrefl

    # Surface emissivity/reflectance arrays must be initialised *before every call to RTTOV*
    # Negative values will cause RTTOV to supply emissivity/BRDF values
    surfemisrefl[:,:,:] = -1.
    surfemisrefl[0,:,:] = emisToRttov
    # where ocean, recalulate rttov emis with winds zero
    surfemisrefl[0,indexOcean,:] = -1
    
    # in case we activate the Lambertian option
    # Default: surface is treated as fully specular reflector --> specularity = 1 (DoLambertian == false)
    # fully Lambertian --> specularity = 0
    if (cfg.DoLambertian):
        surfemisrefl[3,indexIce,:]   = 0
        surfemisrefl[3,indexOcean,:] = 1

    # Check rttov options
    if (cfg.logRTTOVOptions):
        check_rttov_options(mwSatRttov)

    try:
        LOG.info("****************************************************************")        
        LOG.info("RTTOV runDirect started (compute BT for 'no atm')")
        LOG.info("Tot number of profiles: {}".format(nprofiles))
        LOG.info("Tot number of threads: {}".format(cfg.Nthreads))
        LOG.info("Do Lambertian: {}".format(cfg.DoLambertian))
        mwSatRttov.runDirect()
        LOG.info("RTTOV runDirect finished")
        LOG.info("****************************************************************")
    except pyrttov.RttovError as e:
        LOG.error("Error running RTTOV direct model: {!s}".format(e))
        raise Exception("Error running RTTOV direct model: {!s}".format(e))

    allSimulatedBT_no_atm = mwSatRttov.Bt[:, :]

    return allSimulatedBT_no_atm

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
pi = const.pi

# if (sicThr > 0) identify ocean/seaice surface such as:

# where fgSic > sicThr:
#    Emis = (1-sic)*EmisOcean + sic*EmisIce
# otherwise
#    Emis = EmisOcean
sicThr = cfg.sicFgThr

#rttov version
rttovVer = cfg.rttovVer

# jsut for testing, possibility to inflate ice emis
inflateIceEmis = False
# % 
inflatePer = 10

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
# interface to calculate surface emissvity
########################################################
def calc_emis(nwp_rttov_lats,
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
              useEmisAtlas,
              satAzimuthAngle=None):

    # get sensor dependent rtCoef&channels
    rtCoefFile, rtHydFile, chsToKeepList, strChsList, emisIceChsList, constSatZenAngle = cfg.get_mw_sensor_cfg(strSatId)

    #########################
    # Subset channels
    #########################
    chsToProcessList = []
    chsStringList    = []
    chsEmisIceList   = []

    for nCh in range(0, len(chsToKeepList)):

        if (chsToKeepList[nCh]):

            chsToProcessList.append(nCh+1)
            chsStringList.append(strChsList[nCh])

            checkEmisIce = emisIceChsList[nCh]

            if (checkEmisIce < 0):
                LOG.warning("Constant value of ice emis not found for channel {}; use constant value of 0.90".format(nCh+1))
                chsEmisIceList.append(0.90)
            else:
                chsEmisIceList.append(checkEmisIce)


    # Estimate emisvity
    
    # where fgSic > sicThr:
    #    Emis = (1-sic)*EmisOcean + sic*EmisIce
    # otherwise
    #    Emis = EmisOcean
    
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

    # we do not have a sat zen angle at every obs location, so we use a const value
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

    # first RTTOV run: all ocean points to get an estimate of ocean emisivvity from RTTOV
    rttovSurfType = np.asarray([1]*nprofiles)

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

    if rttovVer == 'v13.2':
        surfemisrefl = np.zeros((5,nprofiles,nchans), dtype=np.float64)
    else:
        surfemisrefl = np.zeros((4,nprofiles,nchans), dtype=np.float64)

    mwSatRttov.SurfEmisRefl = surfemisrefl

    # Surface emissivity/reflectance arrays must be initialised *before every call to RTTOV*
    # Negative values will cause RTTOV to supply emissivity/BRDF values
    surfemisrefl[:,:,:] = -1.

    # in case we activate the Lambertian option
    # Default: surface is treated as fully specular reflector --> specularity = 1 (DoLambertian == false)
    # fully Lambertian --> specularity = 0
    if (cfg.DoLambertian):
        surfemisrefl[3,indexOcean,:] = 1
        
    # Check rttov options
    if (cfg.logRTTOVOptions):
        check_rttov_options(mwSatRttov)

    try:
        LOG.info("****************************************************************")        
        LOG.info("RTTOV runDirect started (get ocean surface emissvity from RTTOV)")
        LOG.info("MW ocean emis model (0=TESSEM2 else FASTEM version): {}".format(cfg.FastemVersion))
        LOG.info("Tot number of profiles: {}".format(nprofiles))
        LOG.info("Tot number of threads: {}".format(cfg.Nthreads))
        mwSatRttov.runDirect()
        LOG.info("RTTOV runDirect finished")
        LOG.info("****************************************************************")
    except pyrttov.RttovError as e:
        LOG.error("Error running RTTOV direct model: {!s}".format(e))
        raise Exception("Error running RTTOV direct model: {!s}".format(e))

    # get emis computed by RTTOV
    emisOceanRttov = surfemisrefl[0, :, :]

    if (sicThr > 0):
        LOG.info("Estimate surface emis as (1-sic)*EmisOcean + sic*EmisIce where sic > {}, otherwise emis = EmisOcean".format(sicThr))
    else:
        LOG.info("Estimate surface emis as (1-sic)*EmisOcean + sic*EmisIce")

    # Associciate again surface type to the profile
    rttovSurfType[indexIce] = 2

    surftype = np.array([rttovSurfType,
                         allOnes], dtype=np.int32)

    myProfiles.SurfType = surftype.transpose()
    mwSatRttov.Profiles = myProfiles

    ################################
    # Emissivity over ice
    ################################
    if (useEmisAtlas):

        LOG.info("Emis over ice: getting values from Atlas")
        # TELSEM2 atlas does not require an Rttov object to initialise
        mwAtlas = pyrttov.Atlas()
        mwAtlas.AtlasPath = rttovAtlasEmisDir
        # loadMwEmisAtlas(month, inst=None, atlas_id=-1)
        # - the inst argument can be a loaded Rttov or RttovScatt instance: this is required for the CNRM atlas, but is ignored by TELSEM2
        # - atlas_id = 1 (default); CNRM MW atlas: atlas_id = 2
        mwAtlas.loadMwEmisAtlas(int(swathDate.month))
        mwAtlas.IncSea  = False
        mwAtlas.IncLand = False
        mwAtlas.IncSeaIce = True

        # Call emissivity atlas
        try:
            emisIceToRttov = mwAtlas.getEmisBrdf(mwSatRttov)
        except pyrttov.RttovError as e:
            LOG.error("Error calling atlas: {!s}\n".format(e))
            raise Exception("Error calling atlas: {!s}\n".format(e))

        constEmisIce = np.asarray(chsEmisIceList[:])

        # check only 1 freq (if atlas does not have values returns -1 for all freq)
        indexAtlasNotAv = (emisIceToRttov[:,0] == -1) & (rttovSurfType == 2)

        # Not av emis from atlas over ice --> use estimates from fastem
        if ((indexAtlasNotAv == True).sum() != 0):
            emisIceToRttov[indexAtlasNotAv,:] = constEmisIce[:]
            LOG.info("Ice points where atlas value is not availabe {}: in this case use constant ice emis".format((indexAtlasNotAv == True).sum()))
            LOG.info("Constant ice emis: {}".format(constEmisIce))
        
        if inflateIceEmis:
            
            inflateEmis = emisIceToRttov.copy()
            indexLT = np.logical_and(inflateEmis > 0, inflateEmis < 0.9)
            inflateEmis[indexLT] = inflateEmis[indexLT] + (inflateEmis[indexLT]*inflatePer)/100. 
            emisIceToRttov[indexLT] = inflateEmis[indexLT]
    else:

        LOG.info("Emis over ice: using constant value")

        constEmisIce = np.asarray(chsEmisIceList[:])
        LOG.info("Constant ice emis: {}".format(constEmisIce))

        emisIceToRttov = np.zeros((nprofiles,nchans), dtype=np.float64)
        emisIceToRttov[:,:] = constEmisIce[:]

    # Finally: Emis = (1-sic)*EmisOcean + sic*EmisIce
    emisToRttov = np.zeros((nprofiles,nchans), dtype=np.float64)

    for nCh in range(0,nchans):
        emisToRttov[indexIce,nCh] = (1-nwp_rttov_fg_sic[indexIce])*emisOceanRttov[indexIce,nCh] + nwp_rttov_fg_sic[indexIce]*emisIceToRttov[indexIce,nCh]
        emisToRttov[indexOcean,nCh] = emisOceanRttov[indexOcean,nCh]

    return  emisToRttov

##########################
# FASTEM MW LAND/ICE EMIS 
# FB comment: Fastem land/ice emis not recommended - kept just in case
##########################
def fastem_mw_seaice_emis (satZenithAngle, satHeight, ff_cwn, pol_id):
    """
    This function emulates the RTTOV mw fastem land/ice surface emissvity model (rttov_calcemis_mw.F90)
    
    satZenithAngle: sat zenith angle (single value)
    satHeight: sat height in km
    ff_cwn: list of freq in Ghz
    pol_id: list of pol (numeric values must match those provided in the RTTOV coef file, e.g. 'rtcoef_gcom-w_1_amsr2.dat')
    """
    # fastemCoef can be found in the RTTOV user guide
    # CompactIFastem = [2.0, 1700000.0, 49000000.0, 0., 0.]
    # MYIFastem     = [1.5, 85000.0,   4700000.0,  0., 0.]
    fastemCoef = [2.0, 1700000.0, 49000000.0, 0., 0.]
    
    # Earth Radius km
    earthradius = 6371.0

    pol_v = np.array( [ [0.5, 0.0, 0.0],
                        [0.0, 0.0, 1.0],
                        [0.0, 1.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],])

    pol_h = np.array( [ [0.5, 0.0, 0.0],
                        [0.0, 1.0, 0.0],
                        [0.0, 0.0, 1.0],
                        [0.0, 0.0, 0.0],
                        [1.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],
                        [0.0, 0.0, 0.0],])

    pol_s3 = np.array( [ [0.0, 0.0],
                         [0.0, 0.0],
                         [0.0, 0.0],
                         [0.0, 0.0],
                         [0.0, 0.0],
                         [1.0, 0.0],
                         [0.0, 1.0],])


    angleRad  = np.deg2rad(satZenithAngle)
    sinzen_sq = np.sin(angleRad)*np.sin(angleRad)
    coszen    = np.cos(angleRad)
    coszen_sq = np.cos(angleRad)*np.cos(angleRad)


    ratoe =  (earthradius + satHeight) / earthradius
    sinzen = np.sin(angleRad)
    sinview = sinzen/ratoe

    sinview_sq = sinview * sinview
    cosview_sq = 1.0 - sinview_sq

    nchannels = len(ff_cwn)
    emissivity = np.zeros(nchannels)
    reflectivity = np.zeros(nchannels)

    perm_static         = fastemCoef[0]
    perm_infinite       = fastemCoef[1]
    freqr               = fastemCoef[2]
    small_rough         = fastemCoef[3]
    large_rough         = fastemCoef[4]

    emissstokes   = np.zeros((nchannels,4))
    reflectstokes = np.zeros((nchannels,4))

    for nCh in range(0,nchannels):

        freq_ghz = ff_cwn[nCh] * speedl * 1.0E-09
        polId = int(pol_id[nCh])

        #!Simple Debye + Fresnel model gives reflectivities
        fen = freq_ghz / freqr
        fen_sq              = fen * fen
        den1 = 1.0 + fen_sq
        perm_Real           = (perm_static + perm_infinite * fen_sq) / den1
        perm_imag           = fen * (perm_static - perm_infinite) / den1
        permittivity        = np.complex(perm_Real, perm_imag)
        perm1               = np.sqrt(permittivity - sinzen_sq)
        perm2               = permittivity * coszen
        rhth = (coszen - perm1) / (coszen + perm1)
        rvth = (perm2 - perm1) / (perm2 + perm1)
        fresnel_v_Real      = np.real(rvth)
        fresnel_v_imag      = np.imag(rvth)
        fresnel_v           = fresnel_v_Real * fresnel_v_Real + fresnel_v_imag * fresnel_v_imag
        fresnel_h_Real      = np.real(rhth)
        fresnel_h_imag      = np.imag(rhth)
        fresnel_h           = fresnel_h_Real * fresnel_h_Real + fresnel_h_imag * fresnel_h_imag

        #!Small scale roughness correction
        delta               = 4.0 * pi * ff_cwn[nCh] * 0.1 * small_rough
        delta2              = delta * delta
        small_rough_cor     = np.exp( - delta2 * coszen_sq)

        #Large scale roughness correction
        qdepol              = 0.35 - 0.35 * np.exp( - 0.60 * freq_ghz * large_rough * large_rough)
        emissfactor_v       = 1.0 - fresnel_v * small_rough_cor
        emissfactor_h       = 1.0 - fresnel_h * small_rough_cor
        emissfactor         = emissfactor_h - emissfactor_v
        emissstokes[nCh, 0]   = emissfactor_v + qdepol * emissfactor
        emissstokes[nCh, 1]   = emissfactor_h - qdepol * emissfactor
        emissstokes[nCh, 2]   = 0.0
        emissstokes[nCh, 3]   = 0.0
        reflectstokes[nCh, :] = 1.0 - emissstokes[nCh, :]

        # Now calc channel emissivity after mixing v and h pol
        emissfactor_v = pol_v[polId,0] + pol_v[polId,1] * sinview_sq + pol_v[polId,2] * cosview_sq
        emissfactor_h = pol_h[polId,0] + pol_h[polId,1] * sinview_sq + pol_h[polId,2] * cosview_sq

        emissivity[nCh] = emissstokes[nCh,0] * emissfactor_v + emissstokes[nCh,1] * emissfactor_h + emissstokes[nCh,2] * pol_s3[polId,0] + emissstokes[nCh,3] * pol_s3[polId,1]

        reflectivity[nCh] = reflectstokes[nCh,0] * emissfactor_v + reflectstokes[nCh,1] * emissfactor_h + reflectstokes[nCh,2] * pol_s3[polId,0] + reflectstokes[nCh,3] * pol_s3[polId,1]

    return emissivity
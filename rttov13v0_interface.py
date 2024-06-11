
import mw_rttov_cfg as cfg
from calc_emis_for_rttov import calc_emis
from rttov13v0_calc_bt import calc_bt_rscatt, calc_bt, calc_bt_no_atm

import logging

LOG = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='[%(levelname)s: %(asctime)s: %(name)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
#rttov version
rttovVer = cfg.rttovVer

def run_rttov(nwp_rttov_lats,
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
              hPaFullLevels,
              hPaHalfLevels,
              strSatId,
              satZenAngleIn,
              centralSwathTime,
              rttovCoefDir,
              rttovHydDir,
              rttovAtlasEmisDir,
              rttovScatt,
              useEmisAtlas,
              calbBTnoAtm,
              satAzimuthAngleIn=None):
    '''
    Main function to call and run rttov. Implementaion is based on RTTOV V13.0
    '''
    
    LOG.info('RTTOV VERSION from mw_rttov_cfg.py is {0}'.format(rttovVer))
    
    # we do not have a sat zen angle at every obs location, so we use a const value
    if len(satZenAngleIn) == 1 and satZenAngleIn[0] == 0:
        # get sensor dependent rtCoef&channels
        dummy1, dummy2, dummy3, dummy4, dummy5, constSatZenAngle = cfg.get_mw_sensor_cfg(strSatId)
        LOG.info('Constant sat zenith angle of {0} is used at every obs location'.format(constSatZenAngle))
        satZenAngle = [constSatZenAngle]
    else:
        LOG.info('--> Sat zenith angle varies at every obs location, mean values is {}'.format(satZenAngleIn.mean()))
        satZenAngle = satZenAngleIn
        
    if satAzimuthAngleIn is not None:
        LOG.info('--> Sat Azimuth angle is provided and it varies at every obs location, mean values is {}'.format(satAzimuthAngleIn.mean()))
        
        
    LOG.info("Calculating surface emissvity for RTTOV")
    
    surf_emis = calc_emis(nwp_rttov_lats,
                          nwp_rttov_lons,
                          nwp_rttov_fg_sic,
                          nwp_rttov_skTemp,
                          nwp_rttov_2mTemp,
                          nwp_rttov_u10,
                          nwp_rttov_v10,
                          nwp_rttov_surfPre,
                          nwp_rttov_T,
                          nwp_rttov_Q,
                          hPaFullLevels,
                          strSatId,
                          centralSwathTime,
                          satZenAngle,
                          rttovCoefDir,
                          rttovAtlasEmisDir,
                          useEmisAtlas,
                          satAzimuthAngle=satAzimuthAngleIn)
    
    if (rttovScatt):

        LOG.info("Simulating BT using RTTOV-SCATT")

        Tb_nwp, strChsList, emis_param_dict = calc_bt_rscatt(nwp_rttov_lats,
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
                                                             hPaFullLevels,
                                                             hPaHalfLevels,
                                                             strSatId,
                                                             centralSwathTime,
                                                             satZenAngle,
                                                             rttovCoefDir,
                                                             rttovHydDir,
                                                             rttovAtlasEmisDir,
                                                             surf_emis,
                                                             satAzimuthAngle=satAzimuthAngleIn)
        
    else:

        LOG.info("Simulating BT using RTTOV (T+Q+CLWC)")

        Tb_nwp, strChsList, emis_param_dict = calc_bt(nwp_rttov_lats,
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
                                                      hPaFullLevels,
                                                      strSatId,
                                                      centralSwathTime,
                                                      satZenAngle,
                                                      rttovCoefDir,
                                                      rttovAtlasEmisDir,
                                                      surf_emis,
                                                      satAzimuthAngle=satAzimuthAngleIn)


    if (calbBTnoAtm):
        LOG.info("Simulating BT for reference atmospheric state ('no atm')")
        
        Tb_nwp_no_atm = calc_bt_no_atm(nwp_rttov_lats,
                                       nwp_rttov_lons,
                                       nwp_rttov_fg_sic,
                                       nwp_rttov_skTemp,
                                       nwp_rttov_2mTemp,
                                       nwp_rttov_u10,
                                       nwp_rttov_v10,
                                       nwp_rttov_surfPre,
                                       nwp_rttov_T,
                                       nwp_rttov_Q,
                                       hPaFullLevels,
                                       strSatId,
                                       centralSwathTime,
                                       satZenAngle,
                                       rttovCoefDir,
                                       rttovAtlasEmisDir,
                                       surf_emis,
                                       satAzimuthAngle=satAzimuthAngleIn)
        
    
    else:
        LOG.info("calbBTmoAtm is False --> We do not simulate BT for reference atmospheric state")
        Tb_nwp_no_atm = []
    
    return Tb_nwp, Tb_nwp_no_atm, surf_emis, strChsList, emis_param_dict 
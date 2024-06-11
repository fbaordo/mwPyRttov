#Author: Fabrizio Baordo, DMI: fab@dmi.dk

from mpl_toolkits.basemap import Basemap
import matplotlib.pyplot as plt
import logging
import numpy as np
from scipy.stats import skew
from os import path, makedirs, remove
from netCDF4 import Dataset
import datetime

#import warnings
#warnings.filterwarnings("ignore", category=RuntimeWarning) 

LOG = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO,
                    format='[%(levelname)s: %(asctime)s: %(name)s] %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

prop_cycle = plt.rcParams['axes.prop_cycle']
colorList = prop_cycle.by_key()['color']


# INFO: in the fucntions below: 
# - the variables lats, lons, obsCh are expected to be 1D array
# - the variables named as '***List' (e.g. simChList) are list containing 1D array for every simulation which must be compared


def diag_departures(lats, lons, obsCh, simChList, strCh, strSatId, expLabelList, strYYYYMMDDHHMM, plotSavePath=None, strPltId=None):

    nExp = len(simChList)
    
    nRow = 2
    nCol = nExp + 1
        
    index_nh = lats > 0
    index_sh = lats < 0

    cmap = plt.cm.seismic
    
    bounds = [-20, -15, -10, -5, -2, 0, 2, 5, 10, 15, 20]
        
    histMin = -50
    histMax = 50 
    histBinSize = 1
    
    myBins = np.arange(histMin, histMax+histBinSize, histBinSize)     

    width  = 6*nCol
    height = 12
    
    fig = plt.figure(figsize=(width,height),dpi=100)

    fig.suptitle('Freq (GHz) and pol: ' +  strCh + ' - ' + strSatId + ' ' + strYYYYMMDDHHMM, fontsize=14, fontweight='bold')

    #first row: plots for NH        

    # Maps
    for nE in range(0,nExp):

        dep = obsCh - simChList[nE]        
        strExp = expLabelList[nE]
        
        fig.add_subplot(nRow, nCol, nE+1)
        
        mylat = lats[index_nh]
        mylon = lons[index_nh]
        myDep = dep[index_nh]

        m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
        m.drawcoastlines()
        m.fillcontinents(color='gray')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))    
                
        x, y = m(mylon, mylat)
        cs = m.scatter(x,y,c=myDep,cmap=cmap, s=10, edgecolors='none', vmin=bounds[0], vmax=bounds[len(bounds)-1])
        plt.colorbar(cs, shrink = .8, orientation="horizontal",boundaries=bounds,ticks=bounds, spacing='proportional')

        dataMin = str('{:.2f}'.format(np.nanmin(myDep))) 
        dataMax  = str('{:.2f}'.format(np.nanmax(myDep)))   
        dataMean = str('{:.2f}'.format(np.nanmean(myDep)))         

        strTitleExp = 'min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
        plt.title('Dep (Obs-Sim) [K] - ' + strExp + '\n' + strTitleExp)
        
    # Histogram
    fig.add_subplot(nRow, nCol, nCol)
    
    for nE in range(0,nExp):

        dep = obsCh - simChList[nE]        
        strExp = expLabelList[nE]
        myDep = dep[index_nh]

        dataMean = str('{:.2f}'.format(np.nanmean(myDep))) 
        dataStd  = str('{:.2f}'.format(np.nanstd(myDep)))   
        dataSkew = str('{:.2f}'.format(skew(myDep[np.isfinite(myDep)])))         

        plt.hist(myDep, bins=myBins, histtype='step',color=colorList[nE],density=False, linewidth=1.2,label=expLabelList[nE])
        
        if (nE==0):
            strTitleExp = 'mean/std/skew: ' + dataMean + '/' + dataStd + '/' + dataSkew
        else:
            strTitle = 'mean/std/skew: ' + dataMean + '/' + dataStd + '/' + dataSkew
            
            strTitleExp = strTitleExp + '\n' + strTitle
            
    #plt.yscale('log')
    plt.ticklabel_format(style='sci',axis='y',scilimits=(3,3))
    plt.grid(True)
    plt.title(strTitleExp)  #, fontsize=8, fontweight='bold')
    #plt.xlabel('Dep (Obs-Sim) [K]')
    plt.ylabel('Counts')
    plt.legend(loc=2)
                       
    #second row: plots for SH
    # Maps
    for nE in range(0,nExp):

        dep = obsCh - simChList[nE]        
        strExp = expLabelList[nE]
        
        fig.add_subplot(nRow, nCol, nCol+nE+1)
        
        mylat = lats[index_sh]
        mylon = lons[index_sh]
        myDep = dep[index_sh]

        m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)                   
        m.drawcoastlines()
        m.fillcontinents(color='gray')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))    
                
        x, y = m(mylon, mylat)
        cs = m.scatter(x,y,c=myDep,cmap=cmap, s=10, edgecolors='none', vmin=bounds[0], vmax=bounds[len(bounds)-1])
        plt.colorbar(cs, shrink = .8, orientation="horizontal",boundaries=bounds,ticks=bounds, spacing='proportional')        

        dataMin = str('{:.2f}'.format(np.nanmin(myDep))) 
        dataMax  = str('{:.2f}'.format(np.nanmax(myDep)))   
        dataMean = str('{:.2f}'.format(np.nanmean(myDep)))         
        
        strTitleExp = 'min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
        plt.title('Dep (Obs-Sim) [K] - ' + strExp + '\n' + strTitleExp)
        
    # Histogram
    fig.add_subplot(nRow, nCol, nCol+nCol)
    
    for nE in range(0,nExp):

        dep = obsCh - simChList[nE]        
        strExp = expLabelList[nE]
        myDep = dep[index_sh]

        dataMean = str('{:.2f}'.format(np.nanmean(myDep))) 
        dataStd  = str('{:.2f}'.format(np.nanstd(myDep)))   
        dataSkew = str('{:.2f}'.format(skew(myDep[np.isfinite(myDep)])))         

        plt.hist(myDep, bins=myBins, histtype='step',color=colorList[nE],density=False, linewidth=1.2,label=expLabelList[nE])
        
        if (nE==0):
            strTitleExp = 'mean/std/skew: ' + dataMean + '/' + dataStd + '/' + dataSkew
        else:
            strTitle = 'mean/std/skew: ' + dataMean + '/' + dataStd + '/' + dataSkew
            
            strTitleExp = strTitleExp + '\n' + strTitle
            
    #plt.yscale('log')
    plt.ticklabel_format(style='sci',axis='y',scilimits=(3,3))
    plt.grid(True)
    plt.title(strTitleExp)  #, fontsize=8, fontweight='bold')
    plt.xlabel('Dep (Obs-Sim) [K]')
    plt.ylabel('Counts')
    #plt.legend(loc=0) 
    
    # save plots....
    if plotSavePath is not None:
        
        if strPltId is not None:
            plotName = 'diag_dep.'+ strSatId + '.' + strCh + '.' + strYYYYMMDDHHMM + '.' + strPltId +'.png' 
        else:
            plotName = 'diag_dep.'+ strSatId + '.' + strCh + '.' + strYYYYMMDDHHMM + '.png' 

        plt.savefig(path.join(plotSavePath,plotName))         
        plt.close()
        LOG.info('Diagnostic plot for departure saved as: {0} '.format(path.join(plotSavePath,plotName)))
            
    return None
    
def diag_atmCorr(lats,lons, obsCh, atmCorrList, strCh, strSatId, expLabelList, strYYYYMMDDHHMM, plotSavePath=None, strPltId=None):
    
    nExp = len(atmCorrList)
    
    nRow = 2
    nCol = nExp + 1
        
    index_nh = lats > 0
    index_sh = lats < 0

    cmap = plt.cm.jet
        
    width  = 6*nCol
    height = 12
    
    fig = plt.figure(figsize=(width,height),dpi=100)

    fig.suptitle('Freq (GHz) and pol: ' +  strCh + ' - ' + strSatId + ' ' + strYYYYMMDDHHMM, fontsize=14, fontweight='bold')

    #first row: plots for NH        

    # Maps of atm corr
    for nE in range(0,nExp):

        btCorr = atmCorrList[nE][index_nh]        
        strExp = expLabelList[nE]

        dataMean = str('{:.2f}'.format(np.nanmean(btCorr))) 
        dataMin  = str('{:.2f}'.format(np.nanmin(btCorr)))   
        dataMax  = str('{:.2f}'.format(np.nanmax(btCorr)))   

        if (nE == 0):
            minAtmCorr = dataMin
            maxAtmCorr = dataMax
        
        fig.add_subplot(nRow, nCol, nE+1+1)
        
        mylat = lats[index_nh]
        mylon = lons[index_nh]

        m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
        m.drawcoastlines()
        m.fillcontinents(color='gray')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))    
                
        x, y = m(mylon, mylat)
        cs = m.scatter(x,y,c=btCorr,cmap=cmap, s=10, edgecolors='none', vmin=minAtmCorr, vmax=maxAtmCorr)
        plt.colorbar(cs, shrink = .8, orientation="horizontal")
        
        strTitleExp = 'min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
        plt.title('Atm Corr [K] - ' + strExp + '\n' + strTitleExp)
        
    # maps of obs
    fig.add_subplot(nRow, nCol, 1)
    
    minBT = np.nanmin(obsCh[index_nh])
    maxBT = np.nanmax(obsCh[index_nh])

    m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=obsCh[index_nh],cmap=cmap, s=10, edgecolors='none', vmin=minBT, vmax=maxBT)
    plt.colorbar(cs, shrink = .8, orientation="horizontal")
    plt.title('Observed BT [K]')
                       
    #second row: plots for SH
    # Maps of atm corr
    for nE in range(0,nExp):

        btCorr = atmCorrList[nE][index_sh]        
        strExp = expLabelList[nE]

        dataMean = str('{:.2f}'.format(np.nanmean(btCorr))) 
        dataMin  = str('{:.2f}'.format(np.nanmin(btCorr)))   
        dataMax  = str('{:.2f}'.format(np.nanmax(btCorr)))   

        if (nE == 0):
            minAtmCorr = dataMin
            maxAtmCorr = dataMax
        
        fig.add_subplot(2, nCol, nCol+nE+1+1)
        
        mylat = lats[index_sh]
        mylon = lons[index_sh]

        m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)
        m.drawcoastlines()
        m.fillcontinents(color='gray')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))    
                
        x, y = m(mylon, mylat)
        cs = m.scatter(x,y,c=btCorr,cmap=cmap, s=10, edgecolors='none', vmin=minAtmCorr, vmax=maxAtmCorr)
        plt.colorbar(cs, shrink = .8, orientation="horizontal")
        
        strTitleExp = 'min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
        plt.title('Atm Corr [K] - ' + strExp + '\n' + strTitleExp)
        
    # maps of obs
    fig.add_subplot(nRow, nCol, nCol+1)
    
    minBT = np.nanmin(obsCh[index_sh])
    maxBT = np.nanmax(obsCh[index_sh])

    m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=obsCh[index_sh],cmap=cmap, s=10, edgecolors='none', vmin=minBT, vmax=maxBT)
    plt.colorbar(cs, shrink = .8, orientation="horizontal")
    plt.title('Observed BT [K]')
    
    # save plots....
    if plotSavePath is not None:
        
        if strPltId is not None:
            plotName = 'diag_atmCorr.'+ strSatId + '.' + strCh + '.' + strYYYYMMDDHHMM + '.' + strPltId +'.png' 
        else:
            plotName = 'diag_atmCorr.'+ strSatId + '.' + strCh + '.' + strYYYYMMDDHHMM + '.png' 

        plt.savefig(path.join(plotSavePath,plotName))         
        plt.close()
        LOG.info('Diagnostic plot for atm corr saved as: {0} '.format(path.join(plotSavePath,plotName)))
    
    return None

def diag_corrBT(lats,lons, obsCh, atmCorrList, strCh, strSatId, expLabelList, strYYYYMMDDHHMM, plotSavePath=None, strPltId=None):
    
    nExp = len(atmCorrList)
    
    nRow = 2
    nCol = nExp + 1
        
    index_nh = lats > 0
    index_sh = lats < 0

    cmap = plt.cm.jet
        
    width  = 6*nCol
    height = 12
    
    fig = plt.figure(figsize=(width,height),dpi=100)

    fig.suptitle('Freq (GHz) and pol: ' +  strCh + ' - ' + strSatId + ' ' + strYYYYMMDDHHMM, fontsize=14, fontweight='bold')

    #first row: plots for NH        
    minBT = np.nanmin(obsCh[index_nh])
    maxBT = np.nanmax(obsCh[index_nh])
    meanBT = np.nanmean(obsCh[index_nh])

    strMinBT = str('{:.2f}'.format(minBT)) 
    strMaxBT  = str('{:.2f}'.format(maxBT))   
    strMeanBT  = str('{:.2f}'.format(meanBT))   
    
    # Maps of atm corr
    for nE in range(0,nExp):

        btCorr = obsCh[index_nh] + atmCorrList[nE][index_nh]        
        strExp = expLabelList[nE]

        dataMean = str('{:.2f}'.format(np.nanmean(btCorr))) 
        dataMin  = str('{:.2f}'.format(np.nanmin(btCorr)))   
        dataMax  = str('{:.2f}'.format(np.nanmax(btCorr)))   
        
        if (nE == 0):
            minAtmCorr = dataMin
            maxAtmCorr = dataMax

        fig.add_subplot(nRow, nCol, nE+1+1)
        
        mylat = lats[index_nh]
        mylon = lons[index_nh]

        m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
        m.drawcoastlines()
        m.fillcontinents(color='gray')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))    
                
        x, y = m(mylon, mylat)
        cs = m.scatter(x,y,c=btCorr,cmap=cmap, s=10, edgecolors='none', vmin=minBT, vmax=maxBT)
        plt.colorbar(cs, shrink = .8, orientation="horizontal")
        
        strTitleExp = 'min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
        plt.title('Corr BT [K] - ' + strExp + '\n' + strTitleExp)
        
    # maps of obs
    fig.add_subplot(nRow, nCol, 1)
    
    m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=obsCh[index_nh],cmap=cmap, s=10, edgecolors='none', vmin=minBT, vmax=maxBT)
    plt.colorbar(cs, shrink = .8, orientation="horizontal")
    
    strTitleBT = 'min/max/mean: ' + strMinBT + '/' + strMaxBT + '/' + strMeanBT
    plt.title('Observed BT [K]'+ '\n' + strTitleBT)    
 
    #second row: plots for SH
    minBT = np.nanmin(obsCh[index_sh])
    maxBT = np.nanmax(obsCh[index_sh])
    meanBT = np.nanmean(obsCh[index_sh])

    strMinBT = str('{:.2f}'.format(minBT)) 
    strMaxBT  = str('{:.2f}'.format(maxBT))   
    strMeanBT  = str('{:.2f}'.format(meanBT))       
    # Maps of atm corr
    for nE in range(0,nExp):

        btCorr = obsCh[index_sh] + atmCorrList[nE][index_sh]        
        strExp = expLabelList[nE]

        dataMean = str('{:.2f}'.format(np.nanmean(btCorr))) 
        dataMin  = str('{:.2f}'.format(np.nanmin(btCorr)))   
        dataMax  = str('{:.2f}'.format(np.nanmax(btCorr)))   

        if (nE == 0):
            minAtmCorr = dataMin
            maxAtmCorr = dataMax
        
        fig.add_subplot(2, nCol, nCol+nE+1+1)
        
        mylat = lats[index_sh]
        mylon = lons[index_sh]

        m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)
        m.drawcoastlines()
        m.fillcontinents(color='gray')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))    
                
        x, y = m(mylon, mylat)
        cs = m.scatter(x,y,c=btCorr,cmap=cmap, s=10, edgecolors='none', vmin=minBT, vmax=maxBT)
        plt.colorbar(cs, shrink = .8, orientation="horizontal")
        
        strTitleExp = 'min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
        plt.title('Corr BT [K] - ' + strExp + '\n' + strTitleExp)
        
    # maps of obs
    fig.add_subplot(nRow, nCol, nCol+1)
    
    minBT = np.nanmin(obsCh[index_sh])
    maxBT = np.nanmax(obsCh[index_sh])

    m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=obsCh[index_sh],cmap=cmap, s=10, edgecolors='none', vmin=minBT, vmax=maxBT)
    plt.colorbar(cs, shrink = .8, orientation="horizontal")
 
    strTitleBT = 'min/max/mean: ' + strMinBT + '/' + strMaxBT + '/' + strMeanBT
    plt.title('Observed BT [K]'+ '\n' + strTitleBT)   
    
    # save plots....
    if plotSavePath is not None:
        
        if strPltId is not None:
            plotName = 'diag_corrBT.'+ strSatId + '.' + strCh + '.' + strYYYYMMDDHHMM + '.' + strPltId +'.png' 
        else:
            plotName = 'diag_corrBT.'+ strSatId + '.' + strCh + '.' + strYYYYMMDDHHMM + '.png' 

        plt.savefig(path.join(plotSavePath,plotName))         
        plt.close()
        LOG.info('Diagnostic plot for corr BT saved as: {0} '.format(path.join(plotSavePath,plotName)))
    
    return None

def diag_simBT(lats,lons, obsCh, simChList, strCh, strSatId, expLabelList, strYYYYMMDDHHMM, plotSavePath=None, strPltId=None):
    
    nExp = len(simChList)
    
    nRow = 2
    nCol = nExp + 1
        
    index_nh = lats > 0
    index_sh = lats < 0

    cmap = plt.cm.jet
        
    width  = 6*nCol
    height = 12
    
    fig = plt.figure(figsize=(width,height),dpi=100)

    fig.suptitle('Freq (GHz) and pol: ' +  strCh + ' - ' + strSatId + ' ' + strYYYYMMDDHHMM, fontsize=14, fontweight='bold')

    #first row: plots for NH        
    minBT = np.nanmin(obsCh[index_nh])
    maxBT = np.nanmax(obsCh[index_nh])
    meanBT = np.nanmean(obsCh[index_nh])

    strMinBT = str('{:.2f}'.format(minBT)) 
    strMaxBT  = str('{:.2f}'.format(maxBT))   
    strMeanBT  = str('{:.2f}'.format(meanBT))   

    # Maps of sim BT
    for nE in range(0,nExp):

        btSim = simChList[nE][index_nh]        
        strExp = expLabelList[nE]

        dataMean = str('{:.2f}'.format(np.nanmean(btSim))) 
        dataMin  = str('{:.2f}'.format(np.nanmin(btSim)))   
        dataMax  = str('{:.2f}'.format(np.nanmax(btSim)))   
        
        fig.add_subplot(nRow, nCol, nE+1+1)
        
        mylat = lats[index_nh]
        mylon = lons[index_nh]

        m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
        m.drawcoastlines()
        m.fillcontinents(color='gray')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))    
                
        x, y = m(mylon, mylat)
        cs = m.scatter(x,y,c=btSim,cmap=cmap, s=10, edgecolors='none', vmin=minBT, vmax=maxBT)
        plt.colorbar(cs, shrink = .8, orientation="horizontal")
        
        strTitleExp = 'min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
        plt.title('Sim BT [K] - ' + strExp + '\n' + strTitleExp)
        
    # maps of obs
    fig.add_subplot(nRow, nCol, 1)
    
    m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=obsCh[index_nh],cmap=cmap, s=10, edgecolors='none', vmin=minBT, vmax=maxBT)
    plt.colorbar(cs, shrink = .8, orientation="horizontal")
    
    strTitleBT = 'min/max/mean: ' + strMinBT + '/' + strMaxBT + '/' + strMeanBT
    plt.title('Observed BT [K]'+ '\n' + strTitleBT)
                       
    #second row: plots for SH
    minBT = np.nanmin(obsCh[index_sh])
    maxBT = np.nanmax(obsCh[index_sh])
    meanBT = np.nanmean(obsCh[index_sh])

    strMinBT = str('{:.2f}'.format(minBT)) 
    strMaxBT  = str('{:.2f}'.format(maxBT))   
    strMeanBT  = str('{:.2f}'.format(meanBT))   
    
    # Maps of sim BT
    for nE in range(0,nExp):

        btSim = simChList[nE][index_sh]        
        strExp = expLabelList[nE]

        dataMean = str('{:.2f}'.format(np.nanmean(btSim))) 
        dataMin  = str('{:.2f}'.format(np.nanmin(btSim)))   
        dataMax  = str('{:.2f}'.format(np.nanmax(btSim)))   
        
        fig.add_subplot(2, nCol, nCol+nE+1+1)
        
        mylat = lats[index_sh]
        mylon = lons[index_sh]

        m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)
        m.drawcoastlines()
        m.fillcontinents(color='gray')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))    
                
        x, y = m(mylon, mylat)
        cs = m.scatter(x,y,c=btSim,cmap=cmap, s=10, edgecolors='none', vmin=minBT, vmax=maxBT)
        plt.colorbar(cs, shrink = .8, orientation="horizontal")
        
        strTitleExp = 'min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
        plt.title('Sim BT [K] - ' + strExp + '\n' + strTitleExp)
        
    # maps of obs
    fig.add_subplot(nRow, nCol, nCol+1)
    
    minBT = np.nanmin(obsCh[index_sh])
    maxBT = np.nanmax(obsCh[index_sh])

    m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=obsCh[index_sh],cmap=cmap, s=10, edgecolors='none', vmin=minBT, vmax=maxBT)
    plt.colorbar(cs, shrink = .8, orientation="horizontal")

    strTitleBT = 'min/max/mean: ' + strMinBT + '/' + strMaxBT + '/' + strMeanBT
    plt.title('Observed BT [K]'+ '\n' + strTitleBT)
    
    # save plots....
    if plotSavePath is not None:
        
        if strPltId is not None:
            plotName = 'diag_simBT.'+ strSatId + '.' + strCh + '.' + strYYYYMMDDHHMM + '.' + strPltId +'.png' 
        else:
            plotName = 'diag_simBT.'+ strSatId + '.' + strCh + '.' + strYYYYMMDDHHMM + '.png' 

        plt.savefig(path.join(plotSavePath,plotName))         
        plt.close()
        LOG.info('Diagnostic plot for simulated BT saved as: {0} '.format(path.join(plotSavePath,plotName)))
    
    return None
    
def diag_emis(lats,lons, obsCh, emisChList, strCh, strSatId, expLabelList, strYYYYMMDDHHMM, plotSavePath=None, strPltId=None):
    
    nExp = len(emisChList)
    
    nRow = 2
    nCol = nExp + 1
        
    index_nh = lats > 0
    index_sh = lats < 0
    
    cmap = plt.cm.jet
        
    width  = 6*nCol
    height = 12
    
    fig = plt.figure(figsize=(width,height),dpi=100)

    fig.suptitle('Freq (GHz) and pol: ' +  strCh + ' - ' + strSatId + ' ' + strYYYYMMDDHHMM, fontsize=14, fontweight='bold')

    #first row: plots for NH        
    minBT = np.nanmin(obsCh[index_nh])
    maxBT = np.nanmax(obsCh[index_nh])
    meanBT = np.nanmean(obsCh[index_nh])

    strMinBT = str('{:.2f}'.format(minBT)) 
    strMaxBT  = str('{:.2f}'.format(maxBT))   
    strMeanBT  = str('{:.2f}'.format(meanBT))   

    # Maps of emis
    for nE in range(0,nExp):

        emis = emisChList[nE][index_nh]        
        strExp = expLabelList[nE]

        dataMean = str('{:.2f}'.format(np.nanmean(emis))) 
        dataMin  = str('{:.2f}'.format(np.nanmin(emis)))   
        dataMax  = str('{:.2f}'.format(np.nanmax(emis)))   

        if nE == 0:
            emisMin = dataMin
            emisMax = dataMax
        
        fig.add_subplot(nRow, nCol, nE+1+1)
        
        mylat = lats[index_nh]
        mylon = lons[index_nh]

        m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
        m.drawcoastlines()
        m.fillcontinents(color='gray')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))    
                
        x, y = m(mylon, mylat)
        cs = m.scatter(x,y,c=emis,cmap=cmap, s=10, edgecolors='none', vmin=emisMin, vmax=emisMax)
        plt.colorbar(cs, shrink = .8, orientation="horizontal")
        
        strTitleExp = 'min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
        plt.title('Emissivity - ' + strExp + '\n' + strTitleExp)
        
    # maps of obs
    fig.add_subplot(nRow, nCol, 1)
    
    m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=obsCh[index_nh],cmap=cmap, s=10, edgecolors='none', vmin=minBT, vmax=maxBT)
    plt.colorbar(cs, shrink = .8, orientation="horizontal")
    
    strTitleBT = 'min/max/mean: ' + strMinBT + '/' + strMaxBT + '/' + strMeanBT
    plt.title('Observed BT [K]'+ '\n' + strTitleBT)
                       
    #second row: plots for SH
    minBT = np.nanmin(obsCh[index_sh])
    maxBT = np.nanmax(obsCh[index_sh])
    meanBT = np.nanmean(obsCh[index_sh])

    strMinBT = str('{:.2f}'.format(minBT)) 
    strMaxBT  = str('{:.2f}'.format(maxBT))   
    strMeanBT  = str('{:.2f}'.format(meanBT))   
    
    # Maps of emis
    for nE in range(0,nExp):

        emis = emisChList[nE][index_sh]        
        strExp = expLabelList[nE]

        dataMean = str('{:.2f}'.format(np.nanmean(emis))) 
        dataMin  = str('{:.2f}'.format(np.nanmin(emis)))   
        dataMax  = str('{:.2f}'.format(np.nanmax(emis)))   

        if nE == 0:
            emisMin = dataMin
            emisMax = dataMax
        
        fig.add_subplot(2, nCol, nCol+nE+1+1)
        
        mylat = lats[index_sh]
        mylon = lons[index_sh]

        m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)
        m.drawcoastlines()
        m.fillcontinents(color='gray')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))    
                
        x, y = m(mylon, mylat)
        cs = m.scatter(x,y,c=emis,cmap=cmap, s=10, edgecolors='none', vmin=emisMin, vmax=emisMax)
        plt.colorbar(cs, shrink = .8, orientation="horizontal")
        
        strTitleExp = 'min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
        plt.title('Emissivity [K] - ' + strExp + '\n' + strTitleExp)
        
    # maps of obs
    fig.add_subplot(nRow, nCol, nCol+1)
    
    minBT = np.nanmin(obsCh[index_sh])
    maxBT = np.nanmax(obsCh[index_sh])

    m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=obsCh[index_sh],cmap=cmap, s=10, edgecolors='none', vmin=minBT, vmax=maxBT)
    plt.colorbar(cs, shrink = .8, orientation="horizontal")

    strTitleBT = 'min/max/mean: ' + strMinBT + '/' + strMaxBT + '/' + strMeanBT
    plt.title('Observed BT [K]'+ '\n' + strTitleBT)
    
    # save plots....
    if plotSavePath is not None:
        
        if strPltId is not None:
            plotName = 'diag_emis.'+ strSatId + '.' + strCh + '.' + strYYYYMMDDHHMM + '.' + strPltId +'.png' 
        else:
            plotName = 'diag_emis.'+ strSatId + '.' + strCh + '.' + strYYYYMMDDHHMM + '.png' 

        plt.savefig(path.join(plotSavePath,plotName))         
        plt.close()
        LOG.info('Diagnostic plot for emissivity BT saved as: {0} '.format(path.join(plotSavePath,plotName)))
    
    return None
    
def diag_lee_sohn_mw_emis(strSatId, strFreq, strYYYYMMDDHHMM, lats, lons, bt_v, bt_h, emis_v, emis_h, sic, plotSavePath=None, strPltId=None):
    
    nRow = 3
    nCol = 2

    index_nh = lats > 50
    index_sh = lats < -50
    
    cmap = plt.cm.jet
        
    ######################
    # FIG NH
    ######################
    width  = 12
    height = 16
        
    mylat = lats[index_nh]
    mylon = lons[index_nh]

    myBt_v = bt_v[index_nh]
    myBt_h = bt_h[index_nh]

    myEmis_v = emis_v[index_nh]
    myEmis_h = emis_h[index_nh]
    
    sic_nh = sic[index_nh]
    
    fig = plt.figure(figsize=(width,height),dpi=100)

    fig.suptitle('Freq (GHz): ' +  strFreq + ' - ' + strSatId + ' ' + strYYYYMMDDHHMM, fontsize=14, fontweight='bold')

    # V pol - BT
    fig.add_subplot(nRow, nCol, 1)

    plt.title('Observed BT V pol')
 
    m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=myBt_v,cmap=cmap, s=10, edgecolors='none', vmin=np.nanmin(myBt_v), vmax=np.nanmax(myBt_v))
    plt.colorbar(cs, shrink = .5, orientation="horizontal")

    # V pol - emis
    fig.add_subplot(nRow, nCol, 2)

    dataMean = str('{:.2f}'.format(np.nanmean(myEmis_v))) 
    dataMin  = str('{:.2f}'.format(np.nanmin(myEmis_v)))   
    dataMax  = str('{:.2f}'.format(np.nanmax(myEmis_v)))   

    strTitle = 'min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
    plt.title('Emis V pol'+ '\n' + strTitle)

    m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=myEmis_v,cmap=cmap, s=10, edgecolors='none', vmin=0.4, vmax=1)
    plt.colorbar(cs, shrink = .5, orientation="horizontal")

    #emis_gt1 = myEmis_v > 1
    #emis_lt0 = myEmis_v < 0
    
    # #black gt1
    # x, y = m(mylon[emis_gt1], mylat[emis_gt1])
    # cs = m.scatter(x,y,c=myEmis_v[emis_gt1],cmap=cmap, s=10, edgecolors='k', vmin=0.4, vmax=1)

    # # white lt0    
    # x, y = m(mylon[emis_lt0], mylat[emis_lt0])
    # cs = m.scatter(x,y,c=myEmis_v[emis_lt0],cmap=cmap, s=10, edgecolors='white', vmin=0.4, vmax=1)
        
    # H pol - BT
    fig.add_subplot(nRow, nCol, 3)

    plt.title('Observed BT H pol')

    m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=myBt_h,cmap=cmap, s=10, edgecolors='none', vmin=np.nanmin(myBt_h), vmax=np.nanmax(myBt_h))
    plt.colorbar(cs, shrink = .5, orientation="horizontal")

    # H pol - emis
    fig.add_subplot(nRow, nCol, 4)

    dataMean = str('{:.2f}'.format(np.nanmean(myEmis_h))) 
    dataMin  = str('{:.2f}'.format(np.nanmin(myEmis_h)))   
    dataMax  = str('{:.2f}'.format(np.nanmax(myEmis_h)))   

    strTitle = 'min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
    plt.title('Emis H pol'+ '\n' + strTitle)
 
    m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=myEmis_h,cmap=cmap, s=10, edgecolors='none', vmin=0.4, vmax=1)
    plt.colorbar(cs, shrink = .5, orientation="horizontal")

    #emis_gt1 = myEmis_h > 1
    #emis_lt0 = myEmis_h < 0
    
    # #black gt1
    # x, y = m(mylon[emis_gt1], mylat[emis_gt1])
    # cs = m.scatter(x,y,c=myEmis_h[emis_gt1],cmap=cmap, s=10, edgecolors='k', vmin=0.4, vmax=1)

    # # white lt0    
    # x, y = m(mylon[emis_lt0], mylat[emis_lt0])
    # cs = m.scatter(x,y,c=myEmis_h[emis_lt0],cmap=cmap, s=10, edgecolors='white', vmin=0.4, vmax=1)

    # sic
    fig.add_subplot(nRow, nCol, 5)

    plt.title('SIC')

    m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=sic_nh,cmap=cmap, s=10, edgecolors='none', vmin=np.nanmin(sic_nh), vmax=np.nanmax(sic_nh))
    plt.colorbar(cs, shrink = .5, orientation="horizontal")

    # save plots....
    if plotSavePath is not None:
        
        if strPltId is not None:
            plotName = 'diag_lee_sohn_emis.NH.'+ strSatId + '.' + strFreq + '.' + strYYYYMMDDHHMM + '.' + strPltId +'.png' 
        else:
            plotName = 'diag_lee_sohn_emis.NH.'+ strSatId + '.' + strFreq + '.' + strYYYYMMDDHHMM + '.png' 

        plt.savefig(path.join(plotSavePath,plotName))         
        plt.close()
        LOG.info('Diagnostic plot for lee sohn emissivity saved as: {0} '.format(path.join(plotSavePath,plotName)))

    ######################
    # FIG SH
    ######################
    width  = 12
    height = 16
        
    mylat = lats[index_sh]
    mylon = lons[index_sh]

    myBt_v = bt_v[index_sh]
    myBt_h = bt_h[index_sh]

    myEmis_v = emis_v[index_sh]
    myEmis_h = emis_h[index_sh]
    
    sic_sh = sic[index_sh]
    
    fig = plt.figure(figsize=(width,height),dpi=100)

    fig.suptitle('Freq (GHz): ' +  strFreq + ' - ' + strSatId + ' ' + strYYYYMMDDHHMM, fontsize=14, fontweight='bold')

    # V pol - BT
    fig.add_subplot(nRow, nCol, 1)

    plt.title('Observed BT V pol')
 
    m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=myBt_v,cmap=cmap, s=10, edgecolors='none', vmin=np.nanmin(myBt_v), vmax=np.nanmax(myBt_v))
    plt.colorbar(cs, shrink = .5, orientation="horizontal")

    # V pol - emis
    fig.add_subplot(nRow, nCol, 2)

    dataMean = str('{:.2f}'.format(np.nanmean(myEmis_v))) 
    dataMin  = str('{:.2f}'.format(np.nanmin(myEmis_v)))   
    dataMax  = str('{:.2f}'.format(np.nanmax(myEmis_v)))   

    strTitle = 'min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
    plt.title('Emis V pol'+ '\n' + strTitle)

    m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=myEmis_v,cmap=cmap, s=10, edgecolors='none', vmin=0.4, vmax=1)
    plt.colorbar(cs, shrink = .5, orientation="horizontal")

    #emis_gt1 = myEmis_v > 1
    #emis_lt0 = myEmis_v < 0
    
    # #black gt1
    # x, y = m(mylon[emis_gt1], mylat[emis_gt1])
    # cs = m.scatter(x,y,c=myEmis_v[emis_gt1],cmap=cmap, s=10, edgecolors='k', vmin=0.4, vmax=1)

    # # white lt0    
    # x, y = m(mylon[emis_lt0], mylat[emis_lt0])
    # cs = m.scatter(x,y,c=myEmis_v[emis_lt0],cmap=cmap, s=10, edgecolors='white', vmin=0.4, vmax=1)
        
    # H pol - BT
    fig.add_subplot(nRow, nCol, 3)

    plt.title('Observed BT H pol')

    m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=myBt_h,cmap=cmap, s=10, edgecolors='none', vmin=np.nanmin(myBt_h), vmax=np.nanmax(myBt_h))
    plt.colorbar(cs, shrink = .5, orientation="horizontal")

    # H pol - emis
    fig.add_subplot(nRow, nCol, 4)

    dataMean = str('{:.2f}'.format(np.nanmean(myEmis_h))) 
    dataMin  = str('{:.2f}'.format(np.nanmin(myEmis_h)))   
    dataMax  = str('{:.2f}'.format(np.nanmax(myEmis_h)))   

    strTitle = 'min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
    plt.title('Emis H pol'+ '\n' + strTitle)
 
    m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=myEmis_h,cmap=cmap, s=10, edgecolors='none', vmin=0.4, vmax=1)
    plt.colorbar(cs, shrink = .5, orientation="horizontal")

    #emis_gt1 = myEmis_h > 1
    #emis_lt0 = myEmis_h < 0
    
    # #black gt1
    # x, y = m(mylon[emis_gt1], mylat[emis_gt1])
    # cs = m.scatter(x,y,c=myEmis_h[emis_gt1],cmap=cmap, s=10, edgecolors='k', vmin=0.4, vmax=1)

    # # white lt0    
    # x, y = m(mylon[emis_lt0], mylat[emis_lt0])
    # cs = m.scatter(x,y,c=myEmis_h[emis_lt0],cmap=cmap, s=10, edgecolors='white', vmin=0.4, vmax=1)

    # sic
    fig.add_subplot(nRow, nCol, 5)

    plt.title('SIC')

    m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)
    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon, mylat)
    cs = m.scatter(x,y,c=sic_sh,cmap=cmap, s=10, edgecolors='none', vmin=np.nanmin(sic_nh), vmax=np.nanmax(sic_nh))
    plt.colorbar(cs, shrink = .5, orientation="horizontal")

    # save plots....
    if plotSavePath is not None:
        
        if strPltId is not None:
            plotName = 'diag_lee_sohn_emis.SH.'+ strSatId + '.' + strFreq + '.' + strYYYYMMDDHHMM + '.' + strPltId +'.png' 
        else:
            plotName = 'diag_lee_sohn_emis.SH.'+ strSatId + '.' + strFreq + '.' + strYYYYMMDDHHMM + '.png' 

        plt.savefig(path.join(plotSavePath,plotName))         
        plt.close()
        LOG.info('Diagnostic plot for lee sohn emissivity saved as: {0} '.format(path.join(plotSavePath,plotName)))
    
    return None

def diag_dep_histograms(strFigSubtitle, strHem, chsList, lats, depList, sic=None, sicThr=None, plotSavePath=None, plotName=None):
    
    # single plot comapring all departures provided in depList
    width  = 12
    height = 14
    
    fig = plt.figure(figsize=(width,height),dpi=100)

    fig.suptitle(strFigSubtitle, fontsize=14, fontweight='bold')

    histMin = -50
    histMax = 50 
    histBinSize = 1    
    myBins = np.arange(histMin, histMax+histBinSize, histBinSize)     
    
    nChs = len(chsList)
    mylat = lats
    
    for nCh in range(0,nChs):

        strCh = chsList[nCh]
        myDep = depList[nCh]
        
        if sic is not None and sicThr is not None:   
            
            if sicThr > 0:
                LOG.info('departures are selected according to sic threshold greater than {0} '.format(sicThr))
                index = np.where(sic > sicThr)
                myDep = myDep[index[0]]
                mylat = lats[index[0]]
            else:
                LOG.info('departures are selected over ocean (where sic = 0)')
                index = np.where(sic == 0)
                myDep = myDep[index[0]]
                mylat = lats[index[0]]

        if strHem == 'NH':
            indexlat = mylat > 0
        else:
            indexlat = mylat < 0                
                
        dataMean = str('{:.2f}'.format(np.nanmean(myDep[indexlat]))) 
        dataStd  = str('{:.2f}'.format(np.nanstd(myDep[indexlat])))   
        dataSkew = str('{:.2f}'.format(skew(myDep[indexlat])))         

        plt.hist(myDep[indexlat], bins=myBins, histtype='step',color=colorList[nCh],density=False, linewidth=1.2,label=strCh)
        
        if (nCh==0):
            strTitleExp = 'mean/std/skew: ' + dataMean + '/' + dataStd + '/' + dataSkew
        else:
            strTitle = 'mean/std/skew: ' + dataMean + '/' + dataStd + '/' + dataSkew
            
            strTitleExp = strTitleExp + '\n' + strTitle
            
    #plt.yscale('log')
    plt.ticklabel_format(style='sci',axis='y',scilimits=(3,3))
    plt.grid(True)
    plt.title(strTitleExp)  #, fontsize=8, fontweight='bold')
    plt.xlabel('Departure (Obs-Sim) [K] - Hem: ' + strHem)
    plt.ylabel('Counts')
    plt.legend(loc=2)

    # save plots....
    if plotSavePath is not None and plotName is not None:        
        fig.patch.set_facecolor('xkcd:white')        
        #plt.rcParams['axes.facecolor'] = 'white'
        plt.savefig(path.join(plotSavePath,plotName))         
        plt.close()
        LOG.info('Histogram of departures diagnostic plot saved as: {0} '.format(path.join(plotSavePath,plotName)))
    
    return None

def diag_dep_maps(strFigSubtitle, strHem, chsList, lats, lons, depList, sic=None, sicThr=None, plotSavePath=None, plotName=None):

    # single plot comapring all departures provided in depList
    nChs = len(chsList)
    
    nRow = 2
    nCol = int(np.ceil(nChs/nRow))
        

    cmap = plt.cm.seismic    
    bounds = [-20, -15, -10, -5, -2, 0, 2, 5, 10, 15, 20]
        
    width  = 6*nCol
    height = 12
      
    fig = plt.figure(figsize=(width,height),dpi=100)

    fig.suptitle(strFigSubtitle, fontsize=14, fontweight='bold')

    nChs = len(chsList)
    mylat = lats
    mylon = lons
    
    for nCh in range(0,nChs):

        strCh = chsList[nCh]
        myDep = depList[nCh]
        
        if sic is not None and sicThr is not None:   
            
            if sicThr > 0:
                LOG.info('departures are selected according to sic threshold greater than {0} '.format(sicThr))
                index = np.where(sic > sicThr)
                myDep = myDep[index[0]]
                mylat = lats[index[0]]
                mylon = lons[index[0]]
            else:
                LOG.info('departures are selected over ocean (where sic = 0)')
                index = np.where(sic == 0)
                myDep = myDep[index[0]]
                mylat = lats[index[0]]
                mylon = lons[index[0]]
                
        if strHem == 'NH':
            indexlat = mylat > 0
        else:
            indexlat = mylat < 0
                
        dataMin = str('{:.2f}'.format(np.nanmin(myDep[indexlat]))) 
        dataMax  = str('{:.2f}'.format(np.nanmax(myDep[indexlat])))   
        dataMean = str('{:.2f}'.format(np.nanmean(myDep[indexlat]))) 

        fig.add_subplot(nRow, nCol, nCh+1)
        
        if strHem == 'NH':
            m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
        else:
            m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)

        m.drawcoastlines()
        m.fillcontinents(color='gray')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))    
                
        x, y = m(mylon[indexlat], mylat[indexlat])
        cs = m.scatter(x,y,c=myDep[indexlat],cmap=cmap, s=10, edgecolors='none', vmin=bounds[0], vmax=bounds[len(bounds)-1])
        plt.colorbar(cs, shrink = .8, orientation="horizontal",boundaries=bounds,ticks=bounds, spacing='proportional')

        strTitle = strCh + '\n Dep min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
        plt.title(strTitle)
 
    # save plots....
    if plotSavePath is not None and plotName is not None:   
        fig.patch.set_facecolor('xkcd:white')
        #plt.rcParams['axes.facecolor'] = 'white'        
        plt.savefig(path.join(plotSavePath,plotName))         
        plt.close()
        LOG.info('Maps of departures diagnostic plot saved as: {0} '.format(path.join(plotSavePath,plotName)))
    
    return None

def diag_obsBT_maps(strFigSubtitle, strHem, chsList, lats, lons, obsList, simBTList=None, sic=None, sicThr=None, plotSavePath=None, plotName=None, plotNameSimBT=None):

    # single plot comapring all departures provided in depList
    nChs = len(chsList)
    
    nRow = 2
    nCol = int(np.ceil(nChs/nRow))
        
    cmap = plt.cm.jet    
        
    width  = 6*nCol
    height = 12
      
    fig = plt.figure(figsize=(width,height),dpi=100)

    fig.suptitle(strFigSubtitle, fontsize=14, fontweight='bold')

    nChs = len(chsList)
    mylat = lats
    mylon = lons
    
    for nCh in range(0,nChs):

        strCh = chsList[nCh]
        myObs = obsList[nCh]
        
        if sic is not None and sicThr is not None:   
            
            if sicThr > 0:
                LOG.info('observations are selected according to sic threshold greater than {0} '.format(sicThr))
                index = np.where(sic > sicThr)
                myObs = myObs[index[0]]
                mylat = lats[index[0]]
                mylon = lons[index[0]]
            else:
                LOG.info('observations are selected over ocean (where sic = 0)')
                index = np.where(sic == 0)
                myObs = myObs[index[0]]
                mylat = lats[index[0]]
                mylon = lons[index[0]]
                
        if strHem == 'NH':
            indexlat = mylat > 0
        else:
            indexlat = mylat < 0
                
        vmin = int(np.nanmin(myObs[indexlat]))        
        vmax = int(np.nanmax(myObs[indexlat]))
        
        dataMin = str('{:.2f}'.format(np.nanmin(myObs[indexlat]))) 
        dataMax  = str('{:.2f}'.format(np.nanmax(myObs[indexlat])))   
        dataMean = str('{:.2f}'.format(np.nanmean(myObs[indexlat]))) 

        fig.add_subplot(nRow, nCol, nCh+1)
        
        if strHem == 'NH':
            m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
        else:
            m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)

        m.drawcoastlines()
        m.fillcontinents(color='gray')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))    
                
        x, y = m(mylon[indexlat], mylat[indexlat])
        cs = m.scatter(x,y,c=myObs[indexlat],cmap=cmap, s=10, edgecolors='none', vmin=vmin, vmax=vmax)
        plt.colorbar(cs, shrink = .8, orientation="horizontal")

        strTitle = strCh + '\n Obs BT min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
        plt.title(strTitle)
 
    # save plots....
    if plotSavePath is not None and plotName is not None:   
        fig.patch.set_facecolor('xkcd:white')
        #plt.rcParams['axes.facecolor'] = 'white'        
        plt.savefig(path.join(plotSavePath,plotName))         
        plt.close()
        LOG.info('Maps of observed BT diagnostic plot saved as: {0} '.format(path.join(plotSavePath,plotName)))

    ###########
    # if simBT provided, geenrate equivalent plot to match scale (min/max) of observed BT
    ###########
    if simBTList is not None:

        fig = plt.figure(figsize=(width,height),dpi=100)
    
        fig.suptitle(strFigSubtitle, fontsize=14, fontweight='bold')
    
        nChs = len(chsList)
        mylat = lats
        mylon = lons
        
        for nCh in range(0,nChs):
    
            strCh = chsList[nCh]
            myObs = obsList[nCh]
            mySim = simBTList[nCh]
            
            if sic is not None and sicThr is not None:   
                
                if sicThr > 0:
                    LOG.info('observations are selected according to sic threshold greater than {0} '.format(sicThr))
                    index = np.where(sic > sicThr)
                    myObs = myObs[index[0]]
                    mylat = lats[index[0]]
                    mylon = lons[index[0]]
                    mySim = mySim[index[0]]                    
                else:
                    LOG.info('observations are selected over ocean (where sic = 0)')
                    index = np.where(sic == 0)
                    myObs = myObs[index[0]]
                    mylat = lats[index[0]]
                    mylon = lons[index[0]]
                    mySim = mySim[index[0]]                    
                    
            if strHem == 'NH':
                indexlat = mylat > 0
            else:
                indexlat = mylat < 0
                    
            vmin = int(np.nanmin(myObs[indexlat]))        
            vmax = int(np.nanmax(myObs[indexlat]))
            
            dataMin = str('{:.2f}'.format(np.nanmin(mySim[indexlat]))) 
            dataMax  = str('{:.2f}'.format(np.nanmax(mySim[indexlat])))   
            dataMean = str('{:.2f}'.format(np.nanmean(mySim[indexlat]))) 
    
            fig.add_subplot(nRow, nCol, nCh+1)
            
            if strHem == 'NH':
                m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
            else:
                m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)
    
            m.drawcoastlines()
            m.fillcontinents(color='gray')
            # draw parallels and meridians.
            m.drawparallels(np.arange(-80.,81.,20.))
            m.drawmeridians(np.arange(-180.,181.,20.))    
                    
            x, y = m(mylon[indexlat], mylat[indexlat])
            cs = m.scatter(x,y,c=mySim[indexlat],cmap=cmap, s=10, edgecolors='none', vmin=vmin, vmax=vmax)
            plt.colorbar(cs, shrink = .8, orientation="horizontal")
    
            strTitle = strCh + '\n Sim BT min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
            plt.title(strTitle)
     
        # save plots....
        if plotSavePath is not None and plotNameSimBT is not None:   
            fig.patch.set_facecolor('xkcd:white')
            #plt.rcParams['axes.facecolor'] = 'white'        
            plt.savefig(path.join(plotSavePath,plotNameSimBT))         
            plt.close()
            LOG.info('Maps of simulated BT diagnostic plot saved as: {0} '.format(path.join(plotSavePath,plotNameSimBT)))        

    return None

def diag_emis_maps(strFigSubtitle, strHem, chsList, lats, lons, emisList, sic=None, sicThr=None, plotSavePath=None, plotName=None):

    # single plot comapring all departures provided in depList
    nChs = len(chsList)
    
    nRow = 2
    nCol = int(np.ceil(nChs/nRow))
        
    cmap = plt.cm.jet    
        
    width  = 6*nCol
    height = 12
      
    fig = plt.figure(figsize=(width,height),dpi=100)

    fig.suptitle(strFigSubtitle, fontsize=14, fontweight='bold')

    nChs = len(chsList)
    mylat = lats
    mylon = lons
    
    for nCh in range(0,nChs):

        strCh = chsList[nCh]
        myEmis = emisList[nCh]
        
        if sic is not None and sicThr is not None:   
            
            if sicThr > 0:
                LOG.info('observations are selected according to sic threshold greater than {0} '.format(sicThr))
                index = np.where(sic > sicThr)
                myEmis = myEmis[index[0]]
                mylat = lats[index[0]]
                mylon = lons[index[0]]
            else:
                LOG.info('observations are selected over ocean (where sic = 0)')
                index = np.where(sic == 0)
                myEmis = myEmis[index[0]]
                mylat = lats[index[0]]
                mylon = lons[index[0]]
                
        if strHem == 'NH':
            indexlat = mylat > 0
        else:
            indexlat = mylat < 0

        vmin = np.nanmin(myEmis[indexlat])        
        vmax = np.nanmax(myEmis[indexlat])
                
        dataMin = str('{:.2f}'.format(np.nanmin(myEmis[indexlat]))) 
        dataMax  = str('{:.2f}'.format(np.nanmax(myEmis[indexlat])))   
        dataMean = str('{:.2f}'.format(np.nanmean(myEmis[indexlat]))) 

        fig.add_subplot(nRow, nCol, nCh+1)
        
        if strHem == 'NH':
            m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
        else:
            m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)

        m.drawcoastlines()
        m.fillcontinents(color='gray')
        # draw parallels and meridians.
        m.drawparallels(np.arange(-80.,81.,20.))
        m.drawmeridians(np.arange(-180.,181.,20.))    
                
        x, y = m(mylon[indexlat], mylat[indexlat])
        cs = m.scatter(x,y,c=myEmis[indexlat],cmap=cmap, s=10, edgecolors='none', vmin=vmin, vmax=vmax)
        plt.colorbar(cs, shrink = .8, orientation="horizontal")

        strTitle = strCh + '\n Emis min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
        plt.title(strTitle)
 
    # save plots....
    if plotSavePath is not None and plotName is not None:   
        fig.patch.set_facecolor('xkcd:white')
        #plt.rcParams['axes.facecolor'] = 'white'        
        plt.savefig(path.join(plotSavePath,plotName))         
        plt.close()
        LOG.info('Maps of emis diagnostic plot saved as: {0} '.format(path.join(plotSavePath,plotName)))
    
    return None

def diag_nwpSurf_map(strFigSubtitle, strHem, strNwpPar, nwpData, lats, lons, sic=None, sicThr=None, plotSavePath=None, plotName=None):
        
    cmap = plt.cm.jet    
        
    width  = 10
    height = 8
      
    fig = plt.figure(figsize=(width,height),dpi=100)

    fig.suptitle(strFigSubtitle, fontsize=14, fontweight='bold')

    mylat = lats
    mylon = lons
    
    
    if sic is not None and sicThr is not None:   
        
        if sicThr > 0:
            LOG.info('data are selected according to sic threshold greater than {0} '.format(sicThr))
            index = np.where(sic > sicThr)
            nwpData = nwpData[index[0]]
            mylat = lats[index[0]]
            mylon = lons[index[0]]
        else:
            LOG.info('data are selected over ocean (where sic = 0)')
            index = np.where(sic == 0)
            nwpData = nwpData[index[0]]
            mylat = lats[index[0]]
            mylon = lons[index[0]]
            
    if strHem == 'NH':
        indexlat = mylat > 0
    else:
        indexlat = mylat < 0
            
    vmin = int(np.nanmin(nwpData[indexlat]))        
    vmax = int(np.nanmax(nwpData[indexlat]))
    
    dataMin = str('{:.2f}'.format(np.nanmin(nwpData[indexlat]))) 
    dataMax  = str('{:.2f}'.format(np.nanmax(nwpData[indexlat])))   
    dataMean = str('{:.2f}'.format(np.nanmean(nwpData[indexlat]))) 
    
    if strHem == 'NH':
        m = Basemap(projection='npstere',boundinglat=40,lon_0=270,resolution='l',round=True)
    else:
        m = Basemap(projection='spstere',boundinglat=-40,lon_0=270,resolution='l',round=True)

    m.drawcoastlines()
    m.fillcontinents(color='gray')
    # draw parallels and meridians.
    m.drawparallels(np.arange(-80.,81.,20.))
    m.drawmeridians(np.arange(-180.,181.,20.))    
            
    x, y = m(mylon[indexlat], mylat[indexlat])
    cs = m.scatter(x,y,c=nwpData[indexlat],cmap=cmap, s=10, edgecolors='none', vmin=vmin, vmax=vmax)
    plt.colorbar(cs, shrink = .8, orientation="horizontal")

    strTitle = strNwpPar + ' min/max/mean: ' + dataMin + '/' + dataMax + '/' + dataMean
    plt.title(strTitle)
 
    # save plots....
    if plotSavePath is not None and plotName is not None:   
        fig.patch.set_facecolor('xkcd:white')
        #plt.rcParams['axes.facecolor'] = 'white'        
        plt.savefig(path.join(plotSavePath,plotName))         
        plt.close()
        LOG.info('Maps of NWP surf diagnostic plot saved as: {0} '.format(path.join(plotSavePath,plotName)))    

    return None

############
# This is not a plotting function, but added the possibility to save a netcdf containing simulated BT and emis
############

def writeNcSimBT(myLats2nc,
                 myLons2nc,
                 simBTarr,
                 strChsList,
                 saveToPath,
                 strSensor,
                 strYYYYMMDDHHMM,    
                 scanStartTime,
                 scanEndTime,
                 strRTtype,
                 emisArr=None,
                 strEmisType=None,
                 strDescr=None):
    
        
    strYYYYMMDDHHMMs = scanStartTime.strftime("%Y%d%m%H%M")
    strYYYYMMDDHHMMe = scanEndTime.strftime("%Y%d%m%H%M")
    
    ncFilename = 'simBT.' + strRTtype + '.' + strSensor + '.' +  strYYYYMMDDHHMM + '.nc'    
    
    if not path.exists(saveToPath):
        makedirs(saveToPath)
        
    myncFile = path.join(saveToPath, ncFilename)
        
    # remove file if exists
    if path.exists(myncFile):
        remove(myncFile)
    
    # open file
    netCDF_file = Dataset(myncFile, 'w', format = 'NETCDF4_CLASSIC')    
            
    netCDF_file.sensor   = strSensor
    netCDF_file.sensor_scan_start_time = strYYYYMMDDHHMMs
    netCDF_file.sensor_scan_end_time   = strYYYYMMDDHHMMe
    if strDescr is not None:
        netCDF_file.description = strDescr
    else:
        netCDF_file.description = "Simulated brightness temperatures [k] over ocean/seaice above/below 40N/40S"
    
    netCDF_file.rt_type  =  strRTtype         
    
    if strEmisType is not None:
        netCDF_file.surface_emis = "Surface emissivity: " + strEmisType         
                        
    netCDF_file.creation_date = str(datetime.datetime.now())[0:19] + ' CET'    


    netCDF_file.createDimension('obs',  len(myLats2nc))
    
    lat = netCDF_file.createVariable('lat', 'f4', ('obs',), zlib=True,least_significant_digit=3)
    lat.long_name = "latitude of satellite swath"
    netCDF_file.variables['lat'][:] = myLats2nc[:]

    lon = netCDF_file.createVariable('lon', 'f4', ('obs',), zlib=True,least_significant_digit=3)
    lon.long_name = "longitude of satellite swath"
    netCDF_file.variables['lon'][:] = myLons2nc[:]

    nChs = len(strChsList)
    
    for nCh in range(0,nChs):

        strCh   = strChsList[nCh]
        
        netCDF_file.createVariable(strCh, 'f4', ('obs', ),zlib=True,least_significant_digit=3)
        netCDF_file.variables[strCh][:] = simBTarr[:,nCh]     
        
        if emisArr is not None:
            
            netCDF_file.createVariable('emis_'+strCh , 'f4', ('obs', ),zlib=True,least_significant_digit=3)
            netCDF_file.variables['emis_'+strCh][:] = emisArr[:,nCh]      
            
        
    netCDF_file.close()

    LOG.info('nc simulated BT file saved as: {0} '.format(myncFile))

    return myncFile

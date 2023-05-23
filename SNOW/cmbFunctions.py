
import numpy as np
import netCDF4 as nc
def readCMB(fname): # reads relevant data from the CMB file
    fh_cmb=nc.Dataset(fname)
    dZ=0 # offset between the CMB and DPR
    qv=fh_cmb["KuKaGMI/vaporDensity"][:,:,:]
    press=fh_cmb["KuKaGMI/airPressure"][:,:,:]
    envNodes=fh_cmb["KuKaGMI/envParamNode"][:,:,:]
    airTemp=fh_cmb["KuKaGMI/airTemperature"][:,:,:]
    skTemp=fh_cmb["KuKaGMI/skinTemperature"][:,:]
    binNodes=fh_cmb["KuKaGMI/phaseBinNodes"][:,:]
    pwc=fh_cmb["KuKaGMI/precipTotWaterCont"][:,:,:]
    sfcEmiss=fh_cmb["KuKaGMI/surfEmissivity"][:,:,:]
    dm=fh_cmb["KuKaGMI/precipTotDm"][:,:,:]
    cldw=fh_cmb["KuKaGMI/cloudLiqWaterCont"][:,:,:]
    sfcBin=fh_cmb["KuKaGMI/Input/surfaceRangeBin"][:,:,:]
    zCorrected=fh_cmb["KuGMI/correctedReflectFactor"][:,:,:]+dZ
    pType=fh_cmb["KuKaGMI/Input/precipitationType"][:,:]
    lon=fh_cmb["KuKaGMI/Longitude"][:,:]
    lat=fh_cmb["KuKaGMI/Latitude"][:,:]
    return qv,press,envNodes,airTemp,skTemp,binNodes,pwc,sfcEmiss,dm,cldw,sfcBin,zCorrected,pType,lon,lat

import sys
sys.path.append('/Users/mgrecu/myPythonPackages/')
import lkTables as lkT
lkTables=lkT.scattTables()
from bisectm import bisectm

def simulateZKa(binNodes,zCorrected,dm,dm_factor,lkTables,pType): # calculates Ka-band reflectivity from attenuated corrected Ku-band reflectivity
    # inputs:
    # binNodes: 2D array of bin nodes that define the snow (0 to 1), mixed phase (1 to 3), and rain (3 to 4)
    # zCorrected: 4D array of Ku-band reflectivity (ns,nr,nbins,2) where 0 is for Ku and Ka is 1
    # dm: 3D array of CMB mean drop size (ns,nr,nbins)
    # lkTables: look-up tables for scattering properties
    # pType: 2D array of precipitation type (ns,nr)
    # output:
    # zKaSim: 3D array of simulated Ka-band reflectivity (ns,nr,nbins)
    ns,nr,nbins=zCorrected.shape[:3]
    zKaSim=np.zeros((ns,nr,nbins),float)-99
    pRate=np.zeros((ns,nr,nbins),float)-99
    pwc=np.zeros((ns,nr,nbins),float)-99
    zkaSfc=np.zeros((ns,nr),float)-99
    #print(pType.max(),ns,nr)
    dr=0.25
    for i in range(zCorrected.shape[0]):
        for j in range(zCorrected.shape[1]):
            if pType[i,j]>0:
                #print(binNodes[i,j,:])
                piaKa=0.0
                ik_set=0
                for k in range(binNodes[i,j,0],min(binNodes[i,j,4],binNodes[i,j,1])+1):  # added range limit to avoid out of bounds
                    if dm[i,j,k]>0:
                        ind=bisectm(lkTables.dms,253,dm_factor*dm[i,j,k])
                        dnw=(zCorrected[i,j,k]-lkTables.zKuS[ind])/10
                        attKa=lkTables.attKaS[ind]*10**dnw
                        piaKa+=attKa*dr
                        zKaSim[i,j,k]=10*dnw+lkTables.zKaS[ind]-piaKa
                        pRate[i,j,k]=lkTables.snowRate[ind]*10**dnw
                        pwc[i,j,k]=lkTables.swc[ind]*10**dnw
                        piaKa+=attKa*dr
                        ik_set=1
                for k in range(binNodes[i,j,1]+1,min(binNodes[i,j,4],binNodes[i,j,3])):  # added range limit to avoid out of bounds
                    if dm[i,j,k]>0:
                        ind1=bisectm(lkTables.dms,253,dm_factor*dm[i,j,k])
                        ind2=bisectm(lkTables.dmr,289,dm[i,j,k])
                        fsnow=((binNodes[i,j,3]-k)/(binNodes[i,j,3]-binNodes[i,j,1]))**0.5
                        zKuMPH=np.log10(fsnow*10**(0.1*lkTables.zKuS[ind1])+(1-fsnow)*10**(0.1*lkTables.zKuR[ind2]))*10
                        dnw=(zCorrected[i,j,k]-zKuMPH)/10
                        zKaMPH=np.log10(fsnow*10**(0.1*lkTables.zKaS[ind1])+(1-fsnow)*10**(0.1*lkTables.zKaR[ind2]))*10

                        attKa=fsnow*lkTables.attKaBB[ind1]*10**dnw+(1-fsnow)*lkTables.attKaR[ind2]*10**dnw
                        piaKa+=attKa*dr
                        zKaSim[i,j,k]=10*dnw+zKaMPH-piaKa
                        piaKa+=attKa*dr
                        ik_set=1
                for k in range(binNodes[i,j,3],binNodes[i,j,4]):
                    if dm[i,j,k]>0:
                        ind=bisectm(lkTables.dmr,289,dm[i,j,k])
                        dnw=(zCorrected[i,j,k]-lkTables.zKuR[ind])/10
                        attKa=lkTables.attKaR[ind]*10**dnw
                        piaKa+=attKa*dr
                        zKaSim[i,j,k]=10*dnw+lkTables.zKaR[ind]-piaKa
                        piaKa+=attKa*dr
                        ik_set=1
                if ik_set==1:
                    zkaSfc[i,j]=zKaSim[i,j,k]
    return zKaSim,zkaSfc,pRate,pwc
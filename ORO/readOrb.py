from netCDF4 import Dataset
from scipy.ndimage import gaussian_filter
import glob
fs=sorted(glob.glob("2A-CS/2A*CONUS*"))
import numpy as np
from numpy import *
def readOrb(orb):
    fh=Dataset(orb)
    #print(fh)
    sfcPrecip=fh['FS/SLV/precipRateNearSurface'][:,:]
    precipRate=fh['FS/SLV/precipRate'][:,:]
    lon=fh['FS/Longitude'][:,:]
    lat=fh['FS/Latitude'][:,:]
    hzero=fh['FS/VER/heightZeroDeg'][:,:]
    pType=fh['FS/CSF/typePrecip'][:,:]
    stormTop=fh['FS/PRE/heightStormTop'][:,:]
    pType=(pType/1e7).astype(np.int)
    bzd=fh['FS/VER/binZeroDeg'][:,:]
    zku=fh['FS/PRE/zFactorMeasured'][:,:,:,0]
    zka=fh['FS/PRE/zFactorMeasured'][:,:,:,0]
    bcf=fh['FS/PRE/binClutterFreeBottom'][:,:]
    dm=fh['FS/SLV/paramDSD'][:,:,:,1]
    zkuc=fh['FS/SLV/zFactorFinal'][:,:,:,0]
    piaF=fh['FS/SLV/piaFinal'][:,:,:]
    return sfcPrecip,hzero,pType,stormTop,bzd,zku,bcf,precipRate,dm,zkuc,piaF,zka,fh
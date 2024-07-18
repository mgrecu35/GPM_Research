

import netCDF4 as nc
import numpy as np
import matplotlib.pyplot as plt
def readDBabse(nclass,fname,n5):
    with nc.Dataset(fname) as fh:
        zkum=fh.variables["zkum%2.2i"%(nclass)][:,:-n5]
        zkam=fh.variables["zkam%2.2i"%(nclass)][:,:-n5]
        pRate=fh.variables["pRate%2.2i"%(nclass)][:,:-n5]
        n52=int(n5/2)
        pRateCMB=fh.variables["pRateCMB%2.2i"%(nclass)][:,:-n52]
        bzd_obs=fh.variables["bzd_obs%2.2i"%(nclass)][:]
    zkum[zkum<0]=0
    FL=168+2*nclass
    return zkum,zkam,pRate,pRateCMB,FL,bzd_obs



ncl=-5

import tqdm

def make_sampled_data(fnameL,npoints=25,ncl=-5):
    yL=[]
    ncL=[]
    irangeL=[]
    z_Ku_L=[]
    h_obs_L=[]
    precip_lcfb_L=[]
    bzd_out_L=[]
    pRateCMBL=[]
    zkumL=[]
    for ncl in tqdm.tqdm(range(-8,16)):
        zkumO,zkamO,pRateO,pRateCMBO,FL_O,bzdo=readDBabse(ncl,fnameL,6)
        zkumO=zkumO[:,:78]
        zkumL.append(zkumO[:,:78])
        pRateO=pRateO[:,:78]   
        pRateCMBO=pRateCMBO[:,:39] 
        h=np.arange(zkumO.shape[1])+2*ncl-70
        pRateCMBO_ResampledL=np.zeros((zkumO.shape[0],78))
        for i in range(zkumO.shape[0]):
            pRateCMBO_ResampledL[i,:]=np.interp(h[0:78],h[0::2][:],pRateCMBO[i,:][:])
            pRateCMBO_ResampledL[i,-1]=2*pRateCMBO_ResampledL[i,-2]-pRateCMBO_ResampledL[i,-3]
        pRateCMBL.append(pRateCMBO_ResampledL)
   
        for it in range(15):
            for i,zkum in enumerate(zkumL[ncl+8]):
                zkum=zkum[:-2]
                k=np.random.randint(0,17)
                nlev_end=k+1
                x_zku=zkum[-npoints-nlev_end:-nlev_end]
                x_zku[x_zku<0]=0
                x_h=h[-npoints-nlev_end:-nlev_end]*0.125
                pRate_lcfb=pRateCMBO[i,-nlev_end]
                if pRate_lcfb<1e-4:
                    if pRateCMBO[i,-nlev_end-1]>1e-4:
                        pRate_lcfb=pRateCMBO[i,-nlev_end-1]
                    else:
                        if pRateCMBO[i,-nlev_end+1]>1e-4:
                            pRate_lcfb=pRateCMBO[i,-nlev_end+1]
                if pRate_lcfb>0:
                    z_Ku_L.append(x_zku)
                    h_obs_L.append(x_h)
                    irangeL.append(k)
                    precip_lcfb_L.append(pRate_lcfb)            
                    ncL.append(ncl)
                    yL.append(pRateCMBL[ncl+8][i,-1])
                    bzd_out_L.append(bzdo[i])
    return z_Ku_L,h_obs_L,precip_lcfb_L,ncL,yL,bzd_out_L,irangeL

import xarray as xr
def save_to_netcdf(z_Ku_L,h_obs_L,precip_lcfb_L,ncL,yL,bzd_out_L,irangeL,fname_out):
    yL=np.array(yL)
    a=np.where(yL>0)
    zKuX=xr.DataArray(np.array(z_Ku_L),dims=["sample","height"])
    h_obsX=xr.DataArray(np.array(h_obs_L),dims=["sample","height"])
    precip_lcfbX=xr.DataArray(np.array(precip_lcfb_L),dims=["sample"])
    ncX=xr.DataArray(np.array(ncL),dims=["sample"])
    yX=xr.DataArray(np.array(yL),dims=["sample"])
    bzd_outX=xr.DataArray(np.array(bzd_out_L),dims=["sample"])
    rangebinX=xr.DataArray(np.array(irangeL),dims=["sample"])
    ds=xr.Dataset({"zKu":zKuX,"h_obs":h_obsX,"precip_lcfb":precip_lcfbX,"nc":ncX,"precip":yX,"bzd":bzd_outX,"rangebin":rangebinX})
    compLev=5
    comp=dict(zlib=True, complevel=compLev)

    encoding = {var: comp for var in ds.data_vars}
    ds.to_netcdf(fname_out,encoding=encoding)

#fnameL="strat_Profiles_Ocean_train.nc"
#z_Ku_L,h_obs_L,precip_lcfb_L,ncL,yL,bzd_out_L,irangeL=make_sampled_data(fnameL,npoints=50,ncl=-5)
#fname_out="strat_Profiles_Ocean_sampled_train.nc"
#save_to_netcdf(z_Ku_L,h_obs_L,precip_lcfb_L,ncL,yL,bzd_out_L,irangeL,fname_out)

#fnameL="strat_Profiles_Ocean_test.nc"
#z_Ku_L,h_obs_L,precip_lcfb_L,ncL,yL,bzd_out_L,irangeL=make_sampled_data(fnameL,npoints=50,ncl=-5)
#fname_out="strat_Profiles_Ocean_sampled_test.nc"
#save_to_netcdf(z_Ku_L,h_obs_L,precip_lcfb_L,ncL,yL,bzd_out_L,irangeL,fname_out)

fnameL="conv_Profiles_Ocean_train.nc"
z_Ku_L,h_obs_L,precip_lcfb_L,ncL,yL,bzd_out_L,irangeL=make_sampled_data(fnameL,npoints=50,ncl=-5)
fname_out="conv_Profiles_Ocean_sampled_train.nc"
save_to_netcdf(z_Ku_L,h_obs_L,precip_lcfb_L,ncL,yL,bzd_out_L,irangeL,fname_out)



fnameL="conv_Profiles_Ocean_test.nc"
z_Ku_L,h_obs_L,precip_lcfb_L,ncL,yL,bzd_out_L,irangeL=make_sampled_data(fnameL,npoints=50,ncl=-5)
fname_out="conv_Profiles_Ocean_sampled_test.nc"
save_to_netcdf(z_Ku_L,h_obs_L,precip_lcfb_L,ncL,yL,bzd_out_L,irangeL,fname_out)

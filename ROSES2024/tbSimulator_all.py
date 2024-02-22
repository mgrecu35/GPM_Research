#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import io_subs
import importlib
importlib.reload(io_subs)
import sim_tb
importlib.reload(sim_tb)
import pyresample
from pyresample import kd_tree, geometry
from pyresample.kd_tree import resample_nearest, resample_gauss, resample_custom
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from pyresample import image, geometry
import lkTables

lookupT=lkTables.scattTables()

import xarray as xr
wf = lambda r: 1 - r/20000.0
import glob

f2BCMB=glob.glob('/Volumes/T7 Shield/SUBSETS/Subsets/2B-CS-CONUS*')
f2BCMB=sorted(f2BCMB)
f1CGMI=glob.glob('/Volumes/T7 Shield/SUBSETS/Subsets/1C-*CONUS*')
f1CGMI=sorted(f1CGMI)
print(len(f2BCMB),len(f1CGMI))


# In[2]:


importlib.reload(sim_tb)
for f in f2BCMB[310:]:
    orbit_cmb=f.split('/')[-1].split('.')[-3]
    for f1 in f1CGMI:
        orbit_gmi=f1.split('/')[-1].split('.')[-3]
        if orbit_cmb==orbit_gmi:
            break
    if orbit_cmb==orbit_gmi:
        print(orbit_cmb,orbit_gmi)
        # Read the CMB and GMI data
        qv,press,envNodes,airTemp,skTemp,binNodes,pwc,sfcEmiss,dm,cldw,sfcBin,zCorrected,pType,lon,lat=io_subs.readCMB(f)
        pType=(pType/1e7).astype(int)
        lat_s1,lon_s1,tb_s1=io_subs.read1CGMI(f1)

# Define the input and output grids
        input_def = geometry.SwathDefinition(lons=lon_s1[:,:], lats=lat_s1[:,:])
        output_def = geometry.SwathDefinition(lons=lon, lats=lat)
# Resample the tb_s1 data to the CMB grid using gaussian resampling

        tb_s1_resampled = resample_custom(input_def, tb_s1[:,:,:], output_def, radius_of_influence=30000, neighbours=10, 
        weight_funcs=[wf for k in range(9)], fill_value=None)

        hFreqs=[1,1,2,2,3,4,4,5,5,6,6,7,7]

        tb_sim2,iwp2,rwp2,wvp=sim_tb.sim_tb(sfcEmiss,skTemp,envNodes,binNodes,pType,pwc,dm,airTemp,press,qv,sfcBin,lookupT)
        nt=pwc.shape[0]

        tbX=xr.DataArray(tb_sim2,coords={'scan':np.arange(nt),'ray':np.arange(49),'freq':hFreqs},dims=['scan','ray','freq'])
        tb_resampledX=xr.DataArray(tb_s1_resampled,coords={'scan':np.arange(nt),'ray':np.arange(49),'freq2':hFreqs[:9]},dims=['scan','ray','freq2'])
        iwpX=xr.DataArray(iwp2,coords={'scan':np.arange(nt),'ray':np.arange(49)},dims=['scan','ray'])
        rwpX=xr.DataArray(rwp2,coords={'scan':np.arange(nt),'ray':np.arange(49)},dims=['scan','ray'])
        sfcEmissX=xr.DataArray(sfcEmiss,coords={'scan':np.arange(nt),'ray':np.arange(49),'freq':hFreqs},dims=['scan','ray','freq'])
        skTempX=xr.DataArray(skTemp,coords={'scan':np.arange(nt),'ray':np.arange(49)},dims=['scan','ray'])
        sfcBinX=xr.DataArray(sfcBin[:,:,0],coords={'scan':np.arange(nt),'ray':np.arange(49)},dims=['scan','ray'])
        pTypeX=xr.DataArray(pType,coords={'scan':np.arange(nt),'ray':np.arange(49)},dims=['scan','ray'])
        wvpX=xr.DataArray(wvp,coords={'scan':np.arange(nt),'ray':np.arange(49)},dims=['scan','ray'])
        lonX=xr.DataArray(lon,coords={'scan':np.arange(nt),'ray':np.arange(49)},dims=['scan','ray'])
        latX=xr.DataArray(lat,coords={'scan':np.arange(nt),'ray':np.arange(49)},dims=['scan','ray'])
        ds=xr.Dataset({'tb':tbX,'tb_resampled':tb_resampledX,'iwp':iwpX,'rwp':rwpX,'sfcEmiss':sfcEmissX,'skTemp':skTempX,'sfcBin':sfcBinX,'pType':pTypeX,'wvp':wvpX,'lon':lonX,'lat':latX})
        complevel=5
        encoding = {var: {'zlib': True, 'complevel': complevel} for var in ds.data_vars}
        ds.to_netcdf('output/CMB_GMI_Tbs_CONUS_'+orbit_cmb+'.nc',encoding=encoding)
        #break


# In[11]:


tb_resampledX=xr.DataArray(tb_s1_resampled,coords={'scan':np.arange(nt),'ray':np.arange(49),'freq2':hFreqs[:9]},dims=['scan','ray','freq2'])


# In[ ]:


f1='2B-CS-CONUS.GPM.DPRGMI.CORRA2022.20180625-S050505-E051345.024557.V07A.HDF5'
f2='1C.GPM.GMI.XCAL2016-C.20180625-S041042-E054316.024557.V07A.HDF5'


ds=xr.Dataset({'tb':tbX,'iwp':iwpX,'rwp':rwpX,'sfcEmiss':sfcEmissX,'skTemp':skTempX,'sfcBin':sfcBinX,'pType':pTypeX,'wvp':wvpX,'tb_resampled':tb_resampledX},coords={'scan':np.arange(nt),'ray':np.arange(49),'freq':hFreqs})
complev=5
ds.to_netcdf('simulated_tb.nc',encoding={'tb':{'zlib':True,'complevel':complev},'iwp':{'zlib':True,'complevel':complev},'rwp':{'zlib':True,'complevel':complev},'sfcEmiss':{'zlib':True,'complevel':complev},'skTemp':{'zlib':True,'complevel':complev},'sfcBin':{'zlib':True,'complevel':complev},'pType':{'zlib':True,'complevel':complev},'wvp':{'zlib':True,'complevel':complev}})


# In[2]:


import rtlib
import importlib
importlib.reload(rtlib)
print(dir(rtlib))


# In[3]:


nt=qv.shape[0]
j=24
freqs=[10.65,18.7,23.8,37.0,89.0,166,183.3+3,183.3+7]
hFreqs=[1,1,2,2,3,4,4,5,5,6,6,7,7]
gFreqs=[1,1,2,2,3,4,4,5,5,6,6,7,8]
print(len(hFreqs))
print(sfcEmiss.shape)
import tqdm
importlib.reload(rtlib)
kext_gases=np.zeros((nt,49,78,8))
temp_layer=np.zeros((nt,49,79))
for i in tqdm.tqdm(range(0,nt)):
    for j in range(0,49):
        n5=envNodes[i,j,:]
        press1d=np.interp(range(10,88),n5,press[i,j,:])
        qv1d=np.interp(range(10,88),n5,qv[i,j,:])
        temp1d=np.interp(range(10,88),n5,airTemp[i,j,:])
        temp_layer1d=np.interp(np.arange(9.5,88.5,1),n5,airTemp[i,j,:])
        qv1d[qv1d<1e-3]=1e-3
        rho1d=press1d/(287.05*temp1d)*1e2
        ireturn=0
        #for k,qv1 in enumerate(qv1d):
        #    rhowv1d=rho1d*qv1d*1e-3
        #    for f in freqs:
        #        kext_wv=rtlib.gasabsr98(f,temp1d[k],rhowv1d,press1d[k]*1e2,ireturn)
        rhowv1d=rho1d*qv1d*1e-3
        kext = rtlib.get_wv_extinction(np.array(freqs),temp1d,rhowv1d,press1d*1e2)
        kext_gases[i,j,:]=kext.T
        temp_layer[i,j,:]=temp_layer1d


a=np.nonzero(pType>0)
print(len(a[0]))
from bisectm import bisectm
kext3d=np.zeros((nt,49,88,8))
salb3d=np.zeros((nt,49,88,8))
asym3d=np.zeros((nt,49,88,8))
iwp=np.zeros((nt,49))
rwp=np.zeros((nt,49))
for i,j in tqdm.tqdm(zip(a[0],a[1])):
    pwc1=pwc[i,j,:]
    a0=np.nonzero(pwc1>0)
    #print(a0[0])
    #print(binNodes[i,j,:])
    a1=np.nonzero(pwc1[0:binNodes[i,j,2]]==pwc1[0:binNodes[i,j,2]])
    if len(a1[0])>0:
        iwp[i,j]=np.sum(pwc1[a1])*0.25
    a2=np.nonzero(pwc1[binNodes[i,j,2]:binNodes[i,j,4]]==pwc1[binNodes[i,j,2]:binNodes[i,j,4]])
    if len(a2[0])>0:
        rwp[i,j]=np.sum(pwc1[a2])*0.25
    for k in a0[0]:
        dm1=dm[i,j,k]
        if k<=binNodes[i,j,3]:
            iphase=0
            if pType[i,j]!=2:
                ind=bisectm(lookupT.dms.data, 253, dm1)
                dns=np.log10(pwc[i,j,k]/lookupT.swc[ind])
                kext_snow=lookupT.kextS[ind,:]*10**dns
                salb_snow=lookupT.salbS[ind,:]
                asym_snow=lookupT.asymS[ind,:]
            else:
                ind=bisectm(lookupT.dmg.data, 272, dm1)
                dng=np.log10(pwc[i,j,k]/lookupT.gwc[ind])
                kext_snow=lookupT.kextG[ind,:]*10**dng
                salb_snow=lookupT.salbG[ind,:]
                asym_snow=lookupT.asymG[ind,:]
        
        if k>binNodes[i,j,1]:
            ind=bisectm(lookupT.dmr.data, 289, dm1)
            #print(ind,pwc[k])
            dnr=np.log10(pwc[i,j,k]/lookupT.rwc[ind])
            kext_rain=lookupT.kextR[ind,:]*10**dnr
            salb_rain=lookupT.salbR[ind,:]
            asym_rain=lookupT.asymR[ind,:]
        if k<=binNodes[i,j,1]:
            kext=kext_snow
            salb=salb_snow
            asym=asym_snow
        else:
            if k<binNodes[i,j,3]:
                f=(k-binNodes[i,j,1])/(binNodes[i,j,3]-binNodes[i,j,1])
                kext=kext_snow*(1-f)+kext_rain*f
                salb=(salb_snow*(1-f)*kext_snow+salb_rain*f*kext_rain)
                asym=(asym_snow*(1-f)*kext_snow*asym_snow+asym_rain*f*kext_rain*asym_rain)/salb
                salb=salb/kext
            else:
                kext=kext_rain
                salb=salb_rain
                asym=asym_rain
        kext3d[i,j,k,:]=kext
        salb3d[i,j,k,:]=salb
        asym3d[i,j,k,:]=asym     
        if k==a0[0][-1]:
            #print(k,binNodes[i,j,-1],sfcBin[i,j])
            kext3d[i,j,k+1:,:]=kext_rain
            salb3d[i,j,k+1:,:]=salb_rain
            asym3d[i,j,k+1:,:]=asym_rain
            if(pwc1[k]==pwc1[k]) and sfcBin[i,j,0]==sfcBin[i,j,0]:
                if k>=binNodes[i,j,2]:
                    rwp[i,j]+=(pwc1[k])*0.25*(sfcBin[i,j,0]-k)
                else:
                    iwp[i,j]+=(pwc1[k])*0.25*(sfcBin[i,j,0]-k)
#print(rwp.max())
#a=np.nonzero(pType>0)
#print(np.corrcoef(rwp[a],iwp[a]))

umu=np.cos(53/180*np.pi)
dz=0.25
print(kext_gases[:,:,::-1,:].shape)
print(kext3d[:,:,88:9:-1,:].shape)
binStart=87-sfcBin[:,:,0]
a=np.nonzero(pType>0)
fisot=2.7
#for i,j in zip(a[0],a[1]):
tb_sim=np.zeros((nt,49,len(hFreqs)))
for i in tqdm.tqdm(range(nt)):
    for j in range(49):
        kext1d_g=kext_gases[i,j,::-1,:]
        kext1d_h=kext3d[i,j,88:9:-1,:]
        asym1d_h=asym3d[i,j,88:9:-1,:]
        salb1d_h=salb3d[i,j,88:9:-1,:]
        temp1d=temp_layer[i,j,::-1]
        tbL=[]
        for ic in range(len(hFreqs)):
            kext1d=kext1d_g[binStart[i,j]:,gFreqs[ic]-1]+kext1d_h[binStart[i,j]:,hFreqs[ic]-1]
            asym1d=asym1d_h[binStart[i,j]:,hFreqs[ic]-1]*kext1d_h[binStart[i,j]:,hFreqs[ic]-1]/kext1d
            salb1d=salb1d_h[binStart[i,j]:,hFreqs[ic]-1]        
            temp1d1=temp1d[binStart[i,j]:]
            nz=kext1d.shape[0]
            height=np.arange(nz+1)*dz
        #break  
            tb=rtlib.radtran(umu,skTemp[i,j],temp1d1,height,kext1d,salb1d,asym1d,fisot,sfcEmiss[i,j,ic],sfcEmiss[i,j,ic])
            tbL.append(tb)
        tb_sim[i,j,:]=tbL
    #break

#tb = rtlib.get_tb(kext_gases[:,:,::-1,:],kext3d[:,:,88:9:-1,:],asym3d[:,:,88:9:-1,:],salb3d[:,:,88:9:-1,:],\
#                88-sfcBin[:,:,0],skTemp,temp_layer[:,:,::-1],sfcEmiss,hFreqs,gFreqs,pType,umu,dz)

#tb = get_tb(kext_gases,kext_hyd,asym_hyd,salb_hyd,sfcbin,sktemp,hfreqs,gfreqs,umu,dz


# In[92]:


import sim_tb
importlib.reload(sim_tb)

import xarray as xr
tb_sim2,iwp2,rwp2,wvp=sim_tb.sim_tb(sfcEmiss,skTemp,envNodes,binNodes,pType,pwc,dm,airTemp,press,qv,sfcBin,lookupT)


tbX=xr.DataArray(tb_sim2,coords={'scan':np.arange(nt),'ray':np.arange(49),'freq':hFreqs},dims=['scan','ray','freq'])
iwpX=xr.DataArray(iwp2,coords={'scan':np.arange(nt),'ray':np.arange(49)},dims=['scan','ray'])
rwpX=xr.DataArray(rwp2,coords={'scan':np.arange(nt),'ray':np.arange(49)},dims=['scan','ray'])
sfcEmissX=xr.DataArray(sfcEmiss,coords={'scan':np.arange(nt),'ray':np.arange(49),'freq':hFreqs},dims=['scan','ray','freq'])
skTempX=xr.DataArray(skTemp,coords={'scan':np.arange(nt),'ray':np.arange(49)},dims=['scan','ray'])
sfcBinX=xr.DataArray(sfcBin[:,:,0],coords={'scan':np.arange(nt),'ray':np.arange(49)},dims=['scan','ray'])
pTypeX=xr.DataArray(pType,coords={'scan':np.arange(nt),'ray':np.arange(49)},dims=['scan','ray'])
wvpX=xr.DataArray(wvp,coords={'scan':np.arange(nt),'ray':np.arange(49)},dims=['scan','ray'])
ds=xr.Dataset({'tb':tbX,'iwp':iwpX,'rwp':rwpX,'sfcEmiss':sfcEmissX,'skTemp':skTempX,'sfcBin':sfcBinX,'pType':pTypeX,'wvp':wvpX})
complev=5
ds.to_netcdf('simulated_tb.nc',encoding={'tb':{'zlib':True,'complevel':complev},'iwp':{'zlib':True,'complevel':complev},'rwp':{'zlib':True,'complevel':complev},'sfcEmiss':{'zlib':True,'complevel':complev},'skTemp':{'zlib':True,'complevel':complev},'sfcBin':{'zlib':True,'complevel':complev},'pType':{'zlib':True,'complevel':complev},'wvp':{'zlib':True,'complevel':complev}})


# In[95]:


plt.pcolormesh(lon,lat,wvp)


# In[68]:


ax=plt.subplot(1,2,1,projection=ccrs.PlateCarree())
ax.coastlines()
plt.pcolormesh(lon,lat,tb_sim[:,:,6],transform=ccrs.PlateCarree(),cmap='jet',vmin=160,vmax=300)
#plt.colorbar()
ax=plt.subplot(1,2,2,projection=ccrs.PlateCarree())
ax.coastlines()
plt.pcolormesh(lon,lat,tb_s1_resampled[:,:,6],transform=ccrs.PlateCarree(),cmap='jet',vmin=160,vmax=300)
#plt.colorbar()


# In[ ]:





# In[8]:


print(kext.shape)


# In[ ]:





# In[66]:


import netCDF4 as nc
with nc.Dataset('Data/2A-ENV.GPM.DPR.V9-20211125.20180625-S041042-E054316.024557.V07A.HDF5') as f:
    qv_env=f['FS/VERENV/waterVapor'][4500:5500,:,:]
    press_env=f['FS/VERENV/airPressure'][4500:5500,:,:]
    lon_env=f['FS/Longitude'][4500:5500,:]
    lat_env=f['FS/Latitude'][4500:5500,:]

print(qv.shape)


# In[67]:


ax=plt.subplot(1,1,1,projection=ccrs.PlateCarree())
plt.pcolormesh(lon_env,lat_env,press_env[:,:,-1],cmap='jet')

ax.coastlines()
plt.colorbar()


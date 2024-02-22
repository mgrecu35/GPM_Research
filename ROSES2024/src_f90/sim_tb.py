import tqdm
#importlib.reload(rtlib)
import numpy as np
import rtlib
from bisectm import bisectm

def sim_tb(sfcEmiss,skTemp,envNodes,binNodes,pType,pwc,dm,airTemp,press,qv,sfcBin,lookupT):
    nt=qv.shape[0]
    freqs=[10.65,18.7,23.8,37.0,89.0,166,183.3+3,183.3+7]
    hFreqs=[1,1,2,2,3,4,4,5,5,6,6,7,7]
    gFreqs=[1,1,2,2,3,4,4,5,5,6,6,7,8]
    #print(len(hFreqs))
    #print(sfcEmiss.shape)

    kext_gases=np.zeros((nt,49,78,8))
    temp_layer=np.zeros((nt,49,79))
    wvp=np.zeros((nt,49))
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
            n1=87-sfcBin[i,j,0]
            if n1>0:
                wvp[i,j]=np.sum(rhowv1d[:-n1])*0.25
            else:
                wvp[i,j]=np.sum(rhowv1d[:])*0.25
    a=np.nonzero(pType>0)
    print(len(a[0]))
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
                kext3d[i,j,k+1:,:]=kext
                salb3d[i,j,k+1:,:]=salb
                asym3d[i,j,k+1:,:]=asym
                if(pwc1[k]==pwc1[k]) and sfcBin[i,j,0]==sfcBin[i,j,0]:
                    if k>=binNodes[i,j,2]:
                        rwp[i,j]+=(pwc1[k])*0.25*(sfcBin[i,j,0]-k)
                    else:
                        iwp[i,j]+=(pwc1[k])*0.25*(sfcBin[i,j,0]-k)


    umu=np.cos(53/180*np.pi)
    dz=0.25
    print(kext_gases[:,:,::-1,:].shape)
    print(kext3d[:,:,88:9:-1,:].shape)
    binStart=87-sfcBin[:,:,0]
    a=np.nonzero(pType>0)
    fisot=2.7

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

    return tb_sim,iwp,rwp,wvp
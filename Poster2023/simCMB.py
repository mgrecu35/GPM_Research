import sys
#sys.path.append('/Users/mgrecu/myPythonPackages/')
import lkTables as lkT
lkTables=lkT.scattTables()
from bisectm import bisectm
import numpy as np
# zms = multiscatterf(kext,salb,g,ztrue,dr,noms,alt,theta,freq,nonorm,[nrange])
# tbout = radtran(umu,btemp,lyrtemp,lyrhgt,kext,salb,asym,fisot,emis,ebar,[nlyr])
#  absair,abswv = gasabsr98(f,tk,rhowv,pa,ireturn)
from numba import jit
print("simCMB is being imported")
def simulateZKa(binNodes,zCorrected,dm,dm_factor,lkTables,pType,kextKa,libSc,bsfc,msflag,dr): # calculates Ka-band reflectivity from attenuated corrected Ku-band reflectivity
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
    zKuSim=np.zeros((ns,nr,nbins),float)-99
    zWSim=np.zeros((ns,nr,nbins),float)-99
    zKa_true=np.zeros((ns,nr,nbins),float)-99
    zkaSfc=np.zeros((ns,nr),float)-99
    piaKa2d=np.zeros((ns,nr),float)
    piaKu2d=np.zeros((ns,nr),float)
    #print(pType.max(),ns,nr)
    #print(ns,nr,nbins)
    nfreqm=8
    kext1D=np.zeros((ns,nr,nbins,nfreqm),float)
    salb1D=np.zeros((ns,nr,nbins,nfreqm),float)
    asym1D=np.zeros((ns,nr,nbins,nfreqm),float)
    kexttot=np.zeros((ns,nr,nbins),float)
    salbtot=np.zeros((ns,nr,nbins),float)
    zms=np.zeros((ns,nr,nbins),float)
    pRate=np.zeros((ns,nr,nbins),float)-99
    f1=0.0
    for i in range(zCorrected.shape[0]):
        for j in range(zCorrected.shape[1]):
            if pType[i,j]>0:
                #print(binNodes[i,j,:])
                piaKa=0.0
                piaKu=0.0
                piaW=0.0
                attKa=0.0
                attKu=0.0
                fh1=np.interp(range(nbins),binNodes[i,j,:],[1,1.0+f1,1+f1,1+f1,1])
                #print(fh1.shape)
                #print(fh1)
                #stop
                for k in range(binNodes[i,j,0],min(binNodes[i,j,4],binNodes[i,j,2])+1):  # added range limit to avoid out of bounds
                    #print(dm[i,j,k])
                    if dm[i,j,k]>0:
                        #print(dm[i,j,k]*dm_factor)
                        ind=bisectm(lkTables.dms.data,253,dm_factor*dm[i,j,k])
                        dnw=(zCorrected[i,j,k]-lkTables.zKuS[ind])/10
                        attKa=lkTables.attKaS[ind]*10**dnw
                        attKu=lkTables.attKuS[ind]*10**dnw
                        attW=lkTables.attWS[ind]*10**dnw
                        kext1D[i,j,k,:]=lkTables.kextS[ind,:]*10**dnw
                        salb1D[i,j,k,:]=lkTables.salbS[ind,:]
                        asym1D[i,j,k,:]=lkTables.asymS[ind,:]
                        piaKa+=attKa*dr
                        piaKu+=attKu*dr
                        piaW+=attW*dr
                        zKaSim[i,j,k]=10*dnw+lkTables.zKaS[ind]-piaKa
                        zKuSim[i,j,k]=10*dnw+lkTables.zKuS[ind]-piaKu
                        zWSim[i,j,k]=10*dnw+lkTables.zWS[ind]-piaW
                        zKa_true[i,j,k]=10*dnw+lkTables.zKaS[ind]
                        pRate[i,j,k]=lkTables.snowRate[ind]*10**dnw
                        piaKa+=attKa*dr
                        piaKu+=attKu*dr
                        piaW+=attW*dr
                for k in range(binNodes[i,j,2]+1,min(binNodes[i,j,4],binNodes[i,j,3])):  # added range limit to avoid out of bounds
                    if dm[i,j,k]>0:
                        ind1=bisectm(lkTables.dms.data,253,dm_factor*dm[i,j,k])
                        ind2=bisectm(lkTables.dmr.data,289,dm[i,j,k]*fh1[k])
                        fsnow=((binNodes[i,j,3]-k)/(binNodes[i,j,3]-binNodes[i,j,1]))**0.5
                        zKuMPH=np.log10(fsnow*10**(0.1*lkTables.zKuS[ind1])+(1-fsnow)*10**(0.1*lkTables.zKuR[ind2]))*10
                        dnw=(zCorrected[i,j,k]-zKuMPH)/10
                        zKaMPH=np.log10(fsnow*10**(0.1*lkTables.zKaS[ind1])+(1-fsnow)*10**(0.1*lkTables.zKaR[ind2]))*10
                        zWMPH=np.log10(fsnow*10**(0.1*lkTables.zWS[ind1])+(1-fsnow)*10**(0.1*lkTables.zWR[ind2]))*10
                        attKa=fsnow*lkTables.attKaBB[ind1]*10**dnw+(1-fsnow)*lkTables.attKaR[ind2]*10**dnw
                        attKu=fsnow*lkTables.attKuBB[ind1]*10**dnw+(1-fsnow)*lkTables.attKuR[ind2]*10**dnw
                        attW=fsnow*lkTables.attWS[ind1]*10**dnw+(1-fsnow)*lkTables.attWR[ind2]*10**dnw
                        kext1D[i,j,k,:]=fsnow*lkTables.kextS[ind1,:]*10**dnw+(1-fsnow)*lkTables.kextR[ind2,:]*10**dnw
                        salb1D[i,j,k,:]=fsnow*lkTables.salbS[ind1,:]+(1-fsnow)*lkTables.salbR[ind2,:]
                        asym1D[i,j,k,:]=fsnow*lkTables.asymS[ind1,:]+(1-fsnow)*lkTables.asymR[ind2,:]
                        piaKa+=attKa*dr
                        piaKu+=attKu*dr
                        piaW+=attW*dr
                        zKaSim[i,j,k]=10*dnw+zKaMPH-piaKa
                        zKuSim[i,j,k]=10*dnw+zKuMPH-piaKu
                        zWSim[i,j,k]=10*dnw+zWMPH-piaW
                        zKa_true[i,j,k]=10*dnw+zKaMPH
                        pRate[i,j,k]=fsnow*lkTables.snowRate[ind1]*10**dnw+(1-fsnow)*lkTables.rainRate[ind2]*10**dnw
                        piaKa+=attKa*dr
                        piaKu+=attKu*dr
                        piaW+=attW*dr
                for k in range(binNodes[i,j,3],binNodes[i,j,4]):
                    if dm[i,j,k]>0:
                        ind=bisectm(lkTables.dmr.data,289,dm[i,j,k])
                        dnw=(zCorrected[i,j,k]-lkTables.zKuR[ind])/10
                        attKa=lkTables.attKaR[ind]*10**dnw
                        attKu=lkTables.attKuR[ind]*10**dnw
                        attW=lkTables.attWR[ind]*10**dnw
                        piaKa+=attKa*dr
                        piaKu+=attKu*dr
                        piaW+=attW*dr
                        zKaSim[i,j,k]=10*dnw+lkTables.zKaR[ind]-piaKa
                        zKuSim[i,j,k]=10*dnw+lkTables.zKuR[ind]-piaKu
                        zWSim[i,j,k]=10*dnw+lkTables.zWR[ind]-piaW
                        zKa_true[i,j,k]=10*dnw+lkTables.zKaR[ind]
                        piaKa+=attKa*dr
                        piaKu+=attKu*dr
                        piaW+=attW*dr
                        kext1D[i,j,k,:]=lkTables.kextR[ind,:]*10**dnw
                        salb1D[i,j,k,:]=lkTables.salbR[ind,:]
                        asym1D[i,j,k,:]=lkTables.asymR[ind,:]
                        pRate[i,j,k]=lkTables.rainRate[ind]*10**dnw
                piaKa2d[i,j]=piaKa+(bsfc[i,j,1]-binNodes[i,j,4])*dr*attKa
                piaKu2d[i,j]=piaKu+(bsfc[i,j,1]-binNodes[i,j,4])*dr*attKu
                zkaSfc[i,j]=zKaSim[i,j,k]
                kexttot[i,j,:]=kext1D[i,j,:,3]+kextKa[:] 
                salbtot[i,j,:]=kext1D[i,j,:,3]*salb1D[i,j,:,3]
                n1=binNodes[i,j,4]
                noms=0
                alt=400
                freqKa=35.5
                nonorm=1
                theta=0.5
                if msflag==1:
                    zms[i,j,:n1] = libSc.multiscatterf(kexttot[i,j,:n1],salbtot[i,j,:n1],asym1D[i,j,:n1,3],\
                            zKa_true[i,j,:n1],dr,noms,alt,theta,freqKa,nonorm)
                #zms = libSc.multiscatterf(kexttot[i,j,:],salbtot[i,j,:],asymtot[i,j,:],\
                #                          ztrue,dr,noms,alt,theta,freq,nonorm,[nrange])
    return zKaSim,zkaSfc,kexttot,salbtot,asym1D,\
        zKa_true,zms,piaKa2d,piaKu2d,zKuSim,zWSim,pRate

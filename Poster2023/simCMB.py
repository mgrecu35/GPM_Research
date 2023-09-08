import sys
#sys.path.append('/Users/mgrecu/myPythonPackages/')
import lkTables as lkT
lkTables=lkT.scattTables()
from bisectm import bisectm
import numpy as np
# zms = multiscatterf(kext,salb,g,ztrue,dr,noms,alt,theta,freq,nonorm,[nrange])
# tbout = radtran(umu,btemp,lyrtemp,lyrhgt,kext,salb,asym,fisot,emis,ebar,[nlyr])
#  absair,abswv = gasabsr98(f,tk,rhowv,pa,ireturn)

def simulateZKa(binNodes,zCorrected,dm,dm_factor,lkTables,pType,kextKa,libSc): # calculates Ka-band reflectivity from attenuated corrected Ku-band reflectivity
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
    zKa_true=np.zeros((ns,nr,nbins),float)-99
    zkaSfc=np.zeros((ns,nr),float)-99
    #print(pType.max(),ns,nr)
    dr=0.25
    #print(ns,nr,nbins)
    nfreqm=8
    kext1D=np.zeros((ns,nr,nbins,nfreqm),float)
    salb1D=np.zeros((ns,nr,nbins,nfreqm),float)
    asym1D=np.zeros((ns,nr,nbins,nfreqm),float)
    kexttot=np.zeros((ns,nr,nbins),float)
    salbtot=np.zeros((ns,nr,nbins),float)
    zms=np.zeros((ns,nr,nbins),float)
    for i in range(zCorrected.shape[0]):
        for j in range(zCorrected.shape[1]):
            if pType[i,j]>0:
                #print(binNodes[i,j,:])
                piaKa=0.0
                for k in range(binNodes[i,j,0],min(binNodes[i,j,4],binNodes[i,j,1])+1):  # added range limit to avoid out of bounds
                    #print(dm[i,j,k])
                    if dm[i,j,k]>0:
                        #print(dm[i,j,k]*dm_factor)
                        ind=bisectm(lkTables.dms.data,253,dm_factor*dm[i,j,k])
                        dnw=(zCorrected[i,j,k]-lkTables.zKuS[ind])/10
                        attKa=lkTables.attKaS[ind]*10**dnw
                        kext1D[i,j,k,:]=lkTables.kextS[ind,:]*10**dnw
                        salb1D[i,j,k,:]=lkTables.salbS[ind,:]
                        asym1D[i,j,k,:]=lkTables.asymS[ind,:]
                        piaKa+=attKa*dr
                        zKaSim[i,j,k]=10*dnw+lkTables.zKaS[ind]-piaKa
                        zKa_true[i,j,k]=10*dnw+lkTables.zKaS[ind]
                        piaKa+=attKa*dr
                for k in range(binNodes[i,j,1]+1,min(binNodes[i,j,4],binNodes[i,j,3])):  # added range limit to avoid out of bounds
                    if dm[i,j,k]>0:
                        ind1=bisectm(lkTables.dms.data,253,dm_factor*dm[i,j,k])
                        ind2=bisectm(lkTables.dmr.data,289,dm[i,j,k])
                        fsnow=((binNodes[i,j,3]-k)/(binNodes[i,j,3]-binNodes[i,j,1]))**0.5
                        zKuMPH=np.log10(fsnow*10**(0.1*lkTables.zKuS[ind1])+(1-fsnow)*10**(0.1*lkTables.zKuR[ind2]))*10
                        dnw=(zCorrected[i,j,k]-zKuMPH)/10
                        zKaMPH=np.log10(fsnow*10**(0.1*lkTables.zKaS[ind1])+(1-fsnow)*10**(0.1*lkTables.zKaR[ind2]))*10
                        attKa=fsnow*lkTables.attKaBB[ind1]*10**dnw+(1-fsnow)*lkTables.attKaR[ind2]*10**dnw
                        kext1D[i,j,k,:]=fsnow*lkTables.kextS[ind1,:]*10**dnw+(1-fsnow)*lkTables.kextR[ind2,:]*10**dnw
                        salb1D[i,j,k,:]=fsnow*lkTables.salbS[ind1,:]+(1-fsnow)*lkTables.salbR[ind2,:]
                        asym1D[i,j,k,:]=fsnow*lkTables.asymS[ind1,:]+(1-fsnow)*lkTables.asymR[ind2,:]
                        piaKa+=attKa*dr
                        zKaSim[i,j,k]=10*dnw+zKaMPH-piaKa
                        zKa_true[i,j,k]=10*dnw+zKaMPH
                        piaKa+=attKa*dr
                for k in range(binNodes[i,j,3],binNodes[i,j,4]):
                    if dm[i,j,k]>0:
                        ind=bisectm(lkTables.dmr.data,289,dm[i,j,k])
                        dnw=(zCorrected[i,j,k]-lkTables.zKuR[ind])/10
                        attKa=lkTables.attKaR[ind]*10**dnw
                        piaKa+=attKa*dr
                        zKaSim[i,j,k]=10*dnw+lkTables.zKaR[ind]-piaKa
                        zKa_true[i,j,k]=10*dnw+lkTables.zKaR[ind]
                        piaKa+=attKa*dr
                        kext1D[i,j,k,:]=lkTables.kextR[ind,:]*10**dnw
                        salb1D[i,j,k,:]=lkTables.salbR[ind,:]
                        asym1D[i,j,k,:]=lkTables.asymR[ind,:]
                zkaSfc[i,j]=zKaSim[i,j,k]
                kexttot[i,j,:]=kext1D[i,j,:,3]+kextKa[:] 
                salbtot[i,j,:]=kext1D[i,j,:,3]*salb1D[i,j,:,3]
                n1=binNodes[i,j,4]
                noms=0
                alt=400
                freqKa=35.5
                nonorm=1
                theta=0.5
                #print(kexttot[i,j,n1-10:n1])
                #print(salbtot[i,j,n1-10:n1])
                #print(asym1D[i,j,n1-10:n1,3])
                zms[i,j,:n1] = libSc.multiscatterf(kexttot[i,j,:n1],salbtot[i,j,:n1],asym1D[i,j,:n1,3],\
                        zKa_true[i,j,:n1],dr,noms,alt,theta,freqKa,nonorm)
                #zms = libSc.multiscatterf(kexttot[i,j,:],salbtot[i,j,:],asymtot[i,j,:],\
                #                          ztrue,dr,noms,alt,theta,freq,nonorm,[nrange])
    return zKaSim,zkaSfc,kexttot,salbtot,asym1D,zKa_true,zms  # added unatted Ka-band reflectivity at surface 

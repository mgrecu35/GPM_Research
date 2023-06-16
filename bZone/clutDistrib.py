import glob
import netCDF4 as nc
fs=glob.glob("bZoneProfs/*clutter*.nc")
import matplotlib.pyplot as plt
clutDepthL=[]
jL=[]
bzdL=[]
bcfL=[]
import numpy as np
for f in sorted(fs):
    fh=nc.Dataset(f)
    bcf=fh["bcFree"][:]
    bsfc=fh["bsfc"][:]
    jRay=fh["jRay"][:]
    pType=fh["pType"][:]
    a=np.nonzero(pType[:,0]==1)
    bzdL.extend(fh["bzd"][a[0]])

    jL.extend(abs(jRay[a[0]]-24))
    clutDepthL.extend(bsfc[a[0],0]-bcf[a])
    bcfL.extend(bcf[a])

jL=np.array(jL)
clutDepthL=np.array(clutDepthL)
bzdL=np.array(bzdL)
bcfL=np.array(bcfL)
a=np.nonzero(clutDepthL==clutDepthL)
clutDepthL=clutDepthL[a]
jL=jL[a]
bzdL=bzdL[a]
bcfL=bcfL[a]
for j in range(13,24):
    a=np.nonzero((jL==j) & (bzdL>165))
    b=np.nonzero((bzdL[a]-bcfL[a])>0)
    #print(
    a2=np.nonzero((jL==j) & (bzdL<165))
    b2=np.nonzero((bzdL[a2]-bcfL[a2])>0)
    print(len(a[0]),clutDepthL[a].mean(),len(b[0])/len(a[0]),\
          len(a2[0]),clutDepthL[a2].mean(),len(b2[0])/len(a2[0]))

plt.hist(clutDepthL)

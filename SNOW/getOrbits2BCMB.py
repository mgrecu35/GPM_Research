import pickle
#orbs=pickle.load(open("orogOrbitsIPHEX.pklz","rb"))
path="https://pmm-gv.gsfc.nasa.gov/pub/gpm-validation/data/gpmgv/orbit_subset/GPM/DPRGMI/2BDPRGMI/V07A/CONUS/%4.4i/%2.2i/%2.2i"


import urllib.request


from bs4 import BeautifulSoup
import os
import datetime
st_date=datetime.datetime(2022,1,1)
for iday in range(0,59):
    cdate=st_date+datetime.timedelta(days=iday)
    yy=cdate.year
    mm=cdate.month
    dd=cdate.day
    #if yy==2018:
    #    break
    loc=path%(yy,mm,dd)
    print(path%(yy,mm,dd))
    try:
        req = urllib.request.Request(path%(yy,mm,dd))
        response = urllib.request.urlopen(req)
        the_page = response.read()
        soup = BeautifulSoup(the_page, 'html.parser')
        #orb=fname[73:90]
        for link in soup.find_all('a'):
            item=link.get('href')
            if '2B-CS' in item:
                fname=loc+'/'+item
                print(fname)
                os.system('curl -O '+fname)
                os.system('mv 2B-CS*HDF5 ../SNOW/data')
                #stop
    except:
        pass
    #break

import matplotlib.pyplot as plt
import numpy as np
import sympy
from astropy.io import fits
from matplotlib.pyplot import MultipleLocator
import scipy.signal as signal
from astropy.io import fits
from astropy.table import Table, vstack
from uncertainties import ufloat
from uncertainties.umath import * 
from uncertainties.umath import log
import math

sightlines=np.load('/home/bwang/data/Everest/Healpix/main/dark/100-sightlines.npy',allow_pickle = True,encoding='latin1')
tbl=Table.read('/home/bwang/data/Everest/Healpix/main/dark/100-catalog.fits')
outputpath='/home/bwang/vi_DESIdata/test_100'

keys=tbl['INDEX'][(tbl['S/N']>1)&(tbl['SPECTYPE']=='QSO')]
sightlines_qso=sightlines[keys]

linelist=['Ly'+'$\\alpha$ 1216','NV 1240','CIV 1549','HeII 1640','CIII 1908','MgII 2799']
linerest=[1215.67,1240,1549,1640,1908,2799]
lineheight=[0.9,0.7,0.9,0.7,0.9,0.7]

for i in range(0,100):
    flux=sightlines_qso[i].flux[0:4880]
    wave=10**(sightlines_qso[i].loglam)[0:4880]
    error=sightlines_qso[i].error[0:4880]
    z=sightlines_qso[i].z_qso
    spectraid=sightlines_qso[i].id
    specclass=sightlines_qso[i].spectype
    objclass=sightlines_qso[i].objtype
    w1=sightlines_qso[i].w1
    w1error=1/np.sqrt(sightlines_qso[i].w1_ivar)
    w2=sightlines_qso[i].w2
    w2error=1/np.sqrt(sightlines_qso[i].w2_ivar)
    #smooth_flux=signal.medfilt(flux,3)-np.max(flux)
    smooth_flux=signal.medfilt(flux,3)-np.max(flux)+np.percentile(flux, 95)
    
    if (w1 > 0.0) & (w1error>0.0):
        w1value=ufloat(w1,w1error)
        aa=-2.5*log10(w1value)
        w1abs=22.5+aa
    else : 
        w1abs = 99

    
    
    if (w2 > 0.0) & (w2error>0.0):
        w2value=ufloat(w2,w2error)
        bb=-2.5*log10(w2value)
        w2abs=22.5+bb
    else:
        w2abs = 99
    
    
    plt.figure(figsize=(12,8))
    plt.plot(wave,flux,c='black',drawstyle='steps-mid',label='flux z=%0.2f'%(z))
    plt.plot(wave,smooth_flux,c='dodgerblue',drawstyle='steps-mid',label='smoothed flux')
    plt.plot(wave,error,c='red',ls='--',label='error')
    
    
    #plt.title('%s(%s)-ID%s, z=%0.2f, W1=%.2f W2=%.2f'%(specclass,objclass,spectraid,z,w1,w2),fontsize=20)
    plt.title('%s(%s)-ID%s-%s, W1=%s W2=%s'%(specclass,objclass,spectraid,i,w1abs,w2abs),fontsize=20)

    for j in range(0,len(linerest)):
        lineloc=(1+z)*(linerest[j])
        if (lineloc > 3600)&(lineloc < 7500):
            plt.axvline(x=lineloc,ls="--",c="blue",linewidth=2)
            plt.text(lineloc+10,np.max(flux)*lineheight[j],linelist[j],fontsize=16,color='green')

    plt.axhline(y=0,ls="--",c="grey",linewidth=2)
    plt.ylabel('Flux[$\mathregular{(10)^-17}$erg/$\mathregular{(cm)^2}$/s/'+'$\AA$'+']',fontsize=20)
    plt.xlabel('Wavelength'+'['+'$\\rm\\AA$'+']',fontsize=20)
    plt.xlim(3600,7500)
    plt.tick_params(labelsize=20)
    plt.legend(fontsize=16)
    plt.savefig(outputpath+'/'+'100-%s.png'%(spectraid),dpi=300)

        #plt.clf()
        #plt.close()
    plt.close()
    
    
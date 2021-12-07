import matplotlib.pyplot as plt
import numpy as np
import sympy
from astropy.io import fits
from matplotlib.pyplot import MultipleLocator
import xlwt
import os
import scipy.signal as signal
from astropy.io import fits
from astropy.table import Table, vstack
from uncertainties import ufloat
from uncertainties.umath import * 
from uncertainties.umath import log
import math
from astropy.modeling import models,fitting
from astropy import modeling
import scipy.integrate as integrate

#linelist_new=['DLA 1216','NV 1238','NV 1242','CII 1334','CII 1335','SiII 1526','CIV 1548','CIV 1550','FeII 1608','HeII 1640','Al 1670','CIII 1908','MgII 2795','MgII 2802']
#linerest_new=[1215.67,1238.8,1242.8,1334.5,1335.7,1526,1548.2,1550.8,1608,1640,1670,1908,2795.5,2802.7]

linelist_new=['CII 1334','SiII 1526','CIV 1549','FeII 1608','HeII 1640','Al 1670','CIII 1908','FeII 2344','FeII 2374','FeII 2382',
'FeII 2586','FeII 2600','MgII 2795','MgII 2802','MgII 2852']
linerest_new=[1334,1526,1549,1608,1640,1670,1908,2344,2374,2382,2586,2600,2795.5,2802.7,2852]

from gfit import fit_cont, fit_gauss, get_data, get_ew


def estimate(wave,flux,error,z_dla,z_qso,linerest):


    xdata, ydata, edata = get_data(wave,flux,error,z_dla,linerest)


    fitted_line=fit_cont(xdata,ydata)

    ydata_n=ydata/(fitted_line(xdata))
    edata_n=edata/(fitted_line(xdata))

    '''
    g,fg=fit_gauss(xdata,ydata_n)

    def func(x):
        return -fg(x)

    area = integrate.quad(func, xdata[0], xdata[-1])
    '''
    ew_ob , error_ob = get_ew(xdata,ydata_n,edata_n)

    ew_rest = ew_ob / (1 + z_qso)
    error_rest = error_ob / (1 + z_qso)

    if ew_rest / error_rest < 3:
        ew_rest = 0

    return ew_rest, error_rest


if __name__ == '__main__':

    qso=fits.open('/Users/samwang/Desktop/DESI_project/DLA_catalog/vi_catalog/Everest_qsocatalog_sightlines_vi.fits')
    dla_cnn=fits.open('/Users/samwang/Desktop/DESI_project/DLA_catalog/vi_catalog/Everest_dlacatalog_714model_20to23_vi.fits')

    qso_tab=Table.read('/Users/samwang/Desktop/DESI_project/DLA_catalog/vi_catalog/Everest_qsocatalog_sightlines_vi.fits')
    cnn_tab=Table.read('/Users/samwang/Desktop/DESI_project/DLA_catalog/vi_catalog/Everest_dlacatalog_714model_20to23_vi.fits')

    sightlines=np.load('/Users/samwang/Desktop/DESI_project/DLA_catalog/Everest_qsosightlines.npy',allow_pickle = True,encoding='latin1')

    idlist=[]
    for qsoid in qso_tab['TARGETID']:
        idlist.append(qsoid)


    dlalist=[]
    for j in range(0,len(cnn_tab)):
        if cnn_tab[j]['DLA'] == 1:
            dlalist.append(j)
    print('%s DLAs'%(len(dlalist)))

    dataset={}

    #metalt = []
    #errort = []

    for i in dlalist:
        z_dla=cnn_tab[i]['Z']
        qsoid=cnn_tab[i]['TARGETID']
        nhi=cnn_tab[i]['NHI']
    
        k=idlist.index(qsoid)
        flux=sightlines[k].flux
        wave=10**(sightlines[k].loglam)
        error=sightlines[k].error
    
        z_qso=sightlines[k].z_qso
        spectraid=sightlines[k].id

        mline=[]
        eline=[]

        for pp in range(0,len(linerest_new)):
            lineloc=(1+z_dla)*(linerest_new[pp])
            if (lineloc > 1216*(1+z_qso))&(lineloc < 9800):
                ew_rest, error_rest = estimate(wave,flux,error,z_dla,z_qso,linerest_new[pp])
                mline.append(ew_rest)
                eline.append(error_rest)
            else : 
                mline.append(-1)
                eline.append(-1)

        #metalt.append(mline)
        #errort.append(eline)

        dataset[qsoid] = {'DLA':1,'NHI':nhi,'z_dla':z_dla,'CII_1334':mline[0],'CII_error':eline[0],
        'SiII_1526':mline[1],'Si_error':eline[1],'CIV_1549':mline[2],'CIV_error':eline[2],
        'FeII_1608':mline[3],'Fe1608_error':eline[3],'HeII_1640':mline[4],'He_error':eline[4],'Al_1670':mline[5],'Al_error':eline[5],
        'CIII_1908':mline[6],'CIII_error':eline[6],'FeII_2344':mline[7],'Fe2344_error':eline[7],'FeII_2374':mline[8],'Fe2374_error':eline[8],
        'FeII_2382':mline[9],'Fe2382_error':eline[9],'FeII_2586':mline[10],'Fe2586_error':eline[10],
        'FeII_2600':mline[11],'Fe2600_error':eline[11],'MgII_2795':mline[12],'Mg2795_error':eline[12],
        'MgII_2802':mline[13],'Mg2802_error':eline[13],'MgII_2852':mline[14],'Mg2852_error':eline[14]}

        print(i)


    np.save('/Users/samwang/Desktop/DESI_project/DLA_catalog/test_EW/code/DLA_EW_catalog_v2.npy',dataset)









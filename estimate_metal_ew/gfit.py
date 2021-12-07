import numpy as np
from astropy.modeling import models,fitting
from astropy import modeling

def fit_cont(xdata,ydata):
    an_amplitude = ydata.min()
    an_mean = xdata[ydata.argmin()]
    
    #print(an_mean)
    #print(np.sum((xdata - an_mean)**2) / (len(xdata) - 1))
    
    an_stddev = np.sqrt(np.sum((xdata - an_mean)**2) / (len(xdata) - 1))
    
    c1 = 5
    c2 = 15
    
    for i in range(0,len(xdata)-1):
        if (xdata[i] <= (an_mean - an_stddev/2)) & (xdata[i+1] >= (an_mean - an_stddev/2)):
            c1 = i
            #print('c1 is %s'%c1)
        if (xdata[i] <= (an_mean + an_stddev/2)) & (xdata[i+1] >= (an_mean + an_stddev/2)):
            c2 = i
            #print('c2 is %s'%c2)
    
    xdata_1=xdata[0:c1]
    xdata_2=xdata[c2:-1]
    xdata_c=np.hstack((xdata_1,xdata_2))
    ydata_1=ydata[0:c1]
    ydata_2=ydata[c2:-1]
    ydata_c=np.hstack((ydata_1,ydata_2))
    
    fit = fitting.LinearLSQFitter()
    line_init = models.Linear1D()
    fitted_line = fit(line_init, xdata_c, ydata_c)
    
    return fitted_line


def fit_gauss(xdata,ydata):
    an_amplitude = ydata_n.min()
    an_mean = xdata[ydata_n.argmin()]
    an_stddev = np.sqrt(np.sum((xdata - an_mean)**2) / (len(xdata) - 1))

    an_disp = 1

    an_slope=(ydata_n[-1]-ydata_n[0])/(xdata.max()-xdata.min())
    g_init = (
         models.Const1D(an_disp)+
         models.Gaussian1D(amplitude=(an_amplitude - an_disp), mean=an_mean, stddev=an_stddev)
     )
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, xdata, ydata_n)
    fg=g[1]
    return g,fg


def get_data(wave,flux,error,z_dla,linerest):
	lineob=(linerest*(1+z_dla))
	for i in range(0,len(wave)-1):
		if (wave[i+1] > lineob) & (wave[i] < lineob):
			d=i+1
	xdata=wave[d-10:d+10]
	ydata=flux[d-10:d+10]
	edata=error[d-10:d+10]

	return xdata,ydata,edata

def get_ew(xdata,ydata_n,edata_n):
    area=0
    error=[]
    for ii in range(0,len(xdata)-1):
        s = (xdata[ii] - xdata[ii+1])*(ydata_n[ii] - 1)
        area = area + s
        error.append(((xdata[ii] - xdata[ii+1])**2)*(edata_n[ii]**2))
    error_ob = np.sqrt(np.sum(error))
    return area,error_ob







import time
import numpy as np
import sys
import shutil
import string
import os

from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord

data = fits.getdata('total_offset.fits',1)
ra_all, dec_all = data['target_ra'],data['target_dec']
magz = data['mag_r']
magj = data['mag_r']
flag = data['survey']
priority = data['priority']
#pos_angle = data['pos_angle']
#delta_ra = data['dra_offset']
#delta_dec = data['ddec_offset']

#f = open('lris_with_offset_pa.txt',"w")
f = open('total_starlist.txt',"w")

def degtorad(d):
    """
      Convert degrees into radians.
      Parameters
      ----------
      d : float or array
          Angle in degrees.
      Returns
      -------
      Angle : float or array
          The angle converted into radians.
    """
    return (d / 180.0) * np.pi


def radtodeg(r):
    """
      Convert radians into degrees.
      Parameters
      ----------
      d : float or array
          Angle in radians.
      Returns
      -------
      Angle : float or array
          The angle converted into degrees.
    """
    return (r / np.pi) * 180.0

def hms2deg(ra, dec):
    rah = np.double(ra[0:2])
    ram = np.double(ra[2:4])
    ras = np.double(ra[4:10])
    radeg = 15.0 * (rah + ram / 60 + ras / 3600)
    sig = dec[0]
    decd = np.double(dec[1:3])
    decm = np.double(dec[3:5])
    decs = np.double(dec[5:11])
    if sig == '+':
        decdeg = decd + decm / 60 + decs / 3600
    else:
        decdeg = 0 - (decd + decm / 60 + decs / 3600)
        decdeg = round(decdeg, 8)
    return radeg, decdeg


def deg2hms(ra, dec):
    rah = np.floor(ra / 15.0)
    ram = np.floor((ra / 15.0 - rah) * 60.0)
    ras = ((ra / 15.0 - rah) * 60.0 - ram) * 60.0
    if dec < 0:
        dec = np.abs(dec)
        decd = 0 - np.floor(dec)
        decm = np.floor((dec + decd) * 60.0)
        decs = ((dec + decd) * 60.0 - decm) * 60.0
    else:
        decd = np.floor(dec)
        decm = np.floor((dec - decd) * 60.0)
        decs = ((dec - decd) * 60.0 - decm) * 60.0
    ra = [rah, ram, ras]
    dec = [decd, decm, decs]
    return ra, dec


def hms2name(ra, dec):
    rahms, decdms = deg2hms(ra, dec)
    rahd = int(rahms[0])
    ramd = int(rahms[1])
    rasd = round(rahms[2], 3)

    rastr = ('{rahs:02d}''{rams:02d}''{rass:6.3f}')
    ras = rastr.format(rahs=rahd, rams=ramd, rass=rasd)
    rass = ras.replace(' ', '0')

    dechd = int(np.abs(decdms[0]))
    decmd = int(decdms[1])
    decsd = round(decdms[2], 3)

    decstr = ('{dechs:02d}''{decms:02d}''{decss:6.3f}')
    decs = decstr.format(dechs=dechd, decms=decmd, decss=decsd)
    decss = decs.replace(' ', '0')
    if dec < 0:
        decss = '-' + decss
    else:
        decss = '+' + decss

    name = 'J' + rass + decss
    namebrief = name[0:5] + name[11:-6]
    return name, namebrief

def hms2name_paper(ra, dec):
    rahms, decdms = deg2hms(ra, dec)
    rahd = int(rahms[0])
    ramd = int(rahms[1])
    rasd = round(rahms[2], 2)

    rastr = ('{rahs:02d}''{rams:02d}''{rass:5.2f}')
    ras = rastr.format(rahs=rahd, rams=ramd, rass=rasd)
    rass = ras.replace(' ', '0')

    dechd = int(np.abs(decdms[0]))
    decmd = int(decdms[1])
    decsd = round(decdms[2], 1)

    decstr = ('{dechs:02d}''{decms:02d}''{decss:4.1f}')
    decs = decstr.format(dechs=dechd, decms=decmd, decss=decsd)
    decss = decs.replace(' ', '0')
    if dec < 0:
        decss = '$-$' + decss
    else:
        decss = '$+$' + decss

    name = 'J' + rass + decss
    return name

def offset(ra1, dec1, ra2, dec2, center=False):
    """
    Compute the offset from object1 to object2 in ra and dec
    in units of arcsecond.
    Parameters
    ----------
    ra1 : float
    Right ascension of first object [deg].
    dec1 : float
    Declination of first object [deg].
    ra2 : float
    Right ascension of second object [deg].
    dec2 : float
    Declination of second object [deg].
    center: if true will also return the center position of these two objects
    Returns
    -------
    delta_ra,delta_dec : float
    center : float
    """

    # Convert into rad
    #rarad1 = degtorad(ra1)
    #rarad2 = degtorad(ra2)
    #dcrad1 = degtorad(dec1)
    #dcrad2 = degtorad(dec2)

    radif = (ra2 - ra1) * np.cos(degtorad((dec1 + dec2) / 2)) * 3600
    dcdif = (dec2 - dec1) * 3600

    if center:
        return np.round(radif, 2), np.round(dcdif, 2), (ra1 + ra2) / 2, (dec1 + dec2) / 2
    else:
        return np.round(radif, 2), np.round(dcdif, 2)

def positionAngle(ra1, dec1, ra2, dec2, plot=False, positive=False):
    """
      Compute the position angle.
      The position angle is measured from the first position
      from North through East. If the `positive` flag is set
      True (default) the result will be given as an angle from
      0 to 360 degrees.
      Parameters
      ----------
      ra1 : float
          Right ascension of first object [deg].
      dec1 : float
          Declination of first object [deg].
      ra2 : float
          Right ascension of second object [deg].
      dec2 : float
          Declination of second object [deg].
      Returns
      -------
      Position angle : float
      The position angle in degrees. the output will be
      given as an angle between 0 and 360 degrees.
    """

    # Convert into rad
    rarad1 = degtorad(ra1)
    rarad2 = degtorad(ra2)
    dcrad1 = degtorad(dec1)
    dcrad2 = degtorad(dec2)

    radif = rarad2 - rarad1

    angle = np.arctan2(np.sin(radif), np.cos(dcrad1) * np.tan(dcrad2) - np.sin(dcrad1) * np.cos(radif))

    result = radtodeg(angle)

    if positive and (result < 0.0):
        result += 360.0

    racen = (ra1 + ra2) / 2.0
    deccen = (dec1 + dec2) / 2.0

    if plot:
        pyplot_rcparams()
        plt.figure()
        ax1 = plt.subplot(111)
        ax1.set_xlim([racen + 1.5*abs(ra1-ra2), racen - 1.5*abs(ra1-ra2)])
        ax1.set_ylim([deccen - 1.5*abs(dec1-dec2), deccen + 1.5*abs(dec1-dec2)])

        ax1.plot(ra1,dec1,'o',color='darkorange',label='Target 1')
        ax1.plot(ra2,dec2,'bo',label='Target 2')

        ax1.annotate('',xy=(0.9,0.9),xycoords='axes fraction',xytext=(0.9, 0.695),
                    arrowprops=dict(arrowstyle="->",edgecolor='k',facecolor='k'),
                    horizontalalignment='left',verticalalignment='top')
        ax1.annotate('',xy=(0.7,0.7),xycoords='axes fraction',xytext=(0.905, 0.7),
                    arrowprops=dict(arrowstyle="->",edgecolor='k',facecolor='k'),
                    horizontalalignment='left',verticalalignment='top')

        ax1.text(ra1,dec1,'Target 1', fontsize=14, color='darkorange',horizontalalignment='left',verticalalignment='bottom')
        ax1.text(ra2,dec2,'Target 2', fontsize=14, color='b',horizontalalignment='left',verticalalignment='bottom')

        ax1.text(0.75,0.75,'E', fontsize=14, color='k',transform=ax1.transAxes)
        ax1.text(0.85,0.85,'N', fontsize=14, color='k',transform=ax1.transAxes)
        ax1.text(0.1,0.1,'PA={:0.2f}'.format(result), fontsize=14, color='k',transform=ax1.transAxes)
        plt.show()

    return np.round(result, 2)

for iobj in range(len(ra_all)):
    #
    ra, dec = ra_all[iobj], dec_all[iobj]
    mag_target = magz[iobj]
    magj_target = magj[iobj]

    name, namebrief = hms2name(ra,dec)
    rahms = name[1:3]+':'+name[3:5]+':'+name[5:11]
    decdms = name[11:14]+':'+name[14:16]+':'+name[16:22]

    ra_star, dec_star = data['offset_ra'][iobj], data['offset_dec'][iobj]
    mag_star = data['mag_r'][iobj]
    name_star, _ = hms2name(ra_star, dec_star)
    rahms_off = name_star[1:3]+':'+name_star[3:5]+':'+name_star[5:11]
    decdms_off = name_star[11:14]+':'+name_star[14:16]+':'+name_star[16:22]

    # calculate offset
    delta_ra, delta_dec = offset(ra_star, dec_star, ra, dec, center=False)
    pa = positionAngle(ra_star, dec_star, ra, dec,plot=False)
    #pa = pos_angle[iobj]
    #if pa > 180:
    #    pa -= 360

    
    # Save everything to file
    print('{:}{:}      {:} {:} 2000.00 '
          'flag={:}'.format(flag[iobj],namebrief, rahms.replace(':',' '), decdms.replace(':',' '),   flag[iobj]))
    print('{:}{:}      {:} {:} 2000.00 '
          'flag={:}'.format(flag[iobj],namebrief, rahms.replace(':',' '), decdms.replace(':',' '),   flag[iobj]), file=f)

    #print('{:}      {:} {:} 2000.00 rotdest={:} rotmode=PA vmag={:0.1f} #jmag={:0.2f}, zmag={:0.2f}, Pri={:01d}, '
    #      'flag={:}'.format(namebrief, rahms.replace(':',' '), decdms.replace(':',' '), round(pa,2), mag_target, magj_target, mag_target, priority[iobj], flag[iobj]))
    #print('{:}      {:} {:} 2000.00 rotdest={:} rotmode=PA vmag={:0.1f} #jmag={:0.2f}, zmag={:0.2f}, Pri={:01d}, '
    #      'flag={:}'.format(namebrief, rahms.replace(':',' '), decdms.replace(':',' '), round(pa,2), mag_target, magj_target, mag_target, priority[iobj], flag[iobj]), file=f)
    #print('{:}_OFF  {:} {:} 2000.00 rotdest={:} rotmode=PA vmag={:0.1f} raoffset={:} d'
    #      'ecoffset={:}'.format(namebrief, rahms_off.replace(':',' '), decdms_off.replace(':',' '), round(pa,2), mag_star, round(delta_ra,2),round(delta_dec,2)),
    #      file=f)
f.close()

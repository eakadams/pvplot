#Functions for making pv plots

from __future__ import print_function

__author__ = "E.A.K. Adams"

"""
Functions for plotting position-velocity slices
Optimized for SHIELD galaxies, plotting
major and minor slices, following work
by Andrew McNichols

WCS and ImageGrid do not play nicely together
"""

import matplotlib.pyplot as plt
import os
import numpy as np
import astropy.units as u
from mpl_toolkits.axes_grid1 import ImageGrid
from astropy.io import fits

"""
PV plots:
1. Plot major & minor pv slices
"""

def plot_major_minor(agc,nmin=7,datadir=None, vmin=-2, vmax=4, rms=1.3,
                     velmax = None, cenvel = None):
    """
    Plot major and minor axes pv slices
    Assumes standard SHIELD naming
    Could adapt to take major axis name and 
    list of names for minor axies to be more general
    AGC number is used to name the output figure automatically.
    Inputs:
        agc (str): AGC number, used for completing standard SHIELD naming.
                    Also for naming output files
        nmin (int): Number of minor axis slices, default = 7
        datadir (str): Location of data. default = None, implies in local directory
        vmin (float): minimum intensity value for plotting
        vmax (float): maximum intensity value for plotting
        rms (float): rms data value in cube. default = 1.3 mJy/bm. Assumed value in mJy
        cenvel (km/s): center velocity in km/s
        velmax (km/s): max extent of velocity in km/s
    Outputs:
        fig: Returns the figure instance
    """
    # initialize a figure instance
    fig = plt.figure(figsize=(8,6))

    #Open a fits file to use the header for axes info
    #starting with minor
    if datadir == None:
        f = fits.open('{0}.min1.fits'.format(agc))
    else:
        f = fits.open('{0}/{1}.min1.fits'.format(datadir,agc))
    hminor = f[0].header
    #wcsminor = wcs.WCS(hminor,naxis=2) #ignore 3rd stokes axis
    # initialize the physical grid for the plot
    #nmin slices
    gridMinor = ImageGrid(fig, 212, (1,nmin), #ngrids=7,
                          cbar_mode="single", cbar_location='right',
                          cbar_size='12%', axes_pad=0.0,
                          #axes_class=(WCSAxes, dict(wcs=wcsminor)),
                          #would love to enable but doesn't play nice w/ ImageGrid :(
                          aspect=False, label_mode='L',share_all=True)

    #and initialize major axis
    if datadir == None:
        f = fits.open('{0}.maj.fits'.format(agc))
    else:
        f = fits.open('{0}/{1}.maj.fits'.format(datadir,agc))
    hmajor = f[0].header
    #wcsmajor = wcs.WCS(hmajor,naxis=2)
    # and its grid and cosmetcs
    gridMajor = ImageGrid(fig, 211, (1,1), #ngrids=1,
                          direction='column', cbar_mode="single",
                          cbar_location='right', cbar_size='2%',
                          axes_pad=0.0, 
                         #axes_class=(WCSAxes, dict(wcs=wcsmajor)),
                          aspect=False, share_all=True, label_mode='L')

    #set contour levels
    #relative to rms, can pick final levels later
    contours = np.array([2,4,8])*rms

    # loop over all the minor axis slices
    # numbering starts from 1
    for i in range(1,nmin+1):
        # open the file and headers
        if datadir == None:
            f = fits.open('{1}.min{2}.fits'.format(datadir,agc,i))
        else:
            f = fits.open('{0}/{1}.min{2}.fits'.format(datadir,agc,i))
        
        #get rid of third axies, which is first. grr.
        #squeeze so don't hve to know where empty axis is
        #stokes axis is annoying :(
        h, d = f[0].header, np.squeeze(f[0].data)
    
        #get axes from grid
        ax = gridMinor[(i-1)]

        # define some properties of the image
        im = ax.imshow(d * 1000, # mJy/bm
                       aspect = 'auto',
                       vmin=vmin, vmax=vmax,
                       origin='lower') #this makes sure y-axis runs the right way
                       #vmax=vmax, vmin=vmin, # call the scale range set above
                       #aspect='auto') 
                    
        #overplot contours
        ax.contour(d*1000,levels=contours,colors='black')
        
        #overplot velocity values
        #check that they are not none
        #require all both to not be none
        if (velmax is not None) and (cenvel is not None):
            #get pixel value for velocities
            pixmin,pixmax,cenpix = get_rot_pix(cenvel,velmax,h)
            #also need offset min & max for plotting
            #can do this straight from header. 
            #start at zero, go to n1=hminor['naxis1']
            noffset = h['naxis1']
            #plot min vel
            ax.plot([0,noffset-1],[pixmin,pixmin],color='black')
            #plot max vel
            ax.plot([0,noffset-1],[pixmax,pixmax],color='black')
            #plot cenvel
            ax.plot([0,noffset-1],[cenpix,cenpix],color='black',linestyle='--')
            
        
        #get tickvalues
        #first offset, then vel
        offset_vals,offsetunit,offset_ticks = get_offset_vals(h,sep=45*u.arcsec)
        print(offset_vals)
        #convert values to strings to be labels
        offset_labels = map(str,offset_vals)
        #and set tick locations / lavels
        ax.set_xticks(offset_ticks)
        ax.set_xticklabels(offset_labels,size=10)
        ax.set_xlabel('Offset',size=10)
    
        vel_vals,velunit,vel_ticks = get_vel_vals(h)
        vel_labels = map(str,vel_vals)
        ax.set_yticks(vel_ticks)
        ax.set_yticklabels(vel_labels,size=12)
        ax.set_ylabel('Velocity ({0})'.format(velunit),size=12)
        
        #plot the velocity extent, if specified
        
    #set color bar for last minor axis plot
    #stolen from previous code
    cb2 = gridMinor[0].cax.colorbar(im)
    cb2.set_label_text('mJy Bm$^{-1}$',
                       size='medium',
                       family='serif')
                
    #plot major axis
    if datadir is None:
        f = fits.open('{1}.maj.fits'.format(datadir,agc))
    else:
        f = fits.open('{0}/{1}.maj.fits'.format(datadir,agc))
    h, d = f[0].header, np.squeeze(f[0].data)
    ax = gridMajor[0]

    # define some properties of the image
    im = ax.imshow(d * 1000, # mJy/bm
                   origin='lower', # put x-axis below plot
                   vmax=vmax, vmin=vmin, # call the scale range set above
                   aspect='auto') 
    #overplot contours
    ax.contour(d*1000,levels=contours,colors='black')
    
    #overplot velocity values
    #check that they are not none
    #require all both to not be none
    if (velmax is not None) and (cenvel is not None):
        #get pixel value for velocities
        pixmin,pixmax,cenpix = get_rot_pix(cenvel,velmax,h)
        #also need offset min & max for plotting
        #can do this straight from header. 
        #start at zero, go to n1=hminor['naxis1']
        noffset = h['naxis1']
        #plot min vel
        ax.plot([0,noffset-1],[pixmin,pixmin],color='black')
        #plot max vel
        ax.plot([0,noffset-1],[pixmax,pixmax],color='black')
        #plot cenvel
        ax.plot([0,noffset-1],[cenpix,cenpix],color='black',linestyle='--')
    
    #get tickvalues
    #first offset, then vel
    offset_vals,offsetunit,offset_ticks = get_offset_vals(h)
    #convert values to strings to be labels
    offset_labels = map(str,offset_vals)
    #and set tick locations / lavels
    ax.set_xticks(offset_ticks)
    ax.set_xticklabels(offset_labels,size=10)
    ax.set_xlabel('Offset ({0})'.format(offsetunit),size=10)
    
    vel_vals,velunit,vel_ticks = get_vel_vals(h)
    vel_labels = map(str,vel_vals)
    ax.set_yticks(vel_ticks)
    ax.set_yticklabels(vel_labels)
    ax.set_ylabel('Velocity ({0})'.format(velunit),size=12)
    
    #add some colorbar specification
    #copied from previous PV slice code
    cb2 = gridMajor[0].cax.colorbar(im)
    cb2.set_label_text('mJy Bm$^{-1}$',
                       size='medium',
                       family='serif')
    
    #save figure
    fig.savefig('{0}_pv_major_minor.pdf'.format(agc))
    
    #return fig instance
    return fig

"""
Functions for getting tick values and locations.
This is because WCS and ImageGrid do not play nicely together
"""

def get_rot_pix(cenvel,velmax,header,
                headerunit = u.m / u.s, 
                naxis='NAXIS2', cdelt='CDELT2', crval = 'CRVAL2', crpix='CRPIX2'):
    """
    Get pixel values for rotation velocity values
    Inputs:
        cenvel: center velocity
        velmax: max extent of velocity
        header: header to get coordiante transforms from
    Outputs:
        pixmin, pixmax, cenpix: Pixel values for min, max, cen
    """
    vmin = cenvel - velmax/2.
    vmax = cenvel + velmax/2.
    pixmin = ( (vmin.to(headerunit) - header[crval]*headerunit) /
             (header[cdelt]*headerunit) + header[crpix] )
    pixmax = ( (vmax.to(headerunit) - header[crval]*headerunit) /
             (header[cdelt]*headerunit) + header[crpix] )
    cenpix = ( (cenvel.to(headerunit) - header[crval]*headerunit) /
             (header[cdelt]*headerunit) + header[crpix] )
    return pixmin,pixmax,cenpix

def get_offset_vals(header,sep=15*u.arcsec,headerunit = u.arcsec, 
                    naxis='NAXIS1', cdelt='CDELT1', crval = 'CRVAL1', crpix='CRPIX1'):
    """
    Take a header object and return an array of offset values, plus the unit.
    Also return the tick values
    Inputs:
        header: header object from astropy.fits
        sep (int, angular unit): separation for values in an angular unit
        headerunit (float, angular unit): Unit for header for offset axis
        naxis (str): naxis in header, default='NAXIS1'
        cdelt (str): cdelt in header, default='CDELT1'
        crval (str): crval in header, default='CRVAL1'
        crpix (str): crpix in header, default='CRPIX1'
    Outputs:
        offset_vals (nparray): Array of offset values (ints)
        offsetunit: Unit for separation values
        offset_ticks: Locations of ticks in pixel units
    """
    total_offset = header[naxis] * header[cdelt] * headerunit
    nticks = np.floor(total_offset / sep)
    #check if odd or not - want odd tick marks
    #if even, add 1 tick (since I did floor before)
    if nticks % 2 == 0:
        nticks = nticks + 1
    #set up array for values
    #map to integers
    offset_float_vals = (np.arange(nticks) - np.floor(nticks/2)) * sep.value
    offset_vals = map(int,offset_float_vals)
    offsetunit = sep.unit

    #also get tick values that go with
    #have to use original array, float messes things up
    offset_ticks = (((offset_float_vals*offsetunit).to(headerunit) - header[crval]*headerunit) 
                    / (header[cdelt]*headerunit) + header[crpix])
    
    return offset_vals,offsetunit, offset_ticks

def get_vel_vals(header,sep=5*u.km/u.s,headerunit = u.m / u.s, 
                    naxis='NAXIS2', cdelt='CDELT2', crval = 'CRVAL2', crpix='CRPIX2'):
    """
    Take a header object and return an array of velocity values
    Also return pixel coordinate locations for placing ticks
    Inputs:
        header: header object from astropy.fits
        sep (int, vel unit): separation for values in a velocity unit
        headerunit (float, vel unit): Unit for header for velocity axis
        naxis (str): naxis in header, default='NAXIS2'
        cdelt (str): cdelt in header, default='CDELT2'
        crval (str): cdelt in header, default='CRVAL2'
        crpix (str): crpix in header, default='CRPIX2'
    Outputs:
        vel_vals (nparray): Array of velocity values (ints)
        velunit: Unit for separation values
        vel_ticks (nparray): Array of tick locations in pixel units
    """
    total_velocity = header[naxis] * header[cdelt] * headerunit
    nticks = np.floor( (total_velocity.to(sep.unit) / sep) )
    #set starting value, go up to nearest sep value
    #implicitly assuming crpix is 1. should probably add a check for this
    extra_over = ((header[crval]*headerunit).to(sep.unit).value % sep.value)
    extra_over = ((header[crval]*headerunit).to(sep.unit) % sep)
    startval = (header[crval]*headerunit).to(sep.unit) + (sep - extra_over)
    #set up array for values
    vel_float_vals = (np.arange(nticks)) * sep.value + startval.value
    vel_vals = map(int,vel_float_vals)
    velunit = sep.unit
    
    #also get tick values that go with
    vel_ticks = (((vel_float_vals*velunit).to(headerunit) - header[crval]*headerunit) 
                    / (header[cdelt]*headerunit) + header[crpix])
    return vel_vals,velunit, vel_ticks
    
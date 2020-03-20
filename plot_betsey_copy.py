#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

import numpy as np
from mpl_toolkits.axes_grid1 import ImageGrid
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import pywcsgrid2 as wcs

"""
BA:
SUPER IMPORTANT NOTE. Changed to saving with fig.savefig() rather than
plt.savefig(). This got rid of blank images (plt interface historically
is troublesome).
TO DO:
Add labels to figures
Get xaxis labels to be better suited to figures
Comment code more thoroughly
"""


# hardcode picture paths
#NOTE: I had to update these directory paths to get this code to work
#This would be first bug to trouble-shoot in future
fitsDir = '/Users/research/Desktop/whatUneed/slices_fits/'
psDir = '/Users/research/Desktop/whatUneed/slices_eps/'

"""BA: 
Now I want to add some extra code so that it will automatically loop 
over all SHIELD sources in one go. 
I would also like to include noise info so that I can do contours 
as multiples of rms values. 
I will assume robust cubes are used throughout and
copy those noise values over.
While it would be best to find the SHIELD galaxy names automatically
from fits directory listings,
I'll go ahead and hard code it so that I can manually link noise values
"""

gal_names = ('102728','123352','124056','191706','198507','198508','198691','200232','205590','223231','223254','229053', '229379', '238890', '718245', '728909', '731921', '739005', '742601', '747826')
#noise_vals=(1.03,1.44,1.12,1.50,1.29,0.875,0.51,0.685,0.875,0.925,0.79,0.79)
noise_vals=(1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.)
vcen=(566,467,407,561,502,519,514,450,494,571,603,425,624,360,651,699,504,433,539,558)

#These are central velocities in LSR frame calculated from get_vlsr.py
#vcen=(354.74155226,  160.27316496,  364.57795957,  204.34150153,  272.02375135, 342.61225936,  337.60248038,  392.03355037,  453.55487116,  259.96797894, 378.9120415,458.60469094)

#these are rotational velocities estimated from mcgaugh2012/new_plot_mcgaugh.py
vrot=(21.,25.,22.,23.,37.,29.,33.,49.,29.,19.,19.,40.,22.,20.,15.,33.,62.,46.,27.,31.)


#now start my iteration on galaxies before using any of andrew's code:
#(will have to indent everything below!)
for gal,noise,vc,vr in zip(gal_names,noise_vals,vcen,vrot):
    #define contour levels here
    #taking 2,4,8 times the noise
    cont_levs=np.append(np.arange(2,6,2),8)*noise



    # initialize a figure instance
    fig = plt.figure(figsize=(8,6)) #x y tuple in inches
    #fig = plt.figure()
    fontScale = 10

    # matplotlib uses matplotlibrc configuration files to customize all
    # kinds of properties, which we call rc settings or rc parameters.
    # clear and then reset the high-level user settings of matplotlib
    plt.rcdefaults()
    params = {'backend': 'eps',
              'axes.labelsize': fontScale,
              'font.size': fontScale,
              'legend.fontsize': fontScale*3/4,
              'xtick.labelsize': fontScale,
              'ytick.labelsize': fontScale,
              #'text.usetex': True
              }
    plt.rcParams.update(params)
    plt.figure(1)
    plt.clf()

    # Open a fits file to use the header for axes info
    f = pyfits.open(fitsDir + gal+'.min1.fits')
    hminor = f[0].header
    # initialize the physical grid for the plot
    gridMinor = ImageGrid(fig, (2,1,2),
                          nrows_ncols=(1,7),
                          ngrids=7,
                          cbar_mode="single",
                          cbar_location='right',
                          cbar_size='12%',
                          axes_pad=0.0,
                          axes_class=(wcs.Axes, dict(header=hminor)),
                          aspect=False,
                          label_mode='L',
                          share_all=True)


    # A second file for the MAJOR axis slice plot
    f = pyfits.open(fitsDir + gal+'.maj.fits')
    hmajor = f[0].header
    # and its grid and cosmetcs
    gridMajor = ImageGrid(fig, (2,1,1),
                          nrows_ncols=(1,1),
                          ngrids=1,
                          direction='column',
                          cbar_mode="single",
                          cbar_location='right',
                          cbar_size='2%',
                          axes_pad=0.0,
                          axes_class=(wcs.Axes, dict(header=hmajor)),
                          aspect=False,
                          share_all=True,
                          label_mode='L',
                          )
    
    #gridMajor.set_aspect(1)
    #gridMinor.set_aspect(25)


    # set the range properties for the colorbar and velo axis
    vmin = 0.0
    vmax = 3
    #can get velrange from header information
    velRange = (0,hminor['naxis2']) #119 max range for 174585, 79 for others
    #10/20 or 59 for the offset sources


    # define value for marking over images
    annoColor='k'

    # loop over all the minor axis slices
    for i in xrange(1,8):
        # open the file and headers
        f = pyfits.open(fitsDir+gal+'.min'+ str(i)+'.fits')
        h, d = f[0].header, f[0].data[:,:]
        ax = gridMinor[i]

        # define some properties of the image
        im = ax.imshow(d * 1000, # mJy/bm
                       #origin='lower', # put x-axis below plot
                       cmap = plt.cm.YlOrBr, # define color bar scale
                       vmax=vmax, vmin=vmin, # call the scale range set above
                       aspect='auto'
                       )    

        #add the contours
        ax.contour(d*1000,levels=cont_levs,colors='blue')


        #add lines highlighting vcen and vrot
        #need to get x-axis extent because it's not explicitly set anywhere
        #can do this from headerinfo
        #just care about pixel values
        #also need to put velocity into pix values - use wcs conversion
        n1=hminor['naxis1']

        #now get vel values in pixels:
        vc_pix = (vc*1000.-hminor['crval2'])/hminor['cdelt2']+hminor['crpix2']
        vr_delt = vr*1000/hminor['cdelt2']

        #and now plot things
        ax.plot([0,n1-1],[vc_pix,vc_pix],'k--',linewidth=1.3)
        ax.plot([0,n1-1],[(vc_pix+vr_delt),(vc_pix+vr_delt)],'k',linewidth=1.3)
        ax.plot([0,n1-1],[(vc_pix-vr_delt),(vc_pix-vr_delt)],'k',linewidth=1.3)    
        

        #need to get velocity labels
        #want to center near velocity of source
        #or it will plot pixel values instead.
        #but this way it figures out where I want to be
        labs=[]
        start=np.floor((vc-vr*2)/10.)*1e4
        end=np.ceil((vc+vr*2)/10.)*1e4
        vals=np.arange(start,end,10e3)
        for el in vals:
            labs.append(str(int(el/1e3)))


        ax.set_ticklabel2_type('manual', 
                               locs=vals,
                               labels=labs)
                               #locs=[320e3,330e3,340e3,350e3,360e3,370e3],
                               #labels=['320','330','340','350','360','370'])


        #now I want to do something similar for x-axis
        #but in this case I think I will pick a set of values to use based
        #on size of xaxis
        

        cdelt1=hminor['cdelt1']
        span = n1*cdelt1

        if span < 70:
            ax.set_ticklabel1_type('manual',
                                   locs=[-15,0,15],
                                   labels=['-15','0','15'])
        if span >= 70 and span < 120:
            ax.set_ticklabel1_type('manual',
                                   locs=[-30,0,30],
                                   labels=['-30','0','30'])
        if span >= 120:
            ax.set_ticklabel1_type('manual',
                                   locs=[-45,0,45],
                                   labels=['-45','0','45'])

        # set and mark the boundaries of the images
        ax.set_ylim(velRange[0],velRange[1])


        ax.set_xlabel('Offset [\'\']',
                      size = 'medium',
                      family='serif')
        ax.set_ylabel('Velocity [km s$^{-1}$]',
                      size = 'medium',
                      family='serif')

        ax.annotate('['+str(i)+']', (0.05,0.9), # annotation text and position
                    color=annoColor, xycoords='axes fraction', # color and units) 
                    size=15)

        # Set tick and lines to white
        ax.axis[:].major_ticks.set_color("k")

        # force the interior lines and ticks to be white
        if i==0:
            ax.axis['right'].line.set_color("k")
        elif i==7:    
            ax.axis['left'].line.set_color('k')
        else:
            ax.axis['right'].line.set_color("k")
            ax.axis['left'].line.set_color("k")


    #ax.set_aspect(aspect=1/2.,
    #              adjustable='box-forced',
    #              anchor='C')


    cb2 = gridMinor[0].cax.colorbar(im)
    cb2.set_label_text('mJy Bm$^{-1}$',
                       size='medium',
                       family='serif')

    """
    Now do the major axis. This is very analagous to above
    EXCEPT that there's no need ot iterate over axes
    Solutions for vels, labels, etc. are identical
    """

    ######### 
    # Major axis definites and plotting here
    # open the file and headers
    f = pyfits.open(fitsDir+gal+'.maj.fits')
    h, d = f[0].header, f[0].data[:,:]
    ax = gridMajor[0]

    # define some properties of the image
    im = ax.imshow(d * 1000, # mJy/bm
                   origin='lower', # put x-axis below plot
                   cmap = plt.cm.YlOrBr, # define color bar scale
                   vmax=vmax, vmin=vmin, # call the scale range set above
                   aspect='auto')    

    #and plot the contours also
    ax.contour(d*1000,levels=cont_levs,colors='blue')

    #and add labels here
    ax.text(5,hmajor['naxis2']-5,'AGC '+gal)


    n1=hmajor['naxis1']

    #now get vel values in pixels:
    vc_pix = (vc*1000.-hmajor['crval2'])/hmajor['cdelt2']+hmajor['crpix2']
    vr_delt = vr*1000/hmajor['cdelt2']
    ax.plot([0,n1-1],[vc_pix,vc_pix],'k--',linewidth=1.3)
    ax.plot([0,n1-1],[(vc_pix+vr_delt),(vc_pix+vr_delt)],'k',linewidth=1.3)
    ax.plot([0,n1-1],[(vc_pix-vr_delt),(vc_pix-vr_delt)],'k',linewidth=1.3)    
        


    labs=[]
    start=np.floor((vc-vr*2)/10.)*1e4
    end=np.ceil((vc+vr*2)/10.)*1e4
    vals=np.arange(start,end,10e3)
    for el in vals:
        labs.append(str(int(el/1e3)))


    ax.set_ticklabel2_type('manual', 
                           locs=vals,
                           labels=labs)

    #ax.set_ticklabel2_type('manual', 
    #                       locs=[320e3,330e3,340e3,350e3,360e3,370e3],
    #                       labels=['320','330','340','350','360','370'])




    cdelt1=hmajor['cdelt1']
    span = n1*cdelt1

    if span < 90:
        ax.set_ticklabel1_type('manual',
                               locs=[-30,-15,0,15,30],
                               labels=['-30','-15','0','15','30'])
    if span >= 90 and span < 130:
        ax.set_ticklabel1_type('manual',
                               locs=[-45,-30,-15,0,15,30,45],
                               labels=['-45','-30','-15','0','15','30','45'])
    if span >= 130:
        ax.set_ticklabel1_type('manual',
                               locs=[-60,-45,-30,-15,0,15,30,45,60],
                               labels=['-60','-45','-30','-15','0','15','30','45','60'])


    #ax.set_ticklabel1_type('manual',
    #                       locs=[-45,-30,-15,0,15,30,45],
    #                       labels=['-45','-30','-15','0','15','30','45'])

    # set and mark the boundaries of the images
    ax.set_ylim(velRange[0],velRange[1])

    ax.set_xlabel('Offset [\'\']',
                  size = 'medium',
                  family='serif')
    ax.set_ylabel('Velocity [km s$^{-1}$]',
                  size = 'medium',
                  family='serif')

    # Set tick and lines to black
    ax.axis[:].major_ticks.set_color("k")
    # force the interior lines and ticks to be black
    ax.axis['right'].line.set_color("k")
    ax.axis['left'].line.set_color("k")

    cb2 = gridMajor[0].cax.colorbar(im)
    cb2.set_label_text('mJy Bm$^{-1}$',
                       size='medium',
                       family='serif')

    #save image
    #plt.show()
    fig.savefig(psDir+'/'+gal+'.slices.eps', dpi=400)#,bbox_inches='tight')
    
    #and close the plot
    plt.close()


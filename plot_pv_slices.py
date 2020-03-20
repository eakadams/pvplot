#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

import numpy as np
from mpl_toolkits.axes_grid1 import ImageGrid
from astropy.io import fits as pyfits
import matplotlib.pyplot as plt
import pywcsgrid2 as wcs

# hardcode picture paths
fitsDir = '/Users/research/Desktop/sliceview_casa_28.7.15/fits/'
psDir = '/Users/research/Desktop/sliceview_casa_28.7.15/eps/'

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
f = pyfits.open(fitsDir + '174585.min.0.fits')
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
f = pyfits.open(fitsDir + '174585.maj.fits')
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
velRange = (0,59) #119 max range for 174585, 79 for others
                  #10/20 or 59 for the offset sources

# define value for marking over images
annoColor='k'

# loop over all the minor axis slices
for i in xrange(0,7):
    # open the file and headers
    f = pyfits.open(fitsDir+'174585.min.'+ str(i)+'.fits')
    h, d = f[0].header, f[0].data[:,:]
    ax = gridMinor[i]

    # define some properties of the image
    im = ax.imshow(d * 1000, # mJy/bm
                   origin='lower', # put x-axis below plot
                   cmap = plt.cm.YlOrBr, # define color bar scale
                   vmax=vmax, vmin=vmin, # call the scale range set above
                   aspect='auto',
                   )    

    ax.set_ticklabel2_type('manual', 
                           locs=[320e3,330e3,340e3,350e3,360e3,370e3],
                           labels=['320','330','340','350','360','370'])

    ax.set_ticklabel1_type('manual',
                           locs=[-30,0,30],
                           labels=['-30','0','30'])

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

######### 
# Major axis definites and plotting here
# open the file and headers
f = pyfits.open(fitsDir+'174585.maj.fits')
h, d = f[0].header, f[0].data[:,:]
ax = gridMajor[0]

# define some properties of the image
im = ax.imshow(d * 1000, # mJy/bm
               origin='lower', # put x-axis below plot
               cmap = plt.cm.YlOrBr, # define color bar scale
               vmax=vmax, vmin=vmin, # call the scale range set above
               aspect='auto')    

ax.set_ticklabel2_type('manual', 
                       locs=[320e3,330e3,340e3,350e3,360e3,370e3],
                       labels=['320','330','340','350','360','370'])

ax.set_ticklabel1_type('manual',
                       locs=[-45,-30,-15,0,15,30,45],
                       labels=['-45','-30','-15','0','15','30','45'])

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

#plt.savefig(psDir+'/174585.slices.eps', dpi=400)#,bbox_inches='tight')
plt.show()



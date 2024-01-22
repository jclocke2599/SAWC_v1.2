###########################################################################
#
# PYTHON 3 FUNCTIONS FOR tools_plotting
#
# CONTRIBUTORS: Katherine E. Lukens             NOAA/NESDIS/STAR, CISESS at U. of Maryland
#               Kevin Garrett                   NOAA/NWS/OSTI
#               Kayo Ide                        U. of Maryland
#               David Santek                    CIMSS at U. of Wisconsin-Madison
#               Brett Hoover                    NOAA/NWS/NCEP/EMC, Lynker Technologies
#               David Huber                     NOAA/NWS/NCEP/EMC, Redline Performance Solutions, LLC
#               Ross N. Hoffman                 NOAA/NESDIS/STAR, CISESS at U. of Maryland
#               Hui Liu                         NOAA/NESDIS/STAR, CISESS at U. of Maryland
#
# Built with the following conda environment:
#
# name: bhoover-obs_match_3d
# channels:
#   - conda-forge
#   - defaults
# dependencies:
#   - python=3
#   - numpy
#   - pandas
#   - pynio
#   - matplotlib
#   - cartopy
#   - jupyter
#   - netCDF4
#   - scikit-learn
#   - dask
#   - geopy
#   - pip
#   - pip:
#
###########################################################################
#
# Import modules
#

import os
import sys
import numpy as np #....................................................... Array module
import datetime as dt #.................................................... Datetime module
import time #.............................................................. Time module
import warnings
import math
import statistics as stats

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize 
from matplotlib.colors import LogNorm
import matplotlib.ticker as mticker
from matplotlib.lines import Line2D
import matplotlib.cm as cm
from matplotlib.cm import get_cmap

from scipy.interpolate import interpn
from scipy.stats import gaussian_kde
from scipy.stats import kde
from scipy.stats import pearsonr
from scipy.stats import ttest_ind

from mpl_toolkits.basemap import Basemap
from mpl_toolkits.mplot3d import Axes3D

import cartopy.crs as ccrs
import cartopy.mpl.ticker as cticker
from cartopy.util import add_cyclic_point

from tools_analysis import var_to_2d
from tools_analysis import var_to_2d_counts_1dset

#
###########################################################################

	# missing value
fill = -999

	# font sizes for plots
fonttitle  = 16
fontaxis   = 15
fontlegend = 12

###########################################################################
#
# PLOTTING functions for collocation program
#

# -------------------------------------------------------------------------
# Mean vs Height/Pressure of Matched Observations
#
def contour2d(var2d,var2dstr,xvar,yvar,xvarstr,yvarstr,tname,units,plottype,level,levunits,regionstr,**kwargs):

	# create plot
    fig = plt.figure(figsize=(10,7.5))
    
    ax1 = fig.add_subplot()
   
    ax1.set_xlabel(xvarstr, fontsize=fontaxis)
    ax1.set_ylabel(yvarstr, fontsize=fontaxis)

    if xvarstr.find('Time') != -1:
	# get unique dates (yyyymmdd format)
      datestr  = xvar
      ndatestr = len(datestr)
      datelabels = datestr
  
      tmp = var2d; del var2d
      var2d = tmp[:,0:ndatestr-1]; del tmp
      
        # set xaxis
      ilast = ndatestr-1
              # x-axis specs
      timestepstr = "days"              # include 00-23 hours
              # x-axis tickmark values and labels
      if ndatestr<10:
        stride  = 1
      elif ndatestr>=10 and ndatestr<30:
        stride = 2
      elif ndatestr>=30 and ndatestr<93:
        stride = 5
      elif ndatestr>=93: 
        stride = 30
      xaxis   = np.arange(0, ilast, 1)
      xvalues = np.arange(0, ilast, stride)
      xlabels = datelabels[0:ilast:stride]

      title = str(regionstr)+' Mean'

      ax1.set_xticks(xvalues)
      ax1.set_xticklabels(xlabels, rotation=45)      

      del xvalues,xlabels
    else:
      if xvarstr.find('Latitude') != -1:
        title = 'Zonal Mean'
      elif xvarstr.find('Longitude') != -1 and yvarstr.find('Latitude')==-1:
        title = 'Meridional Mean'
      else:
        title = ''
      xvalues = np.arange(0,np.size(xvar),30)
      xlabels = xvar[::30]
      ax1.set_xticks(xvalues)
      ax1.set_xticklabels(xlabels, fontsize=fontaxis-2)
      del xvalues,xlabels

    if yvarstr.find('Pressure') != -1:
      nspacing = 4
      yvalues = np.arange(0,np.size(yvar),nspacing)
      ylabels = yvar[::nspacing]
      ax1.set_yticks(yvalues)
      ax1.set_yticklabels(ylabels, fontsize=fontaxis-2)
      del yvalues,ylabels,nspacing
    elif yvarstr.find('Height') != -1:
      nspacing = 2
      yvalues = np.arange(0,np.size(yvar),nspacing)
      ylabels = yvar[::nspacing]
      ax1.set_yticks(yvalues)
      ax1.set_yticklabels(ylabels, fontsize=fontaxis-2)
      del yvalues,ylabels,nspacing

	# set range of plotted colored contours
    warnings.filterwarnings('ignore', category=RuntimeWarning)		# ignore warnings
    varmax = np.nanmax(abs(var2d))
    varmin = np.nanmin(var2d)

    if var2dstr.find('Count') != -1:
      ppcmap = "jet"
      pltmin = 1
      pltmax = 100
        # reassign color range max
      if varmax!=np.nan and varmax>=100000  : pltmax=1000000
      if varmax!=np.nan and varmax<100000: pltmax=10000
      if varmax!=np.nan and varmax<10000: pltmax=1000
      if varmax!=np.nan and varmax<1000: pltmax=100
      if varmax!=np.nan and varmax<100: pltmax=10
      var2d  = np.where(var2d==np.nan,0,var2d)
    else:
      if var2dstr.find('SD') != -1:
        ppcmap = "jet"
        pltmin = 1
        pltmax = 100
      elif var2dstr.find('Diff') != -1:
        ppcmap = "seismic"
        pltmin = -30
        pltmax = 30
      else:
        ppcmap = "jet"
        pltmin = -50
        pltmax = 50
      # reassign color range min/max
      if varmax!=np.nan and varmax>=60: pltmax=100
      if varmax!=np.nan and varmax<60: pltmax=50
      if varmax!=np.nan and varmax<50: pltmax=40
      if varmax!=np.nan and varmax<40: pltmax=30
      if varmax!=np.nan and varmax<30: pltmax=20
      if varmax!=np.nan and varmax<20: pltmax=10
      # reassign pltmin if not SD plot
      if var2dstr.find('SD') == -1: pltmin = -pltmax
      if var2dstr.find('Diff') == -1 and varmin>=0: pltmin=1
    del varmax,varmin

    originstr = 'lower'

    if plottype=="2d":
      if var2dstr.find('Count') != -1:
        if pltmax==10:
          pp = plt.imshow(var2d,cmap=ppcmap,vmin=pltmin,vmax=pltmax, origin=originstr, aspect='auto')
        else:
          pp = plt.imshow(var2d,cmap=ppcmap,norm=LogNorm(vmin=pltmin,vmax=pltmax), origin=originstr, aspect='auto')
      else:
        pp = plt.imshow(var2d, cmap=ppcmap, vmin=pltmin, vmax=pltmax, origin=originstr, aspect='auto')

      if var2dstr.find('Count') != -1:
        ax1.set_title(title+" Number Density of Matched Obs for "+tname, fontsize=fonttitle-2)
      elif var2dstr.find('SD') != -1:
        ax1.set_title(title+" "+var2dstr+" of Diff for "+tname, fontsize=fonttitle-2)
      else:
        ax1.set_title(title+" "+var2dstr+" "+str(level)+" "+levunits+" for "+tname, fontsize=fonttitle-2)

    elif plottype.find("map") != -1:
      x,y = np.meshgrid(xvar,yvar)
       	# create map
      proj  = ccrs.PlateCarree()
      axmap = plt.subplot(1,1,1,projection=proj)

	# add cyclic longitude coordinates
      var2d_cyc, x_cyc = add_cyclic_point(var2d, coord=np.asarray(xvar), axis=1)

      yvararr = np.asarray(yvar)

	# plot
      if var2dstr.find('Count') != -1:
        if pltmax==10:
          pp = plt.pcolormesh(x_cyc,yvararr,var2d_cyc,cmap=ppcmap,vmin=pltmin,vmax=pltmax,transform=proj)
        else:
          pp = plt.pcolormesh(x_cyc,yvararr,var2d_cyc,cmap=ppcmap,norm=LogNorm(vmin=pltmin,vmax=pltmax),transform=proj)
      else:
        pp = plt.pcolormesh(x_cyc,yvararr,var2d_cyc, cmap=ppcmap, vmin=pltmin, vmax=pltmax, transform=proj)

      axmap.coastlines()
     
      if var2dstr.find('Count') != -1:
        axmap.set_title(title+" Number Density of Matched Obs for "+tname, fontsize=fonttitle-2)
      elif var2dstr.find('SD') != -1:
        axmap.set_title(title+" "+var2dstr+" of Diff for "+tname, fontsize=fonttitle-2)
      else: 
        axmap.set_title(title+" "+var2dstr+" at "+str(level)+" "+levunits+" for "+tname, fontsize=fonttitle-2)

      gl = axmap.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.7, linestyle='--')
      gl.xlines = True            # plot longitudes
      gl.ylines = True            # plot latitudes
      gl.top_labels=None
      gl.right_labels=None

        # set GLOBAL domain limits
      latmin = -90.0
      latmax = 90.0
      lonmin = -180.0
      lonmax = 180.0
      clat   = (latmax + latmin) / 2
      clon   = (lonmax + lonmin) / 2
      axmap.set_extent([lonmin,lonmax,latmin,latmax], crs=proj)
      		# define x-axis tickmarks
      gl.xlocator = mticker.FixedLocator([-180,-120,-60,0,60,120,180])
      		# define y-axis tickmarks
      gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])

      gl.xlabel_style = {'size': fontaxis, 'color': 'black'}#'gray'}
      gl.ylabel_style = {'size': fontaxis, 'color': 'black'}#'gray'}

    if yvarstr.find('Pressure') != -1:
      plt.gca().invert_yaxis()

	# Add a colorbar to the bottom of the plot.
    cbarstr = var2dstr
    if var2dstr.find('Count') != -1:
	# is Count
      del cbarstr
      cbarstr = "Number Density of Matched Obs"
    elif var2dstr.find('Count') == -1:
      cbarstr += " ("+str(units)+")"

    if xvarstr=='Time':
      cbar = fig.colorbar(pp, label=cbarstr, orientation='horizontal', pad=0.125)
    elif xvarstr=='Longitude' and yvarstr=='Latitude':
      frac = 0.024
      cbar = fig.colorbar(pp, label=cbarstr, orientation='vertical', fraction=frac, pad=0.04)
    else:   
      frac = 0.046
      cbar = fig.colorbar(pp, label=cbarstr, orientation='vertical', fraction=frac, pad=0.04)
    cbar.ax.tick_params(labelsize=fontlegend)

    return fig

# -------------------------------------------------------------------------
# Mean vs Height/Pressure of Matched Observations
#
def contour2d_orthomap(datein,var2d,var2dstr,xvar,yvar,xvarstr,yvarstr,tname,units,center_lon,center_lat,level,levunits,**kwargs):

	# create plot
    fig = plt.figure(figsize=(10,8.5))
    
    ax1 = fig.add_subplot()
   
    ax1.set_xlabel(xvarstr)
    ax1.set_ylabel(yvarstr)

    if xvarstr.find('Latitude') != -1:
      title = 'Zonal Mean'
      xvalues = np.arange(0,np.size(xvar),30)
      xlabels = xvar[::30]
      ax1.set_xticks(xvalues)
      ax1.set_xticklabels(xlabels)
      del xvalues,xlabels
    elif xvarstr.find('Longitude') != -1:
      title = ''

	# set range of plotted colored contours
    warnings.filterwarnings('ignore', category=RuntimeWarning)                # ignore warnings
    varmax = np.nanmax(abs(var2d))
    varmin = np.nanmin(var2d)

    if var2dstr.find('Count') != -1:
      ppcmap = "jet"
      pltmin = 1
      pltmax = 100
        # reassign color range max
      if varmax!=np.nan and varmax>=100000  : pltmax=1000000
      if varmax!=np.nan and varmax<100000: pltmax=10000
      if varmax!=np.nan and varmax<10000: pltmax=1000
      if varmax!=np.nan and varmax<1000: pltmax=100
      if varmax!=np.nan and varmax<100: pltmax=10
      var2d  = np.where(var2d==np.nan,0,var2d)
    else:
      if var2dstr.find('SD') != -1:
        ppcmap = "jet"
        pltmin = 1
        pltmax = 50
        stride = 5
      elif var2dstr.find('Diff') != -1:
        ppcmap = "seismic"
        pltmin = -30
        pltmax = 30
        stride = 5
      else:
        ppcmap = "jet"
        pltmin = -50
        pltmax = 50
        stride = 5
      # reassign color range max
      if varmax!=np.nan and varmax>=60: pltmax=100
      if varmax!=np.nan and varmax<60: pltmax=50
      if varmax!=np.nan and varmax<50: pltmax=40
      if varmax!=np.nan and varmax<40: pltmax=30
      if varmax!=np.nan and varmax<30: pltmax=20
      if varmax!=np.nan and varmax<20: pltmax=10
      # reassign pltmin if not SD plot 
      if var2dstr.find('SD') == -1: pltmin = -pltmax
      if varmin>=0: pltmin=1
    del varmax,varmin

    originstr = 'lower'

      # create map
    proj     = ccrs.Orthographic(central_longitude=center_lon,central_latitude=center_lat,globe=None)
    geo      = ccrs.Geodetic()
    poleproj = ccrs.RotatedPole(pole_latitude=center_lat,pole_longitude=center_lon)

    axmap    = plt.subplot(1,1,1,projection=proj)
     
      # add coastlines
    axmap.coastlines(resolution='50m')
    axmap.set_global()

      # set longitude gridlines and labels
    nmerid = 30		# spacing for longitude lines
    nzonal = 15		# spacing for latitude lines
    num_merid = int(360/nmerid + 1)
    num_zonal = int(90/nzonal + 1)
    if center_lat>0: ylocsPole=np.linspace(0,90,num_zonal)	# NH
    if center_lat==0: ylocsPole=np.linspace(-90,90,num_zonal)	# NH
    if center_lat<0: ylocsPole=np.linspace(-90,0,num_zonal)	# SH
    xlocsPole = np.linspace(-180,180,num_merid)
    axmap.gridlines(xlocs=xlocsPole, ylocs=ylocsPole, linestyle="--", linewidth=1, color='gray', alpha=0.5)

    gl = axmap.gridlines(draw_labels=False, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xlines = True            # plot longitudes
    gl.ylines = True            # plot latitudes
	# define x-axis tickmarks
    gl.xlocator = mticker.FixedLocator(xlocsPole)#[-180,-120,-60,0,60,120,180])
        # define y-axis tickmarks
    gl.ylocator = mticker.FixedLocator(ylocsPole)#[-90,-60,-30,0,30,60,90])

	# plot on map
    x,y = np.meshgrid(xvar,yvar)
    if var2dstr.find('Count') != -1:
      pp = axmap.scatter(x,y,c=var2d,cmap=ppcmap,norm=LogNorm(vmin=pltmin,vmax=pltmax),marker='s',s=5,transform=geo)
    else:
      pp = plt.pcolormesh(x,y,var2d,cmap=ppcmap,vmin=pltmin,vmax=pltmax,transform=ccrs.RotatedPole())
   
    if var2dstr.find('Count') != -1:
      axmap.set_title(title+" Number Density of Matched Obs for "+tname)
    else: 
      axmap.set_title(title+" "+var2dstr+" at "+str(level)+" "+levunits+" for "+tname)

	# Add a colorbar to the bottom of the plot.
    cbarstr = var2dstr
    if var2dstr.find('Count') == -1:
	# is Count
      del cbarstr
      cbarstr = "Number Density of Matched Obs"
    cbar = fig.colorbar(pp, label=cbarstr, orientation='vertical')#, pad=0.03)
    cbar.ax.tick_params(labelsize=fontlegend)

    return fig

# -------------------------------------------------------------------------
# Rotating Map of Locations of Matched (Collocated) Observations
#
#	Orthographic Projection
#
def contour2d_orthomap_rotate(datein,outpath,outname,var2d,var2dstr,xvar,yvar,xvarstr,yvarstr,tname,units,center_lat,level,levunits,**kwargs):

	# central points
    loninc = 2 #1
    center_lons = [*range(0,360,loninc)]		# should unpack all values from 0 to (360-1)=359 with an increment of 1
    size_lons = np.size(center_lons)

	# create directory to store all images for gif creation
    gifdir 	= "gif_images_2Dortho_"+str(datein)+"/"
    outpath_gif = outpath+gifdir

    os.system('if [ -d '+outpath_gif+' ]; then rm -Rf '+outpath_gif+'; fi')		# remove old gif directory before continuing
    os.system("mkdir -m 775 -p "+outpath_gif)						# make new empty gif directory
    
    tgifname 	= "image_"

	# create ortho plot for each central point
		# create ROTATE images
    for i in range(size_lons):
		# make each plot
      contour2d_orthomap(datein,var2d,var2dstr,xvar,yvar,xvarstr,yvarstr,tname,units,center_lons[i],center_lat,level,levunits)

		# assign number to each plot
      if (i+1) <= 100:
        gifnum = "0"+str(i)
        if (i+1) <= 10:
          gifnum = "00"+str(i)
      else:
        gifnum = str(i)

    		# save plot
      plt.savefig(outpath_gif+tgifname+gifnum+".png")

    del size_lons

	# make gif
    delay_sec = 10 				# time to display each image
    infiles   = outpath_gif+"image_*.png"	# all images to use to create gif
    outfile   = outname+".gif"			# filename for output gif

    		# linux command to create gif
    cmd = "convert -delay "+str(delay_sec)+" -loop 0 "+str(infiles)+" "+str(outfile)
    		# run 'cmd'
    os.system(cmd)

    plt.close("all")

# -------------------------------------------------------------------------
# Density Scatter Plot
#
def density_scatter(tx,ty,x_name,y_name,units,regionstr,ax=None,**kwargs):

        # create plot
    fig = plt.figure(figsize=(10,8.5))

    if ax is None:
      ax = plt.subplot()

    x = tx
    y = ty

        # get colormap 'cmap'
    ppcmap = plt.cm.get_cmap("jet")

        # set plot specs
    if units=="m/s":
      if x_name.find('Aeolus') != -1 or y_name.find('Aeolus') != -1:
        label = regionstr+" HLOS Wind Velocity"
        axismin = -100.0
        axismax = 100.0
        txpos   = -95
        typos   = 95
        diffpos = 10
      else:
        label = regionstr+" Wind Speed"
        axismin = 0.0
        axismax = 100.0
        txpos   = 5
        typos   = 95
        diffpos = 5
    elif units=="hPa":
      label = regionstr+" Pressure"
      axismin = 0.0
      axismax = 1000.0
      txpos   = 980
      typos   = 50
      diffpos = -50
    elif units=="km":
      label = regionstr+" Height"
      axismin = 0.0
      axismax = 20.0
      txpos   = 1
      typos   = 19
      diffpos = 1

        # conform x,y to 2D array
    gridcellsize = 1
    xaxis = list(np.arange(axismin,axismax,gridcellsize))
    yaxis = xaxis
    var2d = var_to_2d_counts_1dset(x_name,y_name,xaxis,yaxis,x,y)

        # set range of plotted colored contours
    warnings.filterwarnings('ignore', category=RuntimeWarning)          # ignore warnings
    varmax = np.nanmax(var2d)
    varmin = np.nanmin(var2d)
        # reassign color range max
    pltmin = 1
    pltmax = 100
    if varmax!=np.nan and varmax>=100000  : pltmax=1000000
    if varmax!=np.nan and varmax<100000: pltmax=10000
    if varmax!=np.nan and varmax<10000: pltmax=1000
    if varmax!=np.nan and varmax<1000: pltmax=100
    if varmax!=np.nan and varmax<100: pltmax=10
    var2d  = np.where(var2d==np.nan,0,var2d)
    del varmax,varmin

	# plot
    originstr = 'lower'
    if pltmax==10:
      pp = plt.imshow(var2d,cmap=ppcmap,vmin=pltmin,vmax=pltmax, extent=[axismin,axismax,axismin,axismax], origin=originstr, aspect='auto')
    else:
      pp = plt.imshow(var2d,cmap=ppcmap,norm=LogNorm(vmin=pltmin,vmax=pltmax), extent=[axismin,axismax,axismin,axismax], origin=originstr, aspect='auto')

        # statistical significance
    signif = 95.0                                      # statistical signifiance in percent (%). Type: float

    tsigmax = (100.0-signif)/100.0                     # max p-value allowed for statistical signiificance
    tsig = ttest_ind(x,y).pvalue                       # get p-value from Student's t-test
    if tsig<=tsigmax:                                  # if p-value < tsigmax (difference is statistically significant at 95% level)
      statsigstr = "Diffs are signif. at "+str(int(signif))+"%"

        # set axis limits
    ax.set_xlim([axismin,axismax])
    ax.set_ylim([axismin,axismax])
    ax.tick_params(labelsize=fontaxis+2)


    ax.axhline(y=0,color="black")             # horizontal line
    ax.axvline(x=0,color="black")             # vertical line
    ax.axline((0,0),slope=1.0,color="black")  # one-to-one line

        # add text inside the plot
    xpos = txpos                                      # value on x-axis where text will begin
    ypos = typos                                      # value on y-axis where text will begin
    ftsz = 12                                         # font size of text_* (see below)

        # add obs count to dataset name for legend
    legendlabel = str(y_name)+" count = "+str(np.size(y))
    ypos = ypos - diffpos
    plt.text(xpos, ypos, legendlabel, fontsize = ftsz)

        # compute and print stats of differences
    diff = y - x
    tcorr = np.corrcoef(x,y)     # correlation
    corr = tcorr[0,1]
    avg  = np.mean(diff)                  # mean
    sd   = np.std(diff)                   # standard deviation
    rmsd = np.sqrt(np.mean(diff**2))      # RMSD
                # add title for difference stats
    text = "Difference Stats ("+str(y_name)+" - DRIVER):"
    ypos = ypos - diffpos
    plt.text(xpos, ypos, text, fontsize = ftsz)
                # add correlation
    textr = "r = "+str(np.round_(corr,decimals=2))
    ypos = ypos - diffpos/2
    plt.text(xpos, ypos, textr, fontsize = ftsz)
                # add mean diff
    textm = "Mean_Diff = "+str(np.round_(avg,decimals=2))
    ypos = ypos - diffpos/2
    plt.text(xpos, ypos, textm, fontsize = ftsz)
                # add SD of diff
    texts = "SD_Diff = "+str(np.round_(sd,decimals=2))
    ypos = ypos - diffpos/2
    plt.text(xpos, ypos, texts, fontsize = ftsz)
                # add RMSD
    texts = "RMSD = "+str(np.round_(rmsd,decimals=2))
    ypos = ypos - diffpos/2
    plt.text(xpos, ypos, texts, fontsize = ftsz)
    if tsig<=tsigmax:
                # add space
      texts = " "
      ypos = ypos - diffpos/2
      plt.text(xpos,ypos,texts,fontsize=ftsz)
                # state whether differences are statistcially significant
      texts = statsigstr
      ypos = ypos - diffpos/2
      plt.text(xpos,ypos,texts,fontsize=ftsz)
      del tsig,tsigmax,signif
    del diff,tcorr,corr,avg,sd,rmsd,text,textr,textm,texts

      # plot title and axis labels
    ax.set_title(label+": "+str(y_name)+" vs DRIVER", fontsize=fonttitle)
    ax.set_xlabel("DRIVER ("+x_name+") ("+str(units)+")", fontsize=fontaxis+2)
    ax.set_ylabel(y_name+" ("+str(units)+")", fontsize=fontaxis+2)

        # Add a colorbar to the bottom of the plot.
    cbarstr = "Number Density of Matched Obs"
    frac = 0.046
    cbar = fig.colorbar(pp, label=cbarstr, orientation='vertical', fraction=frac, pad=0.04)
    cbar.ax.tick_params(labelsize=fontlegend+2)

    return fig

# -------------------------------------------------------------------------
# Histogram of Collocation Differences
#
def hist_diffs(nDuniq_list,Dlat,Dlon,tmatch,x_name,tname,alphaval,match_str,dcolors,units,**kwargs):

	# create plot
    fig = plt.figure(figsize=(9,8.5))
    ax1 = plt.subplot()
	
	# histogram specs
    if match_str.find("Time") != -1:
      label    = match_str+" Diffs"
      hist_dir = "vertical"
      histmin  = -100
      histmax  = 100
      binsize  = 10
      bins   = int((2*int(histmax))/binsize)
      xpos = histmin 
      posdiff = binsize
    elif match_str.find("Pressure") != -1:
      label  = match_str+" Diffs"
      hist_dir = "horizontal"
      histmin  = -50
      histmax  = 50			# units = hPa
      binsize  = 5
      bins   = int((2*int(histmax))/binsize)
      ypos = histmin
      posdiff = -1*binsize
    elif match_str.find("Height") != -1:
      label  = match_str+" Diffs"
      hist_dir = "horizontal"
      histmin  = -2000
      histmax  = 2000			# units = m
      binsize  = 100
      bins   = int((2*int(histmax))/binsize)
      ypos = histmin
      posdiff = binsize
    elif match_str.find("Distances") != -1:
      label    = match_str
      hist_dir = "vertical"
      histmin  = 0
      histmax  = 100			# units = km
      binsize  = 10
      bins     = int((2*int(histmax))/binsize)
      xpos = histmin 
      posdiff = binsize

    if hist_dir == "vertical":
      ylabel = label+" ("+units+")"
      xlabel = "Number of Matched Pairs"
    elif hist_dir == "horizontal":
      xlabel = label+" ("+units+")"
      ylabel = "Number of Matched Pairs"

    trange = (histmin,histmax) 
    w      = 0.6			# width of bars

	# count unique DRIVER points
    xstat = []
    ystat = []
    for i in range(np.size(dcolors)):
      txx = Dlon[i]
      tyy = Dlat[i]
	    
      x   = txx
      y   = tyy
      xstat.append(x)
      ystat.append(y)
      del txx,tyy

 	# get counts per dependent dataset for legend labels
    legendlabel = np.nan * np.ones(np.size(tname), dtype=list)
    for j in range(np.size(tname)):
      avg  = np.mean(tmatch[j])             # mean
      sd   = np.std(tmatch[j])              # standard deviation
      rmsd = np.sqrt(np.mean(tmatch[j]**2))  # RMSD
      text = "Mean_Diff="+str(np.round_(avg,decimals=2))+" SD_Diff="+str(np.round_(sd,decimals=2))+" RMSD="+str(np.round_(rmsd,decimals=2))
      legendlabel[j] = str(tname[j])+" | count="+str(np.size(tmatch[j]))+" "+text
      del avg,sd,rmsd

	# plot
    y,x,_ = ax1.hist(tmatch, bins, trange, color=dcolors[0:np.size(tname)], histtype='bar', rwidth=w, alpha=alphaval, label=legendlabel, orientation=hist_dir)
    ax1.set_ylabel(xlabel, fontsize=fontaxis)
    ax1.set_xlabel(ylabel,  fontsize=fontaxis)

	# add DRIVER counts to dataset name for legend
    if match_str.find("Height") != -1 or match_str.find("Pressure") != -1:
      xpos = y.max() - y.max()/2
      legendloc = 'upper right'
    else:
      ypos = y.max() - y.max()/4
      legendloc = 'upper left'
    Dlegendlabel = str(x_name)+" count = "+str(nDuniq_list)
    ftsz = 10                                                   # font size of text
    t = plt.text(xpos,ypos, Dlegendlabel, fontsize = ftsz)
    t.set_bbox(dict(facecolor='white', alpha=0.75)) 

    	# plot title
    plt.title("Histogram of "+label+" relative to "+x_name+" Obs", fontsize=fonttitle/2)

    	# add legend to plot
	#	'bbox_to_anchor' = location of legend box. (x,y) = (left=0 and right=1, bottom=0 and top=1)
	# 	'center'         = center of bounding box at coords 'bbox_to_anchor'
	# 	'centerleft '    = center of left edge of bounding box at coords 'bbox_to_anchor'
	# 	'center right'   = center of right edge of bounding box at coords 'bbox_to_anchor'
    leg = plt.legend(loc=legendloc, prop={'size': fontlegend})

    for lh in leg.legendHandles: 
      lh.set_alpha(1)

    legendmarkersize = 30
    for lh in range(np.size(leg.legendHandles)):
      leg.legendHandles[lh]._sizes = [legendmarkersize]

    return fig

# -------------------------------------------------------------------------
# Map of Locations of Matched (Collocated) Observations
#
#	Cylindrical Equidistant Projection
#
def map_locations_ce(nDuniq_list,Ds,ss,Dx,tx,Dy,ty,x_name,y_name,region_flag,marksize,imark,alphaval,dcolors,datein,**kwargs):

	# create plot
    fig = plt.figure(figsize=(11,8.5))

	# create map
    proj = ccrs.PlateCarree()
    ax = plt.subplot(1,1,1,projection=proj)

    ax.coastlines()
    
    gl = ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlines = True		# plot longitudes
    gl.ylines = True		# plot latitudes
    gl.top_labels=None
    gl.right_labels=None

    if region_flag==0:
	# set GLOBAL domain limits
      latmin = -90.0
      latmax = 90.0
      lonmin = -180.0
      lonmax = 180.0
      clat   = (latmax + latmin) / 2
      clon   = (lonmax + lonmin) / 2
      ax.set_extent([lonmin,lonmax,latmin,latmax], crs=proj)
      	# define x-axis tickmarks
      gl.xlocator = mticker.FixedLocator([-180,-120,-60,0,60,120,180])
	# define x-axis tickmarks
      gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])
		
    gl.xlabel_style = {'size': 15, 'color': 'gray'}
    gl.ylabel_style = {'size': 15, 'color': 'gray'}

    	# plot DRIVER points
    Dlons = []
    Dlats = []
    Dverts = []
    for i in range(np.size(y_name)):
      tDs = Ds[i]
      tss = ss[i]
      
      tmp_x = Dx[i]
      tmp_y = Dy[i]

      Dlons  = np.append(Dlons ,tmp_x,axis=0)
      Dlats  = np.append(Dlats ,tmp_y,axis=0)
      del tmp_x,tmp_y
      del tDs,tss  

    	# add obs counts to dataset name for legend
    Dlegendlabel = str(x_name)+", count = "+str(nDuniq_list)
    	# plot DRIVER on map
    ax.scatter(Dlons, Dlats, c="black", s=marksize+5, marker="o", alpha=0.75, label=Dlegendlabel)

	# plot DEPENDENT points
    for i in range(np.size(y_name)):
      tname = y_name[i]

      tDs = Ds[i]
      tss = ss[i]
      
      x   = tx[i]
      y   = ty[i]
              # add obs counts to dataset name for legend
      legendlabel  = str(tname)+", count = "+str(np.size(x)) #str(nuniq_list)

              # plot DEPENDENTS on map
      ax.scatter(x, y, c=dcolors[i], s=marksize, marker=imark[i], alpha=alphaval, label=legendlabel)
      del x,y,tDs,tss,legendlabel
   
    plt.title("Locations of Matched Obs for "+str(datein), fontsize=fonttitle)

	# add legend to plot
	#	'bbox_to_anchor' = location of legend box. (x,y) = (left=0 and right=1, bottom=0 and top=1)
	# 	'center'         = center of bounding box at coords 'bbox_to_anchor'
	# 	'centerleft '    = center of left edge of bounding box at coords 'bbox_to_anchor'
	# 	'center right'   = center of right edge of bounding box at coords 'bbox_to_anchor'
    leg = plt.legend(loc='lower right', prop={'size': 8})

    for lh in leg.legendHandles: 
      lh.set_alpha(1)

    legendmarkersize = 30
    for lh in range(np.size(leg.legendHandles)):
      leg.legendHandles[lh]._sizes = [legendmarkersize]

    return fig
    
# -------------------------------------------------------------------------
# Map of Locations of Matched (Collocated) Observations
#
#	Cylindrical Equidistant Projection
#
def map_locations_ce_1dset(ss,tx,ty,y_name,region_flag,marksize,imark,alphaval,dcolors,datein,**kwargs):

	# create plot
    fig = plt.figure(figsize=(11,8.5))

	# create map
    proj = ccrs.PlateCarree()
    ax = plt.subplot(1,1,1,projection=proj)

    ax.coastlines()
    
    gl = ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlines = True		# plot longitudes
    gl.ylines = True		# plot latitudes
    gl.top_labels=None
    gl.right_labels=None

    if region_flag==0:
	# set GLOBAL domain limits
      latmin = -90.0
      latmax = 90.0
      lonmin = -180.0
      lonmax = 180.0
      clat   = (latmax + latmin) / 2
      clon   = (lonmax + lonmin) / 2
      ax.set_extent([lonmin,lonmax,latmin,latmax], crs=proj)
      	# define x-axis tickmarks
      gl.xlocator = mticker.FixedLocator([-180,-120,-60,0,60,120,180])
	# define x-axis tickmarks
      gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])
		
    gl.xlabel_style = {'size': 15, 'color': 'gray'}
    gl.ylabel_style = {'size': 15, 'color': 'gray'}

	# plot DEPENDENT points
    for i in range(np.size(y_name)):
      tname = y_name[i]

      tss = ss[i]
      
      x   = tx[i]
      y   = ty[i]
              # add obs counts to dataset name for legend
      legendlabel  = str(tname)+", count = "+str(np.size(x)) #str(nuniq_list)

              # plot DEPENDENTS on map
      if np.size(y_name)==1:
        ax.scatter(x, y, c=dcolors, s=marksize, marker=imark, alpha=alphaval, label=legendlabel)
      else:
        ax.scatter(x, y, c=dcolors[i], s=marksize, marker=imark[i], alpha=alphaval, label=legendlabel)
      del x,y,tss,legendlabel
   
    plt.title("Locations of Obs for "+str(datein), fontsize=fonttitle)

	# add legend to plot
	#	'bbox_to_anchor' = location of legend box. (x,y) = (left=0 and right=1, bottom=0 and top=1)
	# 	'center'         = center of bounding box at coords 'bbox_to_anchor'
	# 	'centerleft '    = center of left edge of bounding box at coords 'bbox_to_anchor'
	# 	'center right'   = center of right edge of bounding box at coords 'bbox_to_anchor'
    leg = plt.legend(loc='lower right', prop={'size': 8})

    for lh in leg.legendHandles: 
      lh.set_alpha(1)

    legendmarkersize = 30
    for lh in range(np.size(leg.legendHandles)):
      leg.legendHandles[lh]._sizes = [legendmarkersize]

    return fig
    
# -------------------------------------------------------------------------
# Map of Locations of Matched (Collocated) Observations
#
def map_points2d_ce(varstr,tx,ty,ss,xaxis,yaxis,x_name,marksize,alphaval,datein,units,level,levunits,opt,**kwargs):

	# create plot
    fig = plt.figure(figsize=(10,7.5))

	# create map
    proj = ccrs.PlateCarree()
    ax = plt.subplot(1,1,1,projection=proj)

    ax.coastlines()
    
    gl = ax.gridlines(draw_labels=True, linewidth=2, color='gray', alpha=0.5, linestyle='--')
    gl.xlines = True		# plot longitudes
    gl.ylines = True		# plot latitudes
    gl.top_labels=None
    gl.right_labels=None

    region_flag=0
    if region_flag==0:
	# set GLOBAL domain limits
      latmin = -90.0
      latmax = 90.0
      lonmin = -180.0
      lonmax = 180.0
      clat   = (latmax + latmin) / 2
      clon   = (lonmax + lonmin) / 2
      ax.set_extent([lonmin,lonmax,latmin,latmax], crs=proj)
      	# define x-axis tickmarks
      gl.xlocator = mticker.FixedLocator([-180,-120,-60,0,60,120,180])
	# define x-axis tickmarks
      gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])
		
    gl.xlabel_style = {'size': 15, 'color': 'gray'}
    gl.ylabel_style = {'size': 15, 'color': 'gray'}

    warnings.filterwarnings('ignore', category=RuntimeWarning)		# ignore warnings
    varmax = np.nanmax(abs(ss))
    varmin = np.nanmin(ss)

    if varstr.find('SD') != -1:
        ppcmap = "jet"
        pltmin = 1
        pltmax = 100
    elif varstr.find('Diff') != -1:
        ppcmap = "seismic"
        pltmin = -30
        pltmax = 30
    else:
        ppcmap = "jet"
        pltmin = -50
        pltmax = 50
    # reassign color range min/max
    if varmax!=np.nan and varmax>=60: pltmax=100
    if varmax!=np.nan and varmax<60: pltmax=50
    if varmax!=np.nan and varmax<50: pltmax=40
    if varmax!=np.nan and varmax<40: pltmax=30
    if varmax!=np.nan and varmax<30: pltmax=20
    if varmax!=np.nan and varmax<20: pltmax=10
    # reassign pltmin if not SD plot
    if varstr.find('SD') == -1: pltmin = -pltmax
    if varstr.find('Diff') == -1 and varmin>=0: pltmin=1
    del varmax,varmin

	# add obs counts to dataset name for legend
    Dlegendlabel = str(x_name)

    	# plot DRIVER points
    Dlons  = tx
    Dlats  = ty

    	# plot DRIVER on map
    sc = ax.scatter(Dlons, Dlats, c=ss, cmap=ppcmap, vmin=pltmin, vmax=pltmax, s=marksize, marker="o", label=Dlegendlabel, edgecolors='none')

    plt.title(x_name+" "+varstr+" ("+str(units)+") at "+str(level)+" "+str(levunits)+" for "+str(datein), fontsize=9)
    
        # Add a colorbar to the bottom of the plot.
    cbar = plt.colorbar(sc, orientation='horizontal', pad=0.1)
    cbar.set_label(varstr+" ("+str(units)+")")

    return fig

# -------------------------------------------------------------------------
# 2D Contour Map of Matched (Collocated) Observations
#
#	Cylindrical Equidistant Projection
#
def map_contour2d_ce(var2dstr,var,flons,flats,tname,alphaval,datein,units,level,levunits,opt,**kwargs):

	# create plot
    fig = plt.figure(figsize=(12,10))

	# create map
    proj = ccrs.PlateCarree()
    ax = plt.subplot(1,1,1,projection=proj)

    if opt==1:
        # plot difference, so use Blue/Red colors scheme
      ppcmap = plt.cm.seismic
    else:
        # use rainbow color scheme
      ppcmap = plt.cm.jet

        # set range of plotted colored contours
    if var2dstr.find('Diff') != -1:
      pltmin = -30
      pltmax = 30
    else:
      pltmin = -50
      pltmax = 50
    stride = 5
    clevels = np.arange(pltmin,pltmax+stride,stride)

    x,y = np.meshgrid(flons,flats)
    var_masked = np.ma.masked_where(var==np.nan, var)

    pp = plt.pcolormesh(x, y, var, cmap=ppcmap, transform=proj, zorder=3)

    ax.coastlines()
    
    gl = ax.gridlines(draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='--')
    gl.xlines = True		# plot longitudes
    gl.ylines = True		# plot latitudes
    gl.top_labels=None
    gl.right_labels=None

	# set GLOBAL domain limits
    latmin = -90.0
    latmax = 90.0
    lonmin = -180.0
    lonmax = 180.0 #179.0
    clat   = (latmax + latmin) / 2
    clon   = (lonmax + lonmin) / 2
    ax.set_extent([lonmin,lonmax,latmin,latmax], crs=proj)
      # define x-axis tickmarks
    gl.xlocator = mticker.FixedLocator([-180,-120,-60,0,60,120,180])
      # define x-axis tickmarks
    gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])
		
    gl.xlabel_style = {'size': 15, 'color': 'gray'}
    gl.ylabel_style = {'size': 15, 'color': 'gray'}

    plt.title(str(tname)+" "+str(var2dstr)+" at "+str(np.round_(level,decimals=2))+" "+str(levunits)+" for "+str(datein), fontsize=fonttitle)
    
    	# Add a colorbar to the bottom of the plot.
    fig.colorbar(pp, label=str(var2dstr)+" ("+str(units)+")", orientation='horizontal', pad=0.1)
 
    return fig

# -------------------------------------------------------------------------
# Map of Locations of Matched (Collocated) Observations
#
#	Orthographic Projection
#
def map_locations_ortho(nDuniq_list,Ds,ss,Dx,tx,Dy,ty,x_name,y_name,marksize,imark,alphaval,dcolors,datein,center_lon,center_lat,**kwargs):

	# create plot
    fig = plt.figure(figsize=(8.5,8.5))

	# create map
    proj = ccrs.Orthographic(central_longitude=center_lon,central_latitude=center_lat,globe=None)
    geo  = ccrs.Geodetic()
    ax = plt.axes(projection=proj)
    
    	# add coastlines
    ax.coastlines(resolution='50m')
    ax.set_global()
    ax.gridlines()

	# plot DRIVER points
    Dlatlonlist = []
    Dlons = []
    Dlats = []
    for i in range(np.size(y_name)):
      tname = y_name[i]
      
      tDs = Ds[i]
      tss = ss[i]
      
      	# Driver
      tmp_x = Dx[i]
      tmp_y = Dy[i]
      	# Dependent
      idx_x = tx[i]
      idx_y = ty[i]

		# transform DRIVER lat/lon points to orthographic points
      x_D = np.asarray(tmp_x)
      y_D = np.asarray(tmp_y)
      points = proj.transform_points(geo, x_D, y_D)
      del x_D,y_D
      x_D = points[:,0]
      y_D = points[:,1]
      del points

      del tmp_x,tmp_y
      
      if i==0:
      		# add obs counts to dataset name for legend
        Dlegendlabel = x_name+", count = "+str(nDuniq_list)
      		# plot DRIVER on map
        ax.scatter(x_D, y_D, c="black", s=marksize+5, marker="o", alpha=0.75, label=Dlegendlabel)
      else:
      		# plot DRIVER on map
        ax.scatter(x_D, y_D, c="black", s=marksize+5, marker="o", alpha=0.75)
	
      del x_D,y_D
      
      		# transform DEPENDENT lat/lon points to orthographic points
      x = np.asarray(idx_x)
      y = np.asarray(idx_y)
      points  = proj.transform_points(geo, x, y)
      del idx_x,idx_y
			# add obs counts to dataset name for legend
      legendlabel  = str(tname)+", count = "+str(np.size(x)) 
      del x,y

      ttx = points[:,0]
      tty = points[:,1]
      x   = ttx
      y   = tty
              # plot DEPENDENT on map
      ax.scatter(x, y, c=dcolors[i], s=marksize, marker=imark[i], alpha=alphaval, label=legendlabel)
      del x,y,ttx,tty,legendlabel
      
      del tDs,tss

    plt.title("Locations of Matched Obs for "+str(datein), fontsize=fonttitle)

	# add legend to plot
	#	'bbox_to_anchor' = location of legend box. (x,y) = (left=0 and right=1, bottom=0 and top=1)
	# 	'center'         = center of bounding box at coords 'bbox_to_anchor'
	# 	'centerleft '    = center of left edge of bounding box at coords 'bbox_to_anchor'
	# 	'center right'   = center of right edge of bounding box at coords 'bbox_to_anchor'
    leg = plt.legend(loc='center',bbox_to_anchor=(0.8,0))

    for lh in leg.legendHandles: 
      lh.set_alpha(1)

    legendmarkersize = 30
    for lh in range(np.size(leg.legendHandles)):
      leg.legendHandles[lh]._sizes = [legendmarkersize]

    return fig
    
# -------------------------------------------------------------------------
# Rotating Map of Locations of Matched (Collocated) Observations
#
#	Orthographic Projection
#
def map_locations_ortho_rotate(outpath,outname,nDuniq_list,Ds,ss,Dx,tx,Dy,ty,x_name,y_name,marksize,imark,alphaval,dcolors,datein,center_lat,**kwargs):

	# central points
    center_lons = [*range(0,360,10)]		# unpacks all values from 0 to (360-1)=359 with an increment of 1
    size_lons = np.size(center_lons)

	# create directory to store all images for gif creation
    gifdir 	= "gif_images"+str(datein)+"/"
    outpath_gif = outpath+gifdir

    os.system('if [ -d '+outpath_gif+' ]; then rm -Rf '+outpath_gif+'; fi')		# remove old gif directory before continuing
    os.system("mkdir -m 775 -p "+outpath_gif)						# make new empty gif directory
    
    tgifname 	= "image_"

	# create ortho plot for each central point
    		# create ROTATE images
    for i in range(size_lons):
		# make each plot
      map_locations_ortho(nDuniq_list,Ds,ss,Dx,tx,Dy,ty,x_name,y_name,marksize,imark,alphaval,dcolors,datein,center_lons[i],center_lat)

		# assign number to each plot
      if (i+1) <= 100:
        gifnum = "0"+str(i)
        if (i+1) <= 10:
          gifnum = "00"+str(i)
      else:
        gifnum = str(i)

    		# save plot
      plt.savefig(outpath_gif+tgifname+gifnum+".png")

    del size_lons

	# make gif
    delay_sec = 60 #10 				# time to display each image
    infiles   = outpath_gif+"image_*.png"	# all images to use to create gif
    outfile   = outname+".gif"			# filename for output gif

    		# linux command to create gif
    cmd = "convert -delay "+str(delay_sec)+" -loop 0 "+str(infiles)+" "+str(outfile)
    		# run 'cmd'
    os.system(cmd)

    plt.close("all")
    
# -------------------------------------------------------------------------
# 3D Map of Matched (Collocated) Observation Profiles
#	
#	Cylindrical Equidistant Projection
#
def map_3d_profile(nDuniq_list,ss,xx,yy,zz,x_name,tname,dset_plotted,dset_plotted_name,datein,sslabel,zzlabel,**kwargs):

	# longitude check: set range to 0=360
    for iP in range(np.size(xx)):
      if xx[iP] > 180.0:
        xx[iP] = xx[iP] - 360.

	# create plot
    fig = plt.figure(figsize=(11.5,8.5))

	# create map
    proj = "3d"
    ax = fig.gca(projection=proj)
    
    	# Define horizontal extent (lower left, upper right lontitude and lattitude respectively)
    extent = [-180, 180, -90, 90] 
    	# Create a basemap instance that draws the Earth layer
    bm = Basemap(llcrnrlon=extent[0], llcrnrlat=extent[2],
             urcrnrlon=extent[1], urcrnrlat=extent[3],
             projection='cyl', resolution='l', fix_aspect=False, ax=ax)

	# set colormap
    warnings.filterwarnings('ignore', category=RuntimeWarning)		# ignore warnings
    varmax = np.nanmax(abs(ss))
    varmin = np.nanmin(ss)

    if sslabel.find('Diff') != -1:
      ppcmap = "seismic"
      pltmin = -30
      pltmax = 30
        # reassign color range min/max
      if varmax!=np.nan and varmax>=60: pltmax=100
      if varmax!=np.nan and varmax<60: pltmax=50
      if varmax!=np.nan and varmax<50: pltmax=40
      if varmax!=np.nan and varmax<40: pltmax=30
      if varmax!=np.nan and varmax<30: pltmax=20
      if varmax!=np.nan and varmax<20: pltmax=10
      pltmin = -pltmax
    else:
      ppcmap = "jet"
      pltmin = -50
      pltmax = 50
      	# reassign color range max
      if varmax!=np.nan and varmax>=60: pltmax=100; pltmin=-100
      if varmax!=np.nan and varmax<60: pltmax=50; pltmin=-50
      if varmax!=np.nan and varmax<50: pltmax=40; pltmin=-40
      if varmax!=np.nan and varmax<40: pltmax=30; pltmin=-30
      if varmax!=np.nan and varmax<30: pltmax=20; pltmin=-20
      if varmax!=np.nan and varmax<20: pltmax=10; pltmin=-10
      if varmin>=0: pltmin=1
    del varmax,varmin

    	# set z-axis range
    if zzlabel=="Pressure":
      zs_val = 1000
      ax.set_zlim(0., 1000.)
      plt.gca().invert_zaxis()
    elif zzlabel=="Height":
      zs_val = 0
      ax.set_zlim(0.,20.)
    
	# Add Basemap to the figure
    ax.add_collection3d(bm.drawcoastlines(linewidth=0.25),zs=zs_val)
    ax.add_collection3d(bm.drawcountries(linewidth=0.35),zs=zs_val)

	# Add meridian and parallel gridlines
    ax.set_xlabel('Longitude', fontsize=fontaxis-2,  labelpad=20)
    ax.set_ylabel('Latitude', fontsize=fontaxis-2, labelpad=20)
    if zzlabel=="Pressure":
      zzunits = "(hPa)"
    elif zzlabel=="Height":
      zzunits = "(km)"
    ax.set_zlabel(zzlabel+" "+zzunits, fontsize=fontaxis-2, labelpad=20)

    lon_step = 30
    lat_step = 30
    meridians = np.arange(extent[0], extent[1] + lon_step, lon_step)
    parallels = np.arange(extent[2], extent[3] + lat_step, lat_step)
    ax.set_yticks(parallels)
    ax.set_yticklabels(parallels)
    ax.set_xticks(meridians)
    ax.set_xticklabels(meridians)

    	# 3D viewpoint
    azim_pt = 240 	# 230=lower left corner is center; less than 230 is to the right, more than 230 is to the left
    elev_pt = 20  
    ax.view_init(azim=azim_pt, elev=elev_pt)

        # add obs counts to dataset name for legend
    ftsz = 10                                             # font size of text
    if dset_plotted.find("DRV")!=-1:
		# plot count for DRIVER
      Dlegendlabel = "DRIVER ("+str(x_name)+") count = "+str(nDuniq_list)
      ax.text2D(0.05, 0.95, Dlegendlabel, transform=ax.transAxes)
    if dset_plotted.find("DEP")!=-1:  
		# plot count for DEPENDENT
      legendlabel = str(tname)+" count = "+str(np.size(yy))
      ax.text2D(0.05, 0.90, legendlabel, transform=ax.transAxes)

    	# add scatterplot map based on lons 'xx', lats 'yy', heights 'zz', and colors based on 'ss'
    p = ax.scatter(xx, yy, zz, c=ss, cmap=ppcmap, vmin=pltmin, vmax=pltmax, s=5)#, alpha=0.5)

    	# add colorbar to reference intensity of 'ss'
    cbar = fig.colorbar(p, label=sslabel)
    cbar.ax.tick_params(labelsize=fontlegend)
    
    #plt.title(dset_plotted_name+" "+sslabel+" for "+str(datein), fontsize=9.5)
    plt.title(dset_plotted_name+" "+sslabel, fontsize=fonttitle-2)
   
    return fig
    
# -------------------------------------------------------------------------
# Scatter Plot (not density) of Matched Observations
#
def scatter_matches(nDuniq_list,Dlat,Dlon,Dvert,tx,ty,x_name,tname,marksize,alphaval,match_str,acolors,units,imark,regionstr,latmin,latmax,**kwargs):

	# create plot
    fig = plt.figure(figsize=(9,8.5))
    ax = plt.subplot()

	# scatter specs
    if match_str=="HLOS":
        label = regionstr+" HLOS Wind Velocity"
        axismin = -100.0
        axismax = 100.0
        diffpos = 5
    elif match_str=="Wind Speed":
        label = regionstr+" Wind Speed"
        axismin = 0.0
        axismax = 100.0
        diffpos = 5
    elif match_str=="Pressure":
        label = regionstr+" "+match_str
        axismin = 0.0
        axismax = 1000.0
        diffpos = -50
    elif match_str=="Height":
        label = regionstr+" "+match_str
        axismin = 0.0
        axismax = 20.0
        diffpos = 1

	# axis limits
    ax.set_xlim([axismin,axismax])
    ax.set_ylim([axismin,axismax])

    if match_str=="Pressure":
      xpos = axismax-25
      ypos = axismin+25
    elif match_str=="Height":
      xpos = axismin+0.5
      ypos = axismax-0.5
    else:
      xpos = axismin+5
      ypos = axismax-5

    	# get UNIQUE counts for DRIVER dataset (for legend labels)
    Dlatlonlist = []
    xstat = []
    ystat = []
    for i in range(np.size(acolors)):
      txx = tx[i]
      tyy = ty[i]
      tlat = Dlat[i]
      idx = np.where((tlat>=latmin)*(tlat<latmax))
      x   = txx[idx]
      y   = tyy[idx]
      del txx,tyy,tlat,idx

          # adding text inside the plot
      if match_str=="Pressure":
        xpos = xpos
        ypos = ypos + 50
      elif match_str=="Height":
        xpos = xpos
        ypos = ypos - 1
      else:
        xpos = xpos						# value on x-axis where text will begin
        ypos = ypos						# value on y-axis where text will begin

              # compute and print stats of differences
      ftsz = 7
      diff = y - x
      tcorr = np.corrcoef(x,y)	 # correlation
      corr = tcorr[0,1]
      avg  = np.mean(diff)		  # mean
      sd   = np.std(diff) 		  # standard deviation
      rmsd = np.sqrt(np.mean(diff**2))	  # RMSD
      text = "Difference Stats ("+str(tname[i])+"-"+str(x_name)+"):"
      ypos = ypos - diffpos
      plt.text(xpos, ypos, text, fontsize = ftsz)
      textr = "r = "+str(np.round_(corr,decimals=2))
      ypos = ypos - diffpos/2
      plt.text(xpos, ypos, textr, fontsize = ftsz)
      textm = "Mean_Diff = "+str(np.round_(avg,decimals=2))
      ypos = ypos - diffpos/2
      plt.text(xpos, ypos, textm, fontsize = ftsz)
      texts = "SD_Diff = "+str(np.round_(sd,decimals=2))
      ypos = ypos - diffpos/2
      plt.text(xpos, ypos, texts, fontsize = ftsz)
      texts = "RMSD = "+str(np.round_(rmsd,decimals=2))
      ypos = ypos - diffpos/2
      plt.text(xpos, ypos, texts, fontsize = ftsz)
      del diff,tcorr,corr,avg,sd,rmsd,text,textr,textm,texts

	    # legend label
      legendlabel = str(tname[i])+" | count = "+str(np.size(x))

          # plot scatter
      plt.scatter(x, y, c=acolors[i], marker=imark[i], edgecolor="black", s=marksize, alpha=alphaval, label=legendlabel)

      del x,y

	# add obs counts to dataset name for legend
    Dlegendlabel = "DRIVER ("+str(x_name)+") count = "+str(nDuniq_list)
    ftsz = 8                                                   # font size of text
    plt.text(xpos, ypos-diffpos, Dlegendlabel, fontsize = ftsz)

	# plot labels
    plt.axhline(y=0,color="black")  	      	  # horizontal line
    plt.axvline(x=0,color="black")  	      	  # vertical line
    plt.axline((0,0),slope=1.0,color="black")     # one-to-one line

    plt.title(regionstr+" "+label+" Scatterplot: DEPENDENTS vs DRIVER", fontsize=fonttitle)
    plt.xlabel(x_name+" ("+units+")", fontsize=fontaxis)
    plt.ylabel("DEPENDENT dataset(s)", fontsize=fontaxis)

    if match_str=="Pressure":
      plt.gca().invert_xaxis()
      plt.gca().invert_yaxis()

      	# add legend to plot
	#	'bbox_to_anchor' = location of legend box. (x,y) = (left=0 and right=1, bottom=0 and top=1)
	# 	'center'         = center of bounding box at coords 'bbox_to_anchor'
	# 	'centerleft '    = center of left edge of bounding box at coords 'bbox_to_anchor'
	# 	'center right'   = center of right edge of bounding box at coords 'bbox_to_anchor'
    leg = plt.legend(loc='lower right', prop={'size': fontlegend})

    for lh in leg.legendHandles: 
      lh.set_alpha(1)

    legendmarkersize = 30
    for lh in range(np.size(leg.legendHandles)):
      leg.legendHandles[lh]._sizes = [legendmarkersize]

    return fig

# -------------------------------------------------------------------------
# Time Series of Matched Observations
#
def time_series(nDuniq_list,date_uniq,Dyyyymmdd,Dlat,Dlon,Dvert,tx,ty,x_name,tname,match_str,acolors,units,regionstr,latmin,latmax,**kwargs):
	
	# set font sizes for plotting
    ftsz = 9
    legendmarkersize = ftsz
    
	# create plot
    fig = plt.figure(figsize=(12,12))
    
    	# create grid of panels
    gs = fig.add_gridspec(nrows=4, ncols=1)
    
    	# create each panel
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[1,0])
    ax3 = fig.add_subplot(gs[2,0])
    ax4 = fig.add_subplot(gs[3,0])
   
	# panel specs
    if match_str=="HLOS":
      label = "HLOS Wind Velocity"
      labelyaxis = "Wind"
      axismin = -5.0
      axismax = 5.0
      diffpos = 0.2
    elif match_str=="Wind Speed":
      label = match_str
      labelyaxis = "Wind"
      axismin = -10.0
      axismax = 10.0
      diffpos = 5
    elif match_str=="Pressure":
      label = match_str
      labelyaxis = label
      axismin = -15.0
      axismax = 15.0
      diffpos = -5
    elif match_str=="Height":
      label = match_str
      labelyaxis = label
      axismin = -2.0
      axismax = 2.0
      diffpos = 0.5

	# set initial positions of plot text boxes
    if match_str=="Pressure":
      xpos = 0
      ypos = axismin
    elif match_str=="Height":
      xpos = 0
      ypos = axismin
    else:
      xpos = 0
      ypos = axismin

    warnings.filterwarnings('ignore', category=FutureWarning)		# ignore warnings
    datestr  = date_uniq	
    ndatestr = len(datestr)
    datelabels = datestr

	# set xaxis
    warnings.filterwarnings('ignore', category=RuntimeWarning)		# ignore warnings
    idx = np.where(Dyyyymmdd==datestr[ndatestr-1])
    if np.size(idx)<=0:
      ilast = ndatestr-1
    else:
      ilast = ndatestr
    del idx
      		# x-axis specs
    timestepstr = "days"	      # includes 00-23 hours
    	      # x-axis tickmark values and labels
    if ndatestr<10:
      stride  = 1
    elif ndatestr>=10 and ndatestr<30:
      stride = 5
    elif ndatestr>=30 and ndatestr<93:
      stride = 10
    elif ndatestr>=93: 
      stride = 30
    xaxis   = np.arange(0, ilast, 1)
    xvalues = np.arange(0, ilast, stride)
    xlabels = datelabels[0:ilast:stride]

    	# compute stats per DEPENDENT dataset
    xstat = []
    ystat = []
    nobsmax = 0
    for i in range(np.size(acolors)):
	# variable to plot from DRIVER (x) and DEPENDENT (y)
      txx = tx[i]
      tyy = ty[i]
      ttmp_Dyyyymmdd = Dyyyymmdd[i]
      tlat = Dlat[i]
      idx = np.where((tlat>=latmin)*(tlat<latmax))
      x = txx[idx]
      y = tyy[idx]
      tmp_Dyyyymmdd = ttmp_Dyyyymmdd[idx]
      del txx,tyy,tlat,idx,ttmp_Dyyyymmdd

        # adding text inside the plot
      if match_str=="Pressure":
        xpos = xpos
        ypos = ypos + 50
      elif match_str=="Height":
        xpos = xpos
        ypos = ypos - 1
      else:
        xpos = xpos						# value on x-axis where text will begin
        ypos = ypos						# value on y-axis where text will begin

	# loop thru all obs for dataset
      diff = []
      corr = []
      sd   = []
      rmsd = []
      nobs = []
      difftot = []
      sdtot   = []
      rmsdtot = []
      nobstot = []
      for j in range(ilast):
        warnings.filterwarnings('ignore', category=RuntimeWarning)		# ignore warnings
      		# extract vars based on time
        idx = np.where((x!=np.nan)*(y!=np.nan)*(tmp_Dyyyymmdd==datestr[j]))
        if np.size(idx)>0:
          varxdate = x[idx]
          varydate = y[idx]
          tnobs = np.size(varxdate)
          nobs.append(tnobs)
	  	# find absolute max of nobs
          if tnobs>nobsmax:
            nobsmax = tnobs
		#compute stats per timestep
          tdiff  = varydate - varxdate				# y-x = DEPENDENT minus DRIVER
          ttcorr = np.corrcoef(varxdate,varydate)		# correlation
          tcorr  = ttcorr[0,1]
          tavg   = np.mean(tdiff)				# mean
          tsd    = np.std(tdiff) 		  		# standard deviation
          trmsd  = np.sqrt(np.mean(tdiff**2))			# RMSD
		#append stats per timestep
          diff.append(tavg)
          corr.append(tcorr)
          sd.append(tsd)
          rmsd.append(trmsd)
		# stats for total means
          difftot.append(tavg)
          sdtot.append(tsd)
          rmsdtot.append(trmsd)

          del varxdate,varydate,tdiff,ttcorr,tcorr,tavg,tsd,trmsd,tnobs
        else:
          diff.append(np.nan)
          corr.append(np.nan)
          sd.append(np.nan)
          rmsd.append(np.nan)
          nobs.append(np.nan)
        del idx

      	# compute stats
      diffarr = np.asarray(diff)
      sdarr   = np.asarray(sd)
      rmsdarr = np.asarray(rmsd)

      sdifftot = np.asarray(difftot)
      ssdtot   = np.asarray(sdtot)
      srmsdtot = np.asarray(rmsdtot)

      xarr     = np.asarray(x)
      yarr     = np.asarray(y)
      tcorr    = np.corrcoef(xarr,yarr)	 	  # correlation
      corr_all = tcorr[0,1]

	# average stats over entire time period
      avg_all  = np.mean(sdifftot)
      sd_all   = np.mean(ssdtot)
      rmsd_all = np.mean(srmsdtot)

      legendlabel1 = str(tname[i])+" | Mean_Diff="+str(np.round_(avg_all,decimals=2))
      del tcorr,xarr,yarr
    
    	# plot 1: time series: difference
      ax1.plot(xaxis, diff, '-', color=acolors[i], linewidth=1, label=legendlabel1)
 
        # plot 2: correlation coefficient
      legendlabel2 = str(tname[i])+" | r="+str(np.round_(corr_all,decimals=2))
      ax2.plot(xaxis, corr, color=acolors[i], label=legendlabel2)

        # plot 3: standard deviation, RMSD
		# color lines legend
      legendlabel3 = str(tname[i])+" | SD_Diff="+str(np.round_(sd_all,decimals=2))+" RMSD="+str(np.round_(rmsd_all,decimals=2))

		# black lines legend
      custom_lines = [Line2D([0],[0],color="black",label="SD_Diff",linestyle="-"), Line2D([0],[0],color="black",label="RMSD",linestyle="--")]

      ax3.plot(xaxis, sd, color=acolors[i], label=legendlabel3)
      ax3.plot(xaxis, rmsd, color=acolors[i], linestyle='dashed')

      	# plot 4: number of observations
      legendlabel4 = str(tname[i])+" | count="+str(np.size(x))
      ax4.plot(xaxis, nobs, color=acolors[i], label=legendlabel4)		# add legend to bottom plot

      del x,y,sd,rmsd,tmp_Dyyyymmdd

	# set tick marks for all panels
    ax1.set_xticks(xvalues)
    ax2.set_xticks(xvalues)
    ax3.set_xticks(xvalues)
    ax4.set_xticks(xvalues)
        # set tick mark labels for bottom panel only
    ax1.set_xticklabels([])
    ax2.set_xticklabels([])
    ax3.set_xticklabels([])
    ndatestr = 5
    ax4.set_xticklabels(xlabels, rotation=45)
    
    	# add DRIVER obs count to plot 1 legend
    Dlegendlabel = "DRIVER Count = "+str(nDuniq_list)

	# axis limits
    nobsarr = np.asarray(nobs)

    ax1.set_xlim([min(xaxis),max(xaxis)])
    ax1.set_ylim([axismin,axismax])
    ax2.set_xlim([min(xaxis),max(xaxis)])
    ax2.set_ylim([0,1])
    ax3.set_xlim([min(xaxis),max(xaxis)])
    ax3.set_ylim([0,axismax*2])
    ax4.set_xlim([min(xaxis),max(xaxis)])
    ax4.set_ylim([0,nobsmax])

	# plot 1 labels
    ax1.axhline(y=0,color="black")  	      			# horizontal line
    ax1.axvline(x=0,color="black")  	      			# vertical line
    ax1.set_title("Daily Mean Time Series: "+regionstr+" "+label+" Difference (DEPENDENT - DRIVER)", fontsize=fonttitle-2)
    ax1.set_ylabel(labelyaxis+" Diff. ("+units+")", fontsize=10)
    
		# legend
    leg1 = ax1.legend(loc='lower right', prop={'size': ftsz})
    for lh in leg1.legendHandles:
      lh.set_alpha(1)
    for lh in range(np.size(leg1.legendHandles)):
      leg1.legendHandles[lh]._sizes = [legendmarkersize]

        # plot 2 labels
    ax2.axhline(y=0,color="black")                              # horizontal line
    ax2.axvline(x=0,color="black")                              # vertical line
    ax2.set_ylabel("Corr. Coeff. (r)", fontsize=10)

                # legend
    leg2 = ax2.legend(loc='lower right', prop={'size': ftsz})
    for lh in leg2.legendHandles:
      lh.set_alpha(1)
    for lh in range(np.size(leg2.legendHandles)):
      leg2.legendHandles[lh]._sizes = [legendmarkersize]

        # plot 3 labels
    ax3.axhline(y=0,color="black")                              # horizontal line
    ax3.axvline(x=0,color="black")                              # vertical line
    ax3.set_ylabel("RMSD, SD of Diff ("+units+")", fontsize=10)

		# color lines legend
    leg3a = ax3.legend(loc='lower right', prop={'size': ftsz})
    for lh in leg3a.legendHandles:
      lh.set_alpha(1)
    for lh in range(np.size(leg3a.legendHandles)):
      leg3a.legendHandles[lh]._sizes = [legendmarkersize]

                # black lines legend
    ax3b = ax3.twinx()
    leg3b = ax3b.legend(handles=custom_lines, loc='lower left', prop={'size': ftsz})
    for lh in leg3b.legendHandles:
      lh.set_alpha(1)
    for lh in range(np.size(leg3b.legendHandles)):
      leg3b.legendHandles[lh]._sizes = [legendmarkersize] 
    ax3b.tick_params(right=False,labelright=False)
 
	# plot 4 labels
    ax4.axhline(y=0,color="black")  	      			# horizontal line
    ax4.axvline(x=0,color="black")  	      			# vertical line
    ax4.set_ylabel("Count", fontsize=10)
    
		# legend
    leg4 = ax4.legend(loc='upper right', prop={'size': ftsz})
    for lh in leg4.legendHandles:
      lh.set_alpha(1)
    for lh in range(np.size(leg4.legendHandles)):
      leg4.legendHandles[lh]._sizes = [legendmarkersize]

		# DRIVER count legend
    ax4D = ax4.twinx()
    Dlegend = ax4D.plot([],[],' ',label=Dlegendlabel)
    leg4D = ax4D.legend(handles=Dlegend, loc='upper left', prop={'size': ftsz})
    for lh in leg4D.legendHandles:
      lh.set_alpha(1)
    for lh in range(np.size(leg4D.legendHandles)):
      leg4D.legendHandles[lh]._sizes = [legendmarkersize]
    ax4D.tick_params(right=False,labelright=False)

    del diff,nobs,diffarr,sdarr

    return fig

# -------------------------------------------------------------------------
# X,Z Plot Comparing Matched Observation Statistics in Pressure/Height bins
#
def z_distr(vertbinsP,vertbinsH,Dlat,Dlon,tx,ty,tz,x_name,tname,xvarstr,zvarstr,acolors,xunits,zunits,regionstr,latmin,latmax,**kwargs):

	# create plot: PRESSURE y-axis
    fig = plt.figure(figsize=(12,8))
    		# create grid of panels
    gsP = fig.add_gridspec(nrows=1, ncols=2)
        	# create each panel
    axP1 = fig.add_subplot(gsP[0,0])
    axP2 = fig.add_subplot(gsP[0,1])
    
    x = tx
    y = ty
    z = tz
    
    	# set y-axis range bins
    ylabel    = "Pressure (hPa)"
    vertbins  = vertbinsP
    ylabel2   = "Height (km)"
    vertbins2 = vertbinsH

	# Vert bin ranges
	# ...Pressure
    binmins = [vertbins[0]]
    binmaxs = [-999]
    ip=1
    for i in range(np.size(vertbins)-2):
      tbinmins = vertbins[ip] - (vertbins[ip] - vertbins[ip-1])/2.
      tbinmaxs = vertbins[ip] + (vertbins[ip+1] - vertbins[ip])/2.
      binmins.append(tbinmins)
      binmaxs.append(tbinmaxs)
      del tbinmins,tbinmaxs
      ip+=1
    binmins.append(-999)
    binmaxs.append(vertbins[np.size(vertbins)-1])
	# ...Height
    binmins2 = [vertbins2[0]]
    binmaxs2 = [-999]
    ip=1
    for i in range(np.size(vertbins2)-2):
      tbinmins = vertbins2[ip] - (vertbins2[ip] - vertbins2[ip-1])/2.
      tbinmaxs = vertbins2[ip] + (vertbins2[ip+1] - vertbins2[ip])/2.
      binmins2.append(tbinmins)
      binmaxs2.append(tbinmaxs)
      del tbinmins,tbinmaxs
      ip+=1
    binmins2.append(-999)
    binmaxs2.append(vertbins2[np.size(vertbins2)-1])

    	# loop thru all obs for dataset
    diffP    = [[] for i in range(np.size(acolors))]
    avgP     = [[] for i in range(np.size(acolors))]
    sdP      = [[] for i in range(np.size(acolors))]
    nobsP    = [[] for i in range(np.size(acolors))]
    legendP  = [[] for i in range(np.size(acolors))]
    colorP   = [[] for i in range(np.size(acolors))]
    markers_on = [[] for i in range(np.size(acolors))]		# array to indicate statistical significance
    nobsmaxP = 0
    for j in range(np.size(zvarstr)):
      ttmp_x = x[j]
      ttmp_y = y[j]
      ttmp_z = z[j]

      tlat = Dlat[j]      
      idx = np.where((tlat>=latmin)*(tlat<latmax))
      tmp_x = ttmp_x[idx]
      tmp_y = ttmp_y[idx]
      tmp_z = ttmp_z[idx]
      del tlat,idx,ttmp_x,ttmp_y,ttmp_z
     
      legendP[j] = tname[j]
      colorP[j]  = acolors[j] 
      if zvarstr[j]=="Pressure":
        tbins      = vertbins
        bmin       = binmins
        bmax       = binmaxs
      elif zvarstr[j]=="Height":
        tbins      = vertbins2
        bmin       = binmins2
        bmax       = binmaxs2
      else:
        continue

      idx = []
      tidx = np.where((tmp_z<=bmin[0])*(tmp_x!=0.0)*(tmp_y!=0.0))
      idx.append(tidx)
      del tidx
      for ip in range(np.size(tbins)):
        tidx = np.where((tmp_z>bmin[ip])*(tmp_z<=bmax[ip])*(tmp_x!=0.0)*(tmp_y!=0.0))
        idx.append(tidx)
        del tidx
      tidx = np.where((tmp_z>bmax[np.size(tbins)-1])*(tmp_x!=0.0)*(tmp_y!=0.0))
      idx.append(tidx)
      del tidx

      for ip in range(np.size(tbins)):
        if np.size(idx[ip])>0:
          aidx = idx[ip]
          varxdate = tmp_x[aidx]
          varydate = tmp_y[aidx]
              #compute stats per timestep
          tdiff = varydate - varxdate					# y-x = DEPENDENT minus DRIVER
          tavg  = np.mean(tdiff)
          tsd   = np.std(tdiff)						# standard deviation
          tnobs = np.size(varxdate)
              #append stats per timestep
          diffP[j].append(tavg)
          sdP[j].append(tsd)
          nobsP[j].append(tnobs)
              # find absolute max of nobs
          if tnobs>nobsmaxP:
            nobsmaxP = tnobs
	      # statistical significance
          signif = 95.0                                      # statistical signifiance in percent (%). Type: float

          tsigmax = (100.0-signif)/100.0                     # max p-value allowed for statistical signiificance
          tsig = ttest_ind(varxdate,varydate).pvalue         # get p-value from Student's t-test
          if tsig<=tsigmax:                                  # if p-value < tsigmax: difference is statistically significant at 95% level
            markers_on[j].append(ip)		             # append vertical level index that is stat. signif. at 95%

          del varxdate,varydate,tdiff,tavg,tsd,tnobs,aidx,tsig
        else:
          diffP[j].append(np.nan)
          sdP[j].append(np.nan)
          nobsP[j].append(np.nan)
      del idx
      del tmp_x,tmp_y,tmp_z
      del tbins,bmin,bmax
    
    if xunits=="m/s":
      xmin = -20
      xmax = 20
      xstride = 2

    custom_lines = [Line2D([0],[0],color="black",label="Diff",linestyle="-"), Line2D([0],[0],color="black",label="SD_Diff",linestyle="--"), Line2D([0],[0],color="black",label="Stat Signif. "+str(int(signif))+"%",linestyle="None",marker="o")]

	# set axes
    xvalues = list(np.arange(xmin,xmax,xstride))
    xlabels = xvalues

    axP1.set_xticks(xvalues)
    axP1.set_xticklabels(xlabels, fontsize=fontaxis+2)
    axP1.set_ylim([min(vertbins),max(vertbins)])

    axP2.set_xlim([0,nobsmaxP])
    axP2.set_ylim([min(vertbins),max(vertbins)])

    ax2P1 = axP1.twinx()        # allows for 2nd y-axis to be plotted on same plot
    ax2P1.set_ylim([min(vertbins2),max(vertbins2)])

    ax2P2 = axP2.twinx()        # allows for 2nd y-axis to be plotted on same NOBS plot
    ax2P2.set_ylim([min(vertbins2),max(vertbins2)])

    axP1.tick_params(axis='both', which='major', labelsize=fontaxis)
    axP2.tick_params(axis='both', which='major', labelsize=fontaxis)
    ax2P1.tick_params(axis='both', which='major', labelsize=fontaxis)
    ax2P2.tick_params(axis='both', which='major', labelsize=fontaxis)

	# plot lines
    for j in range(np.size(zvarstr)): 
      if np.size(diffP[j])>0:
        if zvarstr[j]=="Pressure":
          if np.size(markers_on[j])>0:
            axP1.plot(diffP[j], vertbins, '-o', color=colorP[j], markevery=markers_on[j])
          else:
            axP1.plot(diffP[j], vertbins, '-o', color=colorP[j])
          axP1.plot(sdP[j],   vertbins, color=colorP[j], linestyle='dashed')
          axP2.plot(nobsP[j], vertbins, color=colorP[j], label=legendP[j])
        elif zvarstr[j]=="Height":
          if np.size(markers_on[j])>0:
            ax2P1.plot(diffP[j], vertbins2, '-o', color=colorP[j], markevery=markers_on[j])
          else:
            ax2P1.plot(diffP[j], vertbins2, '-o', color=colorP[j])
          ax2P1.plot(sdP[j],   vertbins2, color=colorP[j], linestyle='dashed')
          ax2P2.plot(nobsP[j], vertbins2, color=colorP[j], label=legendP[j])

    if xvarstr=="HLOS":
      xvarstr += " Wind Speed"

    	# left panel
    axP1.axhline(y=0,color="black")  	      			# horizontal line
    axP1.axvline(x=0,color="black")  	      			# vertical line
    axP1.set_title(regionstr+" "+str(xvarstr)+" Stats: DEPENDENT - DRIVER", fontsize=fonttitle-fonttitle/4)
    axP1.set_xlabel("Driver "+str(xvarstr)+" ("+str(xunits)+")", fontsize=fontaxis)
    axP1.set_ylabel(ylabel, fontsize=fontaxis)
    ax2P1.set_ylabel(ylabel2, fontsize=fontaxis)
    
	# legend
    legP1 = axP1.legend(handles=custom_lines, loc='lower left', prop={'size': fontlegend})
    legP2 = axP2.legend(loc='upper right', prop={'size': fontlegend-1})
    legP2 = ax2P2.legend(loc='upper center', prop={'size': fontlegend-1})

    for lhP1 in legP1.legendHandles:
      lhP1.set_alpha(1)
    for lhP2 in legP2.legendHandles:
      lhP2.set_alpha(1)

    legendmarkersize = 30
    for lhP1 in range(np.size(legP1.legendHandles)):
      legP1.legendHandles[lhP1]._sizes = [legendmarkersize]
    for lhP2 in range(np.size(legP2.legendHandles)):
      legP2.legendHandles[lhP2]._sizes = [legendmarkersize]

	# right panel
    #axP2.axhline(y=0,color="black")  	      			# horizontal line
    #axP2.axvline(x=0,color="black")  	      			# vertical line
    axP2.set_xlabel("Matched Obs Count", fontsize=fontaxis)
    axP2.set_ylabel(ylabel, fontsize=fontaxis)
    ax2P2.set_ylabel(ylabel2, fontsize=fontaxis)

    axP1.invert_yaxis()
    axP2.invert_yaxis()
  
    fig.tight_layout()
 
    return fig
    
# -------------------------------------------------------------------------
# X,Y Plot showing bias/error vs binned wind speed
#
def stats_windbins(bins,Dlat,tx,ty,x_name,tname,xvarstr,acolors,xunits,yunits,regionstr,latmin,latmax,**kwargs):

	# create plot: PRESSURE y-axis
    fig = plt.figure(figsize=(12,8))
    		# create grid of panels
    gsP = fig.add_gridspec(nrows=1, ncols=2)
        	# create each panel
    axP1 = fig.add_subplot(gsP[0,0])
    axP2 = fig.add_subplot(gsP[0,1])
    
    x = tx
    y = ty
    
    	# set axis labels and range bins
    if xvarstr == "HLOS": xvarstrplot = "HLOS Wind Velocity"
    xlabel = "DRIVER ("+str(x_name)+") "+str(xvarstrplot)+" (m/s)"
    ylabel = "Wind Diff, SD (m/s)"
    xbins = bins

	# DRIVER wind bin ranges
    binmins = [xbins[0]]
    binmaxs = [-999]
    ip=1
    for i in range(np.size(xbins)-2):
      tbinmins = xbins[ip] - (xbins[ip] - xbins[ip-1])/2.
      tbinmaxs = xbins[ip] + (xbins[ip+1] - xbins[ip])/2.
      binmins.append(tbinmins)
      binmaxs.append(tbinmaxs)
      del tbinmins,tbinmaxs
      ip+=1
    binmins.append(-999)
    binmaxs.append(xbins[np.size(xbins)-1])

    	# loop thru all obs for dataset
    diff    = [[] for i in range(np.size(acolors))]
    avgdiff = [[] for i in range(np.size(acolors))]
    sdx     = [[] for i in range(np.size(acolors))]
    sdy     = [[] for i in range(np.size(acolors))]
    sddiff  = [[] for i in range(np.size(acolors))]
    nobs    = [[] for i in range(np.size(acolors))]
    legends = [[] for i in range(np.size(acolors))]
    colors  = [[] for i in range(np.size(acolors))]
    markers_on = [[] for i in range(np.size(acolors))]		# array to indicate statistical significance
    nobsmax = 0
    for j in range(np.size(acolors)):
      ttmp_x = x[j]
      ttmp_y = y[j]

      tlat = Dlat[j]      
      idx = np.where((tlat>=latmin)*(tlat<latmax))
      tmp_x = ttmp_x[idx]
      tmp_y = ttmp_y[idx]
      del tlat,idx,ttmp_x,ttmp_y
      
      legends[j] = tname[j]
      colors[j]  = acolors[j]

      idx = []
      tidx = np.where((tmp_x<=binmins[0])*(tmp_x!=0.0)*(tmp_y!=0.0))
      idx.append(tidx)
      del tidx
      for ip in range(np.size(xbins)):
        tidx = np.where((tmp_x>binmins[ip])*(tmp_x<=binmaxs[ip])*(tmp_x!=0.0)*(tmp_y!=0.0))
        idx.append(tidx)
        del tidx
      tidx = np.where((tmp_x>binmaxs[np.size(xbins)-1])*(tmp_x!=0.0)*(tmp_y!=0.0))
      idx.append(tidx)
      del tidx

      for ip in range(np.size(xbins)):
        if np.size(idx[ip])>0:
          aidx = idx[ip]
          varxdate = tmp_x[aidx]
          varydate = tmp_y[aidx]	  
              #compute stats per x bin
          tdiff   = varydate - varxdate			        # wind diff: y-x = DEPENDENT minus DRIVER
          tavg    = np.mean(tdiff)				# mean wind diff
          tsdx    = np.std(varxdate)				# SD of x winds
          tsdy    = np.std(varydate)				# SD of y winds
          tsddiff = np.std(tdiff)				# SD of wind diff
          tnobs   = np.size(varxdate)
              #append stats per x bin
          diff[j].append(tavg)
          sdx[j].append(tsdx)
          sdy[j].append(tsdy)
          sddiff[j].append(tsddiff)
          nobs[j].append(tnobs)
              # find absolute max of nobs
          if tnobs>nobsmax:
            nobsmax = tnobs
	      # statistical significance
          signif = 95.0                                      # statistical signifiance in percent (%). Type: float

          tsigmax = (100.0-signif)/100.0                     # max p-value allowed for statistical signiificance
          tsig = ttest_ind(varxdate,varydate).pvalue         # get p-value from Student's t-test
          if tsig<=tsigmax:                                  # if p-value < tsigmax: difference is statistically significant at 95% level
            markers_on[j].append(ip)                         # append vertical level index that is stat. signif. at 95%

          del varxdate,varydate,tdiff,tavg,tsdx,tsdy,tsddiff,tnobs,aidx,tsig
        else:
          diff[j].append(np.nan)
          sdx[j].append(np.nan)
          sdy[j].append(np.nan)
          sddiff[j].append(np.nan)
          nobs[j].append(np.nan)
      del idx
      del tmp_x,tmp_y
    
    if xunits=="m/s":
      ymin = -15
      ymax = 15

    custom_lines = [Line2D([0],[0],color="black",label="Diff",linestyle="-"), Line2D([0],[0],color="black",label="SD_Diff",linestyle="--"), Line2D([0],[0],color="black",label="SD_drv",linestyle="-."), Line2D([0],[0],color="black",label="SD_dep",linestyle=":"), Line2D([0],[0],color="black",label="Stat Signif. "+str(int(signif))+"%",linestyle="None",marker="o")]

	# set axis limits
    if xvarstr=="Wind Speed":
      axP1.set_xlim([0,100])
      axP2.set_xlim([0,100])
    else:
      axP1.set_xlim([-100,100])
      axP2.set_xlim([-100,100])

    axP1.set_ylim([ymin,ymax])
    axP2.set_ylim([0,nobsmax])
    axP1.tick_params(axis='both', which='major', labelsize=fontaxis)              # change font size of major tick labels on both axes
    axP2.tick_params(axis='both', which='major', labelsize=fontaxis)

    for j in range(np.size(acolors)): 
      if np.size(diff[j])>0:
     	# plot lines
        if np.size(markers_on[j])>0:
          axP1.plot(xbins, diff[j], '-o', color=colors[j], markevery=markers_on[j])
        else:
          axP1.plot(xbins, diff[j], '-o', color=colors[j])
        axP1.plot(xbins, sdx[j],    color=colors[j], linestyle='dashdot')
        axP1.plot(xbins, sdy[j],    color=colors[j], linestyle='dotted')
        axP1.plot(xbins, sddiff[j], color=colors[j], linestyle='dashed')
        axP2.plot(xbins, nobs[j],   color=colors[j], label=legends[j])

    	# left panel
    axP1.axhline(y=0,color="black")  	      			# horizontal line
    axP1.axvline(x=0,color="black")  	      			# vertical line
    axP1.set_title(regionstr+" "+str(xvarstrplot)+" Stats: DEPENDENT - DRIVER", fontsize=fonttitle-3)
    axP1.set_ylabel(ylabel, fontsize=fontaxis-3)
    axP1.set_xlabel(xlabel, fontsize=fontaxis-3)
    
	# legend
    legP1 = axP1.legend(handles=custom_lines, loc='upper right', prop={'size': fontlegend})
    legP2 = axP2.legend(loc='upper right', prop={'size': fontlegend})

    for lhP1 in legP1.legendHandles:
      lhP1.set_alpha(1)
    for lhP2 in legP2.legendHandles:
      lhP2.set_alpha(1)

    legendmarkersize = 30
    for lhP1 in range(np.size(legP1.legendHandles)):
      legP1.legendHandles[lhP1]._sizes = [legendmarkersize]
    for lhP2 in range(np.size(legP2.legendHandles)):
      legP2.legendHandles[lhP2]._sizes = [legendmarkersize]

	# right panel
    #axP2.axhline(y=0,color="black")  	      			# horizontal line
    axP2.axvline(x=0,color="black")  	      			# vertical line
    axP2.set_ylabel("Matched Obs Count", fontsize=fontaxis-3)
    axP2.set_xlabel(xlabel, fontsize=fontaxis-3)

    fig.tight_layout()

    return fig


# -------------------------------------------------------------------------

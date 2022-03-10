#Allison Towner
#15 November 2021
#Script to create summary figures for SiO Outflow catalog candidates


#import all the necessary stuff
import numpy as np
import astropy
import pylab as pl
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from astropy import visualization
from astropy.visualization import quantity_support
import warnings
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
from spectral_cube import SpectralCube
import regions
from operator import xor
import os
from pathlib import Path
from mpl_toolkits.axes_grid1.inset_locator import TransformedBbox, BboxPatch, BboxConnector
from matplotlib.transforms import Bbox
from spectral_cube.spectral_axis import convert_spectral_axis
from astropy.constants import c as clight
from matplotlib.patches import *
from astropy.coordinates import SkyCoord
from regions import PixCoord, PolygonSkyRegion, PointSkyRegion, Regions
from sio_masterfile import *

#initialize field-specific information - later this will be the basis for a loop
#fields=['G328.25']
fields=targets.keys()

for field in fields:
    candidates=targets[field]['outflows'].keys()
    #candidates=['s2']
    for cand in candidates:
        if field=='W43-MM1':#iff the field is MM1, have to account for the fact that all .crtf entries are in sexagesimal format rather than decimal degrees. Regions.read doesn't handle this correctly - not sure if it's supposed to or not, but it definitely doesn't.
            f=open(targets[field]['outflows'][cand]['ap'],'r')
            f.readline()
            line=f.readline()
            line=line.split()
            stopindex=line.index('coord=J2000,')
            line=line[1:stopindex]
            x=np.array([])
            y=np.array([])
            for l in line:
                #print(l)
                if ':' in l:
                    l=l.strip('[')
                    l=l.strip(',')
                    xsg = l.split(':')
                    xcoo = 15*(float(xsg[0]) + (float(xsg[1])/60.) + (float(xsg[2])/3600.))
                    x=np.append(x,float(xcoo))
                    #print(x)
                else:
                    l=l.strip('],')
                    ysg = l.split('.')
                    ycoo = abs(float(ysg[0])) + (float(ysg[1])/60.) + ((float(ysg[2]) + (float(ysg[3])/1.0E4))/3600.)
                    if float(ysg[0]) < 0.0:
                        ycoo = -1.0*ycoo
                    y=np.append(y,float(ycoo))
                    #print(y)
            #now take x and y and make a region object with them for the region
            vertices = SkyCoord(x,y,unit='deg',frame='fk5')
            #region=Regions(PolygonSkyRegion(vertices=vertices))
            centers=[]
            for i in range(0,1):
                region_sky = PolygonSkyRegion(vertices=vertices)
                centers.append(region_sky)
            region = Regions(centers)
            f.close()

            f=open(targets[field]['outflows'][cand]['pv'],'r')
            f.readline()
            lines=f.readlines()
            x=np.array([])
            y=np.array([])
            for line in lines:
                line=line.split()
                xsg=line[1].strip('[')
                xsg=xsg.strip(',')
                xsg=xsg.split(':')
                xcoo=15*(float(xsg[0]) + (float(xsg[1])/60.) + (float(xsg[2])/3600.))
                x=np.append(x,xcoo)
                ysg=line[2].strip('],')
                ysg=ysg.split('.')
                ycoo=abs(float(ysg[0])) + (float(ysg[1])/60.) + ((float(ysg[2]) + float(ysg[3])/1E4)/3600.)
                if float(ysg[0]) < 0.0:
                    ycoo = -1.0*ycoo
                y=np.append(y,ycoo)            
            #now take x and y and make a region object with the for the pathregion
            #pathregion = [PointSkyRegion(center=SkyCoord(x[i],y[i],unit='deg',frame='fk5') for ra, dec in [x,y]]
            centers=[]
            for i in range(0,len(x)):
                center_sky = SkyCoord(x[i],y[i],unit='deg',frame='fk5')
                region_sky = PointSkyRegion(center=center_sky)
                centers.append(region_sky)
            pathregion = Regions(centers)
            f.close()

        else:
            region = regions.Regions.read(targets[field]['outflows'][cand]['ap'])
            pathregion = regions.Regions.read(targets[field]['outflows'][cand]['pv'])

        #read in all the necessary files
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')

            #Moment maps
            moment0 = SpectralCube.read(field+'/'+field+'_B6_spw1_12M_sio.image.common.rebin.mom0',format='casa_image')
            moment1 = SpectralCube.read(field+'/'+field+'_B6_spw1_12M_sio.image.common.rebin.mom1.masked',format='casa_image')
            moment2 = SpectralCube.read(field+'/'+field+'_B6_spw1_12M_sio.image.common.rebin.mom2.masked',format='casa_image')
            moment8 = SpectralCube.read(field+'/'+field+'_B6_spw1_12M_sio.image.common.rebin.mom8',format='casa_image')
    
            #Line cube
            cube = SpectralCube.read(field+'/'+field+'_B6_spw1_12M_sio.image.common.rebin',format='casa_image')
    
            #Band6 continuum
            continuum = fits.open(targets[field]['B6'])
    
            #PV Diagram (created separately by /orange/adamginsburg/atowner/ALMA_IMF/sio_outflows/make_pvslices.py)
            pvplot = fits.open(field+'/'+field+'_'+cand+'_pvslice.fits') #don't construct the slices on the fly - use the separate make_pvslices script to make them so that they exist as their own separate files
    

            
        coords = WCS(moment0[0].header)
        contcoords = WCS(continuum[0].header)
        pvcoords = WCS(pvplot[0].header)

        coords = coords.celestial #need .celestial here to force the object to be 2-dimensional
        contcoords = contcoords.celestial

        maskcube = SpectralCube.read(field+'/'+field+'_B6_spw1_12M_sio.mask',format='casa_image')
        mask = maskcube[10,:,:] > 0.99
        mom0_masked = moment0.with_mask(mask)

        #make a subcube of the mom0 map using the region initialized in the first cell
        mom0 = moment0.subcube_from_regions(region)
        mom0[0].meta['slice'][0]

        #get the pixel coordinates of the extremes of the region +/- 20 pixels for a bit of padding
        ymin = mom0[0].meta['slice'][0][1][0] - 20
        ymax = mom0[0].meta['slice'][0][1][1] + 20
        xmin = mom0[0].meta['slice'][0][2][0] - 20
        xmax = mom0[0].meta['slice'][0][2][1] + 20
        
        pathx = []
        pathy = []
        for p in pathregion:
            pathx.append(p.to_pixel(coords).center.x)
            pathy.append(p.to_pixel(coords).center.y)
        pymin=int(np.round(np.min(pathy))) - 10
        pymax=int(np.round(np.max(pathy))) + 10
        pxmin=int(np.round(np.min(pathx))) - 10
        pxmax=int(np.round(np.max(pathx))) + 10

        #make subcubes using the min/max pixel coordinates determined above
        mom0 = moment0[0][ymin:ymax,xmin:xmax]
        mom1 = moment1[0][ymin:ymax,xmin:xmax]
        mom2 = moment2[0][ymin:ymax,xmin:xmax]
        mom8 = moment8[0][ymin:ymax,xmin:xmax]
        
        #create the coordinate axes for the moment-map zoom
        zoomcoords = WCS(mom0.header)
        
        #convert the min-max pixel coordinates to world coordinates
        blc = coords.pixel_to_world(xmin,ymin)
        trc = coords.pixel_to_world(xmax,ymax)

        #convert the world coordinates to pixel coordinates *in the continuum image,* which has different beam and pixel sizes than the moment maps
        cxmin,cymin = contcoords.world_to_pixel(blc)
        cxmax,cymax = contcoords.world_to_pixel(trc)
        cxmin = int(np.round(cxmin)) #make sure that the pixel coordinates are integers
        cymin = int(np.round(cymin))
        cxmax = int(np.round(cxmax))
        cymax = int(np.round(cymax))

        #create the coordinate axes for the continuum zoom
        contzoom = contcoords[cymin:cymax,cxmin:cxmax]
        
        #extract the spectrum for our region of interest only
        regioncube = cube.subcube_from_regions(region)
        spectrum = regioncube.sum(axis=(1,2))
        
        #get x-axis in km/s for spectrum plot
        ##quantity_support()
        restfreq = cube.header['RESTFRQ'] * u.Hz
        refpoint = cube.header['CRPIX3']
        refpointvalue = cube.header['CRVAL3']
        cdelt = cube.header['CDELT3']
        cubestart = refpointvalue - refpoint*cdelt
        xaxis = np.arange(cubestart,cubestart+float(len(spectrum))*cdelt,cdelt) * u.Hz
        if len(xaxis) > len(spectrum):
            xaxis = xaxis[:-1]
        if len(spectrum) > len(xaxis):
            spectrum = spectrum[:-1]
        freq_to_vel = u.doppler_radio(restfreq)
        xaxis = xaxis.to(u.km / u.s,equivalencies=freq_to_vel)
        #

        #make x-axis in km/s starting from km/s (not frequency)
        #clight = clight.to(u.km/u.s)
        #cdelt_vel = cube.header['CDELT3']*clight/cube.header['CRVAL3']
        #cubestart = cube.header['CRPIX3']*u.km/u.s + 2*cdelt_vel
        #xaxis = np.arange(cubestart.value,cubestart.value+float(len(spectrum))*cdelt_vel.value,cdelt_vel.value) * u.km/u.s
        #if len(xaxis) > len(spectrum):
        #    xaxis = xaxis[:-1]
        #if len(spectrum) > len(xaxis):
        #    spectrum = spectrum[:-1]
        #spectrum = spectrum[::-1]
        
        #set mom1 vmin and vmax automatically
        vlsr=targets[field]['vlsr']
        vmin=vlsr-30.0
        vmax=vlsr+30.0
        
        #set mom2 min and max dv to be the following absolute values:
        dvmin=5.0
        dvmax=60.0
        
        ###############################################
        #Start making the actual figure
        ###############################################

        fig = plt.figure(figsize=(24,18))
        gs = GridSpec(3,4,figure=fig)
        fig.suptitle(field+' Structure '+cand[1:],fontsize=22)
            
        data=mom0_masked[0].data
        notnan=np.isfinite(data)
        xminnan=np.argmax(np.any(notnan,axis=0))
        yminnan=np.argmax(np.any(notnan,axis=1))
        xmaxnan=notnan.shape[1] - np.argmax(np.any(notnan[::-1,::-1],axis=0))
        ymaxnan=notnan.shape[0] - np.argmax(np.any(notnan[::-1,::-1],axis=1))
        trimmed=mom0_masked[0][yminnan:ymaxnan,xminnan:xmaxnan]
        trimcoords=coords[yminnan:ymaxnan,xminnan:xmaxnan]
        
        fov = fig.add_subplot(gs[0:2,0:2],projection=trimcoords)
        im = fov.imshow(trimmed.data,cmap='viridis',norm=visualization.simple_norm(data,stretch='power',power=0.5,max_percent=99.95)) 
        fov.set_title(f"Moment 0, full FOV",fontsize=18)
        cbar=fig.colorbar(im,ax=fov,shrink=0.65)
        ##cbar.set_ticks([0.0,0.995*np.max(data)])
        cbar.set_label(f"Integrated Intensity (Jy/beam km/s)")
        fov.set_xlabel(f"Right Ascension (ICRS)")
        fov.set_ylabel(f"Declination (ICRS)")
        tick_fontsize=10
        fontsize=12
        ra = fov.coords['ra']
        dec = fov.coords['dec']
        ra.set_major_formatter('hh:mm:ss.s')
        dec.set_major_formatter('dd:mm:ss')
        ra.set_axislabel(f"RA (ICRS)",fontsize=fontsize)
        dec.set_axislabel(f"Dec (ICRS)",fontsize=fontsize,minpad=-0.5)
        ra.ticklabels.set_fontsize(tick_fontsize)
        dec.ticklabels.set_fontsize(tick_fontsize)
        ra.set_ticklabel(exclude_overlapping=True)
        dec.set_ticklabel(exclude_overlapping=True)
        
        pixregion = region[0].to_pixel(trimcoords)
        pixregion.visual['color'] = 'black'
        pixartist = pixregion.as_artist()
        fov.add_artist(pixartist)
        
        m0reg = fig.add_subplot(gs[0,2],projection=zoomcoords)
        im = m0reg.imshow(mom0.data,cmap='viridis',norm=visualization.simple_norm(mom0.data,stretch='sqrt',max_percent=99.95))
        m0reg.set_title(f"Moment 0, zoom view w/ region",fontsize=16)
        cbar=fig.colorbar(im,ax=m0reg,shrink=0.5)
        ##cbar.set_ticks([0.0,0.995*np.max(data)])
        cbar.set_label(f"Integrated Intensity (Jy/beam km/s)")
        
        m0spec = fig.add_subplot(gs[0,3])
        im = m0spec.plot(xaxis,spectrum.value)
        m0spec.set_title(f"Integrated Spectrum of Region (vlsr={vlsr} {u.km/u.s})",fontsize=16)
        m0spec.set_ylabel(f"Integrated Intensity ({spectrum.unit})")
        
        m0pvpath = fig.add_subplot(gs[1,2],projection=zoomcoords)
        im = m0pvpath.imshow(mom0.data,cmap='viridis',norm=visualization.simple_norm(mom0.data,stretch='power',power=0.7,max_percent=99.95))
        m0pvpath.set_title(f"Moment 0, zoom view w/ PV path",fontsize=16)
        cbar=fig.colorbar(im,ax=m0pvpath,shrink=0.5)
        ##cbar.set_ticks([0.0,0.995*np.max(data)])
        cbar.set_label(f"Integrated Intensity (Jy/beam km/s)")
        pathx = []
        pathy = []
        for p in pathregion:
            pathx.append(p.to_pixel(zoomcoords).center.x)
            pathy.append(p.to_pixel(zoomcoords).center.y)
        coords = np.concatenate((pathx,pathy))
        coords = np.reshape(coords,(2,-1))
        coords = np.transpose(coords)
        if '.ra.crtf' in targets[field]['outflows'][cand]['pv']:
            coords = sorted(coords,key = lambda coo: coo[0]) #be very careful. This line sorts the coordinates in order of increasing RA.
        elif '.dec.crtf' in targets[field]['outflows'][cand]['pv']:
            coords = sorted(coords,key = lambda coo: coo[1]) #be very careful. This line sorts the coordinates in order of increasing Dec.
        pathx=[]
        pathy=[]
        for coo in coords:
            pathx.append(coo[0])
            pathy.append(coo[1])
        m0pvpath.plot(pathx,pathy,color='black',linewidth=2)
        m0pvpath.scatter(pathx,pathy,color='black',linewidth=1)
        
        newpvcoords = convert_spectral_axis(pvcoords,'km/s','VRAD',rest_value=restfreq)
        newpvcoords.wcs.cdelt[1] /= 1000.0
        newpvcoords.wcs.crval[1] -= targets[field]['vlsr']
        newpvcoords.wcs.crval[1] /= 1000.0
        newpvcoords.wcs.cunit[1] = 'km/s'
        newpvcoords.wcs.cdelt[0] *= 3600
        newpvcoords.wcs.cunit[0] = u.arcsec
        newpvcoords.wcs.crval[0] = 0
        m0pvplot = fig.add_subplot(gs[1,3],projection=newpvcoords)
        aspect=1.5*len(pvplot[0].data[0])/len(pvplot[0].data[1])
        m0pvplot.set_aspect('2.0')
        im = m0pvplot.imshow(pvplot[0].data,cmap='magma')#,aspect=str(aspect))
        m0pvplot.set_title(f"PV Diagram, vlsr = {vlsr} {u.km/u.s}",fontsize=16)
        m0pvplot.set_xlabel(f"Offset ('')")
        m0pvplot.set_ylabel(f"Velocity Offset ({u.km/u.s})")
        #m0pvplot.set_aspect(len(pvplot[0].data[0])*1.5/len(pvplot[0].data[1]))
        
        m1 = fig.add_subplot(gs[2,0],projection=zoomcoords)
        im = m1.imshow(mom1.data,cmap='seismic',norm=visualization.simple_norm(mom1.data,stretch='linear',min_cut=vmin,max_cut=vmax))
        m1.set_title(f"Moment 1, zoom view w/ region",fontsize=16)
        cbar=fig.colorbar(im,ax=m1,shrink=0.5)
        ##cbar.set_ticks([0.0,0.995*np.max(data)])
        cbar.set_label(f"Velocity (km/s)")
        
        m2 = fig.add_subplot(gs[2,1],projection=zoomcoords)
        im = m2.imshow(mom2.data,cmap='inferno',norm=visualization.simple_norm(mom2.data,stretch='linear',min_cut=dvmin,max_cut=dvmax))
        m2.set_title(f"Moment 2, zoom view w/ region",fontsize=16)
        cbar=fig.colorbar(im,ax=m2,shrink=0.5)
        ##cbar.set_ticks([0.0,0.995*np.max(data)])
        cbar.set_label(f"Velocity Dispersion (km/s)")
        
        m8 = fig.add_subplot(gs[2,2],projection=zoomcoords)
        im = m8.imshow(mom8.data,cmap='inferno',norm=visualization.simple_norm(mom8.data,stretch='sqrt',min_percent=0.0, max_percent=99.95))
        m8.set_title(f"Maximum Intensity, zoom view w/ region",fontsize=16)
        cbar=fig.colorbar(im,ax=m8,shrink=0.5)
        ##cbar.set_ticks([0.0,0.995*np.max(data)])
        cbar.set_label(f"Maximum Intensity (Jy/beam)")
        
        pixregion = region[0].to_pixel(zoomcoords)
        pixregion.visual['color'] = 'black'
        pixartist = pixregion.as_artist()
        pixartist.set_linewidth(3)
        m0reg.add_artist(pixartist)
        pixartist1 = pixregion.as_artist()
        pixartist1.set_linewidth(3)
        m1.add_artist(pixartist1)
        pixartist2 = pixregion.as_artist()
        pixartist2.set_linewidth(3)
        m2.add_artist(pixartist2)
        pixartist8 = pixregion.as_artist()
        pixartist8.set_linewidth(3)
        m8.add_artist(pixartist8)
        
        rectangle = BboxPatch(TransformedBbox(Bbox(np.array([(xmin-xminnan,ymin-yminnan),(xmax-xminnan,ymax-yminnan)])),fov.transData), fill=False, ec='black',linewidth=3)
        fov.add_patch(rectangle)
        
        c = fig.add_subplot(gs[2,3],projection=contzoom)
        im = c.imshow(continuum[0].data.squeeze()[cymin:cymax,cxmin:cxmax],cmap='Greys',norm=visualization.simple_norm(continuum[0].data.squeeze(),stretch='log',min_percent=5.0,max_percent=99.99))
        c.set_title(f"Band6 Continuum, zoom view w/ region",fontsize=16)
        cbar=fig.colorbar(im,ax=c,shrink=0.5)
        ##cbar.set_ticks([0.0,0.995*np.max(data)])
        cbar.set_label(f"Intensity (Jy/beam)")
        #c.axis([10,20,10,20])
        cpixregion = region[0].to_pixel(contzoom)
        cpixregion.visual['color'] = 'black'
        cpixartist = cpixregion.as_artist()
        cpixartist.set_linewidth(3)
        c.add_artist(cpixartist)

        for panel in [m0reg,m0pvpath,m1,m2,m8,c]:
            panel.set_xlabel(f"Right Ascension (ICRS)")
            panel.set_ylabel(f"Declination (ICRS)")
            tick_fontsize=10
            fontsize=12
            ra = panel.coords['ra']
            dec = panel.coords['dec']
            ra.set_major_formatter('hh:mm:ss.s')
            dec.set_major_formatter('dd:mm:ss')
            ra.set_axislabel(f"RA (ICRS)",fontsize=fontsize)
            dec.set_axislabel(f"Dec (ICRS)",fontsize=fontsize,minpad=-0.75)
            ra.ticklabels.set_fontsize(tick_fontsize)
            dec.ticklabels.set_fontsize(tick_fontsize)
            ra.set_ticklabel(exclude_overlapping=True)
            dec.set_ticklabel(exclude_overlapping=True)
            
        pl.tight_layout()
        fig.savefig(field+'/'+field+'_'+cand+'.summary.png',format='png')

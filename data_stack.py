'''
Created on Jun 20, 2011

@author: Mirna Lerotic
'''
from __future__ import division

import numpy as np
import scipy.interpolate
import scipy.ndimage
import h5py 
import datetime 

import x1a_stk
import aps_hdf5
import data_struct

#----------------------------------------------------------------------
class data(x1a_stk.x1astk,aps_hdf5.h5):
    def __init__(self, data_struct):
        x1a_stk.x1astk.__init__(self)
        aps_hdf5.h5.__init__(self)
        
        self.data_struct = data_struct

#----------------------------------------------------------------------   
    def new_data(self):
        self.n_cols = 0
        self.n_rows = 0
        self.n_ev = 0
               
        self.x_dist = 0       
        self.y_dist = 0

        self.ev = 0 
        self.msec = 0
        self.imagestack = 0             
        self.absdata = 0

        self.original_n_cols = 0
        self.original_n_rows = 0
        self.original_n_ev = 0
        self.original_ev = 0
        self.original_absdata = 0
        
        self.i0data =0
        self.evi0 = 0
        
        self.od = 0
        self.od3d = 0
        self.original_absdata = 0
        self.original_ev = 0
        self.original_n_cols = 0
        self.original_n_ev = 0
        self.original_n_rows = 0
        self.original_od3d = 0
        
           
#----------------------------------------------------------------------   
    def read_stk_i0(self, filename):
        x1a_stk.x1astk.read_stk_i0(self,filename)
        self.calculate_optical_density()
        
        self.fill_h5_struct_normalization()
        
        

#---------------------------------------------------------------------- 
    def read_stk(self, filename):    
        self.new_data()  
        x1a_stk.x1astk.read_stk(self, filename)
        
        self.fill_h5_struct_from_stk()
        
        self.scale_bar()
  
#---------------------------------------------------------------------- 
    def read_h5(self, filename):    
        self.new_data()  
        aps_hdf5.h5.read_h5(self, filename, self.data_struct)
        
        if self.data_struct.spectromicroscopy.normalization.white_spectrum.any():
            self.calculate_optical_density()

            
        self.scale_bar()
        
#---------------------------------------------------------------------- 
    def fill_h5_struct_from_stk(self):   
        
        
        now = datetime.datetime.now()
        
        self.data_struct.implements = 'base, exchange, spectromicroscopy'
        self.data_struct.version = '0.9'
        self.data_struct.file_creation_datetime = now.strftime("%Y-%m-%dT%H:%M")
        self.data_struct.comment = 'Converted from .stk',
        
        self.data_struct.add_experimenter()

        self.data_struct.add_sample()
        
        
        self.data_struct.exchange.add_detector(data=self.absdata, 
                                                  signal = 1, 
                                                  axes='x:y', 
                                                  energy=self.ev, 
                                                  energy_units = 'ev')
        
        
        self.data_struct.exchange.detector[0].add_dimscale(key = 'x', units = 'um', data = self.x_dist)
        self.data_struct.exchange.detector[0].add_dimscale(key = 'y', units = 'um', data = self.y_dist)
        
        self.data_struct.exchange.add_white_data()
        self.data_struct.exchange.add_dark_data()
        
        
        
        self.data_struct.spectromicroscopy.add_beamline()
        self.data_struct.spectromicroscopy.add_positions()
        
        
        self.data_struct.spectromicroscopy.add_normalization()
        
        
#---------------------------------------------------------------------- 
    def fill_h5_struct_normalization(self):   
        
        self.data_struct.spectromicroscopy.normalization.white_spectrum=self.i0data 
        self.data_struct.spectromicroscopy.normalization.white_spectrum_energy = self.evi0
        self.data_struct.spectromicroscopy.normalization.white_spectrum_energy_units='eV'
        
        self.data_struct.spectromicroscopy.optical_density = self.od
        
    
    
#----------------------------------------------------------------------   
    def calc_histogram(self):
        #calculate average flux for each pixel
        self.averageflux = np.mean(self.absdata,axis=2)
        self.histogram = self.averageflux
        
        return
        
        
#----------------------------------------------------------------------   
    def i0_from_histogram(self, fluxmin, fluxmax):

        self.evi0hist = self.ev.copy()

        self.i0datahist = np.zeros(self.n_ev)

        i0_indices = np.where((fluxmin<self.averageflux)&(self.averageflux<fluxmax))
        if i0_indices:
            invnumel = 1./self.averageflux[i0_indices].shape[0]
            for ie in range(self.n_ev):  
                thiseng_abs = self.absdata[:,:,ie]
                self.i0datahist[ie] = np.sum(thiseng_abs[i0_indices])*invnumel


        self.evi0 = self.ev.copy()
        self.i0data = self.i0datahist 

        self.calculate_optical_density()   
        
        self.fill_h5_struct_normalization()
        
    
        return    
    
#----------------------------------------------------------------------   
    def set_i0(self, i0data, evdata):

        self.evi0 = evdata
        self.i0data = i0data 

        self.calculate_optical_density()
        
        self.self.fill_h5_struct_normalization()
    
        return  
    
#----------------------------------------------------------------------   
# Normalize the data: calculate optical density matrix D 
    def calculate_optical_density(self):

        n_pixels = self.n_cols*self.n_rows
        self.od = np.empty((self.n_cols, self.n_rows, self.n_ev))
        
        #little hack to deal with rounding errors
        self.evi0[self.evi0.size-1] += 0.001
        
        fi0int = scipy.interpolate.interp1d(self.evi0,self.i0data, kind='cubic', bounds_error=False, fill_value=0.0)      
        i0 = fi0int(self.ev)
        
        #zero out all negative values in the image stack
        negative_indices = np.where(self.absdata <= 0)
        if negative_indices:
            self.absdata[negative_indices] = 0.01
                           

        for i in range(self.n_ev):
            self.od[:,:,i] = - np.log(self.absdata[:,:,i]/i0[i])
        
        #clean up the result
        nan_indices = np.where(np.isfinite(self.od) == False)
        if nan_indices:
            self.od[nan_indices] = 0
            
        self.od3d = self.od.copy()
        

        #Optical density matrix is rearranged into n_pixelsxn_ev
        self.od = np.reshape(self.od, (n_pixels, self.n_ev), order='F')
        self.original_od3d = self.od3d.copy()

        return
    
#----------------------------------------------------------------------   
    def scale_bar(self): 
           
        x_start = np.amin(self.x_dist)
        x_stop = np.amax(self.x_dist)
        bar_microns = 0.2*np.abs(x_stop-x_start)
        
        if bar_microns >= 10.:
            bar_microns = 10.*int(0.5+0.1*int(0.5+bar_microns))
            bar_string = str(int(0.01+bar_microns)).strip()
        elif bar_microns >= 1.:      
            bar_microns = float(int(0.5+bar_microns))
            if bar_microns == 1.:
                bar_string = '1'
            else:
                bar_string = str(int(0.01+bar_microns)).strip()
        else:
            bar_microns = np.maximum(0.1*int(0.5+10*bar_microns),0.1)
            bar_string = str(bar_microns).strip()
            
        self.scale_bar_string = bar_string

        #Matplotlib has flipped scales so I'm using rows instead of cols!
        self.scale_bar_pixels_x = int(0.5+float(self.n_rows)*
                       float(bar_microns)/float(abs(x_stop-x_start)))
        
        self.scale_bar_pixels_y = int(0.01*self.n_rows)
        
        if self.scale_bar_pixels_y < 2:
                self.scale_bar_pixels_y = 2
                       

             
             
    
#----------------------------------------------------------------------   
    def write_xas(self, filename, evdata, data):
        file = open(filename, 'w')
        print>>file, '*********************  X-ray Absorption Data  ********************'
        print>>file, '*'
        print>>file, '* Formula: '
        print>>file, '* Common name: '
        print>>file, '* Edge: '
        print>>file, '* Acquisition mode: '
        print>>file, '* Source and purity: ' 
        print>>file, '* Comments: Stack list ROI ""'
        print>>file, '* Delta eV: '
        print>>file, '* Min eV: '
        print>>file, '* Max eV: '
        print>>file, '* Y axis: '
        print>>file, '* Contact person: '
        print>>file, '* Write date: '
        print>>file, '* Journal: '
        print>>file, '* Authors: '
        print>>file, '* Title: '
        print>>file, '* Volume: '
        print>>file, '* Issue number: '
        print>>file, '* Year: '
        print>>file, '* Pages: '
        print>>file, '* Booktitle: '
        print>>file, '* Editors: '
        print>>file, '* Publisher: '
        print>>file, '* Address: '
        print>>file, '*--------------------------------------------------------------'
        for ie in range(self.n_ev):
            print>>file, '\t%.6f\t%.6f' %(evdata[ie], data[ie])
        
        file.close()
    
        return  
    
#----------------------------------------------------------------------   
#Read x-ray absorption spectrum
    def read_xas(self, filename):
        
        spectrum_common_name = ' '
 
        f = open(str(filename),'r')
        
        elist = []
        ilist = []    
    
        for line in f:
            if line.startswith('*'):
                if 'Common name' in line:
                    spectrum_common_name = line.split(':')[-1].strip()

            else:
                e, i = [float (x) for x in line.split()] 
                elist.append(e)
                ilist.append(i)
                
        spectrum_evdata = np.array(elist)
        spectrum_data = np.array(ilist) 
                
        f.close()
        
        if spectrum_evdata[-1]<spectrum_evdata[0]:
            spectrum_evdata = spectrum_evdata[::-1]
            spectrum_data = spectrum_data[::-1]
        
        return spectrum_evdata, spectrum_data, spectrum_common_name
        
    
#----------------------------------------------------------------------   
#Register images using Fourier Shift Theorem
    def register_images(self, ref_image, image2, maxshift = 5, have_ref_img_fft = False):
        
        if have_ref_img_fft == False:
            self.ref_fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(ref_image)))
        img2_fft = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(image2)))
        
        fr = (self.ref_fft*img2_fft.conjugate())/(np.abs(self.ref_fft)*np.abs(img2_fft))
        fr = np.fft.fftshift(np.fft.ifft2(np.fft.fftshift(fr)))
        fr = np.abs(fr)
        
        shape = ref_image.shape
        
        xc, yc = np.unravel_index(np.argmax(fr), shape)
        
        # Limit the search to 1 pixel border
        if xc == 0:
            xc=1
        if xc == shape[0]-1:
            xc = shape[0]-2
            
        if yc == 0:
            yc=1
        if yc == shape[1]-1:
            yc = shape[1]-2
            
        #Use peak fit to find the shifts 
        xpts = xc+[-1,0,1]
        ypts = fr[xpts,yc]     
        xf, fit = self.peak_fit(xpts, ypts)

        xpts = yc+[-1,0,1]
        ypts = fr[xc,xpts]     
        yf, fit = self.peak_fit(xpts, ypts)   
              
        
        xshift = xf - np.float(shape[0])/2.0
        yshift = yf - np.float(shape[1])/2.0
        
                
        return xshift, yshift
    

#----------------------------------------------------------------------   
#Apply image registration
    def apply_image_registration(self, image, xshift, yshift):
        
        shape = image.shape
        nx = shape[0]
        ny = shape[1]
        
        outofboundariesval = np.sum(image)/float(nx*ny)        
        shifted_img = scipy.ndimage.interpolation.shift(image,[xshift,yshift],
                                                        mode='constant', 
                                                        cval=outofboundariesval)
        
        return shifted_img
    
#----------------------------------------------------------------------   
#Apply image registration
    def crop_registed_images(self, images, xshifts, yshifts):
        
        min_xshift = np.min(xshifts)
        max_xshift = np.max(xshifts)
        
        min_yshift = np.min(yshifts)
        max_yshift = np.max(yshifts)
        
        # if the image is moved to the right (positive) we need to crop the left side 
        xleft = np.ceil(max_xshift)
        if xleft < 0:
            xleft = 0
        # if the image is moved to the left (negative) we need to crop the right side 
        xright = np.floor(self.n_cols+min_xshift)
        if xright>(self.n_cols):
            xright = self.n_cols
        
        ybottom = np.ceil(max_yshift)
        if ybottom <0:
            ybottom = 0
        ytop = np.floor(self.n_rows+min_yshift)
        if ytop > (self.n_rows):
            ytop = self.n_rows
            
            
        cropped_stack = images[xleft:xright, ybottom:ytop, :]
        
        return cropped_stack


#----------------------------------------------------------------------   
#Quadratic peak fit: Fits the 3 data pairs to y=a+bx+cx^2, returning fit=[a,b,c]'
#  and xpeak at position of inflection'
    def peak_fit(self, x, y):
        
        y1y0=y[1]-y[0]
        y2y0=y[2]-y[0]
        x1x0=np.float(x[1]-x[0])
        x2x0=np.float(x[2]-x[0])
        x1x0sq=np.float(x[1]*x[1]-x[0]*x[0])
        x2x0sq=np.float(x[2]*x[2]-x[0]*x[0])
                
        c_num=y2y0*x1x0-y1y0*x2x0
        c_denom=x2x0sq*x1x0-x1x0sq*x2x0
        
        if c_denom == 0:
            print 'Divide by zero error'
            return 

        c=c_num/np.float(c_denom)
        if x1x0 == 0:
            print 'Divide by zero error'
            return

        b=(y1y0-c*x1x0sq)/np.float(x1x0)
        a=y[0]-b*x[0]-c*x[0]*x[0]
  
        fit=[a,b,c]
        if c == 0:
            xpeak=0.
            print 'Cannot find xpeak'
            return
        else:
            #Constrain the fit to be within these three points. 
            xpeak=-b/(2.0*c)
            if xpeak < x[0]:
                xpeak = np.float(x[0])
            if xpeak > x[2]:
                xpeak = np.float(x[2])
        
        return xpeak, fit
        
        
        

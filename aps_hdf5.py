# 
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
# 
#   Copyright (C) 2011 Mirna Lerotic - 2nd Look, Nicholas Schwarz
#   http://2ndlook.co/products.html
#   License: GNU GPL v3
#
#   Mantis is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   any later version.
#
#   Mantis is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details <http://www.gnu.org/licenses/>.

# File format details can be found on 
#https://confluence.aps.anl.gov/display/NX/Data+Exchange+Basics
 

from __future__ import division

import numpy as np
import scipy as scy
import h5py 

import data_struct


#----------------------------------------------------------------------
class h5:
    def __init__(self):
        pass
    
#----------------------------------------------------------------------
    def read_h5(self, filename, data_struct):
        
       
        # Open HDF5 file
        f = h5py.File(filename, 'r') 
    
 
        # Read basic definitions
        ds = f['implements']
        data_struct.implements = ds[...]
        ds = f['version']
        data_struct.version = ds[...]
        ds = f['file_creation_date']
        data_struct.file_creation_datetime = ds[...]
        ds = f['comment']
        data_struct.comment = ds[...]
           

        # base definition
    
        # Experimenter HDF5 group
        # /experimenter
        experimenterGrp = f['experimenter']
        dsname = experimenterGrp['name']
        dsrole = experimenterGrp['role']
        dsaffil = experimenterGrp['affiliation']
        dsaddr = experimenterGrp['address']
        dsphone = experimenterGrp['phone']
        dsem = experimenterGrp['email']
        dsuid = experimenterGrp['facility_user_id']
        data_struct.add_experimenter(name = dsname[...], 
                                     role = dsrole[...], 
                                     affiliation = dsaffil[...], 
                                     address = dsaddr[...], 
                                     phone = dsphone[...], 
                                     email = dsem[...], 
                                     facility_user_id = dsuid[...])
     
        # Sample HDF5 group
        # /sample
        sampleGrp = f["sample"]
        dsname = sampleGrp['name']
        dsdes = sampleGrp['description']
        dspdt = sampleGrp['preparation_date']
        dscf = sampleGrp['chemical_formula']
        data_struct.add_sample(name = dsname[...], 
                               description = dsdes[...], 
                               preparation_datetime = dspdt[...], 
                               chemical_formula = dscf[...])

    
        # Primary beam energy HDF5 group
        # /primary_beam_energy
        ds = f['primary_beam_energy']
        data_struct.primary_beam_energy = ds[...]
        data_struct.primary_beam_energy_units = f['primary_beam_energy'].attrs['units']  
    
    
    
        # exchange definition 
    
        # exchange HDF5 group
        # /exchange
        exchangeGrp = f['exchange']
        dscom = exchangeGrp['comment']
        dscdt = exchangeGrp['data_collection_datetime']
        data_struct.add_exchange(comment = dscom[...],
                                 collection_datetime = dscdt[...])
    
        # /exchange/detector
        detectorGrp = exchangeGrp['detector']      
        dsdata = detectorGrp['data']
        dseng = detectorGrp['energy']
        data_struct.exchange.add_detector( data = dsdata[...], 
                                           signal = dsdata.attrs['signal'], 
                                           units = dsdata.attrs['units'], 
                                           axes = dsdata.attrs['axes'], 
                                           energy = dseng[...], 
                                           energy_units = dseng.attrs['units'])
        
        axes = data_struct.exchange.detector[0].axes
        axes_list =  axes.split(':')
        
        
        for i in axes_list:
            ds = detectorGrp[i]
            data_struct.exchange.detector[0].add_dimscale(key = i, data = ds[...], 
                                                          units = ds.attrs['units'])
        

    
        # /exchange/white_data
        dswd = exchangeGrp['white_data']
        data_struct.exchange.add_white_data(white_data = dswd[...], 
                                            white_data_units = dswd.attrs['units'])
    
    
        # /exchange/dark_data
        dsdd = exchangeGrp['dark_data']
        data_struct.exchange.add_dark_data(dark_data = dsdd[...], 
                                           dark_data_units = dsdd.attrs['units'])
    
 
        # /exchange/instrument
        instrumentGrp = exchangeGrp['instrument']
        dsname = instrumentGrp['name']
        dsdes = instrumentGrp['description']
        data_struct.exchange.add_instrument(name = dsname[...], 
                                            description = dsdes[...])
    
    
        # /exchange/process
        processGrp = exchangeGrp['process']
        dsdt = processGrp['datetime']
        dspr = processGrp['program']
        dsver = processGrp['version']
        dspars = processGrp['parameters']
        data_struct.exchange.add_process(datetime = dsdt[...], 
                                            program = dspr[...], 
                                            version = dsver[...], 
                                            parameters = dspars[...])
   
    
        # spectromicroscopy definition
    
        # Spectromicroscopy HDF5 group
        # /spectromicroscopy
        spectromicroscopyGrp = f['spectromicroscopy']
        
        # /spectromicroscopy/positions
        dspos = spectromicroscopyGrp['positions']
        data_struct.spectromicroscopy.add_positions(positions = dspos[...], 
                                                    positions_units = dspos.attrs['units'], 
                                                    positions_names = dspos.attrs['names'])
        
    
        # /spectromicroscopy/beamline
        beamlineGrp = spectromicroscopyGrp['beamline']
        dsfn = beamlineGrp['facility_name']
        dsbn = beamlineGrp['beamline_name']
        dsrc = beamlineGrp['ring_current']
        dsbe = beamlineGrp['beam_energy']
        dsmono = beamlineGrp['monostripe']
        data_struct.spectromicroscopy.add_beamline(facility_name = dsfn[...], 
                                                   beamline_name = dsbn[...], 
                                                   ring_current = dsrc[...], 
                                                   ring_current_units = dsrc.attrs['units'], 
                                                   beam_energy = dsbe[...], 
                                                   beam_energy_units = dsbe.attrs['units'], 
                                                   monostripe = dsmono[...])
 
    
        # /spectromicroscopy/normalization
        normalizationGrp = spectromicroscopyGrp['normalization']
        dsws = normalizationGrp['white_spectrum']
        dswse = normalizationGrp['white_spectrum_energy']
        data_struct.spectromicroscopy.add_normalization(white_spectrum = dsws[...], 
                                                        white_spectrum_units = dsws.attrs['units'], 
                                                        white_spectrum_energy = dswse[...], 
                                                        white_spectrum_energy_units = dswse.attrs['units'])
        

    
        # Close
        f.close()
        
        self.absdata = data_struct.exchange.detector[0].data
        
        datadim = np.int32(self.absdata.shape)
        
        
        self.n_cols = datadim[0].copy()
        self.n_rows =  datadim[1].copy()
        self.ev = data_struct.exchange.detector[0].energy 
        self.n_ev = np.int32(self.ev.shape[0]).copy()
        
        npixels = self.n_cols*self.n_rows*self.n_ev
        
        
        self.x_dist = data_struct.exchange.detector[0].ds['x'] 
        self.y_dist = data_struct.exchange.detector[0].ds['y']
        
        self.i0data = data_struct.spectromicroscopy.normalization.white_spectrum           
        self.evi0 = data_struct.spectromicroscopy.normalization.white_spectrum_energy
        
        verbose = 0
        
        if verbose == 1: 
            print  'filename ', filename
            print 'File creation date ', data_struct.file_creation_datetime
            print 'Data array shape: ', self.absdata.shape
            print 'n columns ', self.n_cols
            print 'n_rows ', self.n_rows
            print 'n_ev ', self.n_ev
            print 'ev array ', self.ev
            #print 'x dist ', self.x_dist
            #print 'y_dist ', self.y_dist
            print 'type ', type 
            print 'i0 data ', self.i0data
            print 'evi0 ', self.evi0
        

                
        return
    
#----------------------------------------------------------------------
    def write_h5(self, filename, data_struct):
        
        # Open HDF5 file
        f = h5py.File(filename, 'w')
    
 
        # Create basic definitions
        ds = f.create_dataset('implements', data = data_struct.implements)
        ds = f.create_dataset('version', data = data_struct.version)
        ds = f.create_dataset('file_creation_date', data = data_struct.file_creation_datetime)
        ds = f.create_dataset('comment', data = data_struct.comment)
    

        # base definition
    
        # Experimenter HDF5 group
        # /experimenter
        experimenterGrp = f.create_group("experimenter")
        ds = experimenterGrp.create_dataset('name', data = data_struct.experimenter[0].name)
        ds = experimenterGrp.create_dataset('role', data = data_struct.experimenter[0].role)
        ds = experimenterGrp.create_dataset('affiliation', data = data_struct.experimenter[0].affiliation)
        ds = experimenterGrp.create_dataset('address', data = data_struct.experimenter[0].address)
        ds = experimenterGrp.create_dataset('phone', data = data_struct.experimenter[0].phone)
        ds = experimenterGrp.create_dataset('email', data= data_struct.experimenter[0].email)
        ds = experimenterGrp.create_dataset('facility_user_id', data = data_struct.experimenter[0].facility_user_id)
    
    
        # Sample HDF5 group
        # /sample
        sampleGrp = f.create_group("sample")
        ds = sampleGrp.create_dataset('name', data = data_struct.sample[0].name)
        ds = sampleGrp.create_dataset('description', data = data_struct.sample[0].description)
        ds = sampleGrp.create_dataset('preparation_date', data = data_struct.sample[0].preparation_datetime)
        ds = sampleGrp.create_dataset('chemical_formula', data = data_struct.sample[0].chemical_formula)
    
    
        # Primary beam energy HDF5 group
        # /primary_beam_energy
        ds = f.create_dataset('primary_beam_energy', data = data_struct.primary_beam_energy)
        ds.attrs['units'] =  data_struct.primary_beam_energy_units
    
    
    
        # exchange definition 
    
        # exchange HDF5 group
        # /exchange
        exchangeGrp = f.create_group("exchange")
    
        ds = exchangeGrp.create_dataset('comment', data = data_struct.exchange.comment)
    
    
        # /exchange/detector
        detectorGrp = exchangeGrp.create_group('detector')
        ds = detectorGrp.create_dataset('data', data = data_struct.exchange.detector[0].data)
        ds.attrs['signal'] = data_struct.exchange.detector[0].signal
        ds.attrs['units'] = data_struct.exchange.detector[0].units
        ds.attrs['axes'] = data_struct.exchange.detector[0].axes
        
        for k, v in data_struct.exchange.detector[0].ds.iteritems():          
            ds = detectorGrp.create_dataset(k, data = v)
            ds.attrs['units'] = data_struct.exchange.detector[0].ds_units[k]
        
        ds = detectorGrp.create_dataset('energy', data = data_struct.exchange.detector[0].energy)
        ds.attrs['units'] = data_struct.exchange.detector[0].energy_units
    
    
        # /exchange/white_data
        if data_struct.exchange.white_data:        
            ds = exchangeGrp.create_dataset('white_data', data = data_struct.exchange.white_data[0])
            ds.attrs['units'] = data_struct.exchange.white_data_units[0]
    
    
        # /exchange/dark_data
        if data_struct.exchange.dark_data:
            ds = exchangeGrp.create_dataset('dark_data', data = data_struct.exchange.dark_data[0])
            ds.attrs['units'] = data_struct.exchange.dark_data_units[0]
    
    
        # /exchange/data_collection_datetime
        ds = exchangeGrp.create_dataset('data_collection_datetime', data = data_struct.exchange.collection_datetime)
    
    
        # /exchange/instrument
        instrumentGrp = exchangeGrp.create_group('instrument')
        ds = instrumentGrp.create_dataset('name', data = data_struct.exchange.instrument.name)
        ds = instrumentGrp.create_dataset('description', data = data_struct.exchange.instrument.description)
    
    
        # /exchange/process
        processGrp = exchangeGrp.create_group('process')
        ds = processGrp.create_dataset('datetime', data = data_struct.exchange.process.datetime)
        ds = processGrp.create_dataset('program', data = data_struct.exchange.process.program)
        ds = processGrp.create_dataset('version', data = data_struct.exchange.process.version)
        ds = processGrp.create_dataset('parameters', data = data_struct.exchange.process.parameters)
    
    
    
        # spectromicroscopy definition
    
        # Spectromicroscopy HDF5 group
        # /spectromicroscopy
        spectromicroscopyGrp = f.create_group('spectromicroscopy')
        ds = spectromicroscopyGrp.create_dataset('positions', data = data_struct.spectromicroscopy.positions)
        ds.attrs['units'] = data_struct.spectromicroscopy.positions_units
        ds.attrs['names'] = data_struct.spectromicroscopy.positions_names
        
    
        # /spectromicroscopy/beamline
        beamlineGrp = spectromicroscopyGrp.create_group('beamline')
        ds = beamlineGrp.create_dataset('facility_name', data = data_struct.spectromicroscopy.beamline.facility_name)
        ds = beamlineGrp.create_dataset('beamline_name', data = data_struct.spectromicroscopy.beamline.beamline_name)
        ds = beamlineGrp.create_dataset('ring_current', data = data_struct.spectromicroscopy.beamline.ring_current)
        ds.attrs['units'] = data_struct.spectromicroscopy.beamline.ring_current_units
        ds = beamlineGrp.create_dataset('beam_energy', data = data_struct.spectromicroscopy.beamline.beam_energy)
        ds.attrs['units'] = data_struct.spectromicroscopy.beamline.beam_energy_units
        ds = beamlineGrp.create_dataset('monostripe', data = data_struct.spectromicroscopy.beamline.monostripe)
    
    
        # /spectromicroscopy/normalization
        normalizationGrp = spectromicroscopyGrp.create_group('normalization')
        ds = normalizationGrp.create_dataset('white_spectrum', 
                                             data = data_struct.spectromicroscopy.normalization.white_spectrum)
        ds.attrs['units'] = data_struct.spectromicroscopy.normalization.white_spectrum_units
        ds = normalizationGrp.create_dataset('white_spectrum_energy', 
                                             data = data_struct.spectromicroscopy.normalization.white_spectrum_energy)
        ds.attrs['units'] = data_struct.spectromicroscopy.normalization.white_spectrum_energy_units
    
    
        # Close
        f.close()



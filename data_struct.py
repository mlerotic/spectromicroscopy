# 
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
# 
#   Copyright (C) 2011 Mirna Lerotic, 2nd Look
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


#----------------------------------------------------------------------
class experimenter:
    def __init__(self, name = '', role = '', affiliation = '', address = '', phone = '',
                 email = '', facility_user_id = ''):

        self.name = name   
        self.role = role
        self.affiliation = affiliation 
        self.address = address  
        self.phone = phone
        self.email = email
        self.facility_user_id = facility_user_id 
        
        
#----------------------------------------------------------------------
class sample:
    def __init__(self, name = '', description = '', preparation_datetime = '', chemical_formula = ''):
        self.name = name
        self.description = description
#       preparation_date [string - ISO 8601 format]        
        self.preparation_datetime = preparation_datetime
#       chemical_formula [string - abbreviated CIF format]
        self.chemical_formula = chemical_formula
        
#----------------------------------------------------------------------
class detector:
    def __init__(self, data = 0, signal = 1, units = '', axes = '', energy = 0.0, energy_units = ''):

#        data [3-dim raw data, a stack of x-y images at z energies]
#           @signal [1]
#           @units [absorption or fluorescence]
#           @axes [x:y:z]
#        x [1D array of actual pixel positions in x; x dimension scale: 1..n_columns]
#           @units [microns or if not known integers]
#        y[1D array of pixel positions y; y dimension scale : 1..n_rows]
#           @units [none]
#        z [1D array of pixel positions z direction; z dimension scale: 1..n_energies]
#           @units [none]
#        energy [1D array of energy values corresponding to z-axis]
#           @units [keV]
        self.data = data
        self.signal = signal
        self.units = units
        self.axes = axes
        
        self.energy = energy
        self.energy_units = energy_units
        
#       for axes create a dictionary with lookup values 
        axes_list =  self.axes.split(':')
        #Dimension scales dictionary - look into how you could pass them in the argument list
        self.ds = {}
        self.ds_units = {}
                       
    
            
#---------------------------------------------------------------------- 
# Add dimension scales for the detector - the arguments must match the specification in @axes!
    def add_dimscale(self, key = '', data = 0, units = ''):      

        self.ds[key] = data
        self.ds_units[key] = units

            
#----------------------------------------------------------------------
class instrument:
    def __init__(self, name ='', description = ''):
#            name [Argonne APS 2ID beamline] 
        self.name = name
#            description [hard x-ray fluorescence microprobe]
        self.description = description

#----------------------------------------------------------------------
class process:
    def __init__(self, datetime = '', program = '', version = '', parameters = ''):   
#            date [string]
        self.datetime = datetime
#            program [string]
        self.program = program
#            version [string]  
        self.version = version          
#            parameters [string]   
        self.parameters = parameters      

#----------------------------------------------------------------------
class exchange:
    def __init__(self, comment = '',  collection_datetime = ''):
        pass
    
#        comment [Another comment.]  
        self.comment = comment  
#        detector/[HDF5 group]
        self.detector = []
#        white_data []
        self.white_data = []
#            @units 
        self.white_data_units = []
#        dark_data [1-n? dim array]
        self.dark_data = []
#            @units
        self.dark_data_units = []
#        collection_time [string]
        self.collection_datetime = collection_datetime
#        instrument/ [HDF5 group]
        self.instrument = instrument()
#        process/ [HDF5 group]
        self.process = process()
        
#---------------------------------------------------------------------- 
    def add_detector(self, data = 0, signal = 1, units = '', axes = '', energy = 0.0, energy_units = ''):   
        
        new_detector = detector(data=data, signal=signal, axes=axes, units = units, energy = energy,
                                energy_units = energy_units)
        self.detector.append(new_detector)
        
#---------------------------------------------------------------------- 
    def add_white_data(self, white_data = 0.0, white_data_units=''):  

        self.white_data.append(white_data)
        self.white_data_units.append(white_data_units)
        
#---------------------------------------------------------------------- 
    def add_dark_data(self, dark_data = 0, dark_data_units =''):  

        self.dark_data.append(dark_data)
        self.dark_data_units.append(dark_data_units)
        
#---------------------------------------------------------------------- 
    def add_instrument(self, name = '', description = ''):     
        self.instrument = instrument(name=name, description=description)

#---------------------------------------------------------------------- 
    def add_process(self, datetime = '', program = '', version = '', parameters = ''):   
        self.process = process(datetime=datetime, program=program, version=version,
                               parameters=parameters)

 
 #----------------------------------------------------------------------
class beamline:
    def __init__(self, facility_name='', beamline_name ='', ring_current = 0, ring_current_units='',
                 beam_energy = 0, beam_energy_units='', monostripe = ''):
        self.facility_name = facility_name
        self.beamline_name = beamline_name
        self.ring_current = ring_current
        self.ring_current_units = ring_current_units
        self.beam_energy = beam_energy
        self.beam_energy_units = beam_energy_units
        self.monostripe = monostripe
        
 #----------------------------------------------------------------------
class normalization:
    def __init__(self, white_spectrum=0.0, white_spectrum_units = '', white_spectrum_energy=0.0,
                          white_spectrum_energy_units = ''):
        self.white_spectrum = white_spectrum
        self.white_spectrum_units = white_spectrum_units
        self.white_spectrum_energy=white_spectrum_energy
        self.white_spectrum_energy_units=white_spectrum_energy_units
               
#----------------------------------------------------------------------
class spectromicroscopy:
    def __init__(self):
        
        self.positions = []
        self.positions_units = []
        self.positions_names = []
        
        self.beamline = beamline()
        self.normalization = normalization()
        
        self.optical_density = 0
        
#---------------------------------------------------------------------- 
    def add_positions(self,  positions = 0, positions_units = '', positions_names = ''):   

        self.positions.append(positions)
        self.positions_units.append(positions_units)
        self.positions_names.append(positions_names)
        
#---------------------------------------------------------------------- 
    def add_beamline(self,  facility_name='', beamline_name ='', ring_current = 0, ring_current_units='',
                 beam_energy = 0, beam_energy_units='', monostripe = ''):  
        self.beamline = beamline(facility_name=facility_name, beamline_name =beamline_name, 
                                 ring_current = ring_current, ring_current_units=ring_current_units,
                                 beam_energy = beam_energy, beam_energy_units=beam_energy_units, 
                                 monostripe = monostripe)
        
#---------------------------------------------------------------------- 
    def add_normalization(self,  white_spectrum=0.0, white_spectrum_units = '', white_spectrum_energy=0.0,
                          white_spectrum_energy_units = ''):      
        self.normalization = normalization(white_spectrum=white_spectrum, 
                                           white_spectrum_units=white_spectrum_units, 
                                           white_spectrum_energy=white_spectrum_energy, 
                                           white_spectrum_energy_units=white_spectrum_energy_units)     


        
#----------------------------------------------------------------------
class h5:
    def __init__(self, implements = '', version = '', file_creation_datetime = '', comment = '',
                 primary_beam_energy = 0.0, primary_beam_energy_units = ''):
        
        #implements [string] comma separated string that tells the user which entries file contains
        self.implements = implements
#        version [1.01] file version 
        self.version = version 
        #file_creation_date [string - ISO 8601 format]
        self.file_creation_datetime = file_creation_datetime
        #comment [string]
        self.comment = comment        
#        primary_beam_energy [double]
#            @units[string] 
        self.primary_beam_energy = primary_beam_energy
        self.primary_beam_energy_units = primary_beam_energy_units


# Make H5 group initialization into function calls

#        experimenter/ [HDF5 group]
        self.experimenter = []
        #can have experimenter1, experimenter 2....
        
#        sample/ [HDF5 group]
        self.sample = []    

#       exchange/[HDF5 group]
        self.exchange = exchange()
        
#       spectromicroscopy/[HDF5 group]
        self.spectromicroscopy = spectromicroscopy()
        
      #---------------------------------------------------------------------- 
    def delete_data(self):  
                #implements [string] comma separated string that tells the user which entries file contains
        self.implements = ' '
#        version [1.01] file version 
        self.version = ' '
        #file_creation_date [string - ISO 8601 format]
        self.file_creation_datetime = ' '
        #comment [string]
        self.comment = ' '       
#        primary_beam_energy [double]
#            @units[string] 
        self.primary_beam_energy = 0
        self.primary_beam_energy_units = ' '


# Make H5 group initialization into function calls

#        experimenter/ [HDF5 group]
        self.experimenter = []
        #can have experimenter1, experimenter 2....
        
#        sample/ [HDF5 group]
        self.sample = []    

#       exchange/[HDF5 group]
        self.exchange = exchange()
        
#       spectromicroscopy/[HDF5 group]
        self.spectromicroscopy = spectromicroscopy()
        
        
#---------------------------------------------------------------------- 
    def add_experimenter(self,  name = '', role = '', affiliation = '', address = '', phone = '',
                 email = '', facility_user_id = ''):   
        
        newexp = experimenter(name=name, role = role, affiliation = affiliation, address = address, 
                           phone = phone, email = email, facility_user_id=facility_user_id)
        
        self.experimenter.append(newexp)
        
#---------------------------------------------------------------------- 
    def add_sample(self,  name = '', description = '', preparation_datetime = '', chemical_formula = ''):   
        
        new_sample = sample(name=name, description=description, preparation_datetime=preparation_datetime,
                            chemical_formula=chemical_formula)
        self.sample.append(new_sample)
        
#---------------------------------------------------------------------- 
    def add_exchange(self,  comment = '', collection_datetime = ''):   
        
        self.exchange = exchange(comment=comment, collection_datetime=collection_datetime)
        


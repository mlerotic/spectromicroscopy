# 
#   This file is part of Mantis, a Multivariate ANalysis Tool for Spectromicroscopy.
# 
#   Copyright (C) 2011 Mirna Lerotic, 2nd Look
#   http://2ndlookconsulting.com
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


'''
! This is from file COMPOUND.DAT for HENKE.  Entries are:
! compound    formula                    density (g/cc)
water           H2O                                 1.0
protein         H48.6C32.9N8.9O8.9S0.6              1.35
lipid           H62.5C31.5O6.3                      1.0
nucleosome      H42.1C31.9N10.3O13.9P1.6S0.3        1.5
dna             H35.5C30.8N11.7O18.9P3.1            1.7
helium          He                                  0.1663e-3
chromatin       H49.95C24.64N8.66O15.57P1.07S0.03   1.527
air             N78.08O20.95Ar0.93                  1.2047e-3
pmma            C5H8O2                              1.18
nitride         Si3N4                               3.44
graphite        C                                   2.26
nickel          Ni                                  8.876
beryl           Be                                  1.845
copper          Cu                                  8.96
quartz          SiO2                                2.20
aluminum        Al                                  2.70
gold            Au                                  19.3
ice             H2O                                 0.92
carbon          C                                   1.0
polystyrene     C8H8                                1.06
silicon         Si                                  2.33
germanium       Ge                                  5.323
'''
from __future__ import division

import numpy as np
import scipy.interpolate
from xdrlib import *
import string

#-----------------------------------------------------------------------------
class henke:
    def __init__(self):   
        
        self.compound_name = [ 'water' , 'protein', 'lipid', 'nucleosome', 'dna', 'helium', 'chromatin', 
                               'air', 'pmma', 'nitride', 'graphite', 'nickel', 'beryl', 'copper', 
                               'quartz', 'aluminum', 'gold', 'ice', 'carbon', 'polystyrene', 
                               'silicon', 'germanium']
        
        self.compound_forumula = ['H2O', 'H48.6C32.9N8.9O8.9S0.6', 'H62.5C31.5O6.3', 
                                   'H42.1C31.9N10.3O13.9P1.6S0.3', 'H35.5C30.8N11.7O18.9P3.1',
                                   'He' , 'H49.95C24.64N8.66O15.57P1.07S0.03', 
                                   'N78.08O20.95Ar0.93', 'C5H8O2', 'Si3N4', 'C', 'Ni', 'Be', 
                                   'Cu', 'SiO2', 'Al', 'Au', 'H2O', 'C', 'C8H8',
                                   'Si', 'Ge']
        self.compound_density = [ 1.0, 1.35, 1.0, 1.5, 1.7, 1.66E-04, 1.527, 1.20E-03, 1.18,
                                 3.44, 2.26, 8.876, 1.845, 8.96, 2.2, 2.7, 19.3, 0.92, 1, 1.06,
                                 2.33, 5.323 ]



#-----------------------------------------------------------------------------    
    def compound(self, compound_string, density): 
        
#    Type in a chemical formula and the density.  The routine
#    then figures out the atomic weight in g/mole of "molecules"
#    and how many of what Z atoms you have.   

        z_array = []
        
        
        if compound_string in self.compound_name:
            compound_string = self.compound_forumula[self.compound_name.index(compound_string)] 

        if compound_string in self.compound_forumula: 
            z_array = self.zcompound(compound_string, z_array)
            atwt =  self.zatwt(z_array)

        return z_array, atwt
            
#-----------------------------------------------------------------------------               
    def zcompound(self, compound_string, z_array, paren_multiplier=False):
        
        verbose = False
        if verbose: print ('compound_string', compound_string)
        if paren_multiplier == False:
            z_array = np.zeros(92)
            paren_multiplier=1.

      
        max_z_index=93
        
        last_char_index = len(compound_string) - 1
        
        # If we don't start off with a parenthesis, all we have to do
        # is strip off the first element and process it.  We then
        # call the routine over again to handle the next part of
        # the string...
        if (compound_string[0] != '(') :
            # Look to see if the string has an element
            # like "C" or like "He".
            first_char=compound_string[0]
            if len(compound_string) > 1: 
                second_char = compound_string[1]
            else: 
                second_char = ''
            this_element_name = first_char
      
            if ((second_char >= 'a') and (second_char <= 'z')) :
                this_element_name = this_element_name + second_char
                num_start_index = 2
            else:
                this_element_name = this_element_name + ' '
                num_start_index = 1
      
        if verbose: print ('this_element_name:',this_element_name,', num_start_index:', num_start_index)
        
        this_z=0
        if   this_element_name == 'H ': this_z=1
        elif this_element_name == 'He': this_z=2
        elif this_element_name == 'Li': this_z=3
        elif this_element_name == 'Be': this_z=4
        elif this_element_name == 'B ': this_z=5
        elif this_element_name == 'C ': this_z=6
        elif this_element_name == 'N ': this_z=7
        elif this_element_name == 'O ': this_z=8
        elif this_element_name == 'F ': this_z=9
        elif this_element_name == 'Ne': this_z=10
        elif this_element_name == 'Na': this_z=11
        elif this_element_name == 'Mg': this_z=12
        elif this_element_name == 'Al': this_z=13
        elif this_element_name == 'Si': this_z=14
        elif this_element_name == 'P ': this_z=15
        elif this_element_name == 'S ': this_z=16
        elif this_element_name == 'Cl': this_z=17
        elif this_element_name == 'Ar': this_z=18
        elif this_element_name == 'K ': this_z=19
        elif this_element_name == 'Ca': this_z=20
        elif this_element_name == 'Sc': this_z=21
        elif this_element_name == 'Ti': this_z=22
        elif this_element_name == 'V ': this_z=23
        elif this_element_name == 'Cr': this_z=24
        elif this_element_name == 'Mn': this_z=25
        elif this_element_name == 'Fe': this_z=26
        elif this_element_name == 'Co': this_z=27
        elif this_element_name == 'Ni': this_z=28
        elif this_element_name == 'Cu': this_z=29
        elif this_element_name == 'Zn': this_z=30
        elif this_element_name == 'Ga': this_z=31
        elif this_element_name == 'Ge': this_z=32
        elif this_element_name == 'As': this_z=33
        elif this_element_name == 'Se': this_z=34
        elif this_element_name == 'Br': this_z=35
        elif this_element_name == 'Kr': this_z=36
        elif this_element_name == 'Rb': this_z=37
        elif this_element_name == 'Sr': this_z=38
        elif this_element_name == 'Y ': this_z=39
        elif this_element_name == 'Zr': this_z=40
        elif this_element_name == 'Nb': this_z=41
        elif this_element_name == 'Mo': this_z=42
        elif this_element_name == 'Tc': this_z=43
        elif this_element_name == 'Ru': this_z=44
        elif this_element_name == 'Rh': this_z=45
        elif this_element_name == 'Pd': this_z=46
        elif this_element_name == 'Ag': this_z=47
        elif this_element_name == 'Cd': this_z=48
        elif this_element_name == 'In': this_z=49
        elif this_element_name == 'Sn': this_z=50
        elif this_element_name == 'Sb': this_z=51
        elif this_element_name == 'Te': this_z=52
        elif this_element_name == 'I ': this_z=53
        elif this_element_name == 'Xe': this_z=54
        elif this_element_name == 'Cs': this_z=55
        elif this_element_name == 'Ba': this_z=56
        elif this_element_name == 'La': this_z=57
        elif this_element_name == 'Ce': this_z=58
        elif this_element_name == 'Pr': this_z=59
        elif this_element_name == 'Nd': this_z=60
        elif this_element_name == 'Pm': this_z=61
        elif this_element_name == 'Sm': this_z=62
        elif this_element_name == 'Eu': this_z=63
        elif this_element_name == 'Gd': this_z=64
        elif this_element_name == 'Tb': this_z=65
        elif this_element_name == 'Dy': this_z=66
        elif this_element_name == 'Ho': this_z=67
        elif this_element_name == 'Er': this_z=68
        elif this_element_name == 'Tm': this_z=69
        elif this_element_name == 'Yb': this_z=70
        elif this_element_name == 'Lu': this_z=71
        elif this_element_name == 'Hf': this_z=72
        elif this_element_name == 'Ta': this_z=73
        elif this_element_name == 'W ': this_z=74
        elif this_element_name == 'Re': this_z=75
        elif this_element_name == 'Os': this_z=76
        elif this_element_name == 'Ir': this_z=77
        elif this_element_name == 'Pt': this_z=78
        elif this_element_name == 'Au': this_z=79
        elif this_element_name == 'Hg': this_z=80
        elif this_element_name == 'Tl': this_z=81
        elif this_element_name == 'Pb': this_z=82
        elif this_element_name == 'Bi': this_z=83
        elif this_element_name == 'Po': this_z=84
        elif this_element_name == 'At': this_z=85
        elif this_element_name == 'Rn': this_z=86
        elif this_element_name == 'Fr': this_z=87
        elif this_element_name == 'Ra': this_z=88
        elif this_element_name == 'Ac': this_z=89
        elif this_element_name == 'Th': this_z=90
        elif this_element_name == 'Pa': this_z=91
        elif this_element_name == 'U ': this_z=92
        else: this_z=0
        
            
        if (this_z == 0) :
            print ('zcompound is confused: ',compound_string)
            compound_string=''
            z_array = 0
            return
        
        
        
        # Find the next element or parenthesis, as
        # anything before it must be a number.
        postnum_index = num_start_index
        if len(compound_string) > num_start_index+1 :
            test_char = compound_string[postnum_index]
        else :
            test_char = ''
        while ( ((test_char == '0') or (test_char == '1') or \
                 (test_char == '2') or (test_char == '2') or \
                 (test_char == '3') or (test_char == '4') or \
                 (test_char == '5') or (test_char == '6') or \
                 (test_char == '7') or (test_char == '8') or \
                 (test_char == '9') or (test_char == '.')) and \
                 (postnum_index <= last_char_index) ) :
            postnum_index=postnum_index+1
            if (postnum_index <= last_char_index) :
                test_char=compound_string[postnum_index]
            else:
                test_char=''
              
        # is there more?
        if (num_start_index != postnum_index) :
            number_string=compound_string[num_start_index:postnum_index]
            num_multiplier=1.
            if verbose: print( 'Trying to interpret ',number_string,' as a number.')
            if (len(number_string) != 0) :
                num_multiplier = float(number_string)
        else:
            num_multiplier=1.

        # We've handled this element, so pop it into the
        # matrix and continue.
        if (this_z <= max_z_index) :
            z_array[this_z-1] = z_array[this_z-1] + num_multiplier 
        else:
            print ('zcompound: z_array smaller than ',max_z_index)
            z_array = 0
            return

        # And deal with what's left
        remaining_string=compound_string[postnum_index:last_char_index+1]
 
        if len(remaining_string) > 0:
            z_array = self.zcompound(remaining_string,z_array,paren_multiplier=True)
        
        return z_array
    
    
#-----------------------------------------------------------------------------    
    def zatwt(self, z_array): 
        
        maxz=z_array.size
        atwt=0.

        for i in range(maxz):
            if (z_array[i] != 0.):
                if   i+1 == 1: this_atwt=1.00794
                elif i+1 == 2: this_atwt=4.0026
                elif i+1 == 3: this_atwt=6.941
                elif i+1 == 4: this_atwt=9.01218
                elif i+1 == 5: this_atwt=10.81
                elif i+1 == 6: this_atwt=12.011
                elif i+1 == 7: this_atwt=14.0067
                elif i+1 == 8: this_atwt=15.9994
                elif i+1 == 9: this_atwt=18.9984
                elif i+1 == 10: this_atwt=21.179
                elif i+1 == 11: this_atwt=22.98977
                elif i+1 == 12: this_atwt=24.305
                elif i+1 == 13: this_atwt=26.98154
                elif i+1 == 14: this_atwt=28.0855
                elif i+1 == 15: this_atwt=30.97376
                elif i+1 == 16: this_atwt=32.06
                elif i+1 == 17: this_atwt=35.453
                elif i+1 == 18: this_atwt=39.948
                elif i+1 == 19: this_atwt=39.0983
                elif i+1 == 20: this_atwt=40.08
                elif i+1 == 21: this_atwt=44.9559
                elif i+1 == 22: this_atwt=47.88
                elif i+1 == 23: this_atwt=50.9415
                elif i+1 == 24: this_atwt=51.996
                elif i+1 == 25: this_atwt=54.9380
                elif i+1 == 26: this_atwt=55.847
                elif i+1 == 27: this_atwt=58.9332
                elif i+1 == 28: this_atwt=58.69
                elif i+1 == 29: this_atwt=63.546
                elif i+1 == 30: this_atwt=65.38
                elif i+1 == 31: this_atwt=69.72
                elif i+1 == 32: this_atwt=72.59
                elif i+1 == 33: this_atwt=74.9216
                elif i+1 == 34: this_atwt=78.96
                elif i+1 == 35: this_atwt=79.904
                elif i+1 == 36: this_atwt=83.80
                elif i+1 == 37: this_atwt=85.4678
                elif i+1 == 38: this_atwt=87.62
                elif i+1 == 39: this_atwt=88.9059
                elif i+1 == 40: this_atwt=91.22
                elif i+1 == 41: this_atwt=92.9064
                elif i+1 == 42: this_atwt=95.94
                elif i+1 == 43: this_atwt=98.
                elif i+1 == 44: this_atwt=101.07
                elif i+1 == 45: this_atwt=102.9055
                elif i+1 == 46: this_atwt=106.42
                elif i+1 == 47: this_atwt=107.8682
                elif i+1 == 48: this_atwt=112.41
                elif i+1 == 49: this_atwt=114.82
                elif i+1 == 50: this_atwt=118.69
                elif i+1 == 51: this_atwt=121.75
                elif i+1 == 52: this_atwt=127.60
                elif i+1 == 53: this_atwt=126.9054
                elif i+1 == 54: this_atwt=131.29
                elif i+1 == 55: this_atwt=132.9054
                elif i+1 == 56: this_atwt=137.33
                elif i+1 == 57: this_atwt=138.9055
                elif i+1 == 58: this_atwt=140.12
                elif i+1 == 59: this_atwt=140.9077
                elif i+1 == 60: this_atwt=144.24
                elif i+1 == 61: this_atwt=145.
                elif i+1 == 62: this_atwt=150.36
                elif i+1 == 63: this_atwt=151.96
                elif i+1 == 64: this_atwt=157.25
                elif i+1 == 65: this_atwt=158.9254
                elif i+1 == 66: this_atwt=162.5
                elif i+1 == 67: this_atwt=164.9304
                elif i+1 == 68: this_atwt=167.26
                elif i+1 == 69: this_atwt=168.9342
                elif i+1 == 70: this_atwt=173.04
                elif i+1 == 71: this_atwt=174.967
                elif i+1 == 72: this_atwt=178.49
                elif i+1 == 73: this_atwt=180.9479
                elif i+1 == 74: this_atwt=183.85
                elif i+1 == 75: this_atwt=186.207
                elif i+1 == 76: this_atwt=190.2
                elif i+1 == 77: this_atwt=192.22
                elif i+1 == 78: this_atwt=195.08
                elif i+1 == 79: this_atwt=196.9665
                elif i+1 == 80: this_atwt=200.59
                elif i+1 == 81: this_atwt=204.383
                elif i+1 == 82: this_atwt=207.2
                elif i+1 == 83: this_atwt=208.9804
                elif i+1 == 84: this_atwt=209.
                elif i+1 == 85: this_atwt=210.
                elif i+1 == 86: this_atwt=222.
                elif i+1 == 87: this_atwt=223.
                elif i+1 == 88: this_atwt=226.0254
                elif i+1 == 89: this_atwt=227.0278
                elif i+1 == 90: this_atwt=232.0381
                elif i+1 == 91: this_atwt=231.0359
                elif i+1 == 92: this_atwt=238.0289
                else: this_atwt=0.

                atwt=atwt+z_array[i]*this_atwt

        return atwt


#-----------------------------------------------------------------------------    
    def extra(self, ielement = -1): 
        
        energies, f1, f2, n_extra, energies_extra, f1_extra, f2_extra = self.read(ielement, all = False)
        
        if n_extra != 0:
            energies_all=np.concatenate((energies,energies_extra), axis=0)
            f1_all=np.concatenate((f1,f1_extra), axis=0)
            f2_all=np.concatenate((f2,f2_extra), axis=0)
            sort_order=energies_all.argsort()
            energies_all=energies_all[sort_order]
            f1_all=f1_all[sort_order]
            f2_all=f2_all[sort_order]
        else:
            energies_all=energies
            f1_all=f1
            f2_all=f2


        return energies, f1, f2, energies_extra, f1_extra, f2_extra

#-----------------------------------------------------------------------------    
    def read(self, ielement = -1, all = True): 
        
        
        # If we don't specifiy element return all energies 
        if ielement == -1 : 
            all = True
        
        verbose = False
        expected_pos = 0
        
        filename = 'henke.xdr'
        file = open(str(filename),'rb')
#         try:
#             file = open(str(filename),'rb')
#         except:
#             print ('Could not open file ', filename)
#             return -1
        
        if verbose:
            print ('File: ', filename)
    
        buf = file.read()        
        u = Unpacker(buf)
        
        if all:
            n_elements = u.unpack_int()
            n_energies = u.unpack_int()
        
            if verbose:
                print ('n_energies: ', n_energies)
                print ('n_elements: ', n_elements)
                expected_pos = expected_pos+2*4
                print ('Actual, expected file position before reading energies: ' ,u.get_position(), expected_pos)
            
            
            energies = u.unpack_farray(n_energies, u.unpack_float)   
            energies = np.array(energies) 
            if verbose:
                print ('energies: ', energies)
            
            f1 = np.zeros((n_elements, n_energies))
            f2 = np.zeros((n_elements, n_energies))
            this_f1 = np.zeros((n_energies))
            this_f2 = np.zeros((n_energies))

            if verbose:        
                expected_pos = expected_pos+4*n_energies
                print ('Actual, expected file position before reading elements: ',u.get_position(), expected_pos)
        
            for i_element in range(n_elements):
                this_f1 = u.unpack_farray(n_energies, u.unpack_float)
                this_f2 = u.unpack_farray(n_energies, u.unpack_float)
                f1[i_element, :] = this_f1
                f2[i_element, :] = this_f2
            
            if verbose:
                expected_pos =expected_pos+n_elements*n_energies*2*4
                print ('Actual, expected file position before reading n_extra_energies: ', u.get_position(), expected_pos)

            n_extra_energies  = u.unpack_int()
            if verbose:
                print ('n_extra_energies: ', n_extra_energies)
            
            if verbose:        
                expected_pos = expected_pos+4
                print ('Actual, expected file position before reading extras: ',u.get_position(), expected_pos)
            
            
            n_extra = np.zeros((n_elements), dtype = np.int)
            extra_energies = np.zeros((n_elements, n_extra_energies))
            extra_f1 = np.zeros((n_elements, n_extra_energies))
            extra_f2 = np.zeros((n_elements, n_extra_energies))
            this_n_extra = 0
            this_extra_energies = np.zeros((n_extra_energies))
            this_extra_f1 = np.zeros((n_extra_energies))
            this_extra_f2 = np.zeros((n_extra_energies))
        
            for i_element in range(n_elements):
                this_n_extra = u.unpack_int()
                this_extra_energies = u.unpack_farray(n_extra_energies, u.unpack_float)
                this_extra_f1 = u.unpack_farray(n_extra_energies, u.unpack_float)
                this_extra_f2 = u.unpack_farray(n_extra_energies, u.unpack_float)
                n_extra[i_element] = this_n_extra
                extra_energies[i_element, :] = this_extra_energies
                extra_f1[i_element, :] = this_extra_f1
                extra_f2[i_element, :] = this_extra_f2
                
        else:
            n_elements = u.unpack_int()
            n_energies = u.unpack_int()           
            
            energies = u.unpack_farray(n_energies, u.unpack_float)   
            energies = np.array(energies)  
            if verbose:
                print ('energies: ', energies)
            
            byte_offset = 4+4+4*n_energies + 8*ielement*n_energies
            u.set_position(byte_offset)
            
            f1 = u.unpack_farray(n_energies, u.unpack_float)
            f2 = u.unpack_farray(n_energies, u.unpack_float)
            
            byte_offset = 4+4+4*n_energies + 8*n_elements*n_energies
            u.set_position(byte_offset)
            
            n_extra_energies = u.unpack_int()
            if verbose: print ('n_extra_energies ', n_extra_energies)
            
            
            # Now we have the above plus i_element times the quantity:
            #   (2 for n_extra, and n_extra_energies each of three floats)
            byte_offset = long(4)+long(4)+long(4)*n_energies + long(8)*n_elements*n_energies + long(4) + ielement*(4+12*n_extra_energies)
            u.set_position(byte_offset)
            
            n_extra = u.unpack_int()
            this_extra_energies = u.unpack_farray(n_extra_energies, u.unpack_float)
            this_extra_f1 = u.unpack_farray(n_extra_energies, u.unpack_float)
            this_extra_f2 = u.unpack_farray(n_extra_energies, u.unpack_float)
            
            extra_energies = this_extra_energies[0:n_extra]
            extra_f1 = this_extra_f1[0:n_extra]
            extra_f2 = this_extra_f2[0:n_extra]
                                                
        file.close()
        
        
        return energies, f1, f2, n_extra, extra_energies, extra_f1, extra_f2
    
    
    
#-----------------------------------------------------------------------------    
    def array(self, compound_name, density, graze_mrad = 0):

        z_array, atwt = self.compound(compound_name,density)

        maxz = 92
        first_time = 1
        for i in range(maxz):
            if (z_array[i] != 0.):
                energies, this_f1, this_f2, n_extra, extra_energies, extra_f1, extra_f2 = self.read(ielement=i)
            if (first_time == 1) :
                f1 = z_array[i]*this_f1
                f2 = z_array[i]*this_f2
                first_time = 0
            else:
                f1 = f1+z_array[i]*this_f1
                f2 = f2+z_array[i]*this_f2

  
        num_energies = len(energies)
        AVOGADRO=6.02204531e23
        HC_ANGSTROMS=12398.52
        RE=2.817938070e-13            # in cm
  
        if (atwt != 0.0) :
            molecules_per_cc = density * AVOGADRO / atwt
        else:
            molecules_per_cc = 0.0

  
        wavelength_angstroms = HC_ANGSTROMS/energies
        # This constant has wavelength in angstroms and then
        # they are converted to centimeters.
        constant = RE * (1.0e-16 * wavelength_angstroms * wavelength_angstroms) * molecules_per_cc / (2.0 * np.math.pi)
        delta = constant * f1
        beta = constant * f2
        # Alpha is in inverse meters squared
        alpha = 1.e4 * density * AVOGADRO * RE / (2.*np.math.pi*atwt)
        #alpha = alpha[0]

  
        if (graze_mrad == 0.):
            reflect=np.ones((num_energies))
        else:
            theta = 1.0e-3 * graze_mrad
            sinth = np.sin(theta)
            sinth2 = sinth * sinth
            coscot = np.cos(theta)
            coscot = coscot * coscot / sinth
            alpha = 2.0 * delta - delta * delta + beta * beta
            gamma = 2.0 * (1.0 - delta) * beta
            rhosq = 0.5 * (sinth2 - alpha + np.sqrt((sinth2 - alpha)*(sinth2-alpha) + gamma*gamma) )
            rho = np.sqrt(rhosq)
            i_sigma = (4.0 * rhosq * (sinth - rho) * (sinth - rho) + \
                 gamma * gamma) / \
                 (4.0 * rhosq * (sinth + rho) * (sinth + rho) + \
                  gamma * gamma)
            piosig = (4.0 * rhosq * (rho - coscot) * (rho - coscot) + \
                gamma * gamma) / \
                (4.0 * rhosq * (rho + coscot) * (rho + coscot) + \
                 gamma * gamma)
            reflect= 50.0 * i_sigma * (1 + piosig)


        denom = energies*4.*np.math.pi*beta

        zeroes = np.where(denom == 0.)
        nonzeroes = np.where(denom != 0.)  
        denom[zeroes] = 1e-8     
        
        inverse_mu = np.array((len(energies)))

        inverse_mu = 1.239852/denom
        if (len(zeroes) > 0) :
            inverse_mu[zeroes] = np.inf


        return energies, f1, f2, delta, beta, graze_mrad, reflect, inverse_mu, atwt, alpha
  


#-----------------------------------------------------------------------------    
    def dose_calc(self, stack, i_composition, od_spectrum, i0_signal, dosecalc_detector_eff):    
    
        pix_nm_squared = 1.e6*(stack.x_dist[1]-stack.x_dist[0])*(stack.y_dist[1]-stack.y_dist[0])

        dose = 0.
             

        #rho doesn't matter here
        rho = 1.
        henke_energies, f1, f2, delta, beta, graze_mrad, reflect, inverse_mu, atwt, alpha = self.array(i_composition,rho) 

        func_f2_array = scipy.interpolate.interp1d(henke_energies, f2,  bounds_error=False, fill_value=0.0) 
        f2_array = func_f2_array(stack.ev)
        
        # This is the scaling of f2_array as described in pca_gui_man.tex
        i_max_ev = np.argmax(stack.ev)
        if stack.data_dwell == None: 
            stack.data_dwell = 1.
        if (od_spectrum[i_max_ev] != 0.) :
            f2_array = f2_array[0,i_max_ev]*(stack.ev/stack.ev[i_max_ev])*(od_spectrum/od_spectrum[i_max_ev])
            # See the documentation file pca_gui_man.tex. Is i0_signal
            # normalized to 1 sec? Does not seem to be the case! Might
            # have to change it at some point.
            dose  = 6.74e5*(1./(dosecalc_detector_eff* pix_nm_squared*atwt))*np.sum(i0_signal*1e-3*stack.data_dwell *f2_array)


        return dose

# coding: utf-8

"""
This provides functions for dealing with Materials Project data
"""

__author__ = "Mehmet Topsakal"
__email__ = "metokal@gmail.com"
__status__ = "Development"
__date__ = "March 20, 2018"

import os
# import sys
import numpy as np
import pickle
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from . import t4iss_defaults
from .core import mXANES


# This searches for structures in MP database.
def search_MP(mpr, search_pattern, nmax=None):

    if nmax is None:
        nmax = 11

    mpid_list = []
    data = mpr.get_data(search_pattern, data_type="vasp", prop="nsites")
    for i in range(len(data)):
        if data[i]['nsites'] <= nmax:
            mpid_list.append(data[i]['material_id'])

    l = len(mpid_list)
    print("""Found %i structures""" % l)

    return mpid_list


# ICSD 0nly
def search_MP_icsd(mpr, search_pattern, nmax=None):

    if nmax is None:
        nmax = 11

    found0 = []
    data = mpr.get_data(search_pattern, data_type="vasp", prop="nsites")
    for i in range(len(data)):
        if data[i]['nsites'] <= nmax:
            found0.append(data[i]['material_id'])

    mpid_list_icsd = []
    for i in found0:
        data = mpr.get_data(i, data_type="vasp")
        if data[0]['icsd_ids']:
            mpid_list_icsd.append(i)

    l = len(mpid_list_icsd)
    print("""Found %i structures""" % l)

    return mpid_list_icsd




# This downloads XANES data from MP based on mpid
def download_xanes_from_MP(mpr, mpid, absorption_specie, 
                           download_to=None, return_missing=False):

    here = os.getcwd()

    if download_to is None:
        download_to = t4iss_defaults['t4iss_xanes_data']

    try:
        data = mpr.get_data(mpid, data_type="feff", prop="xas")
        if data[0]['xas']:
            species_in_data = []
            for xas_doc in data[0]['xas']:
                data_abs_specie = xas_doc['structure'].species[xas_doc['absorbing_atom']].name
                species_in_data.append(data_abs_specie)
            if absorption_specie in species_in_data:
                available = True
            else:
                available = False
        else:
            available = False            
            
    except Exception as exc:
        print(exc)
        print('error in reading data from MP')
        available = False
 
    
    if available:
        
        if not os.path.isdir(os.path.join(download_to,mpid)):
            os.mkdir(os.path.join(download_to,mpid))
        os.chdir(os.path.join(download_to,mpid))
        
        
        structure_doc = data[0]['xas'][0]['structure']
        structure_doc.to(fmt='poscar', filename='CONTCAR')
        equiv  = list(set(SpacegroupAnalyzer(structure_doc).get_symmetry_dataset()['equivalent_atoms']))
        
        # first write xas if available in xas_doc
        for xas_doc in data[0]['xas']:
            
            if structure_doc[xas_doc['absorbing_atom']].species_string == absorption_specie:
                
                f = 'feff_%03d_%s-K'%(xas_doc['absorbing_atom']+1,absorption_specie)
                if not os.path.isdir(f):
                    os.mkdir(f)
                    os.chdir(f)
                else:
                    os.chdir(f)                
                x, y = xas_doc['spectrum']
                out = np.column_stack(([x, y]))
                np.savetxt('xanes.dat', out)  
                
                os.chdir('..')
                
        # next determine missing ones
        missing = []
        for i in equiv:
            if structure_doc[i].species_string == absorption_specie:
                f = 'feff_%03d_%s-K'%(i+1,absorption_specie)
                if not os.path.isdir(f):
                    missing.append([mpid,i])

        os.chdir(here)
        if return_missing:
            return missing
        
    else:
        if return_missing: 
            
            structure_mp = mpr.get_structure_by_material_id(mpid,final=True)
            structure_mp_sym = SpacegroupAnalyzer(structure_mp).get_symmetrized_structure()
            equiv  = list(set(SpacegroupAnalyzer(structure_mp_sym).get_symmetry_dataset()['equivalent_atoms']))
            
            missing = [] 
            for i in equiv:
                if structure_mp[i].species_string == absorption_specie:
                    f = 'feff_%03d_%s-K'%(i+1,absorption_specie)
                    if not os.path.isfile(os.path.join(download_to,mpid,f,'xanes.dat')):
                        missing.append([mpid,i])
            return missing
    
    

    #structure = mpr.get_structure_by_material_id(i[0],final=True)
    #structure = SpacegroupAnalyzer(structure).get_symmetrized_structure()
##     structure.to(fmt='poscar',filename='CONTCAR') # this doesn't work

    #cf=open('CONTCAR',"w+") 
    #cf.write('symmetrized structure\n')
    #cf.write('1.0\n')
    #cf.write('%11.6f %11.6f %11.6f\n'%(structure_sym.lattice.matrix[0][0],structure_sym.lattice.matrix[0][1],structure_sym.lattice.matrix[0][2]))
    #cf.write('%11.6f %11.6f %11.6f\n'%(structure_sym.lattice.matrix[1][0],structure_sym.lattice.matrix[1][1],structure_sym.lattice.matrix[1][2]))
    #cf.write('%11.6f %11.6f %11.6f\n'%(structure_sym.lattice.matrix[2][0],structure_sym.lattice.matrix[2][1],structure_sym.lattice.matrix[2][2])) 
    
    #for i in structure_sym.types_of_specie:
        #cf.write(i.symbol+' ')
    #cf.write('\n')

    #for i in structure_sym.types_of_specie:
        #cf.write(str(structure.species.count(i))+' ')
    #cf.write('\n')  
    
    #cf.write('Direct\n')
    
    #for i in structure.sites:
        #cf.write('%8.6f %8.6f %8.6f %s \n'%(i.frac_coords[0],i.frac_coords[1],i.frac_coords[2],i.species_string))
   
    #cf.close()











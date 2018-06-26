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
from t4iss.defaults import t4iss_defaults
from t4iss.core import mXANES


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
    

def download_xanes_from_MP_old(mpr, mpid, absorption_specie, download_to=None):
    import numpy as np

    here = os.getcwd()
    corrupt_paths = []

    if download_to is None:
        download_to = t4iss_defaults['t4iss_xanes_data']

    data = mpr.get_data(mpid, data_type="feff", prop="xas")

    # check if data from MP is available
    if data[0]['xas']:
        species_in_data = []
        for xas_doc in data[0]['xas']:
            data_abs_specie = \
                xas_doc['structure'].species[xas_doc['absorbing_atom']].name
            species_in_data.append(data_abs_specie)
        if absorption_specie in species_in_data:
            available = True
        else:
            available = False
    else:
        available = False

    if available:
        data_structure0 = data[0]['xas'][0]["structure"]
        sites = []
        for i in data_structure0:
            sites.append(i.species_string)
        if absorption_specie not in sites:
            print('ERROR: '+absorption_specie+' is NOT found in this'
                  ' structure. Please check....')
        else:
            os.chdir(download_to)
            os.makedirs(mpid, exist_ok=True)
            os.chdir(mpid)
            spectra = []
            for xas_doc in data[0]['xas']:
                data_abs_specie = \
                    xas_doc['structure'].\
                    species[xas_doc['absorbing_atom']].name
                data_structure = xas_doc["structure"]
                data_absorption_atom = xas_doc['absorbing_atom']
                data_edge = xas_doc["edge"]
                x, y = xas_doc['spectrum']
                if data_abs_specie == absorption_specie:
                    finder = SpacegroupAnalyzer(data_structure)
                    structure = finder.get_symmetrized_structure()
                    [sites, indices] = structure.equivalent_sites, \
                        structure.equivalent_indices
                    for m in indices:
                        if m[0] == data_absorption_atom:
                            multiplicity = len(m)
                            ind = m[0]
                            
                    try:     
                        f = 'feff_{:03d}_{}-{}'.format(ind+1,
                                                       absorption_specie,
                                                       data_edge)
                    except:
                        print('XANES data for '+mpid+' '
                              'in Materials Project is corrupt.')
                        path = os.getcwd() + '/WARNING_CORRUPT'
                        f = open(path, 'w+')
                        f.write('data is corrupt')
                        f.close()
                        return []
                        
                    os.makedirs(f, exist_ok=True)
                    os.chdir(f)
                    xanes = mXANES(data=[x, y], structure=[data_structure,
                                   data_absorption_atom], xanesid=mpid,
                                   source='from_MP',
                                   edge=data_edge, multiplicity=multiplicity)
                    pickle.dump(xanes, open('xanes.pkl', 'wb'))
                    out = np.column_stack((xanes.E0, xanes.I0))
                    np.savetxt('xanes.dat', out)
                    spectra.append(xanes)
                    os.chdir('..')


            data_structure.to(fmt='poscar', filename='CONTCAR')
            os.chdir(here)
            print('XANES data was downloaded to '+download_to+'/'+mpid)
            return [data_structure, spectra]
    else:
        print('XANES data for '+mpid+' is not available in Materials Project'
              ' database...')
        path = t4iss_defaults['t4iss_xanes_data'] + '/' + str(mpid)
        try:
            os.makedirs(path)
        except:
            # file exists already, that's ok
            pass # mp-1071163
        f = open(path + '/WARNING_MISSING', 'w+')
        f.write('data is missing')
        f.close()
        return []













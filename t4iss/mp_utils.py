#!/usr/bin/env python3

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
import pickle
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from defaults import t4iss_defaults
from core import mXANES


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
def download_xanes_from_MP(mpr, mpid, absorption_specie, download_to=None):

    here = os.getcwd()

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
                        f = 'feff_{:03d}_{}-{}'.format(ind+1, absorption_specie,
                                                   data_edge)
                    except:
                        print('XANES data for '+mpid+' in Materials Project is corrupt.')
                        os.chdir('..')
                        
                        if mpid[0:2] == 'mp' or mpid[0:2] == 'mv':
                            import shutil
                            shutil.rmtree(mpid)
                        os.chdir(here)
                        return []
                        
                    os.makedirs(f, exist_ok=True)
                    os.chdir(f)
                    xanes = mXANES(data=[x, y], structure=[data_structure,
                                   data_absorption_atom], xanesid=mpid,
                                   source='from_MP',
                                   edge=data_edge, multiplicity=multiplicity)
                    pickle.dump(xanes, open('xanes.pkl', 'wb'))
                    import numpy as np
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
        os.chdir(here)
        return []
















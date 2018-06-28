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

    l_ength = len(mpid_list)
    print("""Found %i structures""" % l_ength)

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

    l_ength = len(mpid_list_icsd)
    print("""Found %i structures""" % l_ength)

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
                data_abs_specie = xas_doc['structure']\
                    .species[xas_doc['absorbing_atom']].name
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

        if not os.path.isdir(os.path.join(download_to, mpid)):
            os.mkdir(os.path.join(download_to, mpid))
        os.chdir(os.path.join(download_to, mpid))

        structure_doc = data[0]['xas'][0]['structure']
        structure_doc.to(fmt='poscar', filename='CONTCAR')
        equiv = list(set(SpacegroupAnalyzer(structure_doc)
                         .get_symmetry_dataset()['equivalent_atoms']))

        # first write xas if available in xas_doc
        for xas_doc in data[0]['xas']:

            if structure_doc[xas_doc['absorbing_atom']].species_string \
               == absorption_specie:

                f = 'feff_%03d_%s-K' % (xas_doc['absorbing_atom'] + 1,
                                        absorption_specie)
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
                f = 'feff_%03d_%s-K' % (i + 1, absorption_specie)
                if not os.path.isdir(f):
                    missing.append([mpid, i])

        os.chdir(here)
        if return_missing:
            return missing

    else:
        if return_missing:

            structure_mp = mpr.get_structure_by_material_id(mpid, final=True)
            structure_mp_sym = SpacegroupAnalyzer(structure_mp) \
                .get_symmetrized_structure()
            equiv = list(set(SpacegroupAnalyzer(structure_mp_sym)
                             .get_symmetry_dataset()['equivalent_atoms']))

            missing = []
            for i in equiv:
                if structure_mp[i].species_string == absorption_specie:
                    f = 'feff_%03d_%s-K' % (i + 1, absorption_specie)
                    if not os.path.isfile(os.path.join(
                                          download_to, mpid, f, 'xanes.dat')):
                        missing.append([mpid, i])
            return missing


def download_xanes_from_MP_old(mpr, mpid, absorption_specie, download_to=None):
    import numpy as np

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
            print('ERROR:' + absorption_specie + 'is NOT found in this'
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

                # x and y are the energy and intensity in two lists
                x, y = xas_doc['spectrum']

                if data_abs_specie == absorption_specie:
                    finder = SpacegroupAnalyzer(data_structure)
                    structure = finder.get_symmetrized_structure()
                    [sites, indices] = structure.equivalent_sites, \
                        structure.equivalent_indices

                    """
                    If one needs to pull structural data about the sites and
                    equilvalent indices from a CONTCAR file, this is possible.

                    This code produces a structure (structure2) that is
                    slightly different from the above data_structure/structure
                    (missing the last two columns) but appears to produce the
                    same output as when using data_structure from the MP.

                    file_name = '/path/to/CONTCAR'
                    import pymatgen as mg
                    structure2 = mg.Structure.from_file(file_name)
                    finder = SpacegroupAnalyzer(structure2)
                    structure = finder.get_symmetrized_structure()
                    [sites, indices] = structure.equivalent_sites, \
                        structure.equivalent_indices
                    """

                    for m in indices:
                        if m[0] == data_absorption_atom:
                            multiplicity = len(m)
                            ind = m[0]

                    try:
                        f = 'feff_{:03d}_{}-{}'.format(ind + 1,
                                                       absorption_specie,
                                                       data_edge)
                    except: # noqa
                        print('XANES data for' + mpid +
                              ' in Materials Project is corrupt.')
                        path = os.getcwd() + '/WARNING_CORRUPT'
                        f = open(path, 'w+')
                        f.write('data is corrupt')
                        f.close()
                        return []

                    os.makedirs(f, exist_ok=True)
                    os.chdir(f)
                    xanes = mXANES(data=[x, y],
                                   structure=[data_structure,
                                              data_absorption_atom],
                                   xanesid=mpid,
                                   source='from_MP',
                                   edge=data_edge, multiplicity=multiplicity)
                    pickle.dump(xanes, open('xanes.pkl', 'wb'))
                    out = np.column_stack((xanes.E0, xanes.I0))
                    np.savetxt('xanes.dat', out)
                    spectra.append(xanes)
                    os.chdir('..')

            data_structure.to(fmt='poscar', filename='CONTCAR')
            os.chdir(here)
            print('XANES data was downloaded to ' + download_to + '/' + mpid)
            return [data_structure, spectra]
    else:
        print('XANES data for' + mpid +
              ' is not available in Materials Project database...')
        path = t4iss_defaults['t4iss_xanes_data'] + '/' + str(mpid)
        try:
            os.makedirs(path)
        except: # noqa
            # file exists already, that's ok
            pass  # mp-1071163
        f = open(path + '/WARNING_MISSING', 'w+')
        f.write('data is missing')
        f.close()
        return []


def screen_data(species, transition, path=None):
    # TODO: docstring

    import pymatgen as mg

    total_cell_counter = 0
    bad_cell_counter = 0

    correct_string = str(species) + "-" + str(transition)
    forbidden = [".DS_Store", "CONTCAR", "Icon\r"]

    if path is None:
        path = t4iss_defaults['t4iss_xanes_data']

    all_cells_in_path = os.listdir(path)
    all_cells_in_path = [xx for xx in all_cells_in_path
                         if (xx not in forbidden and ".pkl" not in xx)]

    for i in range(len(all_cells_in_path)):
        cell = all_cells_in_path[i]
        cell_path = os.path.join(path, cell)
        all_spectra_in_cell = os.listdir(cell_path)

        # exclusion criteria
        if 'WARNING_MISSING' in all_spectra_in_cell \
           or 'WARNING_CORRUPT' in all_spectra_in_cell\
           or 'CONTCAR' not in all_spectra_in_cell:
            all_spectra_in_cell = []  # skip the entire thing
        else:
            contcar_name = os.path.join(cell_path, 'CONTCAR')
            structure = mg.Structure.from_file(contcar_name)
            finder = SpacegroupAnalyzer(structure)
            structure2 = finder.get_symmetrized_structure()
            indices = structure2.equivalent_indices

        all_spectra_in_cell = [xx for xx in all_spectra_in_cell
                               if xx not in forbidden]

        for j in range(len(all_spectra_in_cell)):
            spectra = all_spectra_in_cell[j]
            csp = os.path.join(cell_path, spectra)

            if correct_string not in csp:
                pass

            elif 'xanes.dat' not in os.listdir(csp) \
                 and 'xmu.dat' not in os.listdir(csp):
                # this means something went wrong, there should be an
                # xanes.dat file for every one of the cells,
                # something is missing and this entire cell could be corrupt
                print("Warning, xanes data missing from %s" % csp)
                f = open(os.path.join(cell_path, 'WARNING_MISSING'), 'w+')
                f.write('xanes data is missing from one or more entry')
                f.close()
                break

            elif 'xanes.pkl' not in os.listdir(csp):
                print("screening %s" % csp)
                for m in indices:
                    if str(m[0] + 1).zfill(3) in spectra:
                        multiplicity = len(m)
                        ind = m[0]
                        break

                try:
                    # this file better exist, otherwise problems
                    f = 'feff_{:03d}_{}-{}'.format(ind + 1,
                                                   species,
                                                   transition)
                except: # noqa
                    print(csp)
                    raise RuntimeError("Unknown error during execution")

                if 'xanes.dat' in os.listdir(csp):
                    xanes_data_path = os.path.join(csp, 'xanes.dat')
                    data = np.loadtxt(xanes_data_path)
                else:
                    xanes_data_path = os.path.join(csp, 'xmu.dat')
                    data = np.loadtxt(xanes_data_path, unpack=True,
                                      comments='#', skiprows=0)

                x = np.array(data[:, 0]).squeeze()
                y = np.array(data[:, 1]).squeeze()
                if len(x) < 3:
                    raise RuntimeError("You did something stupid.")
                xanes = mXANES(data=[x, y],
                               structure=[structure, ind],
                               xanesid=cell,
                               source='from_local_feff',
                               edge=transition, multiplicity=multiplicity)
                if xanes.vcn == []:
                    f = open(os.path.join(cell_path, 'WARNING_MISSING'), 'w+')
                    f.write('xanes data is missing from one or more entry')
                    f.close()
                else:
                    pickle.dump(xanes, open(os.path.join(csp, 'xanes.pkl'),
                                            'wb'))
            else:
                print("passing %s" % csp)

    for i in range(len(all_cells_in_path)):
        total_cell_counter += 1
        cell = all_cells_in_path[i]
        cell_path = os.path.join(path, cell)
        all_spectra_in_cell = os.listdir(cell_path)
        if 'WARNING_MISSING' in all_spectra_in_cell \
           or 'WARNING_CORRUPT' in all_spectra_in_cell:
            bad_cell_counter += 1

    print("good/total = %i/%i" % (total_cell_counter - bad_cell_counter,
                                  total_cell_counter))

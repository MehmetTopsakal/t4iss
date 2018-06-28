#!/usr/bin/env python3

"""All extractors."""

__author__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Development"
__date__ = "June 14, 2018"

# critical imports
import sys
import pickle
import os
from copy import deepcopy
from tqdm import tqdm

import logging

# local imports
from t4iss.cesym import get_cesym
from defaults import t4iss_defaults
from t4iss.xanes import read_xanes

# pymatgen imports
import pymatgen as mg
# from pymatgen import Structure
# from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.chemenv.coordination_environments\
    .coordination_geometry_finder import LocalGeometryFinder
from pymatgen.analysis.chemenv.coordination_environments\
    .chemenv_strategies import MultiWeightsChemenvStrategy
# from pymatgen.analysis.chemenv.coordination_environments\
#   .structure_environments import LightStructureEnvironments


strategy = MultiWeightsChemenvStrategy.stats_article_weights_parameters()


def extract_sites_o6_s5_t4(species, transition, path=None,
                           save_as='site_data.pkl', verbose=0):
    """From a directory coentaining structure/spectra data of various crystal
    structures from the Materials Project, extract and store all site data
    corresponding to the O:6, S:5 and T:4 chemical environments.

    ~ Note: this script takes a long time to run, be patient! ~

    :Input:
    - species (string) the atom of focus for the XANES spectrum. For example:
      'Ti'.
    - path (string) override the default path if desired. Should be the
      absolute path.
    - transition (string) for example: 'K'.
    - save_as (string) the output pkl file name.
    - verbose (int) qqual to 0 or 1. 0 for silent running, 1 for a progress
      indicator. 1 by default.
    """

    logger = logging.getLogger()
    logger.setLevel(logging.CRITICAL)

    # set the path to where the extracted data is located
    if path is None:
        path = t4iss_defaults['t4iss_xanes_data']
    # e.g.
    # path = os.path.join(os.environ['HOME'], "datasets", "FEFF_data_June7")

    cell_counter = 0
    total_counter = 0

    t4_data = []
    s5_data = []
    o6_data = []

    forbidden = [".DS_Store", "CONTCAR", "collection.pkl"]
    correct_string = str(species) + "-" + str(transition)

    # for cell (e.g. mp-390) in the desired data directory
    os_list_directory = os.listdir(path)
    N_os_list_directory = len(os_list_directory)

    for ii in tqdm(range(N_os_list_directory)):

        cell = os_list_directory[ii]

        if verbose == 1:
            if ii % 100 == 0:
                print("%i/%i samples analyzed" % (ii, N_os_list_directory))

        # reset the cell counter
        cell_counter = 0
        pass_cell = False

        # define a path if cell isn't = forbidden
        if cell not in forbidden:
            cell_path = os.path.join(path, str(cell))

            # get the structure
            try:
                structure = mg.Structure.from_file(cell_path + "/CONTCAR")
            except (FileNotFoundError, NotADirectoryError):
                # skip if CONTCAR does not exist in the directory
                # or if it isn't a directory at all, there are a few of those
                pass_cell = True

            # skip if directory is empty
            if not pass_cell:  # call this first to avoid NotADirectoryError
                if os.listdir(cell_path) == []:
                    pass_cell = True

            if not pass_cell:
                # define a path to a given spectra/pickle file
                for sample in os.listdir(cell_path):
                    if sample not in forbidden:
                        pass_sample = False

                        # cell sample path
                        csp = os.path.join(cell_path, str(sample))

                        # assert that we're reading the correct species
                        # and correct transition
                        if correct_string in csp:
                            try:
                                with open(csp +
                                          "/xanes.pkl", 'rb') as pickle_file:
                                    content = pickle.load(pickle_file)

                                # ensure we end up with content to analyze
                                # if not pass it
                                try:
                                    save_stdout = sys.stdout
                                    sys.stdout = open('.trash', 'w')
                                    lgf = LocalGeometryFinder()
                                    cesym = get_cesym(lgf,
                                                      structure,
                                                      content.structure[1])
                                    sys.stdout = save_stdout
                                except IndexError:
                                    pass_sample = True

                                if pass_sample:
                                    pass
                                elif cesym[0] == "T:4":
                                    t4_data.append([cell, cell_counter,
                                                    deepcopy(content), cesym])
                                    cell_counter += 1
                                    total_counter += 1
                                elif cesym[0] == "S:5":
                                    s5_data.append([cell, cell_counter,
                                                    deepcopy(content), cesym])
                                    cell_counter += 1
                                    total_counter += 1
                                elif cesym[0] == "O:6":
                                    o6_data.append([cell, cell_counter,
                                                    deepcopy(content), cesym])
                                    cell_counter += 1
                                    total_counter += 1
                            except FileNotFoundError:
                                pass

    all_data = [t4_data, s5_data, o6_data]
    path_save = os.path.join(path, save_as)

    with open(path_save, 'wb') as f:
        pickle.dump(all_data, f)


def extract_average_spectra(species, transition, path=None,
                            save_as='avg_data.pkl', verbose=0):
    """Documentation is identical to extract_sites_o6_s5_t4 except that this
    function extracts the average spectra for all site-types for every
    crystal structure.
    """

    logger = logging.getLogger()
    logger.setLevel(logging.CRITICAL)

    # set the path to where the extracted data is located
    if path is None:
        path = t4iss_defaults['t4iss_xanes_data']

    if verbose == 1:
        print("path is ", path)

    cell_counter = 0
    total_counter = 0
    key_site_counter = 0
    directory_counter = 0

    # initialize an empty dictionary which can be thought of as one hot
    # encoding the key will be of form O:6, T:4, etc. and the value will be
    # the index of the basis vector corresponding to that site - the labels
    # will be generated after this block is run
    key_site = {}

    # empty lists for all the cn data
    all_site_data = []
    all_cell_data = []

    # these files will end up raising errors during execution of the
    # following loops - avoid them
    forbidden = [".DS_Store", "CONTCAR"]

    # assert that only strings with the correct species and transition are
    # added to the training dataset
    correct_string = str(species) + "-" + str(transition)

    # for cell (e.g. mp-390) in the desired data directory
    os_list_directory = os.listdir(path)
    os_list_directory = [xx for xx in os_list_directory
                         if (xx not in forbidden and ".pkl" not in xx)]

    N_os_list_directory = len(os_list_directory)

    for ii in tqdm(range(3)):
        cell = os_list_directory[ii]

        if verbose == 1:
            if ii % 20 == 0.0:
                print("%i/%i" % (ii, N_os_list_directory))

        # reset the cell counter
        cell_counter = 0
        pass_cell = False
        current_cell_data = []

        # define a path if cell isn't = forbidden
        if cell not in forbidden:
            cell_path = os.path.join(path, str(cell))

            try:
                xanes_data = read_xanes(cell_path, absorption_specie=species)
            except (ValueError, FileNotFoundError):
                # print("No %s in %s. Passing." % (species, cell_path))
                pass_cell = True

            # get the structure
            try:
                structure = mg.Structure.from_file(cell_path + "/CONTCAR")
            except (FileNotFoundError, NotADirectoryError):
                # skip if CONTCAR does not exist in the directory
                # or if it isn't a directory at all, there are a few of those
                pass_cell = True

            # skip if directory is empty
            if not pass_cell:  # call this first to avoid NotADirectoryError
                if os.listdir(cell_path) == []:
                    pass_cell = True

            if not pass_cell:
                # define a path to a given spectra/pickle file
                for sample in os.listdir(cell_path):
                    if sample not in forbidden:

                        pass_sample = False

                        # cell sample path
                        csp = os.path.join(cell_path, str(sample))

                        # assert that we're reading the correct species
                        # and correct transition
                        if correct_string in csp:
                            try:
                                with open(csp +
                                          "/xanes.pkl", 'rb') as pickle_file:
                                    content = pickle.load(pickle_file)

                                # ensure we end up with content to analyze,
                                # if not pass it
                                try:
                                    save_stdout = sys.stdout
                                    sys.stdout = open('.trash', 'w')
                                    lgf = LocalGeometryFinder()
                                    cesym = get_cesym(lgf, structure,
                                                      content.structure[1])
                                    sys.stdout = save_stdout
                                except IndexError:
                                    pass_sample = True

                                if pass_sample:
                                    pass
                                else:
                                    # add in the average data if this is the
                                    # first time the
                                    # directory / cell has been looked at
                                    if cell_counter == 0:
                                        all_cell_data.\
                                            append([cell, directory_counter,
                                                    deepcopy(xanes_data[0])])

                                    # if the symmetry label does not exist,
                                    # append it to the dictionary with its
                                    # type and index, and initialize a new
                                    # list in all_cn_data with key_counter as
                                    # its index
                                    if cesym[0] not in key_site.keys():
                                        key_site[cesym[0]] = key_site_counter
                                        all_site_data.append([])
                                        key_site_counter += 1

                                    x = [cell, directory_counter, cell_counter,
                                         deepcopy(content),
                                         cesym, content.multiplicity]
                                    all_site_data[key_site[cesym[0]]].append(x)
                                    current_cell_data.append(x)
                                    cell_counter += 1
                                    total_counter += 1

                            except FileNotFoundError:
                                pass

        if pass_cell or (pass_sample and cell_counter == 0):
            pass
        else:
            # sanity check - can cross-reference with the atom id's using
            # VESTA
            # finder = SpacegroupAnalyzer(structure)
            # struct = finder.get_symmetrized_structure()
            # [sites, indices]  = struct.equivalent_sites,
            # struct.equivalent_indices
            # print(directory_counter, indices)

            # for every cell / directory (e.g. mvc-16746), generate an
            # individual dictionary which labels the proportion of each
            # coordination number present in the crystal structure
            temp_dictionary = {}

            # running counter of all the multiplicity numbers
            total_multiplicity = 0

            for xx in current_cell_data:

                # current multiplicity of xx
                cm = xx[3].multiplicity

                # if the current ce symbol does not exist in the
                # temp_dictionary, create it and initialize its counter equal
                # to the multiplicity of that entry
                if xx[4][0] not in temp_dictionary.keys():
                    temp_dictionary[xx[4][0]] = cm

                # else, it must have been previously initialized, and so we
                # should add the multiplicity to it
                else:
                    temp_dictionary[xx[4][0]] += cm

                # append the total multiplicity for normalizing later
                total_multiplicity += cm

            # next, modify each value in the dictionary: normalize by total
            # multiplicity
            for xx in temp_dictionary:
                temp_dictionary[xx] /= float(total_multiplicity)

            all_cell_data[directory_counter].append(temp_dictionary)
            directory_counter += 1

    path_save = os.path.join(path, save_as)

    with open(path_save, 'wb') as f:
        pickle.dump(all_cell_data, f)

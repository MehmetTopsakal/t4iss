#!/usr/bin/env python3

"""See main()."""

# pylint: disable=C0103
# pylint disable=E0401

import os
import pickle

from pymatgen import MPRester

from yaml import safe_load

from mp_utils import download_xanes_from_MP as dMP
from mp_utils import search_MP

from feff_utils import xanes_collection_feff
from t4iss.defaults import t4iss_defaults
from t4iss.core import mXANES




def main():
    """Reads in all data from the MP corresponding to the search pattern,
    such as Ti-O, absorption species, such as Ti, and maximum number of atoms
    nmax_atoms. Data is deposited in ~/.t4iss/xanes_data."""


    protocol = safe_load(open("protocol.yml"))
    mpr = MPRester(protocol['MP_code'])  # Matt's code
    mpids = search_MP(mpr, search_pattern=protocol['search_pattern'],
                      nmax=protocol['nmax_atoms'])
    ii = 0

    cutoff = len(mpids)
    if protocol['nmax_total'] > 0:
        cutoff = min(cutoff, protocol['nmax_total'])

    """
    while ii < cutoff:
        dMP(mpr, mpids[ii], absorption_specie=protocol['absorption_species'])
        print('%i/%i' % (ii+1, cutoff))
        ii += 1
    


    # generate the collection
    collection = xanes_collection_feff(absorption_specie='Ti')
    collection.build(mpids, skip_missing=True, download_missing=True, mpr=mpr)

    """

    my_collection = []
    dirs = os.listdir(t4iss_defaults['t4iss_xanes_data'])
    for i in dirs:
        my_collection.append(mXANES(data=(i + 'xanes.dat'), 
                                    structure=(i + 'CONTCAR') ))
    exit(0)

    data_name = "collection.pkl"
    path_save = os.path.join(t4iss_defaults['t4iss_xanes_data'], data_name)
    with open(path_save, 'wb') as f:
        pickle.dump(collection, f)

if __name__ == '__main__':
    main()

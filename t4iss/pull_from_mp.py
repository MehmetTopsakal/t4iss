#!/usr/bin/env python3

from pymatgen import MPRester

# pylint disable=E0401
from mp_utils import download_xanes_from_MP, search_MP
from feff_utils import xanes_collection_feff
from defaults import t4iss_defaults



def main():
    mpr = MPRester('HnEctIuK6wbkQGse')  # Matt's code
    mpids = search_MP(mpr, search_pattern='Ti-O', nmax=201)
    download_xanes_from_MP(mpr, mpids, absorption_specie='Ti')


if __name__ == '__main__':
    main()
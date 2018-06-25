#!/usr/bin/env python3

import os.path
import shutil

def print_defaults():
    for key, val in t4iss_defaults.items():
        print("- {} : {}".format(key, val))

here = os.path.dirname(os.path.realpath(__file__))
user_home = os.path.expanduser('~')

# defaults
t4iss_defaults = dict()
t4iss_defaults['t4iss_data'] = os.path.join(user_home, '.t4iss', 'data')
t4iss_defaults['t4iss_scratch'] = os.path.join(user_home, '.t4iss', 'scratch')
t4iss_defaults['mcr_path'] = os.path.join(here, 'mcr')
t4iss_defaults['scripts_path'] = os.path.join(here, 'scripts')
t4iss_defaults['octave_path'] = shutil.which('octave')
t4iss_defaults['matlab_path'] = shutil.which('matlab')

# dictionary to translate the chemical environments
ce_translator = {
    'S:1':'single neighbor',
    'L:2':'linear',
    'A:2':'angular',
    'TL:3':'trigonal plane',
    'TY:3':'trigonal non-coplanar',
    'T:4':'tetrahedral',
    'S:4':'square plane',
    'SS:4':'see-saw',
    'S:5':'square pyramidal',
    'T:5':'trigonal bipyramid',
    'O:6':'octahedral',
    'T:6':'trigonal prism',
    'PB:7':'pentagonal bipyramid',
    'C:8':'cube',
    'SA:8':'square antiprism',
    'DDPN:8':'dodecahedron with triangular faces',
    'HB:8':'hexoganal bipyramid',
    'C:12':'cuboctahedral'
    }

#!/usr/bin/env python3

__author__ = "Matthew Carbone"
__email__ = "x94carbone@gmail.com"
__status__ = "Development"
__date__ = "June 28, 2018"

from yaml import safe_load

from t4iss.extractors import extract_sites_o6_s5_t4, extract_average_spectra
from t4iss.mp_utils import screen_data


def main():
    """Extracts single site data from local."""

    protocol = safe_load(open("protocol.yml"))

    if protocol['extract_type'] == 'site':
        extract_sites_o6_s5_t4(protocol['absorption_species'],
                               protocol['transition'])

    elif protocol['extract_type'] == 'avg':
        extract_average_spectra(protocol['absorption_species'],
                                protocol['transition'])

    elif protocol['extract_type'] == 'screen_only':
        # screen local data to ensure that all files have a pickle file
        # containing critical information
        # default path is t4iss_defaults['t4iss_xanes_data']
        # = ~/.t4iss/xanes_data
        screen_data(protocol['absorption_species'],
                    protocol['transition'])


if __name__ == '__main__':
    main()

#!/usr/bin/env python3

"""See main()."""

# pylint: disable=C0103
# pylint disable=E0401


from yaml import safe_load

from t4iss.extractors import extract_sites_o6_s5_t4

def main():
    """Extracts single site data from local."""

    protocol = safe_load(open("protocol.yml"))

    if protocol['extract_type'] == 'site':
        extract_sites_o6_s5_t4(protocol['absorption_species'],
                               protocol['transition'])


if __name__ == '__main__':
    main()

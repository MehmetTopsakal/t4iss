from __future__ import (absolute_import, division, print_function)

import setuptools

no_git_reqs = []
with open('requirements.txt') as f:
    required = f.read().splitlines()
    for r in required:
        if not (r.startswith('git') or r.startswith('#') or r.strip() == ''):
            no_git_reqs.append(r)

setuptools.setup(
    name='t4iss',
    license="BSD (3-clause)",
    url="https://github.com/MehmetTopsakal/t4iss",
    packages=setuptools.find_packages(),
    package_data={'t4iss' : ['../examples/*.ipynb', 'mcr/matlab_version/*',
                  'mcr/octave_version/*', 'scripts/*']},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
        ],
    install_requires=no_git_reqs,
    )

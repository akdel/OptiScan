#!/usr/bin/env python

from distutils.core import setup

setup(name='OptiScan',
      version='1.0',
      description='Bionano Tiff to Optical Map Signals',
      author='Mehmet Akdel',
      author_email='mehmet.akdel@wur.nl',
      url='https://gitlab.com/akdel/',
      scripts=["pipelines/extract_molecules",
               "pipelines/extract_molecules_saphyr",
               "pipelines/write_bnx",
               "pipelines/write_bnx_saphyr",
               "pipelines/write_signals.py",
               "pipelines/write_signals_saphyr.py"],
      packages=['OptiScan'],
      install_requires=["sqlalchemy", "scipy", "numba",
                        "intervaltree", "matplotlib", "scikit-image", "pillow",
                        "imageio", "dash", 'imagecodecs'])


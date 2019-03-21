#!/usr/bin/env python

from distutils.core import setup

setup(name='OptiScan',
      version='1.0',
      description='Bionano Tiff to Optical Map Signals',
      author='Mehmet Akdel',
      author_email='mehmet.akdel@wur.nl',
      url='https://gitlab.com/akdel/',
      packages=['OptiScan'],
      install_requires=["sqlalchemy", "numba", "numpy", "scipy", "intervaltree", "matplotlib", "scikit-image", "pillow", "imageio", "dash"])
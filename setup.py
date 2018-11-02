#!/usr/bin/env python

from distutils.core import setup

setup(name='Photomap',
      version='1.0',
      description='Bionano Tiff to Optical Map Signals',
      author='Mehmet Akdel',
      author_email='mehmet.akdel@wur.n',
      url='https://gitlab.com/akdel/',
      packages=['OptiScan'],
      install_requires=["sqlalchemy", "numpy", "scipy", "intervaltree", "matplotlib"])
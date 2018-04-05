
#! /usr/bin/env python

from setuptools import setup, find_packages

import glob

##-------------------------------------------------
## Package setup

setup(name='pyxsis',
      version='0.0',
      description='X-ray Python Spectral Interpretation System',
      author='Lia Corrales',
      author_email='lia@astro.wisc.edu',
      url='https://github.com/eblur/pyxsis',
      packages=find_packages(),
      package_data={'pyxsis': ['models/tables/*']},
      license='LGPL-3.0',
      install_requires=[
        'numpy>=1.10',
        'astropy>=2.0.0',
        'clarsach>=0.0.0']
)

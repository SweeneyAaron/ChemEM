#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 14:43:33 2023

@author: aaron.sweeney
"""
from setuptools import setup, find_packages

setup(
    name='ChemEM',
    version='0.0.1',
    packages=find_packages(),
    include_package_data=True,
    install_requires=['BioTEMPy==2.0.0'],
    entry_points={
        'console_scripts': [
            'chemem = ChemEM.chemem:main',
            'chemem.test = ChemEM.chemem:test'
        ],
    },
    url='https://chemem.topf-group.com/',
    license='Your License',
    description='Software for fitting small molecules into Cryo-EM data.'
)



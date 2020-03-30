#!/usr/bin/env python

from distutils.core import setup

setup(name='MAP4',
            version='1.0',
            description='MinHashed AtomPair Fingerprint of Radius 2',
            author='Alice Capecchi',
            author_email='alice.capecchi@outlook.it',
            url='https://github.com/reymond-group/map4',
            packages=['map4'],
            install_requires=['faerun', 'mhfp']
           )

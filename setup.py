#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

from distutils.core import Extension

if sys.argv[-1] == 'publish':
    os.system('python setup.py sdist upload')
    sys.exit()

readme = open('README.rst').read()
history = open('HISTORY.rst').read().replace('.. :changelog:', '')

module1 = Extension('analytical', sources = ['pas/c_source/analytical_c.c'])

setup(
    name='pas',
    version='0.1.0',
    description='Analytical and semi-analytical solutions for petroleum related equations.',
    long_description=readme + '\n\n' + history,
    author='Omar Al-Hinai',
    author_email='ohinai@gmail.com',
    url='https://github.com/ohinai/pas',
    packages=[
        'pas',
    ],
    package_dir={'pas':
                 'pas'},
    include_package_data=True,
    install_requires=[
    ],
    license="BSD",
    zip_safe=False,
    keywords='pas',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        "Programming Language :: Python :: 2",
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
    ],
    test_suite='tests',
    ext_modules = [module1], 
)

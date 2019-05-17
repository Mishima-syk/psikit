#!/usr/bin/env python

#from distutils.core import setup
from setuptools import setup

with open('README.md') as readme_file:
    long_description = readme_file.read()

setup(name='Psikit',
      version='0.1.4',
      description='A thin wrapper library for Psi4 and RDKit',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Kazufumi Ohkawa, Takayuki Serizawa',
      author_email='kerolinq@gmail.com, seritaka@gmail.com',
      url='https://github.com/Mishima-syk/psikit',
      packages=['psikit'],
      install_requires=[],
      license='MIT',
     )

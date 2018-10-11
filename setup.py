#!/usr/bin/env python

from distutils.core import setup

with open('README.md') as readme_file:
    long_description = readme_file.read()

setup(name='Psikit',
      version='0.1.0',
      description='A thin wrapper library for Psi4 and RDKit',
      long_description=long_description,
      author='Kazufumi Ohkawa, Takayuki Serizawa',
      author_email='kerolinq@gmail.com, ######',
      url='https://github.com/Mishima-syk/psikit',
      packages=['psikit'],
      install_requires=['rdkit','psi4',],
      license='MIT',
     )
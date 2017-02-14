#!/usr/bin/env python

# from distutils.core import setup
from setuptools import setup

setup(name="clusterfast",
      version='0.0.1',
      description="Python script for clustering sequences from closely related\
                    organisms",
      author="Anmol M. Kiran\nJennifer E. Cornick",
      author_email="anmol@liv.ac.uk\nj.cornick@liv.ac.uk",
      url="https://github.com/codemeleon/ClusterFast",
      install_requires=['click',
      					'networkx',
      					'pandas',
      					'numpy',
      					'numba',
      					'biopython'],
      # install_requires=
      license='GPLv3',
      # scripts=['bin/fastclust','bin/pblat']
      scripts=['bin/clusterfast'],
      # packages=["fastclustlib"]
      )

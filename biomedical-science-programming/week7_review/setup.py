import setuptools
from setuptools import setup, Extension
try:
    import numpy as np
except ImportError:
    print("Install numpy first. 'pip3 install numpy --user'")

module = Extension('npmodule.utils', sources = ['c/utils.cpp'], include_dirs=[np.get_include()])

setup(name = 'npmodule',
      version = '1.0',
      description = 'This is a numpy demo package',
      packages = setuptools.find_packages(),
      ext_modules = [module])

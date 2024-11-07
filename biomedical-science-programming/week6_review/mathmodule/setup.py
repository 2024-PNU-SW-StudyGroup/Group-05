from setuptools import setup, Extension

module = Extension('mathmodule', sources = ['mathmodule.cpp'])

setup(name = 'mathmodule',
      version = '1.0',
      description = 'This is a demo package',
      ext_modules = [module])

# python setup.py build_ext --inplace
from distutils.core import setup, Extension

module1 = Extension('gf2c', sources = ['gf2.c'])

setup (name = 'gf2c',
       version = '1.0',
        description = 'Galois field package',
        ext_modules = [module1])

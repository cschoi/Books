from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [Extension("somapCython", ["somapCython.pyx"],libraries=["m"])]

setup(
  name = 'Cython somap app',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)

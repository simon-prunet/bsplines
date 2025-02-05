import numpy as nm
import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# alldir=["/softs/gsl/1.15/lib"]
# allinc=["/softs/gsl/1.15/include"] + [nm.get_include()]
# alllib=["gsl","gslcblas","m"]
#alldir=[" "]
allinc=[nm.get_include()]
alllib=["m"]

setup(
  cmdclass = {'build_ext': build_ext},
  ext_modules = [Extension("bsplines", ["bsplines.pyx", "coeff.c", "interpol.c"],
                            #library_dirs=alldir,
                            include_dirs=allinc,
                            libraries=alllib,
                            extra_compile_args=["-fopenmp"],
                            extra_link_args=["-fopenmp"]
                          )
                ]
)


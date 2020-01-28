from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("fraction_id_calc.pyx")
)

## NOTE: THIS IS COMPILED VIA
##        python setup_fic.py build_ext --inplace

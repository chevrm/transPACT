from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("calculate_distance_transATPKS_cy.pyx")
)

## NOTE: THIS IS COMPILED VIA
##        python setup_fic.py build_ext --inplace

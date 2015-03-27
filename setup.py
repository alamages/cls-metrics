from distutils.core import setup
from Cython.Build import cythonize

setup(
  name = 'clustering evaluation metrics using cython/numpy',
  ext_modules = cythonize("clustering_metrics.pyx"),
)

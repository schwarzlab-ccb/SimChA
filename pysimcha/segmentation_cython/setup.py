from setuptools import setup
from Cython.Build import cythonize
import numpy as np

setup(
    name='consistent_segmentation_c',
    ext_modules=cythonize("consistent_segmentation_step.pyx"),
    include_dirs=[np.get_include()],
    zip_safe=False,
)

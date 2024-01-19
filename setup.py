from setuptools import setup, find_packages

import os
from sys import platform

from setuptools.command.install import install
from setuptools.command.build_py import build_py
from os import listdir, path


PACKAGES = [
    'DynaBench',
]

class install_library(install):
    def run(self):
        """Run the installation."""
        install.run(self)


class build_py_modules(build_py):
    def run(self):
        """Run the build."""
        build_py.run(self)

        os.chdir('interfacea')
        os.system('python setup.py build')
        os.system('python setup.py install') 
        os.chdir('../')


here = path.abspath(path.dirname(__file__))

# Collect names of bin/*py scripts
# e.g. 'pdb_intersect=bin.pdb_intersect:main',
binfiles = listdir(path.join(here, 'DynaBench'))
bin_py = [f[:-3] + '=DynaBench.' + f[:-3] + ':main' for f in binfiles
          if f.endswith('dynabench.py')]

setup(
    name='DynaBench',
    version='0.1.0',
    cmdclass={
          "install": install_library,
          "build_py": build_py_modules,
          # "test": run_tests,
      },
    packages=PACKAGES,
    entry_points={  # Optional
        'console_scripts': bin_py,
    },
)







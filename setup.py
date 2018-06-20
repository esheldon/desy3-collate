import os
import distutils
from distutils.core import setup, Extension, Command
import glob

scripts = glob.glob('./bin/*')
scripts = [os.path.basename(f) for f in scripts if f[-1] != '~']
scripts = [os.path.join('bin',s) for s in scripts]

setup(
    name="desy3collate", 
    packages=['desy3collate'],
    scripts=scripts,
    version="1.0.0",
)





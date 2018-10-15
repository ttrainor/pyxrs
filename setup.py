#!/usr/bin/env python3
"""
Notes:
------
Setup file to install pyxrs into site-packages.

* To build and install run as:
    >>python3 setup.py install

  When installed the resulting package should look like:
    python-modules in:  your-site-packages/pyxrs
    scripts in:         your-bin/
    e.g. on linux standard locations might be:
    your-site-packages = /usr/local/lib/python3.5/dist-packages/
    your-bin           = /usr/local/bin

  To determine where your system installs (and searches for) packages try:
    >>python3 -m site
    >>python -c "import site; print(site.getsitepackages())"

* To keep track of where the files are installed (e.g. so you can remove them)
    >>python3 setup.py install  --record=files.txt

* To seperate build and install:
    >>python3 setup.py build
    >>python3 setup.py install

  After running build you'll have a directory with the following layout:
    ./build/lib/pyxrs/...

  To remove the 'build' directory created by setup
    >>python setup.py clean --all

  You can customize where the build directory goes using:
    >>python3 setup.py build --build-base=path-to-build

* To install the package under python's USER_BASE (ie local to the user) run
    >>python3 setup.py install  --user

  This will result in the following layout:
    python modules in:  USER_BASE/lib/pythonx.y/site-packages/   
    scripts in:         USER_BASE/bin
    data in:            USER_BASE

  Run ">>python3 -m site" to see where USER_BASE points

* To install the package in a custom location you can run
    >>python3 setup.py install  --home=<dir>

  This will result in the following layout:
    python modules in:  home/lib/python/   
    scripts in:         home/bin
    data in:            home

* To build a source distribution:
    >>python3 setup.py sdist --dist-dir=<dir>

  The optional --dist-dir argument specifies the location for the distribution file
  There is also an optional --formats argument that allows specifying the output format
  A similiar method is used to make built distributions using the bdist argument

* see: https://docs.python.org/3/distutils/

"""
#########################################################
import sys, os
import distutils
from distutils.core import setup, Extension

#########################################################
### make sure have python3.5
if sys.version_info.major < 3:
    print("Requires python 3 or greater")
    sys.exit()
if sys.version_info.minor < 5:
    print("Requires python 3.5 or greater")
    sys.exit()

### make sure have necessary modules
required_modules = ['numpy'] #,'matplotlib']
for mod in required_modules:
    try:
        __import__(mod)
    except:
        print("Requires %s" % mod)
        sys.exit()

### package dir - mapping directory names used by packages 
# and package_data 
package_dir  = {'pyxrs':'',
#                'pyxrs.doc':'doc',
#                'pyxrs.examples':'examples',
                'pyxrs.geom':'geom',
#                'pyxrs.scripts':'scripts',
                'pyxrs.xtal':'xtal'}

### packages
packages     = ['pyxrs',
#                'pyxrs.examples',
                'pyxrs.geom',
                'pyxrs.xtal']

### package data
package_data = {}
package_data.update({'pyxrs':['Readme.txt']})
#packages.append('pyxrs.doc')
#package_data.update({'pyxrs.doc':['Install.txt', 'License.txt']})

### include examples...
# note that the individual example dirs arent 'packages'
# so disutils complains, but this seems to get the files in the
# right place (same with 'doc' above)
#def add_example(ex):
#    name = 'pyxrs.examples.' + ex 
#    pth  = os.path.join('examples',ex)
#    package_dir.update({name:pth})
#    packages.append(name)
#    package_data.update({name:['*.py','*.cif','*.spt']})
#add_example('cristobalite')
#add_example('forsterite')
#add_example('forsterite_111')
#add_example('hematite')
#add_example('hematite_012')
#add_example('pyrite')

### scripts
#import glob
#xx = os.path.join('scripts','*')
#scripts = glob.glob(xx)

### call the setup command
name    = 'pyxrs'
version = '0.1'
author  = "Trainor"
email   = "tptrainor@alaska.edu"
descr   = "python xtalography and diffraction calculations"
setup( name = name,
       version = version,
       author =  author,
       author_email = email,
       url = "",
       description  = descr,
       package_dir  = package_dir,
       packages     = packages,
       package_data = package_data,
#       data_files   = datafiles,
#       scripts      = scripts
#       python_requires='>=3.5',
#       install_requires=required_modules,
)



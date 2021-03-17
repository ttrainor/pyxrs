
pyxrs
=====

This is a collection of python modules (python3) for crystallography and 
xray scattering computations.  If you run the install script as described 
below two separate packages are installed ('xtal' and 'xtlgeo').  
You can import them via 

.. code::

    >>import xtal
    >>import xtlgeo 


A directory of example scripts is also provided, but not installed. 


Quick install:
--------------

1.  Unpack the distribution to a temporary folder
2.  Run the setup.py file 

.. code::

    >>python3 setup.py install

For more information on installation see doc/Install.txt


Brief description of modules:
-----------------------------

* xtal

  - Crystallography computations including reading and writing structure
    files, transformation of coordinates, generation of surfaces, coordination 
    sphere computations etc..
  - This module also includes modules for reading and writing cif files 
  - The only dependency for this module is numpy 

* xtlgeo

  - Geometry/Goniometer calculations for crystallography including computation of 
    active area.  Currently defines methods for a six-circle (and kappa) geometry.
  - This uses xtal as well as numpy and Matplotlib 

Docs:
-----

* See docs/Install.txt for more installation information
* See docs/License.txt for information on copyright
* See examples for several example application scripts 

Resources:
----------

* Sources of structure data:
   - http://rruff.geo.arizona.edu/AMS/amcsd.php
   - http://www.crystallography.net/cod/
   - http://www.chemspider.com/

* Jmol for structure viewing
   - http://wiki.jmol.org/index.php/Main_Page
   - https://chemapps.stolaf.edu/jmol/docs/

* Webatoms: http://millenia.cars.aps.anl.gov/webatoms/

* Cif file formats:  
   - https://www.iucr.org/resources/cif/spec/version1.1/cifsyntax
   - https://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/index.html
   - http://xray.tamu.edu/pdf/manuals/cifguide.pdf
   - https://pypi.python.org/pypi/PyCifRW/4.1

* OpenBabel:
   - http://openbabel.org/docs/dev/




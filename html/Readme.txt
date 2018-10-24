
The files here can be used to setup a simple web server to do xtal 
computations and display structures using jmol.  

Setup
=====
* Pre-requisites 
  - Apache server configured to run cgi scripts
  - Python3 with numpy and xtal modules installed

* The following file structure is assumed:
  
  server_root/cgi-bin
             /jmol
             /xtal

* To display structure you need Jsmol installed:
  - Download from: https://sourceforge.net/projects/jmol/files/
  - Open the archive and copy the files from the jsmol directory
    so you end up with the following
          server_root/jmol/j2s
                          /java
                          /php
                          /JSmol.min.js
  - Set permisions / security as appropriate...

* If your file structure is different than above you'll need to edit 
  the paths in:
    xtal/xtal.html   --> path to xtal.cgi
    cgi-bin/xtal.cgi --> path to xtal_dir
    cgi-bin/_xtal.py --> path to JSmol.min.js
                     --> j2sPath
                     --> jarPath
                     --> serverURL      


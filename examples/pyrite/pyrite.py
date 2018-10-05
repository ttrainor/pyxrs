"""
Test some xtal calcs using pyrite
"""
##############################
from pyxrs.xtal.unitcell import cif_to_uc, uc_to_xyz
##############################

## Read pyrite cif file
## and display unit cell info
uc = cif_to_uc('pyrite.cif')
uc.write()

## list P1 coordinates
print("P1 cell")
at_list = uc.atom_list()
at_list.write()

## write an xyz file
uc_to_xyz(uc,fname="pyrite.xyz",na=2,nb=2,nc=2)

#######
# coordination calculations
from pyxrs.xtal.coord import coord_calcs
coord = coord_calcs(uc,rmax=2.5)
coord.write()
coord.write(fname="coord.out")

# run jmol
import subprocess
subprocess.run("jmol jmol_script.spt", shell=True)


"""
Test some xtal calcs using pyrite
"""
##############################
from xtal.unitcell import read_cif
##############################

## Read pyrite cif file
## and display unit cell info
uc = read_cif('pyrite.cif')
uc.write()

## list P1 coordinates
print("P1 cell")
at_list = uc.atom_list()
at_list.write_xyz()

## write an xyz file
uc.write_xyz(fname="pyrite.xyz",na=2,nb=2,nc=2)

#######
# coordination calculations
from xtal.coord import coord_calcs
coord = coord_calcs(uc,rmax=2.5)
coord.write()
coord.write(fname="coord.out")

# run jmol
import subprocess
subprocess.run("jmol jmol_script.spt", shell=True)


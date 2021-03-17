"""
Test some xtal calcs using forsterite
"""
##############################
from xtal.unitcell import read_cif, write_cif
##############################

## generate a UnitCell instance 
## display cell contents
print("** Forsterite unit cell")
uc = read_cif('forsterite_1.cif')
uc.show()

## print the P1 cell 
print("P1 cell")
atom_list = uc.atom_list()
atom_list.show()

## output xyz file
uc.write_xyz("forsterite.xyz",cartesian=True,na=3,nb=3,nc=3,long_fmt=False)

## output cif file with assymetric unit (similiar to what we read in)
#write_cif(uc,fname="forsterite_frac.cif",p1_list=False)

## output cif file with P1 unit cell contents (no symmetry)
#write_cif(uc,fname="forsterite_p1.cif",p1_list=True,na=2,nb=2,nc=2)

#######
# coordination calculations
from xtal.coord import coord_calcs
coord = coord_calcs(uc,rmax=2.5)
coord.show(long_fmt=False)
coord.write(fname="coord.out")

# run jmol
import subprocess
subprocess.run("jmol jmol_script.spt", shell=True)


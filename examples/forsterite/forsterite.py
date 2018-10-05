"""
Test some xtal calcs using forsterite
"""
##############################
from pyxrs.xtal.unitcell import cif_to_uc, uc_to_cif, uc_to_xyz 
##############################

## generate a UnitCell instance 
## display cell contents
print("** Forsterite unit cell")
uc = cif_to_uc('forsterite_1.cif')
uc.write()

## print the P1 cell 
print("P1 cell")
atom_list = uc.atom_list()
atom_list.write()

## output xyz file
uc_to_xyz(uc,fname="forsterite.xyz",cartesian=True,na=3,nb=3,nc=3,long_fmt=False)
#uc_to_xyz(uc,fname="forsterite.xyz",cartesian=False)

## output cif file with assymetric unit (similiar to what we read in)
uc_to_cif(uc,fname="forsterite_frac.cif",p1_list=False)

## output cif file with P1 unit cell contents (no symmetry)
uc_to_cif(uc,fname="forsterite_p1.cif",p1_list=True,na=2,nb=2,nc=2)

#######
# coordination calculations
from pyxrs.xtal.coord import coord_calcs
coord = coord_calcs(uc,rmax=2.5)
coord.write(long_fmt=False)
coord.write(fname="coord.out")

# run jmol
import subprocess
subprocess.run("jmol jmol_script.spt", shell=True)


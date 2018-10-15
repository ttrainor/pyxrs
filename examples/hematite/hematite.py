"""
Test some xtal calcs using hematite
"""
##############################
from pyxrs.xtal.unitcell import cif_to_uc, uc_to_cif, uc_to_xyz 
##############################

## generate a UnitCell instance for hematite
## display cell contents
uc = cif_to_uc('hematite.cif')
uc.write(long_fmt=False)

## print the P1 cell 
print("P1 cell")
at_list = uc.atom_list()
at_list.write()

## output xyz file
uc_to_xyz(uc,fname="hematite.xyz",na=2,nb=2,nc=2,long_fmt=True)
#uc_to_xyz(uc,fname="hematite.xyz",cartesian=False)

## output cif file with assymetric unit (similiar to what we read in)
#uc_to_cif(uc,fname="hematite_frac.cif",p1_list=False)

## output cif file with P1 unit cell contents (no symmetry)
#uc_to_cif(uc,fname="hematite_p1.cif",p1_list=True,na=2,nb=2,nc=2)

#######
# coordination calculations
from pyxrs.xtal.coord import coord_calcs
coord = coord_calcs(uc,rmax=2.5)
coord.write()
coord.write(fname="coord.out")

# run jmol
import subprocess
subprocess.run("jmol jmol_script.spt", shell=True)


"""
Test some xtal calcs using hematite
"""
##############################
from xtal.unitcell import read_cif, write_cif
##############################

## generate a UnitCell instance for hematite
## display cell contents
uc = read_cif('hematite.cif')
#uc = read_cif('COD_Fe2O3_1.cif')
#uc = read_cif('COD_Fe2O3_2.cif')

uc.write(long_fmt=False)

## print the P1 cell 
print("P1 cell")
at_list = uc.atom_list()
at_list.write_xyz()

## output xyz file
uc.write_xyz(fname="hematite.xyz",na=2,nb=2,nc=2,long_fmt=True)
#uc.write_xyz(fname="hematite.xyz",cartesian=False)

## output cif file with assymetric unit (similiar to what we read in)
#write_cif(uc,fname="hematite_frac.cif",p1_list=False)

## output cif file with P1 unit cell contents (no symmetry)
#write_cif(uc,fname="hematite_p1.cif",p1_list=True,na=2,nb=2,nc=2)

#######
# coordination calculations
from xtal.coord import coord_calcs
coord = coord_calcs(uc,rmax=2.5)
coord.write()
coord.write(fname="coord.out")

# run jmol
import subprocess
subprocess.run("jmol jmol_script.spt", shell=True)




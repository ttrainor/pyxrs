"""
Follow example 5.2.3 from Int Tables Vol A
"""
##############################
from xtal.unitcell import read_cif
##############################

## generate UnitCell instance for cristobalite_low
## display the unit cell contents and write an xyz file
print("**Cristobalite low")
uc_low = read_cif('cristobalite_low.cif')
uc_low.show()
uc_low.write_xyz("cristobalite_low.xyz",na=3,nb=3,nc=3)

## generate UnitCell instance for cristobalite_high
## display the unit cell contents and write an xyz file
print("**Cristobalite high")
uc_high = read_cif('cristobalite_high.cif')
uc_high.show()
uc_high.write_xyz("cristobalite_high.xyz",na=3,nb=3,nc=3)

####
# transform the cristobalite_low lattice to compare with cristobalite_high
# the below vectors define a new lattice for cristoblite low that should 
# map on the cristobalite high lattice
#
# note the new cell below is a centered cell.  therefore there is an additional
# set of symmetry operations that would be generated by adding the centering 
# vector "1/2,1/2,0" to all of the transformed symmetry operations...  (ie 
# the centered cell has double the number of symmetry ops)
#
# compare these results to example 5.2.3 in Int Tables Vol A
###
print("**Cristobalite low transformed")
Va = [1,1,0]; Vb = [-1,1,0]; Vc = [0,0,1]; shift = [0.25, 0.25, 0]
uc_low_new = uc_low.transform(Va=Va,Vb=Vb,Vc=Vc,shift=shift)
uc_low_new.show()
uc_low_new.write_xyz("cristobalite_low_new.xyz",na=3,nb=3,nc=3)

## print P1 coordinates of the new cell
print("**Cristobalite low transformed P1 cell")
at_list = uc_low_new.atom_list()
at_list.show()



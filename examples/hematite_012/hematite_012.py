"""
Generate a surface cell for hematite (1,-1,2) surface
"""
##############################
import numpy as num
from xtal.unitcell import read_cif, write_cif
from xtal.surface import SurfaceCell
##############################

### get the bulk unit cell
uc = read_cif('hematite.cif')

### this computes the surface cell.  note we pass the hex_to_rhom
### matrix so that the bulk hexagonal cell is transformed to a
### rhombahedral primitive cell before trying to find the surface
### indexing (ie, we want to find the mimimum surface unit cell) 
hex_to_rhom = [ [2/3, 1/3, 1/3], [-1/3, 1/3, 1/3], [-1/3, -2/3, 1/3]]
surf = SurfaceCell(uc,hkl=[1,-1,2],nd=2,term=+1,bulk_trns=hex_to_rhom)

### alternativley compute the surface cell by specifying the in-plane 
### and repeat vectors.  Note you can have the routine find the repeat
### vector, but you need to specify the hex_to_rhom matrix above in 
### order to find the one used below...
#surf = SurfaceCell(uc)
#Va=[1.,1.,0.]
#Vb=[-1/3., 1/3., 1/3.]
#Vr=[-2/3., 2/3., -1/3.]
#surf.set_surf_lattice(hkl=[1,-1,2],Va=Va,Vb=Vb,n=2,Vr=Vr)

### display the surface cell data
surf.write()

### this shows all the possible in-plane and repeat vectors
#print("\n\nSurface Vectors")
#for j in range(len(surf.Vs_lst)): print(surf.Vs_lst[j])
#print("\n\nRepeat Vectors")
#for j in range(len(surf.Vr_lst)): print(surf.Vr_lst[j])

## write surface coordinate files
surf.write_xyz(fname="hematite_012.xyz",cartesian=True,na=2,nb=2,nbulk=2,long_fmt=True)
surf.write_xyz(fname="hematite_012.frac",cartesian=False,na=1,nb=1,nbulk=2,long_fmt=True)

# coordination calculations
from xtal.coord import coord_calcs
#coord = coord_calcs(surf,rmax=2.5)
coord = coord_calcs(surf,rmax=2.5, labels=['O_1:t22','Fe_1:t3','Fe_1:b3'])
coord.write(long_fmt=False)
coord.write(fname="coord.out")

# run jmol
import subprocess
subprocess.run("jmol jmol_script.spt", shell=True)





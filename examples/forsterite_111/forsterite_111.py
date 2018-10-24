"""
Test some xtal calcs using forsterite
"""
##############################
from xtal.unitcell import cif_to_uc, uc_to_cif, uc_to_xyz 
from xtal.surface import SurfaceCell, surface_to_xyz
##############################

## generate a UnitCell instance 
## display cell contents
uc = cif_to_uc('forsterite.cif')
surf = SurfaceCell(uc,hkl=[1,1,1],nd=1,term=0)
surf.write()
surface_to_xyz(surf,fname="forsterite_111.xyz",cartesian=True,na=1,nb=1,nbulk=5,term=-99,long_fmt=False)

#######
# coordination calculations
from pyxrs.xtal.coord import coord_calcs
coord = coord_calcs(surf,rmax=2.5)
coord.write(long_fmt=False)
coord.write(fname="coord.out")

# run jmol
import subprocess
subprocess.run("jmol jmol_script.spt", shell=True)


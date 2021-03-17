"""
Test some xtal calcs using forsterite
"""
##############################
from xtal.unitcell import read_cif, write_cif
from xtal.surface import SurfaceCell
##############################

## generate a UnitCell instance 
## display cell contents
uc = read_cif('forsterite.cif')
surf = SurfaceCell(uc,hkl=[1,1,1],nd=1,term=0)
surf.write()
surf.write_xyz(fname="forsterite_111.xyz",cartesian=True,na=1,nb=1,nbulk=5,term=-99,long_fmt=False)

#######
# coordination calculations
from xtal.coord import coord_calcs
coord = coord_calcs(surf,rmax=2.5)
coord.write(long_fmt=False)
coord.write(fname="coord.out")

# run jmol
import subprocess
subprocess.run("jmol jmol_script.spt", shell=True)


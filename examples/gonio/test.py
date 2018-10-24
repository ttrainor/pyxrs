"""
Test psic geometry calcs
"""
#######################################################################
import numpy as num
from matplotlib import pyplot as plt
plt.ion()

from specfile import spec_scan
from psic_calcs import psic_from_spec_G, calc_active_area

#######################################################################
## file name and scan number
spfile = 'test2.spc'
sc_num = 3

## get scan and psic instance
scan = spec_scan(spfile,sc_num,geo='PSIC_APS_S13')
print(repr(scan))

psic = psic_from_spec_G(scan['G'],calcUB=True)
psic.set_angles(phi=scan['phi'][0],chi=scan['chi'][0],eta=scan['eta'][0],
                mu=scan['mu'][0],nu=scan['nu'][0],delta=scan['del'][0])
print(repr(psic))

## sample parameters
beam_slits = {'horz':.6,'vert':.8}
#det_slits  = None
det_slits  = {'horz':1,'vert':2}
#sample     = 10.
sample = {}
sample['polygon'] = [[1.,1.], [.5,1.5], [-1.,1.], [-1.,-1.],[0.,.5],[1.,-1.]]
sample['angles']  = {'phi':148.01,'chi':-14.7}

## do an area correction of scan data
npts = scan.dims[0]
Ic = num.zeros(npts)
ca = num.zeros(npts)
for j in range(npts):
    psic.set_angles(phi=scan['phi'][j],chi=scan['chi'][j],eta=scan['eta'][j],
                    mu=scan['mu'][j],nu=scan['nu'][j],delta=scan['del'][j])
    ca[j] = calc_active_area(psic, beam_slits, det_slits, sample, plot=False)
Ic = scan['IROI'] * ca / scan['io']
plt.figure(1)
plt.plot(scan['L'], num.log10(Ic))
plt.plot(scan['L'], num.log10(ca))

## make an active area plot
plt.figure(2)
j = 15
psic.set_angles(phi=scan['phi'][j],chi=scan['chi'][j],eta=scan['eta'][j],
                mu=scan['mu'][j],nu=scan['nu'][j],delta=scan['del'][j])
ca[j] = calc_active_area(psic, beam_slits, det_slits, sample, plot=True)

##
input('enter to continue')




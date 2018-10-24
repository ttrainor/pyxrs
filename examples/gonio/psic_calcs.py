"""
For psic geometry

Authors/Modifications:
----------------------
* Tom Trainor (tptrainor@alaska.edu)

"""
#######################################################################
import numpy as num
from xtal.lattice import Lattice
from xtlgeo.gonio_psic import Psic
from xtlgeo import gonio_psic, active_area

##########################################################################
def psic_from_spec_G(G, calcUB=True):
    """
    Given spec G array, return a psic instance

    * If calcUB is True recompute the orientation 
      matrix from the or information in G
    * If calcUB is False then use the UB array from 
      the spec G array.
    """
    psic = Psic()
    (cell,lam,or0,or1,UB,n) = _spec_psic_G(G)
    psic.n   = n
    psic.or0 = or0
    psic.or1 = or1
    psic.lattice = Lattice(*cell)
    psic.lam  = lam
    if calcUB==True:
        psic._calc_UB()
    else:
        psic._set_UB(UB)
    return psic

def _spec_psic_G(G):
    """
    Parse lattice and OR data
    from the spec G array for psic geometry

    Notes:
    ------
    If parsing angles from the G array:
      angles = G[x:y]
    where x:y depend on whether you are
    parsing out or0 or or1.

    We then assume the following:
      del = angles[0]; eta = angles[1]
      chi = angles[2]; phi = angles[3]
      nu  = angles[4]; mu  = angles[5]
    """
    #azimuthal reference vector, n (hkl)
    n = num.array(G[3:6],dtype=float)
    #lattice params a,b,c,alp,bet,gam
    cell = G[22:28]
    cell = num.array(cell,dtype=float)
    #lambda
    lam = float(G[66])
    #
    def _or_angles(angles):
        return {'phi':angles[3],'chi':angles[2],'eta':angles[1],'mu':angles[5],
                'delta':angles[0],'nu':angles[4]}
    # or0
    or0 = {}
    or0['h']   = num.array(G[34:37],dtype=float)
    or0.update(_or_angles(num.array(G[40:46],dtype=float)))
    or0['lam'] = float(G[52])
    # or1
    or1 = {}
    or1['h']   = num.array(G[37:40],dtype=float)
    or1.update(_or_angles(num.array(G[46:52],dtype=float)))
    or1['lam'] = float(G[53])
    # UB
    UB = num.array(G[54:63], dtype=float)
    UB.shape = (3,3)
    UB = UB / (2*num.pi)
    return (cell,lam,or0,or1,UB,n)

def calc_kappa_angles(phi,chi,eta):
    """
    kappa angles
    To calc kappa, keta and kphi use:
    kap_alp = 50.031;
    keta    = eta - asin(-tan(chi/2)/tan(kap_alp))
    kphi    = phi - asin(-tan(chi/2)/tan(kap_alp))
    kappa   = asin(sin(chi/2)/sin(kap_alp))
    """
    kap_alp = 50.031
    keta    = eta - arcsind(-tand(chi/2.)/tand(kap_alp))
    kphi    = phi - arcsind(-tand(chi/2.)/tand(kap_alp))
    kappa   = asind(sind(chi/2.)/sind(kap_alp))

def calc_active_area(gonio, beam_slits, det_slits, sample, plot=False):
    """
    Compute active area correction (c_a = A_beam/A_int**2)
    
    Use to correct scattering data for area effects including:
        Detector view and spilloff (A_int/A_beam)
        Normalization to unit surface area (1/A_beam)
    Therefore, 
        Ic = Im * ca = Im/A_ratio 
        A_ratio = A_int/(A_beam**2) 
    where
        A_int = intersection area (area of the beam on sample
                that is viewed by detector)
        A_beam = total beam area
    """
    alpha = gonio.pangles['alpha']
    beta  = gonio.pangles['beta']
    if plot == True:
        print('Alpha = ', alpha, ', Beta = ', beta)
    if alpha < 0.0:
        print('alpha is less than 0.0')
        return 0.0
    elif beta < 0.0:
        print('beta is less than 0.0')
        return 0.0
    # get beam vectors
    bh = beam_slits['horz']
    bv = beam_slits['vert']
    beam = gonio_psic.beam_vectors(h=bh,v=bv)
    # get det vectors
    if det_slits == None:
        det = None
    else:
        dh = det_slits['horz']
        dv = det_slits['vert']
        det  = gonio_psic.det_vectors(h=dh,v=dv,
                                      nu=gonio.angles['nu'],
                                      delta=gonio.angles['delta'])
    # get sample poly
    if type(sample) == dict:
        sample_dia    = sample.get('dia',0.)
        sample_vecs   = sample.get('polygon',None)
        sample_angles = sample.get('angles',{})
        #
        if sample_vecs != None and sample_dia <= 0.:
            sample = gonio_psic.sample_vectors(sample_vecs,
                                               angles=sample_angles,
                                               gonio=gonio)
        elif sample_dia > 0.:
            sample = sample_dia
        else:
            sample = None
    else:
        sample = sample
    # compute active_area
    (A_beam,A_int) = active_area.active_area(gonio.nm,ki=gonio.ki,
                                             kr=gonio.kr,beam=beam,det=det,
                                             sample=sample,plot=plot)
    if A_int == 0.:
        ca = 0.
    else:
        ca = A_beam/(A_int**2)
    return ca

#######################################################################
#######################################################################
def test_1():
    """
    G parsing from specfile.py
    """
    G = [0.0, 0.0, 1.0, 0.0039915744589999998, 0.00075650941450000001, 1.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 50.0, 0.0, 0.0, 1.0, 4.0, 4.0, 5.0, 4.0, 0.0, 0.0, 8.0939999999999994,
         4.9880000000000004, 6.0709999999999997, 90.0, 90.0, 90.0, 0.77627690969999996,
         1.2596602459999999, 1.034950635, 90.0, 90.0, 90.0, 0.0, 0.0, 4.0, 2.0, 0.0,
         2.4809999999999999, -0.00089999999999999998, 0.00080000000000000004, -0.1244,
         175.2192, 31.965, 16.27375, 15.090299999999999, 7.5477999999999996, 11.9002,
         -96.048199999999994, 17.594249999999999, 0.35625000000000001, 0.83580100000000002,
         0.83580100000000002, -0.23926186720000001, 1.1983306439999999, -0.0026481987470000001,
         -0.73847260000000003, -0.38826052760000002, -0.0050577973450000001, -0.004221190899,
         0.001168841303, 1.034934888, -0.00011346681140000001, 0.0002801762336, 5.9400141580000003,
         0.83580100000000002, 23.955432519999999, 24.314067479999999, 0.28474999779999999,
         48.269500000000001, 0.075104988760000005, 0.17931763710000001, -0.00040199175790000001,
         -0.00014478299110000001, 0.48309999999999997, 108.00069999999999, 2.0, 0.0, 0.0, 0.0, 0.0,
         12.0, 0.0, 0.0, 2.0802999999999998, 123.1461, 0.0, 0.0, 0.0, 0.0, -180.0, -180.0, -180.0,
         -180.0, -180.0, -180.0, -180.0, -180.0, -180.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0]
    # 
    psic = psic_from_spec_G(G)
    psic.set_angles(phi=178.1354,chi=-0.1344,eta=0.0002,
                    mu=24.4195,nu=48.2695,delta=-0.0003)
    print("#########################################")
    print("h calc is:", psic.h)
    print("should be: [-0.000113467 0.000280176 5.94001]")
    print("----")
    print("alpha calc is:", psic.pangles['alpha'])
    print("should be:", G[67])
    print("----")
    print("beta calc is:", psic.pangles['beta'])
    print("should be:", G[68])
    print("----")
    print("omega calc is:", psic.pangles['omega'])
    print("should be:", G[69])
    print("----")
    print("tth calc is:", psic.pangles['tth'])
    print("should be:", G[70])
    print("----")
    print("psi calc is:", psic.pangles['psi'])
    print("should be:", G[71])
    print("----")
    print("tau calc is:", psic.pangles['tau'])
    print("should be:", G[72])
    print("----")
    print("qaz calc is:", psic.pangles['qaz'])
    print("should be:", G[73])
    print("----")
    print("naz calc is:", psic.pangles['naz'])
    print("should be:", G[74])
    print("----")
    print("sigma_az calc is:", psic.pangles['sigma_az'])
    print("should be:", G[75])
    print("----")
    print("tau_az calc is:", psic.pangles['tau_az'])
    print("should be:", G[76])
    print("#########################################")
    return psic

def test_2(spfile='test1.spc', sc_num=1, pnt=7):
    """
    read a spec file and parse spec G
    """
    scan = spec_scan(spfile,sc_num,geo='PSIC_APS_S13')
    print(repr(scan))
    #
    psic = psic_from_spec_G(scan['G'],calcUB=False)
    psic.set_angles(phi=scan['phi'][pnt],chi=scan['chi'][pnt],eta=scan['eta'][pnt],
                    mu=scan['mu'][pnt],nu=scan['nu'][pnt],delta=scan['del'][pnt])
    print(repr(psic))
    return scan, psic

##########################################################################
if __name__ == "__main__":
    """
    test 
    """
    from specfile import spec_scan
    #psic = test_1()
    scan, psic = test_2()





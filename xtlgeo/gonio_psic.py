"""
Gonio calcs for 6 circle psic geometry

Authors / Modifications:
------------------------
Tom Trainor (tptrainor@alaska.edu)
Frank Heberling (Frank.Heberling@ine.fzk.de)

Notes:
------
The following describes the calculations used to define the orientation
matrix of a xtal mounted on a goniometer and to convert between
motor angles and hkl and visa-versa.  Additional functions are also
provided for computing beam and detector slit aperature vectors
and sample position vectors etc..  

Define the matrix B which transforms the indicies of the 
reciprocal lattice vector, h=[h,k,l], to a cartesian basis,
hc=[xc,yc,zc], where the cartesian system has been chosen so that:
    x is parrallel to ar, 
    y is in the plane of ar and br, 
    z is perp to the plane of ar and br
and (ar,br,cr) are the reciprocal lattice basis vectors
(see Busing and Levy)

Therefore:  
    hc = B*h
Note that hc (and the cartesian basis) is defined with respect to
the crystal lattice and therefore is independant of the gonio angles.  

U is the matrix that defines how the sample is mounted on the
diffractometer.  For example we define the lab frame coordinate
system for the psic geometry such that:
    x is vertical (pointing to the ceiling of the hutch, perpendicular to beam path)
    y is directed along the incident beam path
    z make the system right handed and lies in the horizontal scattering plane
    (i.e. z is parallel to the phi axis)
Therefore with all gonio angles set to zero the lab xyz axis are coincident
with the goniometer axes in a well defined way (depending on the axis definitions).
Note that the yz plane is the horizontal scattering plane, and yx is the vertical
scattering plane.

With the instrument settings all at zero angles (phi frame), and with the
sample oriented in an arbitrary maner, the matrix U is used to calculate
the lab frame indicies (hphi) of a given reciprocal lattice vector, h=[h,k,l],
according to:
    hphi = U*hc = U*B*h

Therefore, U is a simple rotation matrix which gives
the indicies of h in the phi frame accounting for how 
the sample is mounted. ie this matrix (or its transpose)
would take the cartesian reciprocal basis vectors (hc) and
rotate them to be coincident with the lab frame (phi-frame)
basis vectors

The algorithm to determine U from a set of reflections and angles
is given by Busing and Levy.

We can then orient the hphi vector in the lab frame applying the rotation
matricies of the various gonio axes.  

For psic geom from H. You's paper, the order of sample
axis rotations is:
1. phi rotation (matrix = P)
2. chi rotation (matrix = X)
3. eta rotation (matrix = H)
4. mu rotation  (matrix = M)

To calculate the lab frame coords of h --> hm:
    hm = M*H*X*P*hphi = M*H*X*P*U*B*h
Or letting Z = M*H*X*P, 
    hm = Z*hphi

Therefore, hm gives the indicies of the recip lattice vector, h,
in the lab frame cartesian basis after the sample is rotated.

The diffraction condition specifies that:
    Q = (2*pi)*h 

The lab frame coordinates of Q can be calc from:
    Qm = kr - ki
where ki and kr in the lab frame are given (for psic) by :
    ki = (2*pi/lam)*[ 0, 1, 0 ]
    kr = (2*pi/lam)*[ sin(del), cos(nu)*cos(del),sin(nu)*cos(del) ]
    Qm = kr - ki;
or
        |      sin(del)      |
Qm = k* |cos(del)*cos(nu) - 1|
        |  cos(del)*sin(nu)  |

Therefore the diffraction condition in the lab frame is
    Qm = (2*pi)*hm = (2*pi)*Z*hphi = (2*pi)*Z*U*B*h

Now given an orientation matrix and set of gonio
angles we can then solve for hphi
    hphi = inv(Z) * Qm / (2*pi)

The reciprocal lattice indicies (h) are then calc from
    h = inv(UB)*hphi
This gives the hkl values of the vector that is in the
diffraction condition for a given set of angles.

References:
-----------
1. H. You, J. Appl. Cryst. (1999) 32, 614-623
2. Busy and Leving

"""
##########################################################################
import numpy as num
import types
import copy

from xtal._common import cosd, sind, tand
from xtal._common import arccosd, arcsind, arctand
from xtal._common import cartesian_mag, cartesian_angle
from xtal.lattice import Lattice

##########################################################################
class Psic:
    """
    Orientation calculations for Psic geometry.

    The default (dummy) orientation matrix is set up 
    assuming the sample is mounted such that (001) plane
    is perpendicular to the eta and phi rot axes
    (ie c-axis is parrallel to the eta and phi rot axes)
    and that the b-axis is parallel to the nu and mu rot axes
    (ie parrallel to the lab frame Z)
    """
    ###################################################
    def __init__(self,a=10.,b=10.,c=10.,alpha=90.,beta=90.,gamma=90.,lam=1.0):
        """
        Initialize

        Parameters:
        -----------
        * a,b,c in angstroms 
        * alpha, beta, gamma in degrees
        * lam wavelength in angstroms
        """
        # set lattice and lambda
        self.lattice = Lattice(a,b,c,alpha,beta,gamma)
        self.lam     = lam
        
        # gonio angles 
        self.angles={'phi':0.0,'chi':0.0,'eta':0.0,'mu':0.0,
                     'nu':0.0,'delta':0.0}

        # psuedo angles
        self.pangles = {}
        self.calc_psuedo = True

        # n reference vector (HKL)
        # ie surface normal vector for psuedo angles
        self.n = num.array([0.,0.,1.], dtype=float)
        
        # Z and calculated h
        self.Z  = []
        self.Q  = []
        self.ki = []
        self.kr = []
        self.h = [0.,0.,0.]

        # dummy primary reflection
        tth = self.lattice.tth([0.,0.,1.],lam=self.lam)
        self.or0={'h':num.array([0.,0.,1.]),
                  'phi':0.0,'chi':0.0,'eta':0.0,'mu':tth/2.,
                  'nu':tth,'delta':0.0,'lam':self.lam}
        
        # dummy secondary reflection
        tth = self.lattice.tth([0.,1.,0.],lam=self.lam)
        self.or1={'h':num.array([0.,1.,0.]),
                  'phi':0.0,'chi':0.0,'eta':tth/2.,'mu':0.0,
                  'nu':0.0,'delta':tth,'lam':self.lam}

        # Compute OR matricies
        self.U = []
        self.B = []
        self.UB = []
        self._calc_UB()

    ###################################################
    def __repr__(self,):
        """ display """
        lout = self.lattice.__repr__()
        lout = "%sPrimary:\n   h=%3.2f,k=%3.2f," % (lout,self.or0['h'][0],self.or0['h'][1])
        lout = "%sl=%3.2f, lam=%6.6f\n" % (lout,self.or0['h'][2],self.or0['lam'])
        lout = "%s   phi=%6.3f,chi=%6.3f," % (lout,self.or0['phi'],self.or0['chi'])
        lout = "%seta=%6.3f,mu=%6.3f," % (lout,self.or0['eta'],self.or0['mu'])
        lout = "%snu=%6.3f,delta=%6.3f\n" % (lout,self.or0['nu'],self.or0['delta'])
        #
        lout = "%sSecondary:\n   h=%3.2f,k=%3.2f," % (lout,self.or1['h'][0],self.or1['h'][1])
        lout = "%sl=%3.2f, lam=%6.6f\n" % (lout,self.or1['h'][2],self.or1['lam'])
        lout = "%s   phi=%6.3f,chi=%6.3f," % (lout,self.or1['phi'],self.or1['chi'])
        lout = "%seta=%6.3f,mu=%6.3f," % (lout,self.or1['eta'],self.or1['mu'])
        lout = "%snu=%6.3f,delta=%6.3f\n" % (lout,self.or1['nu'],self.or1['delta'])
        #
        lout = "%sSetting:" % (lout)
        lout = "%s   h=%3.2f,k=%3.2f,l=%3.2f\n" % (lout,self.h[0],self.h[1],self.h[2])
        lout = "%s   phi=%6.3f,chi=%6.3f," % (lout,self.angles['phi'],self.angles['chi'])
        lout = "%seta=%6.3f,mu=%6.3f," % (lout,self.angles['eta'],self.angles['mu'])
        lout = "%snu=%6.3f,delta=%6.3f\n" % (lout,self.angles['nu'],self.angles['delta'])
        #
        if self.calc_psuedo:
            lout = "%s   TTH=%6.3f," % (lout,self.pangles['tth'])
            lout = "%sSIGMA_AZ=%6.3f," % (lout,self.pangles['sigma_az'])
            lout = "%sTAU_AZ=%6.3f," % (lout,self.pangles['tau_az'])
            lout = "%sN_AZ=%6.3f," % (lout,self.pangles['naz'])
            lout = "%sALPHA=%6.3f," % (lout,self.pangles['alpha'])
            lout = "%sBETA=%6.3f\n" % (lout,self.pangles['beta'])
            lout = "%s   TAU=%6.3f," % (lout,self.pangles['tau'])
            lout = "%sPSI=%6.3f," % (lout,self.pangles['psi'])
            lout = "%sQ_AZ=%6.3f," % (lout,self.pangles['qaz'])
            lout = "%sOMEGA=%6.3f," % (lout,self.pangles['omega'])
        #
        return lout
    
    ###################################################
    def set_lat(self,a=None,b=None,c=None,alpha=None,
                beta=None,gamma=None,lam=None):
        """
        Update lattice parameters and lambda

        Parameters:
        -----------
        * a,b,c in angstroms 
        * alpha, beta, gamma in degrees,
        * lambda in angstroms
        """
        self.lattice.update(a=a,b=b,c=c,alpha=alpha,
                            beta=beta,gamma=gamma)
        self.lam = lam
        self._calc_UB()

    ################################################### 
    def set_or0(self,h=None,phi=None,chi=None,eta=None,
                mu=None,nu=None,delta=None,lam=None):
        """
        Set / adjust the primary orientation reflection

        Parameters:
        -----------
        * h is the hkl array of the reflection
        * the rest of the parameters are motor angles
          in degrees,
        * lam is the wavelength in angstroms
          If lam = None, then lambda defined for the lattice
          is used.
        """
        if h!=None:     self.or0['h'] = num.array(h,dtype=float)
        if phi!=None:   self.or0['phi']=float(phi)
        if chi!=None:   self.or0['chi']=float(chi)
        if eta!=None:   self.or0['eta']=float(eta)
        if mu!=None:    self.or0['mu']=float(mu)
        if nu!=None:    self.or0['nu']=float(nu)
        if delta!=None: self.or0['delta']=float(delta)
        if lam!= None:  self.or0['lam']=float(lam)
        self._calc_UB()

    def set_or1(self,h=None,phi=None,chi=None,eta=None,
                mu=None,nu=None,delta=None,lam=None):
        """
        Set / adjust the secondary orientation reflection

        Parameters:
        -----------
        * h is the hkl array of the reflection
        * the rest of the parameters are motor angles
          in degrees,
        * lam is the wavelength in angstroms
          If lam = None, then lambda defined for the lattice
          is used.
        """
        if h!=None:     self.or1['h'] = num.array(h,dtype=float)
        if phi!=None:   self.or1['phi']=float(phi)
        if chi!=None:   self.or1['chi']=float(chi)
        if eta!=None:   self.or1['eta']=float(eta)
        if mu!=None:    self.or1['mu']=float(mu)
        if nu!=None:    self.or1['nu']=float(nu)
        if delta!=None: self.or1['delta']=float(delta)
        if lam!= None:  self.or1['lam']=float(lam)
        self._calc_UB()

    def swap_or(self,):
        """
        Swap the primary and secondary reflection
        """
        tmp = copy.copy(self.or0)
        self.or0 = copy.copy(self.or1)
        self.or1 = tmp
        self._calc_UB()

    ################################################### 
    def _calc_UB(self,):
        """
        Calculate the orientation matrix, U,
        from the primary and secondary
        reflectons and given lattice
        """
        # check to make sure the OR's arent parallel
        angle = self.lattice.angle(self.or0['h'],self.or1['h'],recip=True)
        angle = num.fabs(num.round(angle,decimals=1))
        if angle == 0:
            print("**Error, orientation reflections are parallel\n")
            return 

        # use these, note they are used below on vectors
        # defined in the cartesian lab frame basis
        cross = num.cross
        norm  = num.linalg.norm

        #Calculate the B matrix
        (a,b,c,alp,bet,gam)       = self.lattice.cell()
        (ar,br,cr,alpr,betr,gamr) = self.lattice.rcell()
        B = num.array([[ar,  br*cosd(gamr),     cr*cosd(betr)        ],
                       [0.,  br*sind(gamr),  -cr*sind(betr)*cosd(alp)],
                       [0.,      0.,                 1./c            ]])
        self.B = B
        
        # calc Z and Q for the OR reflections
        Z1 = calc_Z(self.or0['phi'],self.or0['chi'],self.or0['eta'],self.or0['mu'])
        Q1 = calc_Q(self.or0['nu'],self.or0['delta'],self.or0['lam'])
        #
        Z2 = calc_Z(self.or1['phi'],self.or1['chi'],self.or1['eta'],self.or1['mu'])
        Q2 = calc_Q(self.or1['nu'],self.or1['delta'],self.or1['lam'])

        # calc the phi frame coords for diffraction vectors
        # note divide out 2pi since the diffraction condition
        # is 2pi*h = Q
        vphi_1 = num.dot(num.linalg.inv(Z1), Q1/(2.*num.pi) )
        vphi_2 = num.dot(num.linalg.inv(Z2), Q2/(2.*num.pi) )
        
        #calc cartesian coords of h vectors
        hc_1 = num.dot(self.B, self.or0['h'])
        hc_2 = num.dot(self.B, self.or1['h'])

        #So at this point the following should be true:
        #     vphi_1 = U*hc_1
        #     vphi_2 = U*hc_2
        # and we could use these relations to solve for U.
        # But U solved directly from above is likely not to be orthogonal
        # since the angles btwn (vphi_1 and vphi_2) and (hc_1 and hc_2) are 
        # not exactly the same due to experimental errors.....
        # Therefore, get an orthogonal solution for U from the below treatment
        
        #define the following normalized vectors from hc vectors 
        tc_1 = hc_1 / norm(hc_1)
        tc_3 = cross(tc_1, hc_2) / norm(cross(tc_1, hc_2))
        tc_2 = cross(tc_3, tc_1) / norm(cross(tc_3, tc_1))

        #define tphi vectors from vphi vectors
        tphi_1 = vphi_1 /  norm(vphi_1)
        tphi_3 = cross(tphi_1,vphi_2) / norm(cross(tphi_1,vphi_2))
        tphi_2 = cross(tphi_3,tphi_1) / norm(cross(tphi_3,tphi_1))

        #define the following matrices
        Tc   = num.transpose(num.array([tc_1,tc_2,tc_3]))
        Tphi = num.transpose(num.array([tphi_1,tphi_2,tphi_3]))
        
        # calc orientation matrix U
        # note either of the below work since Tc is orthogonal
        #self.U = num.dot(Tphi, Tc.transpose())
        self.U = num.dot(Tphi, num.linalg.inv(Tc))

        # calc UB
        self.UB = num.dot(self.U,self.B)

        #update h and psuedo angles...
        self.set_angles()

    def _set_UB(self,UB):
        """
        Set the orientation matrix
        """
        self.UB = UB
        self.set_angles()

    ###################################################
    def set_angles(self,phi=None,chi=None,eta=None,
                   mu=None,nu=None,delta=None):
        """
        Set goniometer angles (all in degrees)
        """
        if phi!=None:   self.angles['phi']=float(phi)
        if chi!=None:   self.angles['chi']=float(chi)
        if eta!=None:   self.angles['eta']=float(eta)
        if mu!=None:    self.angles['mu']=float(mu)
        if nu!=None:    self.angles['nu']=float(nu)
        if delta!=None: self.angles['delta']=float(delta)
        # update h, also calc Z etc..
        self._calc_h()
        # update psuedo
        self._update_psuedo()
        
    def _calc_h(self,):
        """
        Calculate the hkl values of the vector that is in the
        diffraction condition for the given set of angles.  

        Notes:
        ------
        Solve for hphi using Z and lab frame Q:
           hphi = inv(Z) * Q / (2*pi)
        then calc h from
           h = inv(UB)*hphi
        """
        self.Z = calc_Z(phi=self.angles['phi'],chi=self.angles['chi'],
                        eta=self.angles['eta'],mu=self.angles['mu'])
        (Q,ki,kr) = calc_Q(self.angles['nu'],
                           self.angles['delta'],
                           self.lam,ret_k=True)
        self.Q=Q
        self.ki=ki
        self.kr=kr
        
        hphi = num.dot(num.linalg.inv(self.Z), self.Q) / (2.*num.pi) 
        h    = num.dot(num.linalg.inv(self.UB), hphi)
        self.h = h
        
    ###################################################
    def set_n(self,n=[0,0,1]):
        """
        Set n, the reference vector used for psuedo angles.
        The vector n is given in hkl values.  see calc_n
        to determine n from chi and phi settings
        """
        self.n = num.array(n,dtype=float)
        self._update_psuedo()

    def calc_n(self,fchi=0.0,fphi=0.0):
        """
        Calculate the hkl values of a reference vector given
        the chi and phi settings that align this
        vector with the eta axis.

        Notes:
        ------
        This algorith is used, for example, 
        to compute the surface normal from the (flat) chi and
        (flat) phi angles that leave an optical reflection in
        a fixed position during an eta rotation 

        Note the vector is normalized such that
        the largest component is unity,
        ie n_hkl isn't a unit vector!
        """
        # polar angles
        sig_az = -fchi
        tau_az = -fphi

        # this block converts the chi and phi values to correctly 
        # defined polar coordinates, ie 0<= sig_az <= 180deg .....
        if sig_az < 0.:
            sig_az = -1.*sig_az
            if tau_az < 0.:
                tau_az = 180. + tau_az
            elif tau_az > 0.:
                tau_az = tau_az - 180.

        # n in the unrotated lab frame (ie phi frame):
        # this is a unit vector!
        n_phi = num.array([ sind(sig_az)*cosd(tau_az),
                           -sind(sig_az)*sind(tau_az), 
                                  cosd(sig_az)        ])
        # n in HKL
        n_hkl = num.dot(num.linalg.inv(self.UB), n_phi)
        n_hkl = n_hkl / num.max(num.abs(n_hkl))
        
        # note if l-component is negative, then its
        # pointing into the surface (ie assume positive L
        # is the positive direction away from the surface)
        # careful here!!
        if n_hkl[2] < 0.:
            n_hkl = -1.*n_hkl

        # set n which triggers recalc of
        # all the psuedo angles
        self.set_n(n_hkl)

    ################################################### 
    ## Pseudo angles
    ###################################################
    def _update_psuedo(self):
        """
        Compute psuedo angles
        
        Note:
        -----
        use this to compute psuedo angles rather than
        individual calls.  ie some psuedo angles depend on others
        so its important that the calcs are executed in the correct
        order.  Also important is that _calc_h is called before this...
        """
        self.pangles = {}
        if self.calc_psuedo == True:
            self._calc_tth()
            self._calc_nm()
            self._calc_sigma_az()
            self._calc_tau_az()
            self._calc_naz()
            self._calc_alpha()
            self._calc_beta()
            self._calc_tau()
            self._calc_psi()
            self._calc_qaz()
            self._calc_omega()
    
    def _calc_tth(self):
        """
        Calculate 2Theta, the scattering angle

        Notes:
        ------
        This should be the same as:
          (ki,kr) = calc_kvecs(nu,delta,lambda)
           tth = cartesian_angle(ki,kr)

        You can also get this given h, the reciprocal lattice
        vector that is in the diffraction condition.  E.g.
          h   = self.calc_h()
          tth = self.lattice.tth(h,self.lam)
        """
        nu    = self.angles['nu']
        delta = self.angles['delta']
        tth   = arccosd(cosd(delta)*cosd(nu))
        self.pangles['tth'] = tth

    def _calc_nm(self):
        """
        Calculate the rotated cartesian lab indicies
        of the reference vector n = nm.  Note nm is
        normalized.  

        Notes:
        ------
        The reference vector n is given in recip
        lattice indicies (hkl)
        """
        # calc n in the rotated lab frame and make a unit vector
        n  = self.n
        Z  = self.Z
        UB = self.UB
        nm = num.dot(num.dot(Z,UB),n)
        nm = nm / cartesian_mag(nm)
        self.nm = nm
    
    def _calc_sigma_az(self):
        """
        sigma_az = angle between the z-axis and n
        in the phi frame
        """
        # calc n in the lab frame (unrotated) and make a unit vector
        n_phi = num.dot(self.UB,self.n)
        n_phi = n_phi / cartesian_mag(n_phi)
        
        # note result of acosd is between 0 and pi
        # get correct sign from the sign of the x-component
        #sigma_az = num.sign(n_phi[0])*arccosd(n_phi[2])
        sigma_az = arccosd(n_phi[2])
        self.pangles['sigma_az'] = sigma_az

    def _calc_tau_az(self):
        """
        tau_az = angle between the projection of n in the
        xy-plane and the x-axis in the phi frame
        """
        # calc n in the lab frame (unrotated) and make a unit vector
        n_phi = num.dot(self.UB,self.n)
        n_phi = n_phi / cartesian_mag(n_phi)

        tau_az = num.arctan2(-n_phi[1], n_phi[0])
        tau_az = tau_az*180./num.pi
        self.pangles['tau_az'] = tau_az

    def _calc_naz(self):
        """
        calc naz, this is the angle btwn the reference vector n 
        and the yz plane at the given angle settings
        """
        # get norm reference vector in cartesian lab frame
        nm  = self.nm
        naz = num.arctan2( nm[0], nm[2] )
        naz = num.degrees(naz)
        self.pangles['naz'] = naz

    def _calc_alpha(self):
        """
        Calc alpha, ie incidence angle or angle btwn 
        -1*k_in (which is parallel to lab-y) and the
        plane perp to the reference vector n.
        """
        nm = self.nm
        ki = num.array([0.,-1.,0.])
        alpha = arcsind(num.dot(nm,ki))
        self.pangles['alpha'] = alpha

    def _calc_beta(self):
        """
        Calc beta, ie exit angle, or angle btwn k_r and the
        plane perp to the reference vector n

        Notes:
        ------
        beta = arcsind(2*sind(tth/2)*cosd(tau)-sind(alpha))
        """
        # calc normalized kr
        #delta = self.angles['delta']
        #nu    = self.angles['nu']
        #kr = num.array([sind(delta),
        #                cosd(nu)*cosd(delta),
        #                sind(nu)*cosd(delta)])
        nm = self.nm
        kr = self.kr / cartesian_mag(self.kr)
        beta = arcsind(num.dot(nm, kr))
        self.pangles['beta'] = beta

    def _calc_tau(self):
        """
        Calc tau, this is the angle btwn n and the scattering-plane
        defined by ki and kr.  ie the angle between n and Q

        Notes:
        ------
        Can also calc from:
         tau = acos( cosd(alpha) * cosd(tth/2) * cosd(naz - qaz) ...
                    + sind(alpha) * sind(tth/2) ) 
        """
        tau = cartesian_angle(self.Q, self.nm)
        self.pangles['tau'] = tau

    def _calc_psi(self):
        """
        calc psi, this is the azmuthal angle of n wrt Q. 
        ie for tau != 0, psi is the rotation of n about Q

        Notes:
        ------
        Note this must be calc after tth, tau, and alpha!
        """
        tau   = self.pangles['tau']
        tth   = self.pangles['tth']
        alpha = self.pangles['alpha']
        #beta = self.calc_beta()
        #xx = (-cosd(tau)*sind(tth/2.) + sind(beta))
        xx    = (cosd(tau)*sind(tth/2.) - sind(alpha))
        denom = (sind(tau)*cosd(tth/2.))
        if denom == 0: 
            self.pangles['psi'] = 0.
            return
        xx    = xx / denom
        psi = arccosd( xx )
        self.pangles['psi'] = psi

    def _calc_qaz(self):
        """
        Calc qaz, the angle btwn Q and the yz plane 
        """
        nu    = self.angles['nu']
        delta = self.angles['delta']
        qaz = num.arctan2(sind(delta), cosd(delta)*sind(nu) )
        qaz = num.degrees(qaz)
        self.pangles['qaz'] = qaz

    def _calc_omega(self):
        """
        calc omega, this is the angle between Q and the plane
        which is perpendicular to the axis of the chi circle.

        Notes:
        ------
        For nu=mu=0 this is the same as the four circle def:
        omega = 0.5*TTH - TH, where TTH is the detector motor (=del)
        and TH is the sample circle (=eta).  Therefore, for 
        mu=nu=0 and del=0.5*eta, omega = 0, which means that Q
        is in the plane perpendicular to the chi axis.

        Note check sign of results??? 
        """
        phi=self.angles['phi']
        chi=self.angles['chi']
        eta=self.angles['eta']
        mu=self.angles['mu']
        H = num.array([[ cosd(eta), sind(eta), 0.],
                       [-sind(eta), cosd(eta), 0.],
                       [   0.,         0.,     1.]],float)
        M  = num.array([[  1.,         0.,     0.      ],
                        [  0.,      cosd(mu), -sind(mu)],
                        [  0.,      sind(mu), cosd(mu)]],float)
        # check the mult order here!!!!
        # T = num.dot(H.transpose(),M.transpose())
        T     = num.dot(M.transpose(),H.transpose())
        Qpp   = num.dot(T,self.Q)
        #omega = -1.*cartesian_angle([Qpp[0], 0, Qpp[2]],Qpp)
        omega = cartesian_angle([Qpp[0], 0, Qpp[2]],Qpp)
        self.pangles['omega'] = omega

##########################################################################
def calc_Z(phi=0.0,chi=0.0,eta=0.0,mu=0.0):
    """
    Calculate the psic goniometer rotation matrix Z
    for the 4 sample angles. Angles are in degrees

    Notes:
    ------
    Z is the matrix that rotates a vector defined in the phi frame
    ie a vector defined with all angles zero => vphi.  After rotation
    the lab frame coordinates of the vector => vm are given by:
         vm = Z*vphi
    """
    P = num.array([[ cosd(phi), sind(phi), 0.],
                   [-sind(phi), cosd(phi), 0.],
                   [  0.,          0.,     1.]],float)
    X = num.array([[ cosd(chi), 0., sind(chi)],
                   [   0.,      1.,    0.],
                   [-sind(chi), 0., cosd(chi)]],float)
    H = num.array([[ cosd(eta), sind(eta), 0.],
                   [-sind(eta), cosd(eta), 0.],
                   [   0.,         0.,     1.]],float)
    M  = num.array([[  1.,         0.,     0.      ],
                    [  0.,      cosd(mu), -sind(mu)],
                    [  0.,      sind(mu), cosd(mu)]],float)
    Z = num.dot(num.dot(num.dot(M,H),X),P)
    return Z

def calc_Q(nu=0.0,delta=0.0,lam=1.0,ret_k=False):
    """
    Calculate psic Q in the cartesian lab frame.
    nu and delta are in degrees, lam is in angstroms

    if ret_k == True return tuple -> (Q,ki,kr)
    """
    (ki,kr) = calc_kvecs(nu=nu,delta=delta,lam=lam)
    Q = kr - ki
    if ret_k == True:
        return (Q,ki,kr)
    else:
        return Q

def calc_kvecs(nu=0.0,delta=0.0,lam=1.0):
    """
    Calculate psic ki, kr in the cartesian lab frame.
    nu and delta are in degrees, lam is in angstroms
    """
    k  = (2.* num.pi / lam)
    ki = k * num.array([0.,1.,0.], dtype=float)
    kr = k * num.array([sind(delta),
                        cosd(nu)*cosd(delta),
                        sind(nu)*cosd(delta)], dtype=float)
    return (ki,kr)

def calc_D(nu=0.0,delta=0.0):
    """
    Calculate the detector rotation matrix.
    Angles are in degrees

    Notes:
    ------
    D is the matrix that rotates a vector defined in the phi frame
    ie a vector defined with all angles zero => vphi.  After rotation
    the lab frame coordinates of the vector => vm are given by:
         vm = D*vphi
    For example 
                            |0|   
         kr_phi = (2pi/lam) |1|
                            |0|
    Since kr is defined by the detector rotation, the lab frame 
    coordinates of the kr vector after detector rotation are
         kr_m = D*kr_phi
    
    """
    D1 = num.array([[cosd(delta),  sind(delta),  0.], 
                    [-sind(delta), cosd(delta),  0.],
                    [     0.     ,     0.     ,  1.]])
          
    D2 = num.array([[    1.,     0.   ,      0.  ],
                    [    0.,  cosd(nu), -sind(nu)], 
                    [    0.,  sind(nu),  cosd(nu)]])
          
    D = num.dot(D2,D1)
    return (D)

def beam_vectors(h=1.0,v=1.0):
    """
    Compute the beam apperature vectors in lab frame

    Parameters:
    -----------
    * h = beam horz width (total slit width in lab-z,
      or the horizontal scattering plane)
    * v = beam vert hieght (total slit width in lab-x,
      or the vertical scattering plane)

    Notes:
    ------
    The slit settings, defined wrt psic phi-frame
    Assume these are centered on the origin
    """
    # beam vectors, [x,y,z], in lab frame
    bh = num.array([   0.,  0., 0.5*h])
    bv = num.array([0.5*v,  0.,    0.])

    # corners of beam apperature
    a =  bv + bh
    b =  bv - bh
    c = -bv - bh
    d = -bv + bh
    beam = [a,b,c,d]

    return beam

def det_vectors(h=1.0,v=1.0,nu=0.0,delta=0.0):
    """
    Compute detector apperature vectors in lab frame

    Parameters:
    -----------
    * h = detector horz width (total slit width in lab-z,
      or the horizontal scattering plane)
    * v = detector vert hieght (total slit width in lab-x,
      or the vertical scattering plane)

    Notes:
    ------
    The slit settings, defined wrt psic phi-frame
    Assume these are centered on the origin, then rotated
    by del and nu
    """
    # detector vectors, [x,y,z] in lab frame
    # note rotation of the vectors...
    dh = num.array([   0.,  0.,  0.5*h])
    dv = num.array([0.5*v,  0.,  0.   ]) 
    D  = calc_D(nu=nu,delta=delta)
    dh = num.dot(D,dh)
    dv = num.dot(D,dv)

    # corners of detector apperature 
    e =  dv + dh
    f =  dv - dh
    g = -dv - dh
    h = -dv + dh
    det = [e,f,g,h]

    return det

def sample_vectors(sample,angles={},gonio=None):
    """
    Parameters:
    -----------
    * sample = [[x,y,z],[x,y,z],[x,y,z],....]
      is a list of vectors that describe the shape of
      the sample.  They should be given in general lab
      frame coordinates.

    * angles = {'phi':0.,'chi':0.,'eta':0.,'mu':0.}
      are the instrument angles at which the sample
      vectors were determined.

    Notes:
    ------
    The lab frame coordinate systems is defined such that:
    x is vertical (pointing to the ceiling of the hutch, perpendicular to beam path)
    y is directed along the incident beam path
    z make the system right handed and lies in the horizontal scattering plane
    (i.e. z is parallel to the phi axis)

    The center (0,0,0) of the lab frame is the rotation center of the instrument.

    If the sample vectors are given at the flat phi and chi values and with
    the correct sample hieght (sample Z set so the sample surface is on the
    rotation center), then the z values of the sample vectors will be zero.
    If 2D vectors are passed we therefore assume these are [x,y,0].  If this
    is the case then make sure:
        angles = {'phi':flatphi,'chi':flatchi,'eta':0.,'mu':0.}

    Note that the sample_poly that is returned is a list of 3D vectors.
    If gonio == None these are defined in the lab phi frame.
    If a gonio instance is passed then they will be rotated to the m-frame

    The easiest way to determine the sample coordinate vectors is to take a picture
    of the sample with a camera mounted such that is looks directly down the omega
    axis and the gonio angles set at the sample flat phi and chi values and
    eta = mu = 0. Then find the sample rotation center and measure the position
    of each corner (in mm) with up being the +x direction, and downstream
    being the +y direction.  

    """
    if sample == None: return None
    if len(sample) < 3:
        print("Sample polygon must be 3 or more points")
        return None
    # If angles are provided then we need to compute the phi
    # frame vectors first. i.e. assume p's are defined in lab
    # frame at the given set of angles. therefore we need to 
    # unrotate to get back to phi frame vectors
    # Otherwise we assume the vectors passed are phi frame
    # (in that case they should be 3D vectors)
    if angles == None: angles = {}
    if len(angles) > 0:
        # calc sample rotation matrix
        Z = calc_Z(**angles)
        Zinv = num.linalg.inv(Z)
        polygon_phi = []
        # If p's have only two components then we assume they are 
        # xy pairs therefore we can add a third value of zero for z
        for p in sample:
            if len(p) == 2: p = [p[0],p[1],0.]
            p_phi = num.dot(Zinv,p)
            polygon_phi.append(p_phi)
    else:
        polygon_phi = sample

    # If gonio is passed then rotate the
    # vectors into the m-frame.  Otherwise
    # just return the phi frame vectors
    polygon = []
    if gonio != None:
        for p in polygon_phi:
            # If p's have only two components then we assume they are 
            # xy pairs therefore we can add a third value of zero for z
            if len(p) == 2: p = [p[0],p[1],0.]
            p_m = num.dot(gonio.Z,p)
            polygon.append(p_m)
    else:
        polygon = polygon_phi

    return polygon

##########################################################################
##########################################################################
def test1():
    # create a new psic instance
    psic = Psic(5.6,5.6,13,90,90,120,lam=1.3756)

    # diffr. angles:
    psic.set_angles(phi=65.78,chi=37.1005,eta=6.6400,
                    mu=0.0,nu=0.0,delta=46.9587)
    #calc tth, d and magnitude of Q
    print("\nh=",  psic.h)
    print("tth =", psic.pangles['tth'])
    print("tth =", psic.lattice.tth(psic.h,lam=psic.lam))

    # reference vector/surface normal in [h,k,l]
    # n = [0, 0, 1]
    n = [-0.0348357, -0.00243595, 1]
    psic.set_n(n)
    #calc miscut
    print("\nmiscut=",psic.lattice.angle([0,0,1],n,recip=True))

    # test surf norm
    n_phi = num.dot(psic.UB,n)
    n_phi = n_phi / num.fabs(cartesian_mag(n_phi))
    print("\nnphi: ", n_phi)

    print("\nsigma_az=",psic.pangles['sigma_az'])

    print("\ntau_az=",psic.pangles['tau_az'])

    psic.calc_n(-psic.pangles['sigma_az'],-psic.pangles['tau_az'])
    print("calc n from sigma and tau az: ",psic.n) 
    
    #calc miscut
    print("\nmiscut=",psic.lattice.angle([0,0,1],psic.n,recip=True))
    
    return psic

##########################################################################
if __name__ == "__main__":
    """
    test 
    """
    psic = test1()
    


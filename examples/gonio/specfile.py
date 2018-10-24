"""
Read spec files

Authors/Modifications:
----------------------
* Matt Newville (newville@cars.uchicago.edu)
* Tom Trainor (tptrainor@alaska.edu)

"""
#######################################################################
import numpy as num
import os, copy

#######################################################################
def spec_scan(spec,sc_num,geo='PSIC_APS_S13'):
    """
    Return a ScanData instance from a specfile / scan number
    
    * spec can be a specfile instance or string file name
    * sc_num is the scan number
    * geo is a geometry label that helps the reader parse gonio
      angles depending on the particular geometry
      used for data collection.  

    Notes:
    ------
    The vector A is appended to the state info.  This
    is the set of essential gonio angles at the start of the scan
    These are defined based on the 'geo' label and are geometry/
    beamline specific

    Example:
    --------
    >>s = spec_scan('datafile.spc',12)
    >>plot(s['phi'],s['bicron'])
    
    """
    # Define to make sure we get things sorted correclty. 
    if geo=='PSIC_APS_S13':
        POSITIONER_KEYS = ['phi','chi','eta','mu','nu','del']
    else:
        POSITIONER_KEYS =[]
    
    # get the spec scan data
    if type(spec) == str:
        spec = SpecFile(spec)
        if spec == None: return None
    d = spec.scan_dict(sc_num)
    # parse positioner and scaler vals
    # note if a positioner or scaler was listed in the data
    # array then we append the array to the positioners/scalers.
    # otherwise the positioner value will be what was given
    # in d['P'] which should be a single value
    scalers = {}
    positioners = copy.copy(d['P'])
    for key in d['data'].keys():
        if key in positioners.keys():
            positioners[key] = num.array(d['data'][key])
        elif key in POSITIONER_KEYS:
            positioners[key] = num.array(d['data'][key])
        else:
            scalers[key] = num.array(d['data'][key])
    name  = d['file'] + ' Scan ' + str(int(sc_num))
    dims  = d['nrow']
    paxis = d['labels'][0]
    pdet  = d['labels'][-1]
    # Grab state info.  
    try:
        A = []
        if geo=='PSIC_APS_S13':
            A.append(d['P'].get('TwoTheta'))
            A.append(d['P'].get('theta'))
            A.append(d['P'].get('chi'))
            A.append(d['P'].get('phi'))
            A.append(d['P'].get('Nu'))
            A.append(d['P'].get('Psi'))
            A.append(d['P'].get('Omega'))
            A.append(d['P'].get('Kappa'))
            A.append(d['P'].get('Phi'))
    except:
        A = []
    state = {'G':d.get('G'),'Q':d.get('Q'),'A':A,
             'ATTEN':d.get('ATTEN'), 'ENERGY':d.get('ENERGY')}
    # create scan data object
    data = ScanData(name=name,dims=[dims],scalers=scalers,
                    positioners=positioners,primary_axis=[paxis],
                    primary_det=pdet,state=state)
    return data

#######################################################################
class ScanData:
    """
    Simple container that holds scan data.

    Attributes:
    -----------
    * name = scan name
    * dims ->[npts1] for a 1-d scan
             [npts1,npts2] for a 2-d scan
             where npts1 is the outer loop
             and npts2 is for the inner loop
    * scalers = {'io':[],'i1':[]}  -> each array has dims = dims
    * positioners  = {'E':[]}      -> these may be single values,
                                      or same size as scalers
    * primary_axis = ['E']         -> if multi-dim, list of outer, inner
    * primary_det  = 'i1'          -> use for default plotting
    * state = {}                   -> dictionary of additional state information
    """
    def __init__(self,name='',dims=[],scalers={},positioners={},
                 primary_axis=[],primary_det=None,state={}):
        """
        Initialize

        Parameters:
        -----------
        * name=''
        * dims=[]
        * scalers={}
        * positioners={}
        * primary_axis=[]
        * primary_det=None
        * state={}
        """
        self.name         = name
        self.dims         = dims
        self.primary_axis = primary_axis
        self.primary_det  = primary_det
        self.scalers      = scalers
        self.positioners  = positioners
        self.state        = state

    def __getitem__(self,arg):
        """
        Get items.  Note that order matters:
        spectra, images, scalers, positioners, state
        """
        for key in self.scalers.keys():
            if key == arg:
                return self.scalers[key]
        for key in self.positioners.keys():
            if key == arg:
                return self.positioners[key]
        for key in self.state.keys():
            if key == arg:
                return self.state[key]
        return None

    def __repr__(self):
        lout = "* Scan Name = %s\n" % self.name
        lout = lout + "* Scan Dimensions = %s\n" % str(self.dims)

        lout = lout + "* Scalers:\n"
        ct = 0
        for sc in self.scalers.keys():
            lout = lout + "%12s" % sc
            ct = ct + 1
            if ct > 7:
                lout = lout + '\n'
                ct = 0
        lout = lout + "\n* Positioners:\n"
        ct = 0
        for p in self.positioners.keys():
            lout = lout + "%12s" % p
            ct = ct + 1
            if ct > 7:
                lout = lout + '\n'
                ct = 0
        if len(self.state) > 0:
            lout = lout + "\n* Additional State Variables:\n"  
            ct = 0
            for p in self.state.keys():
                lout = lout + "%12s" % p
                ct = ct + 1
                if ct > 7:
                    lout = lout + '\n'
                    ct = 0
        lout = lout + "\n* Primary scan axis = %s\n"  % str(self.primary_axis)
        lout = lout + "* Primary detector  = %s\n"  % self.primary_det
        return lout

    def get_scaler(self,label=None):
        """
        return scaler
        """
        if label == None:
            label = self.primary_det[0]
        return self.scalers.get(label)

    def get_positioner(self,label=None):
        """
        return positioner
        """
        if label == None:
            label = self.primary_axis[0]
        return self.positioners.get(label)

#######################################################################
class SpecFile:
    """
    A spec file
    """
    def __init__(self, fname):
        """
        Initialize

        Parameters:
        -----------
        * fname is the specfile name (including full path)
        """
        self.path, self.fname = os.path.split(fname)
        self.max_scan = 0
        self.min_scan = 0
        self._mtime    = 0
        self._lines    = []
        self._summary  = []
        self._ok       = False
        self.read()

    def __repr__(self):
        """ display """
        self.read()
        lout = "Spec file: %s" % self.fname
        lout = "%s\nPath: %s" % (lout, os.path.join(self.path))
        lout = "%s\nFirst scan number: %i" % (lout,self.min_scan)
        lout = "%s\nLast scan number:  %i" % (lout,self.max_scan)
        lout = "%s\nLast scan: %s" % (lout, self._summary[self.max_scan-1]['date'])
        return lout

    def read(self):
        """
        Read the specfile

        This will re-read the file if its time stamp
        has changed since the last read
        """
        try:
            fname = os.path.join(self.path, self.fname)
            if os.path.getmtime(fname) != self._mtime:
                #print("Reading spec file %s" % fname)
                f  = open(fname)
                self._mtime = os.path.getmtime(fname)
                self._lines = f.readlines()
                f.close()
                self._summarize()
                self._ok = True
        except IOError:
            print('**Error reading file ', fname)
            self._ok = False

    def _summarize(self):
        """
        summarize
        """
        lineno = 0
        (mnames,cmnd,date,xtime,Gvals,q,Pvals,atten,energy,lab) = (None,None,None,None,None,None,None,None,None,None)
        (index, ncols, n_sline) = (0,0,0)
        for i in self._lines:
            lineno = lineno + 1
            i  = i[:-1]
            # get motor names: they should be at the top of the file
            # but they can be reset anywhere in the file
            if (i[0:2] == '#O'):
                if i[2] == '0': mnames = ''
                mnames = mnames + i[3:]
            # get scan number
            elif (i[0:3] == '#S '):
                v     = i[3:].split()
                index = int(v[0])
                cmnd  = i[4+len(v[0]):]
                n_sline= lineno
            elif (i[0:3] == '#D '):
                date = i[3:]
            elif (i[0:3] == '#T '):
                xtime = i[3:]
            elif (i[0:2] == '#G'):
                if i[2] == '0': Gvals = ''
                Gvals = Gvals + i[3:]
            elif (i[0:3] == '#Q '):
                q = i[3:]
            elif (i[0:2] == '#P'):
                if i[2] == '0': Pvals = ''
                Pvals = Pvals + i[3:]
            elif (i[0:3] == '#N '):
                ncols = int(i[3:])
            elif (i[0:3] == '#AT'):
                atten = i[6:]
            elif (i[0:3] == '#EN'):
                energy = i[8:]
            elif (i[0:3] == '#L '):
                lab = i[3:]
                ## count how many lines of 'data' we have
                ## and see if the scan was aborted
                xx = self._lines[lineno:]
                nl_dat = 0
                aborted = False
                for ii in xx:
                    if (ii[0:3] == '#S '):
                        break
                    elif (ii[0:1] ==  '#'):
                        if ii.find('aborted') > -1:
                            aborted = True
                    elif (len(ii)  > 3):
                        nl_dat = nl_dat + 1
                ## append all the info...
                self._summary.append({'index':index,
                                     'nl_start':n_sline,
                                     'cmd':cmnd,
                                     'date':date,
                                     'time':xtime,
                                     'G':Gvals,
                                     'Q':q,
                                     'mot_names':mnames,
                                     'P':Pvals,
                                     'ncols':ncols,
                                     'labels':lab,
                                     'atten':atten,
                                     'energy':energy,
                                     'lineno':lineno,
                                     'nl_dat':nl_dat,
                                     'aborted':aborted})
                (cmnd,date,xtime,Gvals,q,Pvals,atten,energy,lab) = (None,None,None,None,None,None,None,None,None)
                (index, ncols, n_sline) = (0,0,0)

        self.min_scan = self._summary[0]['index']
        self.max_scan = self._summary[0]['index']
        for i in self._summary:
            k = i['index']
            if (k > self.max_scan): self.max_scan = k
            if (k < self.min_scan): self.min_scan = k

    def scan_min(self):
        """
        get the minimum scan number
        """
        self.read()
        return self.min_scan

    def scan_max(self):
        """
        get the max scan number
        """
        self.read()
        return self.max_scan

    def nscans(self):
        """
        get the number of scans
        """
        self.read()
        return len(self._summary)

    def _check_range(self,i):
        """
        check if scan number is in range
        """
        self.read()
        j = True
        if ((i > self.max_scan) or (i < self.min_scan)): j = False
        return j

    #def scan_info(self):
    #    self.read()
    #    return self._summary

    #def get_summary(self,sc_num):
    def scan_info(self,sc_num):
        """
        return the scan info in a dictionary
        """
        self.read()
        for s in self._summary:
            if (sc_num == s['index']):
                return s
        return None
    
    def scan_data(self, sc_num):
        """
        return just the column data from the scan 
        """
        self.read()
        s = self.scan_info(sc_num)
        if (s == None): return None
        dat = []
        nl  = s['lineno']
        for i in (self._lines[nl:]):
            if (i[0:3] == '#S '):
                break
            elif (i[0:1] ==  '#'):
                pass
            elif (len(i)  > 3):
                q = i.split()
                #dat.append(map(float,q))
                dat.append([float(f) for f in q])
        return dat

    def scan_dict(self, sc_num):
        """
        return scan information and data in a dictionary 
        """
        self.read()
        sc_dict = {'file':self.fname,
                   'index':sc_num,
                   'cmd':'',
                   'date':'',
                   'G':[],
                   'Q':[],
                   'P':{},
                   'ATTEN':[],
                   'ENERGY':[],
                   'labels':[],
                   'ncol':0,
                   'nrow':0,
                   'data':{}
                   }
        s = self.scan_info(sc_num)
        if (s == None): return sc_dict
        dat = self.scan_data(sc_num)
        # parse the various data into the dict
        sc_dict['cmd']  = s['cmd']
        sc_dict['date'] = s['date']
        if s['G'] != None: sc_dict['G']    = [float(f) for f in s['G'].split()]
        if s['Q'] != None: sc_dict['Q']    = [float(f) for f in s['Q'].split()]
        if s['atten'] != None: sc_dict['ATTEN'] = [int(f) for f in s['atten'].split()]
        if s['energy'] != None: sc_dict['ENERGY'] = [float(f) for f in s['energy'].split()]
        # get the motor positions
        p_dict = {}
        m_names = s['mot_names'].split()
        p_vals  = s['P'].split()
        if len(m_names) != len(p_vals):
            print("Mismatch in Motor Names and Motor Values")
        else:
            for j in range(len(m_names)):
                p_dict.update({m_names[j]:float(p_vals[j])})
        sc_dict['P'] = p_dict
        #
        lbls = s['labels'].split()
        sc_dict['labels'] = lbls
        ncol = len(lbls)
        nrow = len(dat)
        sc_dict['ncol']   = ncol
        sc_dict['nrow']   = nrow
        # data
        data_dict = {}
        for j in range(ncol):
            xx = []
            for k in range(nrow):
                xx.append(dat[k][j])
            data_dict.update({lbls[j]:xx})
        sc_dict['data'] = data_dict
        # all done
        return sc_dict

    def list_scans(self):
        """
        return a list of scans in the file 
        """
        self.read()
        sc_list = [] 
        for s in self._summary:
            line = "%s:%4.4i  (%s)\n   %s" % (self.fname,
                                              s['index'],
                                              s['date'],
                                              s['cmd'])
            sc_list.append(line)
        return sc_list



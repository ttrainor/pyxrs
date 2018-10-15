#!/usr/bin/env python3
"""
xtal structure calcs
"""
##########################################################
import os, sys, traceback
import string
import cgi
import cgitb

from _xtal import *

from pyxrs.xtal.unitcell import cif_to_uc, uc_to_cif, uc_to_xyz 
from pyxrs.xtal.surface import SurfaceCell, surface_to_xyz
from pyxrs.xtal.coord import coord_calcs
##########################################################

xtal_dir = '../xtal'
log_dir = xtal_dir + '/log'
tmp_dir = xtal_dir + '/tmp'
cgitb.enable(logdir=log_dir, format="text")
emsg = ""
test = True

def error_exit(emsg):
    emsg = ERR_HTML % emsg
    print(emsg)
    sys.exit()

###########################################################

## get fields from form 
form = cgi.FieldStorage()

## write the cif file to tmp dir
cif_file = form['cif_file']
if cif_file.filename:
    cif_fname = os.path.basename(cif_file.filename)
    fn = xtal_dir + '/tmp/' + cif_fname
    f = open(fn, 'wb')
    f.write(cif_file.file.read())
    f.close()
    emsg = emsg + "The file '" + fn + "' was uploaded </br>\n"
else:
    emsg = emsg + "Error no cif file</br>\n"
    error_exit(emsg)

## parameter file
if test == True:
    param_file = open(tmp_dir + '/params.temp','w')
    param_file.write("filename = %s\n" % cif_fname)
    for key in form.keys():
        if key != "cif_file":
            t = "key   %s = %s\n" % (key, form[key].value)
            param_file.write(t)
    param_file.write("\n\n")
    param_file.close()
    emsg = emsg + "Param file written</br>\n"

## parse out parameters
try:
    na            = int(form['na'].value)
    nb            = int(form['nb'].value)
    nc            = int(form['nc'].value)
    H             = form['H'].value.strip()
    K             = form['K'].value.strip()
    L             = form['L'].value.strip()
    nd            = form['nd'].value
    term          = form['term'].value
    bulk_trns     = form['bulk_trns'].value.strip()
    show_bulk_hkl = eval(form['show_bulk_hkl'].value.strip())
    long_fmt      = eval(form['long_fmt'].value)
    rmax          = float(form['rmax'].value)
    labels        = form['labels'].value.strip()
except:
    emsg = emsg + "Error getting form parameters</br>\n"
    emsg = emsg + traceback.format_exc()
    error_exit(emsg)

## do surface calcs?
try:
    if H and K and L:
        H = int(H); K = int(K); L = int(L)
        HKL = [H, K, L]
        nd = int(nd)
        term = int(term)
        if bulk_trns:
            bulk_trns = eval(bulk_trns)
        else:
            bulk_trns = None
    else:
        H = K = L = 0
        HKL = None
except:
    emsg = emsg + "Error getting surface parameters</br>\n"
    emsg = emsg + traceback.format_exc()
    error_exit(emsg)

## coord calc params
try:
    if rmax <= 0:
        rmax = 0
    elif rmax > 10.:
        rmax = 6
    if len(labels) == 0:
        labels = None
    else:
        if labels[0] == '[':
            labels = eval(labels)
        else:
            labels = [lbl.strip() for lbl in labels.split(',')]
except:
    emsg = emsg + "Error getting cood calc parameters</br>\n"
    emsg = emsg + traceback.format_exc()
    error_exit(emsg)

## get unit cell
try:
    uc = cif_to_uc(fn)
except:
    emsg = emsg + "Error reading the cif file</br>\n"
    emsg = emsg + traceback.format_exc()
    error_exit(msg)

## get structure info
try:
    xtal_out = ""
    coord_out = ""
    xyz_file = tmp_dir + "/structure.xyz"
    frac_file = tmp_dir + "/structure.frac"
    bulk_cif = tmp_dir + "/bulk.cif"
    coord_file = tmp_dir + "/coord.out" 
    if HKL is None:
        # bulk
        xtal_out = uc._write(long_fmt=long_fmt)
        xtal_out = xtal_out + "\nP1 cell\n"
        atom_list = uc.atom_list()
        xtal_out = xtal_out + atom_list._write(long_fmt=long_fmt)
        # output xyz file
        uc_to_xyz(uc,fname=xyz_file,cartesian=True,na=na,nb=nb,nc=nc,long_fmt=long_fmt)
        uc_to_xyz(uc,fname=frac_file,cartesian=False,na=na,nb=nb,nc=nc,long_fmt=long_fmt)
        coord = coord_calcs(uc,rmax=rmax,labels=labels)
    else:
        # surf
        surf = SurfaceCell(uc,hkl=HKL,nd=nd,term=term,bulk_trns=bulk_trns)
        xtal_out = surf._write(long_fmt=long_fmt)
        surface_to_xyz(surf,fname=xyz_file,cartesian=True,na=na,nb=nb,nbulk=nc,long_fmt=long_fmt)
        surface_to_xyz(surf,fname=frac_file,cartesian=False,na=na,nb=nb,nbulk=nc,long_fmt=long_fmt)
        if show_bulk_hkl == True:
            uc_to_cif(uc,bulk_cif,p1_list=True,na=3,nb=3,nc=3)
        coord = coord_calcs(surf,rmax=rmax,labels=labels)
    #
    coord_out = coord._write(long_fmt=long_fmt)
    coord.write(fname=coord_file)
    emsg = emsg + "Calcs complete</br>\n"
except:
    emsg = emsg + "Error doing calculations </br>\n"
    emsg = emsg + traceback.format_exc()
    error_exit(msg)

## format output page
try:
    subs = {'FNAME':cif_fname,
            'XTAL':xtal_out,
            'COORD':coord_out,
            'H':str(H), 'K':str(K), 'L':str(L)}
    #
    html = START_HEAD
    html = html + JMOL_1_SCRIPT
    if HKL is not None and show_bulk_hkl == True:
        html = html + string.Template(JMOL_2_SCRIPT).substitute(subs)
    html = html + END_HEAD
    html = html + string.Template(START_HTML).substitute(subs)
    html = html + string.Template(XTAL_2_OUT).substitute(subs)
    html = html + string.Template(COORD_OUT).substitute(subs)
    html = html + STRUCT_1_OUT
    if HKL is not None and show_bulk_hkl == True:
        html = html + string.Template(STRUCT_2_OUT).substitute(subs)
    html = html + END
    #
    f = open(tmp_dir + "/output.html",'w')
    f.write(html)
    f.close()
    print("Location: %s/output.html\n\n" % tmp_dir)
except:
    emsg = emsg + "Error generating output page</br>\n"
    emsg = emsg + traceback.format_exc()
    error_exit(emsg)

## if get here display error page
error_exit(emsg)




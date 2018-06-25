# coding: utf-8

"""
This provides classes and functions for dealing with FEFF i/o
"""

__author__ = "Mehmet Topsakal"
__email__ = "metokal@gmail.com"
__status__ = "Development"
__date__ = "March 20, 2018"



import numpy as np
from numpy import linalg as LA
import os,sys

from pymatgen.core.periodic_table import Element

import numpy as np
import os,sys,shutil,subprocess,pickle
from os.path import join

from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import pymatgen as mg

from scipy import interpolate

from . import t4iss_defaults
from .core import mXANES



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# this reads POSCAR
def readposcar(file):
    
    lattice = []; positions_dir = []; labels = []
    
    with open(file, mode='r') as f:
        
        f.readline()
        f.readline()
        
        for i in range(3):
            l =  f.readline()
            l = l.split()
            l = [float(x) for x in l]
            lattice.append(l)
        lattice = np.array(lattice)
        
        f.readline()
        natoms = f.readline().split(); natoms = [int(x) for x in natoms]; natoms = np.array(natoms)
        mode = f.readline().split()
        
        labels = []
        for i in range(sum(natoms)):
            p = f.readline()
            l = p.split()[3]
            p = p.split()[0:3]
            p = [float(x) for x in p]
            positions_dir.append(p)
            labels.append(l)
                              
        positions = []; pnew = [] 
        # TODO: this for loop can be vectorized
        for p in range(len(positions_dir)):
            pnew.append(positions_dir[p][0]*lattice[0][0]+positions_dir[p][1]*lattice[1][0] + positions_dir[p][2]*lattice[2][0] ); 
            pnew.append(positions_dir[p][0]*lattice[0][1]+positions_dir[p][1]*lattice[1][1] + positions_dir[p][2]*lattice[2][1] );
            pnew.append(positions_dir[p][0]*lattice[0][2]+positions_dir[p][1]*lattice[1][2] + positions_dir[p][2]*lattice[2][2] );
            positions.append(pnew); pnew = []   
        positions = np.array(positions) 
        positions = positions.reshape(len(positions),3,order='F').copy()
        

    return labels, natoms, lattice, positions, positions_dir


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
def make_333_supercell(labels, natoms, lattice, positions):
    
    supercell = []; pnew=[]; p=positions; l=lattice; ls = []
    
    # TODO: this for loop can be vectorized
    for s in range(len(p)):

        pnew.append(p[s][0]); pnew.append(p[s][1]); pnew.append(p[s][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[]

        pnew.append(p[s][0]+l[0][0]); pnew.append(p[s][1]+l[0][1]); pnew.append(p[s][2]+l[0][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # +x
        pnew.append(p[s][0]-l[0][0]); pnew.append(p[s][1]-l[0][1]); pnew.append(p[s][2]-l[0][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # -x
        pnew.append(p[s][0]+l[1][0]); pnew.append(p[s][1]+l[1][1]); pnew.append(p[s][2]+l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # +y
        pnew.append(p[s][0]-l[1][0]); pnew.append(p[s][1]-l[1][1]); pnew.append(p[s][2]-l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # -y
        pnew.append(p[s][0]+l[2][0]); pnew.append(p[s][1]+l[2][1]); pnew.append(p[s][2]+l[2][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # +z
        pnew.append(p[s][0]-l[2][0]); pnew.append(p[s][1]-l[2][1]); pnew.append(p[s][2]-l[2][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # -z

        pnew.append(p[s][0]+l[0][0]+l[1][0]); pnew.append(p[s][1]+l[0][1]+l[1][1]); pnew.append(p[s][2]+l[0][2]+l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # +x+y
        pnew.append(p[s][0]+l[0][0]-l[1][0]); pnew.append(p[s][1]+l[0][1]-l[1][1]); pnew.append(p[s][2]+l[0][2]-l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # +x-y
        pnew.append(p[s][0]-l[0][0]+l[1][0]); pnew.append(p[s][1]-l[0][1]+l[1][1]); pnew.append(p[s][2]-l[0][2]+l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # -x+y
        pnew.append(p[s][0]-l[0][0]-l[1][0]); pnew.append(p[s][1]-l[0][1]-l[1][1]); pnew.append(p[s][2]-l[0][2]-l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # -x-y
        pnew.append(p[s][0]+l[2][0]+l[0][0]); pnew.append(p[s][1]+l[2][1]+l[0][1]); pnew.append(p[s][2]+l[2][2]+l[0][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # +z+x
        pnew.append(p[s][0]+l[2][0]-l[0][0]); pnew.append(p[s][1]+l[2][1]-l[0][1]); pnew.append(p[s][2]+l[2][2]-l[0][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # +z-x
        pnew.append(p[s][0]-l[2][0]+l[0][0]); pnew.append(p[s][1]-l[2][1]+l[0][1]); pnew.append(p[s][2]-l[2][2]+l[0][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # -z+x
        pnew.append(p[s][0]-l[2][0]-l[0][0]); pnew.append(p[s][1]-l[2][1]-l[0][1]); pnew.append(p[s][2]-l[2][2]-l[0][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # -z-x
        pnew.append(p[s][0]+l[2][0]+l[1][0]); pnew.append(p[s][1]+l[2][1]+l[1][1]); pnew.append(p[s][2]+l[2][2]+l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # +z+y
        pnew.append(p[s][0]+l[2][0]-l[1][0]); pnew.append(p[s][1]+l[2][1]-l[1][1]); pnew.append(p[s][2]+l[2][2]-l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # +z-y
        pnew.append(p[s][0]-l[2][0]+l[1][0]); pnew.append(p[s][1]-l[2][1]+l[1][1]); pnew.append(p[s][2]-l[2][2]+l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # -z+y
        pnew.append(p[s][0]-l[2][0]-l[1][0]); pnew.append(p[s][1]-l[2][1]-l[1][1]); pnew.append(p[s][2]-l[2][2]-l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # -z-y

        pnew.append(p[s][0]+l[2][0]+l[0][0]+l[1][0]); pnew.append(p[s][1]+l[2][1]+l[0][1]+l[1][1]); pnew.append(p[s][2]+l[2][2]+l[0][2]+l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # +z+x+y
        pnew.append(p[s][0]+l[2][0]+l[0][0]-l[1][0]); pnew.append(p[s][1]+l[2][1]+l[0][1]-l[1][1]); pnew.append(p[s][2]+l[2][2]+l[0][2]-l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # +z+x-y
        pnew.append(p[s][0]+l[2][0]-l[0][0]-l[1][0]); pnew.append(p[s][1]+l[2][1]-l[0][1]-l[1][1]); pnew.append(p[s][2]+l[2][2]-l[0][2]-l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # +z-x-y
        pnew.append(p[s][0]+l[2][0]-l[0][0]+l[1][0]); pnew.append(p[s][1]+l[2][1]-l[0][1]+l[1][1]); pnew.append(p[s][2]+l[2][2]-l[0][2]+l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # +z-x+y
        pnew.append(p[s][0]-l[2][0]+l[0][0]+l[1][0]); pnew.append(p[s][1]-l[2][1]+l[0][1]+l[1][1]); pnew.append(p[s][2]-l[2][2]+l[0][2]+l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # -z+x+y
        pnew.append(p[s][0]-l[2][0]+l[0][0]-l[1][0]); pnew.append(p[s][1]-l[2][1]+l[0][1]-l[1][1]); pnew.append(p[s][2]-l[2][2]+l[0][2]-l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # -z+x-y
        pnew.append(p[s][0]-l[2][0]-l[0][0]-l[1][0]); pnew.append(p[s][1]-l[2][1]-l[0][1]-l[1][1]); pnew.append(p[s][2]-l[2][2]-l[0][2]-l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # -z-x-y
        pnew.append(p[s][0]-l[2][0]-l[0][0]+l[1][0]); pnew.append(p[s][1]-l[2][1]-l[0][1]+l[1][1]); pnew.append(p[s][2]-l[2][2]-l[0][2]+l[1][2]); ls.append(labels[s]); supercell.append(pnew); pnew=[] # -z-x+y

    positions_supercell = np.array(supercell) ;
    positions_supercell = positions_supercell.reshape(len(p)*27,3,order='F').copy()
    
    lattice_supercell = lattice*3
    natoms_supercell = natoms*3
    labels_supercell = ls
    
    return labels_supercell, natoms_supercell, lattice_supercell, positions_supercell



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# this writes feff.inp
def feff_write_inp(structure,dmax=10.1,rFMS=8.1,rSCF=6.1,corehole='RPA',destination=None): 
    
    poscar_path = structure['poscar']
    cai = structure['isite']
    
    [labels, natoms, lattice, positions, positions_dir] = readposcar(poscar_path)
    
    labels_short = []
    for i in labels:
        if i not in labels_short:
            labels_short.append(i)
    
    shifts = positions[cai] 
    
    # move selected to origin
    positions_shifted = []
    for i in range(len(positions)):
        p = [positions[i][0]-shifts[0],positions[i][1]-shifts[1],positions[i][2]-shifts[2]]
        positions_shifted.append(p)
    positions_shifted  

    labels_supercell, natoms_supercell, lattice_supercell, positions_supercell = make_333_supercell(labels, natoms, lattice, positions_shifted)
        
    ds = []
    for i in positions_supercell:
        n = LA.norm(i)
        ds.append(n)    
    #   
    if len(positions_shifted) < 50:
        labels_supercell, natoms_supercell, lattice_supercell, positions_supercell = make_333_supercell(labels_supercell, natoms_supercell, lattice_supercell, positions_supercell)
   
    atoms = []
    ds = []
    for i,s in enumerate(positions_supercell):
        n = LA.norm(s)
        if n <= dmax:
            atoms.append([s[0],s[1],s[2],labels_supercell[i],n])
                     
    def getKey4(item): return item[4]
    
    atoms = sorted(atoms, key=getKey4)
    
    for i,a in enumerate(atoms):
        s = a[3]
        ind = labels_short.index(s)
        atoms[i][3] = ind
 
    if destination is None:
        f=open('feff.inp',"w+")  
    else:
        f=open(os.path.join(destination,'feff.inp'),"w+")
        
    f.write("""TITLE             
                                                   
EDGE      K
S02       0
COREHOLE  %(corehole)s                                  
CONTROL   1 1 1 1 1 1                               
                                                   
XANES 4 0.05 0.1
#ABSOLUTE                               
                                                   
FMS       %(rFMS)2.1f                                
EXCHANGE  0                                 
SCF       %(rSCF)2.1f  0 30 0.2 3                                   
RPATH     -1 
    
""" % vars())
    
    Sca =  labels[cai]      
    el = Element(Sca); d = el.data; Zca = d['Atomic no']    
    
    f.write("""POTENTIALS
*   ipot   Z      element   l_scmt   l_fms   stoichiometry
    0      %(Zca)i     %(Sca)s        -1       -1      0.001 """ % vars())     
    
    for i in range(len(labels_short)):
        n = (i+1)
        s = labels_short[i]
        el = Element(s); d = el.data; z = d['Atomic no'] 
        st = natoms[i]
        f.write("""
    %(n)i      %(z)i     %(s)s        -1       -1      %(st)i """ % vars())         
        
    f.write("\n \n")  
    
    f.write("ATOMS\n")
    f.write("       0.000000     0.000000     0.000000     0    0.0\n")    
    for i in atoms[1:]:      
        f.write('  %13.6f%13.6f%13.6f   %3d   %6.3f\n' % (i[0], i[1], i[2], i[3]+1, i[4]))
    f.write("END\n")    
    f.close()       
    
    for_dist_plot = []
    for i in atoms:
        for_dist_plot.append([i[4],i[3]])
    
    return 








#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# this reads FEFF xmu.dat
def feff_read_xmu(xmudat='xmu.dat'): 
    xmu = np.loadtxt(xmudat, unpack=True, comments='#', usecols=(0,3), skiprows=0)
    return xmu









          
            
 
            
            

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
class xanes_collection_feff:

    def __init__(self, absorption_specie, xanes_data=None):

        
        if xanes_data is None:
            self.xanes_data = t4iss_defaults['t4iss_xanes_data']  
        else:
            self.xanes_data = xanes_data     
        
        self.absorption_specie = absorption_specie        
        self.collection = []
        self.ids = []
        
        self.ave_spectra = []
        self.site_spectra = []
        
        
    def build(self, mpid_list, skip_missing=True,download_missing=False,mpr=None):
        
        if download_missing:
            if not mpr:
                download_missing=False
                print('mpr key is needed for download_missing=True')
                        
        collection0 = []
        for i in mpid_list:
            read = read_xanes_feff(mpid=i,absorption_specie=self.absorption_specie,skip_missing=skip_missing,
                                   download_missing=download_missing,mpr=mpr,xanes_data=self.xanes_data)
            collection0.append(read)
        collection = []
        for i in collection0:
            if i[0]:
                collection.append(i)
        self.collection = collection
        ids = []
        for i in self.collection:
            ids.append(i[0].xanesid)
        self.ids = ids
        
        ave_spectra = []
        for a in self.collection:
            ave_spectra.append(a[0])
        self.ave_spectra = ave_spectra
            
        site_spectra = []
        for s in self.collection:
            for i in s[1]:
                site_spectra.append(i)
        self.site_spectra = site_spectra


                














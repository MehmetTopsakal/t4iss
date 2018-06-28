#!/usr/bin/env python3

# -*- coding: iso-8859-1 -*-
# coding: utf-8

# pylint: disable=C0103
# pylint: disable=C0200
# pylint: disable=R0912
# pylint: disable=R0902
# pylint: disable=R0903
# pylint: disable=R0913
# pylint: disable=R0914
# pylint: disable=R0915

"""This provides core functions needed by t4iss."""

__author__ = "Mehmet Topsakal"
__email__ = "metokal@gmail.com"
__status__ = "Development"
__date__ = "March 20, 2018"

import os
# from os.path import join

# import sys
# import shutil
# import subprocess
import pickle

import numpy as np

from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import pymatgen as mg

from mpl_toolkits.axes_grid1.anchored_artists import AnchoredDrawingArea
from matplotlib.patches import Circle
from matplotlib.offsetbox import AnchoredText

# from scipy.ndimage.filters import gaussian_filter1d
from scipy import interpolate
from scipy import signal
from scipy.stats import rankdata

# from t4iss.defaults import t4iss_defaults
from t4iss.PT import get_c


class mXANES:
    # TODO: docstring

    def __init__(self, data=None, data_loadfrom=None, srange=None,
                 structure=None, vcncutoff=5.0, ca=None, Eonset=None,
                 xanesid=None, source=None, edge=None, multiplicity=1):

        super(mXANES, self).__init__()

        # these are class-variables that were prior not defined in __init__
        self.env = None
        self.rank = None

        if data is None:
            if data_loadfrom:
                data = np.loadtxt(data_loadfrom, unpack=True, comments='#',
                                  skiprows=0)
                self.E0 = np.array(data[0])
                self.I0 = np.array(data[1])
            else:
                self.E0 = np.array([])
                self.I0 = np.array([])
        else:
            self.E0 = np.array(data[0])
            self.I0 = np.array(data[1])

        if srange:
            sel = (self.E0 >= srange[0]) & (self.E0 <= srange[1])
            self.E0 = self.E0[sel]
            self.I0 = self.I0[sel]

        # Energy offset
        if Eonset is None:
            self.Eonset = self.E0[0]
        else:
            self.Eonset = Eonset

        # XANES id
        if xanesid is None:
            self.xanesid = 'not_set'
        else:
            self.xanesid = xanesid

        # XANES edge
        if edge is None:
            self.edge = 'not_set'
        else:
            self.edge = edge

        # source
        if source is None:
            self.source = 'not_set'
        else:
            self.source = source

        # source
        if ca is None:
            self.ca = 'not_set'
        else:
            self.ca = ca

        # structure
        if structure:
            try:
                self.ca = structure[0][structure[1]].species_string
                nnfinder = VoronoiNN(cutoff=vcncutoff, allow_pathological=True)
                vcn = nnfinder.get_cn(structure[0], structure[1],
                                      use_weights=True)
                self.vcn = vcn
                self.structure = [structure[0], structure[1], vcn]
            except Exception as exc:
                print("  Failed to get nearest neighbors.")
                self.structure = [structure[0], structure[1], []]
                self.vcn = []
        else:
            self.structure = [[], [], []]
            self.vcn = []

        self.peaks = []
        self.multiplicity = multiplicity
        
        self.E = None
        self.I = None
        
                                            
    def Interpolate(self, iterprange, stepsize=0.1):
        # TODO: docstring

        # left padding
        if self.E[0] > iterprange[0]:
            npts = int((self.E[0]-iterprange[0])/stepsize)+1
            x_patch = np.linspace(iterprange[0], self.E[0]-stepsize, npts)
            y_patch = np.empty(len(x_patch))
            y_patch.fill(self.I[0])
            self.E = np.concatenate((x_patch, self.E.T), axis=0)
            self.I = np.concatenate((y_patch, self.I.T), axis=0)

        # right padding
        if self.E[-1] < iterprange[1]:
            npts = int((iterprange[1]-self.E[-1])/stepsize)+2
            x_patch = np.linspace(self.E[-1], iterprange[1], npts)
            y_patch = np.empty(len(x_patch))
            y_patch.fill(self.I[-1])
            self.E = np.concatenate((self.E.T, x_patch), axis=0)
            self.I = np.concatenate((self.I.T, y_patch), axis=0)

        f = interpolate.interp1d(self.E, self.I, kind='linear')
        self.E = np.linspace(iterprange[0], iterprange[1],
                             int((iterprange[1]-iterprange[0])/stepsize)+1)
        self.I = f(self.E) 

    def FindPeaks(self, xin=None, yin=None, srangep=None):
        # TODO: docstring

        if (xin is None) or (yin is None):
            xsearch = self.E0
            ysearch = self.I0
        else:
            xsearch = xin
            ysearch = yin                       
        if srangep:
            sel = (xsearch >= srangep[0]) & (xsearch <= srangep[1])
            xsearch = xsearch[sel]
            ysearch = ysearchƒ[sel]
        ipeaks = argrelextrema(ysearch, np.greater)[0]
        peaks = []
        for i in ipeaks:
            peaks.append([xsearch[i], ysearch[i]])
        self.peaks = peaks        
        return peaks
  
    def yscale_by(self, yscale):
        # TODO: docstring

        self.I = self.I*yscale
        
    def normalize_to(self, nstr):
        # TODO: docstring

        if nstr == 'max':            
            self.I = self.I/max(self.I)
        elif nstr == 'tail':
            self.I = self.I/self.I[-1]      
        else:
            self.I = self.I  

            
    def broaden0(self, g_sigma=None, g_fwhm=None, l_gamma=None, l_fwhm=None,
                 lvl=None):
        # TODO: docstring
         
        if g_sigma:
            self.E0, self.I0 = mconv(self.E0, self.I0).Gaussian(sigma=g_sigma)
        elif g_fwhm:
            self.E0, self.I0 = mconv(self.E0, self.I0).Gaussian(fwhm=g_fwhm)
            
        if l_gamma:
            self.E0, self.I0 = mconv(self.E0, self.I0).Lorentzian(gamma=l_gamma)            
        elif l_fwhm:
            self.E0, self.I0 = mconv(self.E0, self.I0).Lorentzian(fwhm=l_fwhm)      
                          
        if lvl:            
            if len(lvl) == 1:
                self.E0, self.I0 = mconv(self.E0, self.I0).\
                                   LorentzianVL(A=lvl[0], B=None, offset=None)
            elif len(lvl) == 2:
                self.E0, self.I0 = mconv(self.E0, self.I0).\
                                   LorentzianVL(A=lvl[0], B=lvl[1],
                                                offset=None)               
            else:
                self.E0, self.I0 = mconv(self.E0, self.I0).\
                                   LorentzianVL(A=lvl[0], B=lvl[1],
                                                offset=lvl[2])  

 
    def transform(self, irange=None, x0shift=False,
                  y0shift=False, normalize='max', xshift=None, std=False):
        # TODO: docstring

        self.E = self.E0.copy()
        self.I = self.I0.copy()        
        if x0shift:
            self.E = self.E -self.Eonset                        
        if xshift:
            self.E = self.E + xshift            
        if irange:
            self.Interpolate(irange)                  
        if y0shift:
            self.I = self.I -self.I[0]              
        if normalize == 'max':            
            self.I = self.I/max(self.I)
        elif normalize == 'tail':
            self.I = self.I/self.I[-1]
        elif normalize == 'none':
            self.I = self.I         
        else:
            self.I = self.I/max(self.I) 
            
        if std:
            self.I = (self.I-np.mean(self.I))/np.std(self.I)

        self.rank = rankdata(self.I, method='average') / len(self.I)
            

    def get_nn(self, nnradius=10.01, axin=None, yshift=0, ms=9, lp='k-', ts=0,
               text_off=False, labels_off=False, atoms_off=False):
        # TODO: docstring

        if self.structure:
            env = self.structure[0].\
                  get_sites_in_sphere(self.structure[0][self.structure[1]].\
                  coords, nnradius)

            def getKey(item):
                # TODO: docstring

                return item[1]

            self.env = sorted(env, key=getKey)
            ss = []
            for i in env:
                ss.append(i[0].specie.name)
            self.species = list(set(ss))
            
            if not atoms_off:
                atext = []
                cs = []
                c = 0
                s = 10
                for i in self.species:
                    ada = AnchoredDrawingArea(s*3, (len(self.species))*s, 0,
                                              0, loc=2, pad=0., frameon=False)
                    cs.append(Circle((s*2, c*s), 3, fc=get_c(i)))
                    atext.append(i)
                    c += 1
                
                at = AnchoredText('\n'.join(atext[::-1]), loc=2, 
                                  prop=dict(size=s), frameon=False,
                                  bbox_to_anchor=(0., 1.),
                                  bbox_transform=axin.transAxes)
                for i in cs:
                    ada.drawing_area.add_artist(i)
                ada.drawing_area.add_artist(at)               
                axin.add_artist(ada)

            if axin:    
                if self.vcn:
                    astr = 'mult.={:d}\nvcn={:04.2f}'.\
                           format(self.multiplicity, self.vcn)
                else:
                    astr = 'mult.={:d}\nvcn=na'.format(self.multiplicity)
                    
                if not text_off:
                    axin.annotate(astr, (-3, ts+yshift-0.1), fontsize=10,
                                  weight='bold')

                ss = []
                ds = []
                for i in self.env:
                    ss.append(i[0].specie.name)
                    ds.append(i[1])
                ds = np.array(ds, np.float)
                
                axin.plot(yshift+ds[0:21], lp)
                
                for i, d in enumerate(ds[0:21]):
                    c = get_c(ss[i])
                    axin.plot(i, yshift+d, 'o', color=c, ms=ms, alpha=0.8)

                axin.set_xticks(list(range(1, 21)))

                if text_off:
                    axin.set_xlim([-0.5, 21])
                else:
                    axin.set_xlim([-3.5, 21])
                    
                if not labels_off:
                    axin.set_xlabel('Neighbour index #');
                    axin.set_ylabel('Distance to absorbing atom ($\AA$)')
                else:
                    axin.set_xticklabels([])
        else:
            env = []
            if axin:
                print('nn info is not available.')
        
        return env


class mconv:
    # TODO: improve docstring
    '''Does convolution with Gaussian or Lorentzian windows.'''
     
    def __init__(self, datax, datay):

        # define self-variables that weren't initialized prior in __init__
        self.sigma = None
         
        self.datax = datax
        self.datay = datay
 
        self.NPoints = len(datax)
        self.X = np.array(datax)
 
    def w_gaussian(self, M, wsigma, dx):
        # TODO: improve docstring
        ''' see : https://github.com/scipy/scipy/blob/v0.19.0/scipy/signal
                  /windows.py#L1159-L1219
        M should be odd number.'''

        if wsigma <= 0:
            wsigma = 0.00001

        wsigma = wsigma/dx
        n = np.arange(0, M) - (M - 1.0) / 2.0
        wsigma2 = 2 * wsigma * wsigma
        wg = np.exp(-n**2/wsigma2)
        return wg
 
    def w_lorentzian(self, M, wgamma, dx):
        # TODO: docstring

        wgamma = wgamma/dx
        if wgamma <= 0:
            wgamma = 0.00001

        n = np.arange(0, M)-(M-1.0)/2.0
        wl = 1.0/(((2*n)/wgamma)**2 + 1)
        return wl             
            
         
    def Gaussian(self, sigma=None, fwhm=None, saveto=None):
        # TODO: doctstring

        if fwhm and sigma:
            print('ignoring input sigma')
        elif sigma:
            fwhm = sigma * np.sqrt(8 * np.log(2))
        else: 
            raise ValueError('sigma/fwhm was not set....')  

        self.sigma = fwhm/np.sqrt(8 * np.log(2))   
         
        M = 101
        diff = [self.datax[i+1] - self.datax[i] \
                for i in range(len(self.datax)-1)]
        dx= np.sum(diff)/len(self.datax) 
        win = self.w_gaussian(M, self.sigma, dx)
        out = signal.convolve(self.datay, win, mode='same')/sum(win)
        
        if saveto:
            if fmt is None:
                fmt="%18.6e %18.6e"
            of = np.column_stack((self.X, out))
            np.savetxt(str(saveto), of, delimiter=" ", fmt=fmt)
             
        return [self.X, out]
             
 
    def Lorentzian(self, gamma=None, fwhm=None, saveto=None, M=None):
        # TODO: docstring
        # in Lorentzian fwhm is equal to gamma
         
        if fwhm:
            if gamma:
                print('in Lorentzian fwhm is equal to gamma')
        elif gamma:
            fwhm = gamma
        else: 
            raise ValueError('sigma/fwhm was not set....')   
             
        self.gamma = fwhm 
         
        if M is None:
            M = 1001
        diff = [self.datax[i+1] - self.datax[i] \
                for i in range(len(self.datax)-1)]
        dx= np.sum(diff)/len(self.datax)
        win = self.w_lorentzian(M, self.gamma, dx)
        out = signal.convolve(self.datay, win, mode='same')/sum(win)
 
        if saveto:
            if fmt is None:
                fmt = "%18.6e %18.6e"
            of = np.column_stack((self.X, out))
            np.savetxt(str(saveto), of, delimiter=" ", fmt=fmt)
             
        return [self.X, out]    
     
     
    def LorentzianVL(self, saveto=None, fmt=None, M=None, A=None, B=None,
                     offset=None):
        # TODO: docstring
        # gamma = A(x-offset) + B
         
        if A is None:
            A = 0.1
        if B is None:
            B = 0
        if offset is None:
            offset = self.datax[0]
         
        gammas = []    
        for i, d in enumerate(self.datax):
            g = max(0, A*(d-offset)) + B
            gammas.append(g)
   
        if M is None:
            M = 1001
        diff = [self.datax[i+1] - self.datax[i] \
                for i in range(len(self.datax)-1)]
        dx = np.sum(diff)/len(self.datax) 
        out = np.zeros(self.NPoints)
         
                 
        for i, [gg, __, __] in enumerate(zip(gammas, self.datax, self.datay)):
            win = self.w_lorentzian(M, gg, dx)
            c = signal.convolve(self.datay, win, mode='same')/sum(win)
            out[i] = c[i]
                         
        if saveto:
            if fmt is None:
                fmt = "%18.6e %18.6e"
            of = np.column_stack((self.X, out))
            np.savetxt(str(saveto), of, delimiter=" ", fmt=fmt)
             
        #return [self.X, out], [self.datax, gammas]  
        return [self.X, out]
 
     
     
    def LorentzianVE2(self, saveto=None, fmt=None, M=None, offset=None):
        # TODO: docstring
        # gamma = (x-offset)^2 + B   
         
        A = 0.1
        B = 0
        if offset is None:
            offset = self.datax[0]  
         
        gammas = []    
        for i, d in enumerate(self.datax):
            g = max(0, A*(d-offset)) + B
            gammas.append(g)
   
        if M is None:
            M = 1001  
        diff = [self.datax[i+1] - self.datax[i] \
                for i in range(len(self.datax)-1)]
        dx = np.sum(diff)/len(self.datax)
        out = np.zeros(self.NPoints)
         
                 
        for i, [gg, __, __] in enumerate(zip(gammas, self.datax, self.datay)):
            win = self.w_lorentzian(M, gg, dx)
            c = signal.convolve(self.datay, win, mode='same')/sum(win)
            out[i] = c[i]
                         
        if saveto:
            if fmt is None:
                fmt = "%18.6e %18.6e"
            of = np.column_stack((self.X, out))
            np.savetxt(str(saveto), of, delimiter=" ", fmt=fmt)
             
        #return [self.X, out], [self.datax, gammas]   
        return [self.X, out]


class xanes_collection:
    # TODO: docstring

    def __init__(self, absorption_specie, collection_name='not_set'):
        
        self.absorption_specie = absorption_specie
        self.collection_name = collection_name        
        self.collection = []
        self.ids = []
        
        self.ave_spectra = []
        self.site_spectra = []
        
        
    def build(self, xanes_list):
        # TODO: docstring
        
        collection = []
        for i in xanes_list:
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
        

def read_xanes(path, absorption_specie, order='eof', skip_missing=False,
               symprec=0.01, ang_tol=5):
    # TODO: docstring

    here = os.getcwd() 

    if os.path.isdir(path):
        os.chdir(path)
        struct = mg.Structure.from_file("CONTCAR")
        finder = SpacegroupAnalyzer(struct,symprec,ang_tol)
        struct = finder.get_symmetrized_structure()
        [sites, indices]  = struct.equivalent_sites, struct.equivalent_indices

        s = []
        for i in struct:
            s.append(i.species_string)
        if absorption_specie not in s:
            if skip_missing:
                os.chdir(here)                
                print('This structure does not have' + absorption_specie
                      + 'in it. Skipping...')
                return [[], []] 
            else:
                os.chdir(here)                
                raise ValueError('This structure does not have'
                                 + absorption_specie
                                 + 'element in it. Please check...')

        xanes = []
        for i, s in enumerate(sites):
            if s[0].species_string is absorption_specie:   
                
                fe = 'exciting_{:03d}_{}-K'.format(indices[i][0]+1,
                                                   absorption_specie)
                fo = 'ocean_{:03d}_{}-K'.format(indices[i][0]+1,
                                                absorption_specie)
                ff = 'feff_{:03d}_{}-K'.format(indices[i][0]+1,
                                               absorption_specie)
                
                if order is 'eo':
                    if os.path.isdir(fe):
                        os.chdir(fe)  
                        pload = pickle.load(open('xanes.pkl', 'rb'))
                        xanes.append(pload)
                        os.chdir('..')
                    elif os.path.isdir(fo):
                        os.chdir(fo)  
                        pload = pickle.load(open('xanes.pkl', 'rb'))
                        xanes.append(pload)
                        os.chdir('..')                         
                    else:
                        if skip_missing:
                            print(fe + 'is not available in' 
                                  + path + '. Skipping...')
                            os.chdir(here)
                            return [[], []] 
                        else:
                            os.chdir(here)                        
                            raise FileNotFoundError(fe + 'is not available in'
                                                    + path)                 

                if order is 'oe':
                    if os.path.isdir(fo):
                        os.chdir(fo)  
                        pload = pickle.load(open('xanes.pkl', 'rb'))
                        xanes.append(pload)
                        os.chdir('..')
                    elif os.path.isdir(fe):
                        os.chdir(fe)  
                        pload = pickle.load(open('xanes.pkl', 'rb'))
                        xanes.append(pload)
                        os.chdir('..')                         
                    else:
                        if skip_missing:
                            print(fo + 'is not available in' + path 
                                  + '. Skipping...')
                            os.chdir(here)
                            return [[], []] 
                        else:
                            os.chdir(here)                        
                            raise FileNotFoundError(fo + 'is not available in'
                                                    + path)  
                                                        
                if order is 'eof':
                    if os.path.isdir(fe):
                        os.chdir(fe)  
                        pload = pickle.load(open('xanes.pkl', 'rb'))
                        xanes.append(pload)
                        os.chdir('..')
                    elif os.path.isdir(fo):
                        os.chdir(fo)  
                        pload = pickle.load(open('xanes.pkl', 'rb'))
                        xanes.append(pload)
                        os.chdir('..')      
                    elif os.path.isdir(ff):
                        os.chdir(ff)  
                        pload = pickle.load(open('xanes.pkl', 'rb'))
                        xanes.append(pload)
                        os.chdir('..')                           
                    else:
                        if skip_missing:
                            print(fe + 'is not available in' + path
                                  + '. Skipping...')
                            os.chdir(here)
                            return [[], []] 
                        else:
                            os.chdir(here)                        
                            raise FileNotFoundError(fe + 'is not available in'
                                                    + path)                                      

                if order is 'f':
                    if os.path.isdir(ff):
                        os.chdir(ff)  
                        pload = pickle.load(open('xanes.pkl', 'rb'))
                        xanes.append(pload)
                        os.chdir('..')                    
                    else:
                        if skip_missing:
                            print(ff + 'is not available in' + path
                                  + '. Skipping...')
                            os.chdir(here)
                            return [[], []] 
                        else:
                            os.chdir(here)                        
                            raise FileNotFoundError(ff + 'is not available in'
                                                    + path)                              

                if order is 'o':
                    if os.path.isdir(fo):
                        os.chdir(fo)  
                        pload = pickle.load(open('xanes.pkl', 'rb'))
                        xanes.append(pload)
                        os.chdir('..')                    
                    else:
                        if skip_missing:
                            print(fo + 'is not available in' + path
                                  + '. Skipping...')
                            os.chdir(here)
                            return [[], []] 
                        else:
                            os.chdir(here)                        
                            raise FileNotFoundError(fo + 'is not available in'
                                                    + path)

                if order is 'e':
                    if os.path.isdir(fe):
                        os.chdir(fe)  
                        pload = pickle.load(open('xanes.pkl', 'rb'))
                        xanes.append(pload)
                        os.chdir('..')                    
                    else:
                        if skip_missing:
                            print(fe + 'is not available in' + path
                                  + '. Skipping...')
                            os.chdir(here)
                            return [[], []] 
                        else:
                            os.chdir(here)                        
                            raise FileNotFoundError(fe + 'is not available in'
                                                    + path)
                        

        # get E limits
        minmaxs = []
        for i in xanes:
            minmaxs.append([min(i.E0), max(i.E0)])
        minmaxs = np.array(minmaxs)    
        irange = [max(minmaxs[:, 0]), min(minmaxs[:, 1])]   
        e_int = np.linspace(irange[0], irange[1],
                            int((irange[1] - irange[0]) / 0.1) + 1)

        ave_xanes = e_int * 0
        ave_vcn = 0
        site_xanes = []
        counter = 0
        for i in xanes:
            i.transform(irange=irange, normalize=None, y0shift=None)
            ave_xanes += i.I * i.multiplicity
            site_onset = i.Eonset
            site_xanes.append(mXANES(data=[i.E, i.I], structure=i.structure,
                                     xanesid=i.xanesid, source=i.source,
                                     Eonset=site_onset,
                                     edge=i.edge,
                                     multiplicity=i.multiplicity))
            if i.vcn:
                ave_vcn += i.vcn*i.multiplicity
            else:
                try:
                    nnfinder = VoronoiNN(cutoff=10, allow_pathological=True)
                    print(i.structure)
                    # vcn = nnfinder.get_cn(i.structure[0], i.structure[1],
                    #                       use_weights=True)
                    ave_vcn += i.vcn*i.multiplicity
                except Exception as exc:
                    # print(exc)
                    # print('warning: vcn info is missing for site. '
                    #       'ave_vnc is not correct. ')
                    pass
            counter += i.multiplicity
        ave_xanes = ave_xanes/counter
        ave_vcn = ave_vcn/counter
        ave_xanes = mXANES(data=[e_int, ave_xanes], Eonset=site_onset,
                           structure=[struct, -1, ave_vcn], xanesid=i.xanesid,
                           source=i.source, edge=i.edge)
        os.chdir(here)
        return [ave_xanes, site_xanes]
    
    else:
        if skip_missing:
            os.chdir(here)             
            print(path + 'is not available. Skipping...')
            return [[], []]
        else:
            os.chdir(here)            
            raise FileNotFoundError(path+' is not available.')
        
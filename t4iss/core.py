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


class Dataset:
    """
    A numpy array-type of data structure which represents the entire
    training set (and labels) of a machine learning problem. Alternatively,
    it can also be initialized as a list or a variety of kinds of other objects
    (as time goes on more features will be implemented).

    : __init__ arguments :
    * data (variable type, can be list, numpy array, etc...) The machine
      learning data, not necessarily sorted in a straightforward way.
    * spec (are these average spectra or something else) The data generated by
      the materials project or by ab initio quantum chemical calculations can
      represent the average or site spectra of a particular crystal structure.
      This property provides a bit more information about what kind of spectra
      they are.
        ~ Options ~
        'site' - These are spectra of individual sites.
        'avg'  - These are weighted average spectra taken across the entire
                 crystal structure (closest to what experimentalists will 
                 observe).
        'mix'  - There are both site and average spectra in data.
    * name (string) The name of the data, for human use in helping to remember
      what the data represents.
    * verbose (int) How much output do you want?
        ~ Options ~
        0 - silent running
        1 - output all information available
    * standardized_grid_classwide (np array) If merging two datasets, one of
      which has a standardized grid already defined for it, then specifying
      the standardized grid ahead of time will interpolate all input data to
      this given grid
    * sanity_cutoff_N (None-type or int) If set, any data with more data
      points than this cutoff are discarded
    * sanity_cutoff_min (None-type or float) If set, any data with a minimum
      value less than this are discarded
    * sanity_cutoff_max (None-type or float) If set, any data with a maximum
      value greater than this are discarded
    """
  
    def __init__(self, data, spec='avg', name=None, verbose=1,
                 standardized_grid_classwide=None,
                 sanity_cutoff_N=None, sanity_cutoff_min=None,
                 sanity_cutoff_max=None, check_has_dictionary=True):

        self.data = data
        self.type = type(data)
        self.spec = spec
        self.verbose = verbose
    
        # set the traindev/test split for prep with use in keras
        self.traindev_X = None
        self.traindev_y = None
        self.test_X = None
        self.test_y = None
    
        # is the current data data_X and data_y shuffled?
        self.is_shuffled = False

        # data and labels
        self.data_X = None
        self.data_y = None

        # track which directory the spectra originated from
        self.tracker = None

        # save the post-numpy array converted X and y values
        self.original_X = None
        self.original_y = None
        self.original_tracker = None

        # name the dataset so you don't forget what it represents
        # obviously optional
        self.name = name

        # total number of training examples
        self.m = None

        # total number of features
        self.n = None

        # are every one of the input feature lengths the same?
        self.n_all_same = True

        # are every one of the E0 grids the same upon initial input?
        self.grid_all_same = True

        if standardized_grid_classwide is None:
            self.standardized_grid = None
        else:
            self.standardized_grid = standardized_grid_classwide

        # total number of labels (classes)
        # also the total number of unique coordination environments
        self.c = 0

        # dictionary of basis vectors
        self.dictionary = {}

        if self.type == list and self.spec == 'avg':

            # determine information about the number of features
            self.n = len(self.data[0][2].E0)

            # default sanity checks on ridiculous spectra
            if check_has_dictionary:
                self.data = list(filter(lambda x: len(x) == 4, self.data))

            if sanity_cutoff_N is not None:
                sCN = sanity_cutoff_N
                self.data = list(filter(lambda x: len(x[2].E0) < sCN,
                                        self.data))

            if sanity_cutoff_min is not None:
                scMin = sanity_cutoff_min
                self.data = list(filter(lambda x: min(x[2].E0) > scMin,
                                        self.data))

            if sanity_cutoff_max is not None:
                scMax = sanity_cutoff_max
                self.data = list(filter(lambda x: max(x[2].E0) < scMax,
                                        self.data))

            self.m = len(self.data)

            # determine whether or not each training input has the same number of
            # features
            for i, cell in enumerate(self.data):
                if self.n != len(cell[2].E0):
                    self.n_all_same = False

                if self.grid_all_same:
                    try:
                        np.testing.assert_array_equal(self.data[i-1][2].E0, cell)
                    except AssertionError:
                        self.grid_all_same = False
              
                # fatal error: if the number of x-points and number of y-points on a
                # single spectra are not equal, terminate
                np.testing.assert_equal(len(cell[2].E0), len(cell[2].I0))

            # if by some chance all the E0 grids are identical, choose the first one
            # and define that as the standardized grid
            if self.grid_all_same:
                self.standardized_grid = self.data[0][2].E0
          
            # determine the total number of classes
            for cell in self.data:
                if len(cell) == 4:
                    for k in cell[3].keys():
                        if k not in self.dictionary:
                            self.dictionary[k] = self.c
                            self.c += 1


        elif self.type == list and self.spec == 'site':
            total_m = 0
            for ce in range(len(self.data)):
                total_m += len(self.data[ce])
            self.m = total_m
            self.n = len(self.data[0][0][2].E0)

            for ce in range(len(self.data)):
                for i, spectra in enumerate(self.data[ce]):
                    if self.n != len(spectra[2].E0):
                        self.n_all_same = False

                    # assert fatal error
                    np.testing.assert_equal(len(spectra[2].E0), len(spectra[2].I0))
            
            # assume grids are non-identical
            self.grid_all_same = False

            # determine the total number of classes
            self.c = len(self.data)

            # critical assertion: code is not setup for more or less than
            # 3 classes for site runs (issue!)
            np.testing.assert_equal(self.c, 3)

        else:
            raise TypeError("Input data/spec type combination not yet supported.")

        # if the standardized grid was defined, that by construction must 
        # represent the number of features
        if self.standardized_grid is not None:
            self.n = len(self.standardized_grid)
  
  
    def sanity_check(self):
        """Execute all assertions on the current self.data_X and self.data_y.
        Recommended to run before feeding the data into a ML algorithm."""
        
        if self.data_X is None or self.data_y is None:
            raise ValueError("data_X or data_y not yet defined. Must convert to "
                             "numpy arrays first.")
        
        failed_tests = []
        
        # check shapes, ensure they match properly
        try:
            np.testing.assert_equal(self.data_X.shape[0], self.data_y.shape[0])
        except AssertionError:
            failed_tests.append("data_X, data_y have different n.o. train examples")
        
        # check the class vectors, they must have norm unity
        y_unity = np.sum(self.data_y, axis=1)
        for val in y_unity:
            if not 0.999 < val < 1.001:
                failed_tests.append("class vectors do not sum to unity")
                break
        
        # ensure no negative values in the training data
        try:
            np.testing.assert_equal(True, np.all(self.data_X >= 0))
        except AssertionError:
            failed_tests.append("some data in self.data_X is negative")
         
        if failed_tests == []:
            print("All tests passed.")
        else:
            print("\nWARNING: CHECKS FAILED")
            for test in failed_tests:
                print(test)
    

    def generate_standardized_grid(self):
        """Inputs data and uses splines to interpolate it all onto a standard
        grid such that the feature entries match up."""
        
        if self.spec not in ['avg', 'site']:
            raise TypeError("Method not yet defined for this type.")
        
        elif self.grid_all_same:
            if self.verbose == 1:
                print("Grid is already standardized, nothing will be changed.")

        elif self.standardized_grid is not None:
            # this is the special case in which a grid is provided, likely
            # from another data set, the current data should be interpolated to
            # this grid
            if self.verbose == 1:
                print("Grid has been provided, nothing will be changed.")
        
        elif self.spec == 'avg':
            all_n = []
            all_max = []
            all_min = []
       
            for cell in self.data:
                spectra_E0 = cell[2].E0
                all_n.append(len(spectra_E0))
                all_max.append(max(spectra_E0))
                all_min.append(min(spectra_E0))

            # these will be the values used in the splines interpolation
            avg_n = int(np.ceil(np.mean(all_n)))
            u = np.ceil(max(all_max))
            l = np.floor(min(all_min))

            # define the standardized grid
            self.standardized_grid = np.linspace(l, u, avg_n, endpoint=True)
          
            # redefine the number of features
            self.n = avg_n
          
            if self.verbose == 1:
                print("Grid has been standardized.")

        else:

            # execute a similar procedure to that in the 'avg' case
            all_n = []
            all_max = []
            all_min = []

            for ce in range(len(self.data)):
                for spectra in self.data[ce]:
                    spectra_E0 = spectra[2].E0
                    all_n.append(len(spectra_E0))
                    all_max.append(max(spectra_E0))
                    all_min.append(min(spectra_E0))

            avg_n = int(np.ceil(np.mean(all_n)))
            u = np.ceil(max(all_max))
            l = np.floor(min(all_min))

            self.standardized_grid = np.linspace(l, u, avg_n, endpoint=True)
            self.n = avg_n
            if self.verbose == 1:
                print("Grid has been standardized.")


    def convert_to_numpy_arrays(self):
        """Takes a list of data of the form
        [directory name, index, t4iss.core.mXANES object, {% of coord env}]
        and converts it to two numpy arrays, one of the data, the other the labels.
        """

        from scipy.interpolate import InterpolatedUnivariateSpline as ius
        
        if self.type == np.ndarray:
            raise TypeError("Data is already of type np.ndarray.")

        # assert that every training example has the same number of features,
        # or that if that isn't the case that a standardized grid has been
        # computed for spline interpolation
        if self.standardized_grid is None:
            raise RuntimeError("All training examples do not have the same number "
                               "of features and grid has not been standardized. "
                               "Standardize the grid before calling this method.")

        if self.spec not in ['avg', 'site']:
            raise TypeError("Method not yet defined for this type.")
        
        # initialize the X and y tensors which represent the training data and
        # the labels respectively
        X = np.empty((self.m, self.n))
        y = np.zeros((self.m, self.c))
        
        # initialize the tracker to the same shape
        tracker = np.array(range(self.m), dtype='a10').reshape(self.m, 1)
        
        # first, set the X matrix by running through all training examples
        # then, set the y matrix by analyzing the classes through the dictionary
        if self.spec == 'avg':
            for mm, cell in enumerate(self.data):
                # X
                if self.grid_all_same:
                    X[mm, :] = np.array(cell[2].I0).reshape(1, self.n)
                else:
                    s = ius(np.squeeze(cell[2].E0), np.squeeze(cell[2].I0), k=1)
                    E1 = np.reshape(s(np.squeeze(self.standardized_grid)), (1, self.n))
                    X[mm, :] = E1.squeeze()

                # y
                temp_dictionary = cell[3]
                for k in temp_dictionary.keys():
                    prop = temp_dictionary[k]
                    basis_index = self.dictionary[k]
                    y[mm, basis_index] = prop
              
                # tracker
                tracker[mm] = cell[0]

        else:

            # temporary dictionary sites
            tds = {'O:6' : 0, 'S:5' : 1, 'T:4' : 2}
            self.dictionary = tds

            # temporary counter sites
            tcs = 0

            for ce in range(len(self.data)):
                for spectra in self.data[ce]:

                    # X 
                    if self.grid_all_same:
                        X[tcs, :] = np.array(spectra[2].I0).reshape(1, self.n)
                    else:
                        s = ius(np.squeeze(spectra[2].E0), np.squeeze(spectra[2].I0), k=1)
                        E1 = np.reshape(s(np.squeeze(self.standardized_grid)), (1, self.n))
                        X[tcs, :] = E1.squeeze()

                    # y (assume only 3 classes for now)
                    # this is a potential long term issue!
                    y[tcs, tds[spectra[3][0]]] = 1
                    tracker[tcs] = spectra[0]
                    tcs += 1

        self.data_X = X
        self.data_y = y
        self.original_X = X
        self.original_y = y
        self.tracker = tracker
        self.original_tracker = tracker

    
    def remove_outliers(self, remove_negative=True, remove_large=False,
                        large_cutoff=3.0, advanced=False,
                        advanced_parameter=100.0):
    
        if self.data_X is None or self.data_y is None:
            raise ValueError("data_X or data_y not yet defined. Must convert to "
                             "numpy arrays first.")
        
        X = self.data_X
        y = self.data_y
        tracker = self.tracker
        
        initial_m = X.shape[0]
        
        # there are some outliers that are mostly negative or zero
        # we want to eliminate those since they are definitely unphysical
        X_bool = [X >= 0.0][0]
        X_bool_s = np.sum(X_bool, axis=1)
        X_lt_cutoff = X_bool_s < self.n/2
        X = np.delete(X, np.where(X_lt_cutoff)[0], axis=0)
        y = np.delete(y, np.where(X_lt_cutoff)[0], axis=0)

        # there are also some spectra that are so large they are likely errors
        # in the electronic structure calculations
        if remove_large:
            X_bool = [X > large_cutoff][0]
            X_bool_s = np.sum(X_bool, axis=1)
            X_cutoff = X_bool_s != 0
            X = np.delete(X, np.where(X_cutoff)[0], axis=0)
            y = np.delete(y, np.where(X_cutoff)[0], axis=0)
            tracker = np.delete(tracker, np.where(X_cutoff)[0], axis=0)
        
        if advanced:
            # compute the average spectrum
            mean_spectrum = np.mean(X, axis=0)

            # use broadcasting to shift all spectra by the average one
            X_shift = X - mean_spectrum

            # element-wise square
            X_shift *= X_shift

            # and sum
            X_sum = np.sum(X_shift, axis=1, keepdims=True)

            worst_val = max(X_sum)

            X_shift = np.delete(X_shift, np.where(X_sum == max(X_sum))[0], axis=0)
            X = np.delete(X, np.where(X_sum == max(X_sum))[0], axis=0)
            y = np.delete(y, np.where(X_sum == max(X_sum))[0], axis=0)
            tracker = np.delete(tracker, np.where(X_sum == max(X_sum))[0], axis=0)

            X_sum = np.delete(X_sum, np.where(X_sum == max(X_sum))[0], axis=0)

            next_worst_val = max(X_sum)

            while worst_val - next_worst_val > advanced_parameter:
                worst_val = next_worst_val
                X_shift = np.delete(X_shift, np.where(X_sum == max(X_sum))[0], axis=0)
                X = np.delete(X, np.where(X_sum == max(X_sum))[0], axis=0)
                y = np.delete(y, np.where(X_sum == max(X_sum))[0], axis=0)
                tracker = np.delete(tracker, np.where(X_sum == max(X_sum))[0], axis=0)
                X_sum = np.delete(X_sum, np.where(X_sum == max(X_sum))[0], axis=0)
                next_worst_val = max(X_sum)


        # ensure the extrapolation didn't make some values negative
        if remove_negative:
            X[X < 0.0] = 0.0
        
        final_m = X.shape[0]

        if self.verbose == 1:
            print("Outlier removal:")
            if final_m != initial_m:
                print("  initial -> final number of training examples: %i->%i"
                      % (initial_m, final_m))
            else:
                print("  No outliers.")
        print("  Negatives removed is set to %a." % remove_negative)
        print("  Remove large is set to %a." % remove_large)
        print("  Advanced removal is set to %a." % advanced)
        
        self.data_X = X
        self.data_y = y
        self.m = final_m
        self.tracker = tracker
    
    
    def keep_only(self, ce=['O:6', 'S:5', 'T:4']):
        """Allows further trimming of the dataset such that only certain labels
        are kept. This is useful since some of the chemical environments are only
        represented by a handful of spectra, which could mess up the training.
        """
        
        if self.data_X is None or self.data_y is None:
            raise ValueError("data_X or data_y not yet defined. Must convert to "
                             "numpy arrays first.")

        initial_m = self.m
        
        X = self.data_X
        y = self.data_y
        tracker = self.tracker
        
        # iterate through the entire data set, keep only the spectra that
        # correspond to the labels in the ce variable
        # when y_multiplier is element-wise multiplied onto a label, it will
        # zero out any labels we don't want to keep anymore
        
        y_multiplier = np.zeros(self.data_y[0].shape)

        # these will be the values we keep
        new_dictionary = {}
        for i in ce:
            y_multiplier[self.dictionary[i]] = 1

        # ensures the correct order
        for key, value in self.dictionary.items():
          if key in ce:
            new_dictionary[key] = value

        y_multiplier_invert = y_multiplier == 0
        y_temp = self.data_y * y_multiplier_invert
        y_bool = [y_temp != 0][0]
        y_bool_s = np.sum(y_bool, axis=1)
        y_cutoff = y_bool_s != 0

        X = np.delete(X, np.where(y_cutoff)[0], axis=0)
        y = np.delete(y, np.where(y_cutoff)[0], axis=0)
        tracker = np.delete(tracker, np.where(y_cutoff)[0], axis=0)

        # also trim the class vector
        y = np.delete(y, np.where(y_multiplier_invert)[0], axis=1)
        
        final_m = y.shape[0]
        
        if self.verbose == 1:
            print("Keep only execution:")
            if final_m == initial_m:
                print("  No data removed.")
            else:
                print("  initial -> final number of training examples: %i->%i"
                      % (initial_m, final_m))

        self.data_X = X
        self.data_y = y
        self.m = final_m
        self.c = y.shape[1]
        self.tracker = tracker
        self.dictionary = new_dictionary

    
    def emergency_revert(self):
        """Revert back to the post-numpy-transform data. Useful for if you
        really screwed up your data during the process of trimming it down."""
        if self.data_X is None or self.data_y is None:
            raise ValueError("data_X or data_y not yet defined. Must convert to "
                             "numpy arrays first.")
        self.data_X = self.original_X
        self.data_y = self.original_y
        self.traindev_X = None
        self.traindev_y = None
        self.test_X = None
        self.test_y = None
        self.tracker = self.original_tracker
        
        if self.verbose == 1:
            print("Data reverted: stored in self.data_X and self.data_y.")
    
    
    def augment(self, energy_shift=True, energy_shift_n=1,
                intensity_stretch=False, only_mixed=False):
        """Data augmentation suite."""
        
        from scipy.interpolate import InterpolatedUnivariateSpline as ius
        
        if self.data_X is None or self.data_y is None:
            raise ValueError("data_X or data_y not yet defined. Must convert to "
                             "numpy arrays first.")
          
        if self.standardized_grid is None:
            raise ValueError("Standardized grid must be defined before running "
                             "data augmentation.")
        
        X = self.data_X
        y = self.data_y
        sg = self.standardized_grid
        tracker = self.tracker
        
        initial_m = self.m

        if energy_shift:
            X_to_append = np.empty((X.shape[0]*energy_shift_n*2, X.shape[1]))
            y_to_append = np.empty((y.shape[0]*energy_shift_n*2, y.shape[1]))
            tracker_to_append = \
              np.array(range(X.shape[0]*energy_shift_n*2), 
                       dtype='a10').reshape(X.shape[0]*energy_shift_n*2, 1)
            counter = 0
            for i in range(self.m):
                current_y = y[i, :]
                current_X = X[i, :]
                current_tracker = tracker[i, :]
                for j in range(energy_shift_n):
                  
                    # value in the units of energy here
                    shift_value = j + 1
                    grid_R = sg + shift_value
                    grid_L = sg - shift_value
                  
                    # treat these new grids as the actual ones and interpolate back
                    # to the standardized grid
                  
                    sR = ius(np.squeeze(grid_R), np.squeeze(current_X), k=1)
                    ER = np.reshape(sR(np.squeeze(sg)), (1, self.n))
                    sL = ius(np.squeeze(grid_L), np.squeeze(current_X), k=1)
                    EL = np.reshape(sL(np.squeeze(sg)), (1, self.n))
                  
                    X_to_append[counter, :] = np.squeeze(ER)
                    y_to_append[counter, :] = current_y
                    tracker_to_append[counter, :] = current_tracker
                    counter += 1
            
                    X_to_append[counter, :] = np.squeeze(EL)
                    y_to_append[counter, :] = current_y
                    tracker_to_append[counter, :] = current_tracker
                    counter += 1
        
            # ensure no negative values
            X_to_append[X_to_append < 0.0] = 0.0
          
            # concatenate
            X = np.concatenate((X, X_to_append))
            y = np.concatenate((y, y_to_append))
            tracker = np.concatenate((tracker, tracker_to_append))
          
            self.data_X = X
            self.data_y = y
            final_m = X.shape[0]
            self.m = X.shape[0]
            self.tracker = tracker
            self.is_shuffled = False
          
            if self.verbose == 1:
                print("Augmentation - energy shift:")
                if final_m == initial_m:
                    print("  No data augmented.")
                else:
                    print("  initial -> final number of training examples: %i->%i"
                          % (initial_m, final_m))

        if intensity_stretch:
            initial_m = self.m
            X = self.data_X
            y = self.data_y
            tracker = self.tracker

            all_counter = 0
            if only_mixed:
                for i in range(self.m):
                    if not np.any(y[i, :] == 1):
                        all_counter += 1
            else:
                all_counter = initial_m

            X_to_append = np.empty((all_counter, self.n))
            y_to_append = np.zeros((all_counter, self.c))
            tracker_to_append = \
                np.array(range(all_counter),
                         dtype='a10').reshape(all_counter, 1)

            if only_mixed:
                counter = 0
                for i in range(self.m):
                    if not np.any(y[i, :] == 1):
                        avg_val = np.mean(X[1, :])
                        val = (X[i, :] - avg_val) * 1.1 + avg_val
                        X_to_append[counter, :] = val
                        y_to_append[counter, :] = y[i, :]
                        tracker_to_append[counter, :] = tracker[i, :]
                        counter += 1
            else:
                for i in range(self.m):
                    avg_val = np.mean(X[1, :])
                    val = (X[i, :] - avg_val) * 1.1 + avg_val
                    X_to_append[i, :] = val
                    y_to_append[i, :] = y[i, :]
                    tracker_to_append[i, :] = tracker[i, :]

            X_to_append[X_to_append < 0] = 0

            # concatenate
            X = np.concatenate((X, X_to_append))
            y = np.concatenate((y, y_to_append))
            tracker = np.concatenate((tracker, tracker_to_append))

            self.data_X = X
            self.data_y = y
            self.m = X.shape[0]
            final_m = X.shape[0]
            self.tracker = tracker
            self.is_shuffled = False

            if self.verbose == 1:
                print("Augmentation - intensity stretch:")
                print("  Mixed only is %a" % only_mixed)
                if final_m == initial_m:
                    print("  No data augmented.")
                else:
                    print("  initial -> final number of training examples: %i->%i"
                          % (initial_m, final_m))
          
    def shuffle(self, seed):
    
        if self.data_X is None or self.data_y is None:
            raise ValueError("data_X or data_y not yet defined. Must convert to "
                             "numpy arrays first.")
        
        # set the seed and shuffle X
        np.random.seed(seed)
        np.random.shuffle(self.data_X)
        
        # set the SAME SEED and shuffle y, it is critical the seeds are the same
        np.random.seed(seed)
        np.random.shuffle(self.data_y)
        
        # same thing for the tracker
        np.random.seed(seed)
        np.random.shuffle(self.tracker)
        
        self.is_shuffled = True
        
        if self.verbose == 1:
            print("Data has been shuffled with seed %i" % seed)
    
  
    def get_distribution_unique(self):
    
        if self.data_X is None or self.data_y is None:
            raise ValueError("data_X or data_y not yet defined. Must convert to "
                             "numpy arrays first.")

        X = self.data_X
        y = self.data_y
        y_counter = np.zeros(y[0, :].shape)
        for i in range(y.shape[0]):
            for j in range(y.shape[1]):
                if y[i][j] == 1:
                    y_counter[j] += 1
                    break
        return y_counter
  
  
    def trim_by_minimum(self):
    
        if self.data_X is None or self.data_y is None:
            raise ValueError("data_X or data_y not yet defined. Must convert to "
                             "numpy arrays first.")
        
        if not self.is_shuffled:
            raise ValueError("Data must be shuffled to call this method.")
        
        X = self.data_X
        y = self.data_y
        tracker = self.tracker
        
        initial_m = X.shape[0]
        
        y_counter = self.get_distribution_unique()
        min_val = int(min(y_counter))
        for i in range(len(y_counter)):
            if y_counter[i] == min_val:
                min_index = i
                break
            
        mix_counter = 0
        for i in range(self.m):
            if not np.any(self.data_y[i, :] == 1):
                mix_counter += 1
        
        X_new = np.empty((min_val*len(y_counter) + mix_counter, self.n))
        y_new = np.empty((min_val*len(y_counter) + mix_counter, self.c))
        tracker_new = np.array(range(min_val*len(y_counter) + mix_counter), 
                               dtype='a10').reshape(min_val*len(y_counter) 
                                                    + mix_counter, 1)

        y_counter_now = [0 for i in range(len(y_counter))]
        now_counter = 0
        
        for i in range(self.m):
            if not np.any(self.data_y[i, :] == 1):
                X_new[now_counter, :] = X[i, :]
                y_new[now_counter, :] = y[i, :]
                tracker_new[now_counter, :] = tracker[i, :]
                now_counter += 1
            else:
                for j in range(self.c):
                    if y[i, j] == 1 and y_counter_now[j] < min_val:
                        X_new[now_counter, :] = X[i, :]
                        y_new[now_counter, :] = y[i, :]
                        tracker_new[now_counter, :] = tracker[i, :]
                        y_counter_now[j] += 1
                        now_counter += 1
                        break
        
        self.data_X = X_new
        self.data_y = y_new
        final_m = X_new.shape[0]
        self.m = X_new.shape[0]
        self.tracker = tracker_new
        
        if self.verbose == 1:
            print("Trim by minimum:")
            if final_m == initial_m:
                print("  No data changed.")
            else:
                print("  initial -> final number of training examples: %i->%i"
                      % (initial_m, final_m))


    def scale(self, method='1'):
        """Normalizes all spectra."""

        if self.data_X is None or self.data_y is None:
            raise ValueError("data_X or data_y not yet defined. Must convert to "
                             "numpy arrays first.")

        X = self.data_X

        if method == '1':
            # scale all spectra such that the maxima are 1.0
            X_max = np.max(X, axis=1, keepdims=True)
            X /= X_max
            self.data_X = X
        else:
            raise RuntimeError("Method not supported.")


    def mass_augment(self, maximum_new_data=1000, randint_min=0,
                     randint_max=6, replace=False):
        """Mass data synthesis of fake data from site-spectra. Should be run
        after keeping only the minimum."""

        if self.data_X is None or self.data_y is None:
            raise ValueError("data_X or data_y not yet defined. Must convert to "
                             "numpy arrays first.")

        if self.spec == 'avg':
            raise RuntimeError("Average data should not be mass augmented.")

        distribution = self.get_distribution_unique().astype(int)

        if np.any(distribution/max(distribution) != 1.0):
            raise RuntimeError("Must run after keep_only.")

        X_tensor = np.empty((distribution[0], self.n, self.c))
        X = self.data_X
        y = self.data_y
        tracker = self.tracker

        X_to_append = np.empty((maximum_new_data, self.n))
        y_to_append = np.empty((maximum_new_data, self.c))
        null_tracker = np.zeros((maximum_new_data, 1))

        initial_m = self.m

        index_counter = np.zeros(distribution.shape).astype(int)

        for i in range(self.m):
            for j in range(self.c):
                if y[i, j] == 1:
                    index = j
                    break
            X_tensor[index_counter[index], :, index] = X[i, :]
            index_counter[index] += 1

        total_counter = 0
        while total_counter < maximum_new_data:

            # for every class, generate a random integer representing the
            # amount of contribution from that class
            cvec = []
            ctotal = 0
            for __ in range(self.c):
                cvec.append(np.random.randint(randint_min, randint_max))
            #cvec = np.array(cvec)
            
            ctotal = np.sum(cvec)

            # ensure we have a mixed result
            if ctotal > 1 and not np.any((np.array(cvec)/ctotal) == 1.0):

                # current counter for the running average
                ccc = 0

                # generate the normalized label
                y_to_append[total_counter, :] = cvec/ctotal

                running_average = np.zeros((ctotal, self.n))

                # for every class type
                for cc_index, cc in enumerate(cvec):

                    # get that many random numbers that represent which
                    # spectra of class ii to get
                    which_X = np.random.randint(0, distribution[cc_index], 
                                                size=cc)

                    for j in which_X:
                        running_average[ccc, :] = X_tensor[j, :, cc_index]
                        ccc += 1

                # sum
                running_average = np.sum(running_average, axis=0, keepdims=True)

                # scale it
                running_average /= ctotal

                # append it
                X_to_append[total_counter, :] = running_average

                total_counter += 1

        if not replace:
            X = np.concatenate((X, X_to_append))
            y = np.concatenate((y, y_to_append))
            tracker = np.concatenate((tracker, null_tracker))

            self.tracker = tracker
            self.data_X = X
            self.data_y = y
            self.m = self.data_X.shape[0]
        else:
            self.data_X = X_to_append
            self.data_y = y_to_append
            self.m = X_to_append.shape[0]
            self.tracker = null_tracker

        if self.verbose == 1:
            print("Mass augment finished: %i -> %i" % (initial_m, self.m))
            print("  Replace is %a" % replace)


    def combine(self, new_data):
        """Combines new_data with the current Dataset class."""

        if self.data_X is None or self.data_y is None:
            raise ValueError("data_X or data_y not yet defined. Must convert to "
                             "numpy arrays first.")

        if not type(new_data).__name__ == 'Dataset':
            raise RuntimeError("New data is not of type core.Dataset.")

        if self.c != new_data.c:
            raise RuntimeError("Class mismatch: %i != %i" 
                               % (self.c, new_data.c))

        if self.n != new_data.n:
            raise RuntimeError("Number of features mismatch: %i != %i" 
                               % (self.n, new_data.n))

        X = self.data_X
        y = self.data_y
        tracker = self.tracker

        # stack them
        X = np.concatenate((X, new_data.data_X))
        y = np.concatenate((y, new_data.data_y))
        tracker = np.concatenate((tracker, new_data.tracker))
        previous_m = self.m
        self.m += new_data.m

        if self.verbose == 1:
            print("Data stacked: data expanded %i -> %i"
                  % (previous_m, self.m))

        self.data_X = X
        self.data_y = y
        self.tracker = tracker


    def print_info(self, n_each_class=False, n_mixed=True, n_unique=True, 
                   suppress=True):
    
        X = self.data_X
        y = self.data_y
        
        if type(X) == np.ndarray:
          
            # rough indicator of how many chemical environments exist in the data
            if n_each_class:
                counter = np.zeros(y.shape[1])
                for i in range(y.shape[0]):
                    for j in range(y.shape[1]):
                        if y[i][j] != 0:
                            counter[j] += 1
                total = np.sum(counter)
                if total == 0:
                    raise ValueError("Total = 0, which means something went wrong "
                                     "with an earlier attempt to trim the data.")
                np.set_printoptions(suppress=suppress)
                print(np.round(counter))
          
            # print the proportion of the data that contains multiple kinds of 
            # chemical environments
            n_mixed_counter = 0
            if n_mixed:
                for i in range(y.shape[0]):
                    for j in range(y.shape[1]):
                        if y[i][j] == 0 or y[i][j] == 1:
                            pass
                        else:
                            # this feature vector is mixed
                            n_mixed_counter += 1
                            break  # the inner loop
                print("The ratio of mixed/total is %i/%i" % (n_mixed_counter, self.m))
          
            if n_unique:
                print("The distribution of unique spectra is ", 
                      self.get_distribution_unique().astype(int))
              
        else:
            print(type(X))
            raise TypeError("Type not yet supported.")

        
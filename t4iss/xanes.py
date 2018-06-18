
import os
import pickle
import numpy as np
from pymatgen.analysis.local_env import VoronoiNN
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
import pymatgen as mg
from .core import mXANES

def read_xanes(path, absorption_specie, order='eof', skip_missing=False, 
               symprec=0.01,ang_tol=5):
               
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
                print('This structure does not have '
                      +absorption_specie+' in it. Skipping...')
                return [[],[]] 
            else:
                os.chdir(here)                
                raise ValueError('This structure does not have '
                                 +absorption_specie+
                                 ' element in it. Please check...')

        xanes = []
        for i,s in enumerate(sites):
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
                            print(fe+' is not available in '+path+'. Skipping...')
                            os.chdir(here)
                            return [[],[]] 
                        else:
                            os.chdir(here)                        
                            raise FileNotFoundError(fe+' is not available in '+path)                 

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
                            print(fo+' is not available in '+path+'. Skipping...')
                            os.chdir(here)
                            return [[],[]] 
                        else:
                            os.chdir(here)                        
                            raise FileNotFoundError(fo+' is not available in '+path)  
                                                        
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
                            print(fe+' is not available in '+path+'. Skipping...')
                            os.chdir(here)
                            return [[],[]] 
                        else:
                            os.chdir(here)                        
                            raise FileNotFoundError(fe+' is not available in '+path)                                      

                if order is 'f':                                    
                    if os.path.isdir(ff):
                        os.chdir(ff)  
                        pload = pickle.load(open('xanes.pkl', 'rb'))
                        xanes.append(pload)
                        os.chdir('..')                    
                    else:
                        if skip_missing:
                            print(ff+' is not available in '+path+'. Skipping...')
                            os.chdir(here)
                            return [[],[]] 
                        else:
                            os.chdir(here)                        
                            raise FileNotFoundError(ff+' is not available in '+path)                              

                if order is 'o':                                    
                    if os.path.isdir(fo):
                        os.chdir(fo)  
                        pload = pickle.load(open('xanes.pkl', 'rb'))
                        xanes.append(pload)
                        os.chdir('..')                    
                    else:
                        if skip_missing:
                            print(fo+' is not available in '+path+'. Skipping...')
                            os.chdir(here)
                            return [[],[]] 
                        else:
                            os.chdir(here)                        
                            raise FileNotFoundError(fo+' is not available in '+path)                                

                if order is 'e':                                    
                    if os.path.isdir(fe):
                        os.chdir(fe)  
                        pload = pickle.load(open('xanes.pkl', 'rb'))
                        xanes.append(pload)
                        os.chdir('..')                    
                    else:
                        if skip_missing:
                            print(fe+' is not available in '+path+'. Skipping...')
                            os.chdir(here)
                            return [[],[]] 
                        else:
                            os.chdir(here)                        
                            raise FileNotFoundError(fe+' is not available in '+path)                             
                        

        # get E limits
        minmaxs = []
        for i in xanes:
            minmaxs.append([min(i.E0),max(i.E0)])
        minmaxs = np.array(minmaxs)    
        irange = [max(minmaxs[:,0]),min(minmaxs[:,1])]   
        e_int = np.linspace(irange[0],irange[1], int((irange[1]-irange[0])/0.1)+1  )

        ave_xanes = e_int*0
        ave_vcn = 0
        site_xanes = []
        counter = 0
        
        # iterate over all spectra in the xanes list, which comes from a
        # single directory (e.g. mp-390)
        for i in xanes:
            i.transform(irange=irange,normalize=None,y0shift=None)
            ave_xanes += i.I*i.multiplicity
            site_onset = i.Eonset
            site_xanes.append(mXANES(data=[i.E,i.I],structure=i.structure,xanesid=i.xanesid,source=i.source,
                                     Eonset=site_onset, edge=i.edge,multiplicity=i.multiplicity))
            
            # if there already exists a vcn for this structure, mutiply it by
            # the number of structures that are symmetrically equivalent
            # divide out the total number of symmetrically equivalent
            # structures after
            if i.vcn:
                ave_vcn += i.vcn*i.multiplicity 

            # else, call VoronoiNN.get_cn to get the coordination number for
            # that particular site
            else:
                try:
                    nnfinder = VoronoiNN(cutoff=10,allow_pathological=True)
                    print(i.structure)
                    vcn = nnfinder.get_cn(i.structure[0], i.structure[1], use_weights=True)
                    ave_vcn += i.vcn*i.multiplicity 
                except Exception as exc:
                    print(exc)
                    print('warning: vcn info is missing for site. ave_vnc is not correct. ')
            counter += i.multiplicity
        
        # execute the averaging of the multiplicities
        ave_xanes = ave_xanes/counter
        ave_vcn = ave_vcn/counter
        ave_xanes = mXANES(data=[e_int,ave_xanes],Eonset=site_onset,structure=[struct,-1,ave_vcn],xanesid=i.xanesid,source=i.source,edge=i.edge)   
        os.chdir(here)
        return [ave_xanes,site_xanes]
    
    else:
        if skip_missing:
            os.chdir(here)             
            print(path+' is not available. Skipping...')
            return [[],[]]       
        else:
            os.chdir(here)            
            raise FileNotFoundError(path+' is not available.')

def main():
    pass

if __name__ == '__main__':
    main()
        
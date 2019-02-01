'''
This module contains the FluxGeneratorClass, which generates synthetic light curves from the charateristic spectral function. (details see Ryan+ 2014)

'''

import numpy as np
import h5py as h5
from InterpolatorClass import InterpolatorClass

class FluxGeneratorClass:
    '''
    FluxGeneratorClass generates synthetic light curves from the charateristic spectral function Table.
    '''
    
    ### Attributes
    TableInfo = None
    ''' Dictionary which contains Table information. '''
    
    def __init__(self, Table, LogTable=True, LogAxis=['tau']):
        '''
        Initialize FluxGeneratorClass.  
        
        Args:
            Table (str): directory to boosted fireball table
            LogTable (bool): whether Table is measured in log scale
            LogAxis (list of str): whether certain axis is measured in log scale
        '''
        self._Interpolator = InterpolatorClass(Table, LogTable=LogTable, LogAxis=LogAxis)
        self.TableInfo = self._Interpolator.GetTableInfo()

    ### Public Function
    def GetTaus(self, Times, P):
        '''
        Args:
            Times (Array): observational time in second.
            P (dictionary): contains all paramters {z, dL, E, n, p, epse, epsb, xiN, Eta0, GammaB,theta_obs}
        Return:
            float: scaled time
        '''

        f0 = (P['n']/P['E'])**(1./3.)/(1.+P['z'])
        Taus = f0 * Times 
        return Taus


    def GetTransformedValue(self, Taus, P):
        '''
        Apply the scale relations. sees Ryan+ 2014.

        Args:
            Taus (Array): scaled time.
            P (dictionary): contains all paramters {z, dL, E, n, p, epse, epsb, xiN, Eta0, GammaB,theta_obs} 
        Return:
            float, float, float: spectral function values corresponding to parameter dict P. 
        '''

        Position = np.array([[tau, P['Eta0'], P['GammaB'], P['theta_obs']]  for tau in Taus])
        f_peak, f_nu_c, f_nu_m = self._Interpolator.GetValue(Position)

        if np.isnan(f_peak[0]):
            return f_peak, f_nu_c, f_nu_m

        else:    
            f1 = (1+P['z'])/(P['dL']*P['dL'])*(P['p']-1)/(3.*P['p']-1.)*P['E']*P['n']**0.5*P['epsb']**0.5*P['xiN']
            f2 = 1./(1+P['z'])*P['E']**(-2./3.)*P['n']**(-5./6.)*P['epsb']**(-3./2.)
            f3 = 1./(1+P['z'])*(P['p']-2.)**2/(P['p']-1.)**2*P['n']**(1./2.)*P['epse']**2*P['epsb']**(1./2.)*P['xiN']**(-2)

            F_peak = f1*f_peak
            nu_c = f2*f_nu_c
            nu_m = f3*f_nu_m

            return F_peak, nu_c, nu_m
    
    
    def GetSpectral(self, Times, Freqs, P):
        '''
        Get synthetic light curve through interpolation. Spectral is constructed as power laws. (Sari+ 1998)

        Args:
            Times (Array): observational time in second
            Freqs (Array): frequencies. The length should be the same as Times
            P (dictionary): contains all paramters {z, dL, E, n, p, epse, epsb, xiN, Eta0, GammaB,theta_obs} 
        Return:
            Array: synthetic light curve
        '''
        Taus = self.GetTaus(Times, P)
        F_peak, nu_c, nu_m = self.GetTransformedValue(Taus, P)

        if np.isnan(F_peak[0]):
            return F_peak
        else:        
            IdxSlow = (nu_m < nu_c)
            IdxSlow1 = IdxSlow * (Freqs < nu_m)
            IdxSlow2 = IdxSlow * (Freqs >= nu_m) * (Freqs < nu_c)
            IdxSlow3 = IdxSlow * (Freqs >= nu_c)
            IdxFast1 = (~IdxSlow) * (Freqs < nu_c)
            IdxFast2 = (~IdxSlow) * (Freqs >= nu_c) * (Freqs < nu_m)
            IdxFast3 = (~IdxSlow) * (Freqs >= nu_m)

            inu_m = 1.0/nu_m
            inu_c = 1.0/nu_c

            p = P['p']
            output = np.zeros(len(Times))
            output[IdxSlow1] = np.power(Freqs[IdxSlow1] * inu_m[IdxSlow1], 1.0/3.0)
            output[IdxSlow2] = np.power(Freqs[IdxSlow2] * inu_m[IdxSlow2], 0.5-0.5*p)
            output[IdxSlow3] = np.power(nu_c[IdxSlow3] * inu_m[IdxSlow3], 0.5-0.5*p) * np.power(Freqs[IdxSlow3] * inu_c[IdxSlow3], -0.5*p)

            output[IdxFast1] = np.power(Freqs[IdxFast1] * inu_c[IdxFast1], 1.0/3.0)
            output[IdxFast2] = np.power(Freqs[IdxFast2] * inu_c[IdxFast2], -0.5)
            output[IdxFast3] = np.power(nu_m[IdxFast3] * inu_c[IdxFast3], -0.5) * np.power(Freqs[IdxFast3] * inu_m[IdxFast3], -0.5*p)

            Spectral = F_peak * output

            return Spectral

    def GetIntegratedFlux(self, Times, Freqs, P):
        '''
        Get synthetic light curve through interpolation. Spectral is constructed as power laws. (Sari+ 1998)

        Args:
            Times (Array): observational time in second
            Freqs (Array): frequencies. The length should be the same as Times
            P (dictionary): contains all paramters {z, dL, E, n, p, epse, epsb, xiN, Eta0, GammaB,theta_obs} 
        Return:
            Array: synthetic light curve (erg cm^-2 s^-1)
        '''

        Taus = self.GetTaus(Times, P)
        F_peak, nu_c, nu_m = self.GetTransformedValue(Taus, P)

        if np.isnan(F_peak[0]):
            return F_peak
        else:        
            LFreqs = Freqs[:,0]
            HFreqs = Freqs[:,1]

            IdxSlow = (nu_m < nu_c)
            IdxSlow1 = IdxSlow * (LFreqs < nu_m)
            IdxSlow2 = IdxSlow * (HFreqs >= nu_m) * (LFreqs < nu_c)
            IdxSlow3 = IdxSlow * (HFreqs >= nu_c)
            IdxFast1 = (~IdxSlow) * (LFreqs < nu_c)
            IdxFast2 = (~IdxSlow) * (HFreqs >= nu_c) * (LFreqs < nu_m)
            IdxFast3 = (~IdxSlow) * (HFreqs >= nu_m)

            inu_m = 1.0/nu_m
            inu_c = 1.0/nu_c

            p = P['p']
            output = np.zeros(len(Times))

            output[IdxSlow1] = 0.75 * nu_m[IdxSlow1] * (np.power(np.minimum(HFreqs[IdxSlow1],nu_m[IdxSlow1])*inu_m[IdxSlow1], 4.0/3.0) - np.power(LFreqs[IdxSlow1]*inu_m[IdxSlow1], 4.0/3.0))
            output[IdxSlow2] = 2.0/(3.0-p) * nu_m[IdxSlow2] * (np.power(np.minimum(HFreqs[IdxSlow2],nu_c[IdxSlow2])*inu_m[IdxSlow2], 1.5-0.5*p) - np.power(np.maximum(LFreqs[IdxSlow2],nu_m[IdxSlow2])*inu_m[IdxSlow2], 1.5-0.5*p))
            output[IdxSlow3] = 2.0/(2.0-p) * nu_c[IdxSlow3] * np.power(nu_c[IdxSlow3]*inu_m[IdxSlow3], 0.5-0.5*p) * (np.power(HFreqs[IdxSlow3]*inu_c[IdxSlow3], 1.0-0.5*p) - np.power(np.maximum(LFreqs[IdxSlow3],nu_c[IdxSlow3])*inu_c[IdxSlow3], 1.0-0.5*p))

            output[IdxFast1] = 0.75 * nu_c[IdxFast1] * (np.power(np.minimum(HFreqs[IdxFast1],nu_c[IdxFast1])*inu_c[IdxFast1], 4.0/3.0) - np.power(LFreqs[IdxFast1]*inu_c[IdxFast1], 4.0/3.0))
            output[IdxFast2] = 2.0 * nu_c[IdxFast2] * (np.sqrt(np.minimum(HFreqs[IdxFast2],nu_m[IdxFast2])*inu_c[IdxFast2]) - np.sqrt(np.maximum(LFreqs[IdxFast2],nu_c[IdxFast2])*inu_c[IdxFast2]))
            output[IdxFast3] = 2.0/(2.0-p) * nu_m[IdxFast3] / np.sqrt(nu_m[IdxFast3]*inu_c[IdxFast3]) * (np.power(HFreqs[IdxFast3]*inu_m[IdxFast3], 1.0-0.5*p) - np.power(np.maximum(LFreqs[IdxFast3],nu_m[IdxFast3])*inu_m[IdxFast3], 1.0-0.5*p))

            Integrated = (1.0e-26) * F_peak * output

        return Integrated
        
        

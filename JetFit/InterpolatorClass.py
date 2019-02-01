'''
This module contains the InterpolatorClass, which performs interpolation in the boosted fireball Table.

'''

import numpy as np
import h5py as h5

class InterpolatorClass:
    '''
    InterpolatorClass performs interpolation in the boosted fireball Table. 
    '''
    ### Private: h5 Table
    _Table = {}
    _Axis = {}
    
    ### Private: interpolation function
    _f_peak = None
    _f_nu_c = None
    _f_nu_m = None
    
    def __init__(self, Table, LogTable=True, LogAxis=['tau']):
        '''
        Initialize InterpolatorClass.
        
        Args:
            Table (str): directory to boosted fireball table
            LogTable (bool): whether Table is measured in log scale
            LogAxis (list of str): whether certain axis is measured in log scale
        '''
        self._LoadTable(Table)
        self._SetScale(LogTable=True,LogAxis=['tau'])
        self._GetInterpolator()
    
    ### Private Function
    def _LoadTable(self, Table):
        '''
        Load boosted fireball Table.

        Args:
            Table (str): directory to boosted fireball table
        '''
        Data = h5.File(Table, 'r')
        for key in Data.keys():
            if key in ['f_peak','f_nu_c','f_nu_m']:
                self._Table[key] = Data[key][...]
            else:
                self._Axis[key] = Data[key][...]

        ### Due to some reasons, we need to convert bytes to string. 
        self._Axis['Axis'] = np.array([x.decode("utf-8") for x in self._Axis['Axis']])
        Data.close()
    
    
    def _SetScale(self, LogTable=True,LogAxis=['tau']):
        '''
        Set proper scales to table and axis.

        Args:
            LogTable (bool): whether Table is measured in log scale
            LogAxis (list of str): whether certain axis is measured in log scale
        '''

        self.Info = self._Axis.copy()

        if 'LogAxis' not in self._Axis.keys():
            self._Axis['LogAxis'] = LogAxis
            for key in LogAxis:
                if key not in self._Axis.keys():
                    raise ValueError('could not find %s in Axis' %(key))
                else:
                    self._Axis[key] = np.log(self._Axis[key])

        if 'LogTable' not in self._Table.keys():
            for key in self._Table.keys():
                temp = np.ma.log(self._Table[key])
                self._Table[key] = temp.filled(-np.inf)
            self._Table['LogTable'] = True
    
    
    def _GetInterpolator(self):
        '''
        Use scipy.interpolate.RegularGridInterpolator to perform interpolation.
        '''
        from scipy.interpolate import RegularGridInterpolator

        Axes = [self._Axis[key] for key in self._Axis['Axis']]
        self._f_peak = RegularGridInterpolator(Axes, self._Table['f_peak'])
        self._f_nu_c = RegularGridInterpolator(Axes, self._Table['f_nu_c'])
        self._f_nu_m = RegularGridInterpolator(Axes, self._Table['f_nu_m'])
    
    ### Public Function
    def GetTableInfo(self):
        '''
        Get Table Information. 

        Return:
            dict
        '''
        return self.Info
    
    
    def GetValue(self, Position):
        '''
        Get the characteristic function values at the Position. 

        Args:
            Position (Array): (tau, Eta0, GammaB, theta_obs) (linear scale)
        Return:
            float, float, float: characteristic function values

        '''

        ScaledPosition = Position.copy()

        ### convert linear scale to log scale
        for key in self._Axis['LogAxis']:
            idx = np.where(self._Axis['Axis'] == key)[0][0]
            ScaledPosition[:,idx] = np.log(ScaledPosition[:,idx])

        ### When lorentz factor is low and observation time is ealry, there is no detection, which is represented by nans. 
        try:
            if self._Table['LogTable']:
                f_peak = np.exp(self._f_peak(ScaledPosition))
                f_nu_c = np.exp(self._f_nu_c(ScaledPosition))
                f_nu_m = np.exp(self._f_nu_m(ScaledPosition))
            else:
                np.seterr(all='ignore')
                f_peak = self._f_peak(ScaledPosition)
                f_nu_c = self._f_nu_c(ScaledPosition)
                f_nu_m = self._f_nu_m(ScaledPosition)
                np.seterr(all='raise')
            return f_peak, f_nu_c, f_nu_m
        except:
            Nans = [np.nan for x in range(len(Position))]
            f_peak, f_nu_c, f_nu_m = Nans, Nans, Nans
            return f_peak, f_nu_c, f_nu_m
    
    
    
    

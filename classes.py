# Holding all classes data for ETA
import numpy as np
import astropy.units as u
import astropy.constants as c


class Source:
    """
    Defining Source object: it has a flux one can observe from it at a certain
    wavelength, used to derive the total photon flux observed later on.
    """
    
    def __init__(self, wavelengths, flux_density_lambda = None, T_rad = None, 
                 BB=True):
        """
        wavelength: numpy ndarray of wavelength values, should have units.
        flux density_lambda: Flux density in wavelength units (cgs)
        T_rad: radiative temperature of source, used to calculate blackbody's
        flux density when BB is True. Does nothing otherwise. Default is None.
        BB: Boolean to state if source is a blackbody. Default is True.  
        
        """
        self.wavelengths = wavelengths.to(u.AA)
        self.__isBB__ = BB
        if BB = True:
            B = (2 * c.h * c.c**2/(wavelengths**5)) * (1/(np.exp(c.h * c.c/(wavelengths*c.k_B*T))-1))
            self.F_lambda = B.to(u.erg * u.s**-1 * u.cm**-2 * u.AA**-1)
        else:
            self.F_lambda = flux_density_lambda
        self.frequencies = self.wavelengths.to(u.Hz, equivalencies=u.spectral())   #what does equivalencies do?
        self.F_nu = ((wavelengths**2/c.c)*self.F_lambda).to(u.erg * u.s**-1 * u.cm**-2 * u.Hz**-1)
        
        
    
class Atmosphere:
    
    def __init__(self, a_lambda):
        
        self.a_lambda = a_lambda
        
class Telescope:
    
    def __init__(self, t_lambda):
        
        self.t_lambda = t_lambda
        elif name == '' :
            self.diameter=diameter
            self.mirrors=mirrors
        
   def area(self):'''
        gives telescope area for diameter (self.diameter) in sq cm'''
            return (np.pi*(self.diameter/2)**2*(1-self.eps**2).to(u.cm**2)
                    
   def throughput(self,wavelength): '''
        gives throughput (efficiency) for the telescope '''
             t=np.ones(len(wavelength))
        for mir in self.mirrors: 
            if type(mir) is float: 
                    t *= mir
            else:
                    tmp = Mirror(mir)
                    t *= tmp.reflectivity(wavelength)
                    
        return t
                    
class Mirror: '''
    class for a mirror with a coating name, which quantifies relfectivity, by changing the constant m_lambda value '''
    def __int__(self,type,m_lambda=0.9): 
        self.type= type
        self.const= m_lambda
    def reflectivity(self,wavelength): 
        self.type == 'const' : 

return np.ones(len(wavelength)) + self.const
        
                    
class Instrument:'''
    describing the instrumental uncertainties'''
    
    def __init__(self, name='', i_lambda):
        
        self.i_lambda = i_lambda
        self.name=name
        self.detector = Detector(efficiency=1.)
                    
    def throughput(self,wavelength): 
        if self.name == '' :
            return np.ones(len(wavelength))*self.i_lambda
        else: 
            raise ValueError('need to implement instrument self.name == "" ')
                                
                    
 
class Filter:
    
    def __init__(self, b_lambda):
        
        self.b_lambda = b_lambda
  
def filter(self,wavelength,filter='',cent=5500*u.AA, wid=850*u.AA, trans=0.8): 
                    '''I'm not sure what this trans function is- is this the same as f_lambda in this case? 
                    the uncertainty constant associated with the filter? what do we set it to here, 
                    if we define it in the if loop? '''
    if filter=='':
        out=np.zeros(len(wavelength))
        out[np.where((wave>cent-wid/2)&(wave<cent+wid/2))] = trans
        return out

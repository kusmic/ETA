# Holding all classes data for ETA
import numpy as np
import astropy.units as u
import astropy.constants as c

VEGA_T = 9600 * u.K
VEGA_V = 5500 * u.AA
VEGA_BBl = 3.63e-9 * u.erg * u.s**-1 * u.cm**-2 * u.AA**-1
VEGA_BBn = ((VEGA_V**2/c.c)*VEGA_BBl).to(u.erg * u.s**-1 * u.cm**-2 * u.Hz**-1)
VEGA_wl_ENG = (c.h * c.c / VEGA_V).to(u.erg)
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
        if BB == True:
            exponent = c.h * c.c / (wavelengths * c.k_B * T_rad)
            exponent = exponent.decompose()
            exponential = np.exp( exponent.value )
            B = (2 * c.h * c.c**2/(wavelengths**5)) * (1/ (exponential-1))
            self.F_lambda = B.to(u.erg * u.s**-1 * u.cm**-2 * u.AA**-1)
        else:
            self.F_lambda = flux_density_lambda.to(u.erg * u.s**-1 * u.cm**-2 * u.AA**-1)/VEGA_BBl.to(u.erg * u.s**-1 * u.cm**-2 * u.AA**-1)
        self.frequencies = self.wavelengths.to(u.Hz, equivalencies=u.spectral())
        self.F_nu = ((wavelengths**2/c.c)*(self.F_lambda)).to(u.erg * u.s**-1 * u.cm**-2 * u.Hz**-1)
        
        self.f_nu = self.F_nu / VEGA_wl_ENG
        self.f_lambda = self.F_lambda / VEGA_wl_ENG
    
class Atmosphere:
    
    def __init__(self, a_lambda):
        
        self.a_lambda = a_lambda
        
class Telescope:
    
    def __init__(self, t_lambda, diameter=None, mirrors=None, name=None, eps=0):
        
        self.t_lambda = t_lambda
        if name == None:
            self.diameter=diameter
            self.mirrors=mirrors
            self.eps=eps
            
    def area(self):
       ''' gives telescope area for diameter (self.diameter) in sq cm'''
       return (np.pi*(self.diameter/2)**2*(1-self.eps**2).to(u.cm**2))
                    
    def throughput(self,wavelength): 
       ''' gives throughput (efficiency) for the telescope '''
       t=np.ones(len(wavelength))
       for mir in self.mirrors: 
           if type(mir) is float: 
               t *= mir
           else:
               tmp = Mirror(mir)
               t *= tmp.reflectivity(wavelength)
                    
       return t
                    
class Mirror: 
    '''
    class for a mirror with a coating name, which quantifies relfectivity, by 
    changing the constant m_lambda value '''
    def __init__(self,typename="const",m_lambda=0.9): 
        self.typename = typename
        self.const = m_lambda
    def reflectivity(self,wavelength): 
        if self.typename == 'const':
            return np.ones(len(wavelength))*self.const
        
                    
class Instrument:
    '''
    describing the instrumental uncertainties'''
    
    def __init__(self, i_lambda, sigma_rn = None, pixel_scale = None, name=''):
        
        self.i_lambda = i_lambda
        self.name=name
        self.sigma_rn = sigma_rn
                    
    def throughput(self,wavelength): 
        if self.name == '' :
            return np.ones(len(wavelength))*self.i_lambda
        else: 
            raise ValueError('need to implement instrument self.name == "" ')
                    
 
class Filter:
    
    def __init__(self, b_lambda):
        
        self.b_lambda = b_lambda
  

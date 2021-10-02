# Holding all classes data for ETA
import numpy as np
import astropy.units as u
import astropy.constants as c
from astropy.modeling.models import BlackBody1D

# parameters to normalize to Vega
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
    
    def __init__(self, wavelengths, mag=None, bands=None, flux_density_lambda=None, 
                 T_rad = None, BB=True, Vega_norm=False):
        """
        -wavelength (float 1darray): numpy ndarray of wavelength values, should have units.
        -flux density_lambda (float 1darray): Flux density in wavelength units
        (cgs)
        -T_rad (float): radiative temperature of source, used to calculate 
        blackbody's flux density when BB is True. Does nothing otherwise. 
        Default is None.
        -BB (bool): Boolean to state if source is a blackbody. Default is True. 
        -mag (float 1darray, float): apparent magnitude of source if flux 
        density to not be used. 
        -bands (string 1darray, string): Bands classifying mag
        -Vega_norm (bool): Boolean to set whether to normalize flux density to
        observed Vega values. Default is False.
        
        """
        self.wavelengths = wavelengths.to(u.AA)
        self.frequencies = self.wavelengths.to(u.Hz, equivalencies=u.spectral())
        # Samir found out you can call class method before its definition
        if mag==None:
            self.f_lambda, self.F_lambda, self.f_nu, self.F_nu = self.__setFluxDensities__(
                    self.wavelengths, self.frequencies, BB, T_rad, Vega_norm)
        else:
            self.mag = mag
            self.bands = bands
        self.__isBB__ = BB
        self.__isVegaNorm__ = Vega_norm        
        
    def __setFluxDensities__(self, lam, nu, isBB=False, T_rad=None, VegaNorm=False):
        if isBB == True:
            bb_model = BlackBody1D(temperature=T_rad)
            # energy flux density in frequency
            F_nu = bb_model(self.wavelengths)
            # energy flux density in wavelength
            F_lambda = ((self.frequencies**2/c.c)*F_nu).to(u.erg * u.s**-1 * u.cm**-2 * u.AA**-1,
                       equivalencies=u.spectral())
                
        else:
            F_lambda = self.flux_density_lambda.to(u.erg * u.s**-1 * u.cm**-2 * u.AA**-1)
            F_nu = ((self.wavelengths**2/c.c)*(F_lambda)).to(u.erg * u.s**-1 * u.cm**-2 * u.Hz**-1)
            
        if VegaNorm == True:
            index_5500 = np.where(self.wavelengths.value==5500)
            flux_5500 = self.F_lambda[index_5500]
            F_lambda = self.F_lambda / flux_5500
            F_lambda = self.F_lambda * VEGA_BBl 
            F_nu = ((self.wavelengths**2/c.c)*(F_lambda)).to(u.erg * u.s**-1 * u.cm**-2 * u.Hz**-1)
        
        # photon flux density in frequency
        f_nu = F_nu / ((c.h*self.frequencies).to(u.erg))
        # photon flux density in wavelength
        f_lambda = F_lambda / ((c.c * c.h / self.wavelengths).to(u.erg))
            
        return(f_lambda, F_lambda, f_nu, F_nu)
    
class Sky:
    
    def __init__(self, sky_absorp=None, sky_emiss=None, a_lambda=None):
        
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
       return (np.pi*(self.diameter/2)**2*(1-self.eps**2)).to(u.cm**2)
                    
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
    
    def __init__(self, i_lambda, datatype="imaging", dispersion=None, 
                 sigma_rn = None, pixel_scale = None, name=''):
        
        self.i_lambda = i_lambda
        self.name=name

        if datatype == "imaging":
            self.sigma_rn = sigma_rn
        elif datatype == "spectra":
            self.dispersion = dispersion
        self.pixel_scale = pixel_scale
                    
    def throughput(self,wavelength): 
        if self.name == '' :
            return np.ones(len(wavelength))*self.i_lambda
        else: 
            raise ValueError('need to implement instrument self.name == "" ')
                    
 
class Filter:
    
    def __init__(self, b_lambda=None), banddata=None, bandname=None):
        
        self.b_lambda = b_lambda
        
class Observation:
    
    def __init__(self, source, sky, telescope, mirror, instrument, band):
        pass
    
    def __set_source__(self, source, sky, telescope, mirror, instrument, band):
        pass
  

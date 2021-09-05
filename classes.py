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
        self.frequencies = self.wavelengths.to(u.Hz, equivalencies=u.spectral())
        self.F_nu = ((wavelengths**2/c.c)*self.F_lambda).to(u.erg * u.s**-1 * u.cm**-2 * u.Hz**-1)
        
    
class Atmosphere:
    
    def __init__(self, a_lambda):
        
        self.a_lambda = a_lambda
        
class Telescope:
    
    def __init__(self, t_lambda):
        
        self.t_lambda = t_lambda
        
class Instrument:
    
    def __init__(self, i_lambda):
        
        self.i_lambda = i_lambda
        
class Filter:
    
    def __init__(self, b_lambda):
        
        self.b_lambda = b_lambda
        

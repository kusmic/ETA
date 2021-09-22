import numpy as np
import classes as cl
import astropy.units as u
from routines import get_SNR, get_exposure_time

if __name__ == "__main__":
    dlam = 1
    lam_arr = np.arange(4000, 9000, dlam) * u.AA
    test_source = cl.Source(lam_arr, BB=True, T_rad = 6000*u.K)
    test_atmos = cl.Atmosphere(np.ones(lam_arr.size)*0.5)
    part_thruput = np.power(0.3, 1/3)
    test_tele = cl.Telescope(part_thruput, diameter=0.6*u.m, mirrors=[part_thruput])
    test_mirror = cl.Mirror(np.ones(lam_arr.size)*part_thruput)
    test_inst = cl.Instrument(np.ones(lam_arr.size)*part_thruput, sigma_rn=10, 
                              pixel_scale=20 / u.arcsec)
    test_band = cl.Filter(np.ones(lam_arr.size)*0.8)
    test_BgSB = 1000 / u.arcsec**2 / u.s
    test_exposure_t = 100 * u.s
    test_aper_radius = 3600 * 5 * u.arcsec
    test_SNR = 10
    
    obs_SNR = get_SNR(test_exposure_t, test_BgSB, test_aper_radius, test_source,
                      test_atmos, test_tele, test_mirror, test_inst, test_band)
    obs_exp_t = get_exposure_time(test_SNR, test_BgSB, test_aper_radius, test_source,
                      test_atmos, test_tele, test_mirror, test_inst, test_band)
    print("SNR: ", obs_SNR)
    print("Exposure time: ", obs_exp_t)

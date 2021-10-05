import numpy as np
import classes as cl
import astropy.units as u

if __name__ == "__main__":
    dlam = 1
    lam_arr = np.arange(4000, 9000, dlam) * u.AA
    test_BgSB = 10 / u.arcsec**2 / u.s

    test_source = cl.Source(mag=5, bands="V")
    test_atmos = cl.Sky(a_lambda=np.ones(lam_arr.size)*0.5, B=test_BgSB)
    part_thruput = np.power(0.3, 1/3)
    test_tele = cl.Telescope(part_thruput, diameter=0.6*u.m, mirrors=[part_thruput])
    test_mirror = cl.Mirror(np.ones(lam_arr.size)*part_thruput)
    test_inst = cl.Instrument(np.ones(lam_arr.size)*part_thruput, sigma_rn=10, 
                              pixel_scale=2 / u.arcsec, zpt=-22)
    test_band = cl.Filter(np.ones(lam_arr.size)*0.8)
    test_exposure_t = 100 * u.s
    test_aper_radius = 2 * u.arcsec
    test_SNR = 100
    
    test_obs = cl.Observation(test_source, test_atmos, test_tele, test_mirror,
                              test_inst, test_band, test_aper_radius)
    obs_SNR = test_obs.get_SNR(test_exposure_t)
    obs_exp_t = test_obs.get_exposure_time(test_SNR)
    print("SNR: ", obs_SNR)
    print("Exposure time: ", obs_exp_t)

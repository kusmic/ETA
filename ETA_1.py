import numpy as np
import classes as cl
import astropy.units as u

dlam = 1
lam_arr = np.arange(4000, 10000, dlam) * u.AA
test_source = cl.Source(lam_arr, BB=True, T_rad = 5500*u.K)
test_atmos = cl.Atmosphere(np.ones(lam_arr.size)*0.8)
part_thruput = np.power(0.5, 1/3)
test_tele = cl.Telescope(part_thruput, diameter=0.6*u.m, mirrors=[part_thruput])
test_mirror = cl.Mirror(np.ones(lam_arr.size)*part_thruput)
test_inst = cl.Instrument(np.ones(lam_arr.size)*part_thruput)
test_band = cl.Filter(np.ones(lam_arr.size)*0.8)

f = test_source.F_lambda.value * 1005 * u.s**-1 * u.cm**-2 * test_atmos.a_lambda * test_tele.t_lambda * test_mirror.const * test_inst.i_lambda * test_band.b_lambda
S_t = test_tele.diameter * np.trapz(f, lam_arr).decompose()
print(S_t)
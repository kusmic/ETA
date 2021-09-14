#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 14 14:00:51 2021

@author: skusmic

This holds all functions to do our work
For now, has functions to get SNR from exposure time, andget exposure time for
SNR
"""

import numpy as np
import classes

def get_SNR(t, BG, aperture_radius, source, atmosphere, telescope, mirror, 
            instrument, band, limit="none"):
    """
    Gets the SNR from given exposure time t and each class.
    BG: background's surface brightness rate [photons/cm2/s]
    aperture_radius: radius of seeing area
    limit: is this statistically limited? It can be signal-, background-, or 
    readout-limited. Each denoted as "signal", "bg", and "readout". Otherwise,
    it is "none" and all three are taken into consideration. Default:"none"
    
    Out:
        SNR - float denoting the signal to noise ratio.
    """
    f = source.f_lambda * atmosphere.a_lambda * telescope.t_lambda * mirror.const * instrument.i_lambda * band.b_lambda
    S = t * telescope.area() * np.trapz(f, source.wavelengths)
    aper_area = aperture_radius**2 * np.pi
    N_pix = 2*(2*aperature_radius * instrument.pixel_scale)
    Nsig = N_pix * instrument.sigma_rn**2
    B = aper_area * BG * t
    if limit == "signal":
        SNR = S / np.sqrt(S)
    elif limit == "bg":
        SNR = S / np.sqrt(B)
    elif limit == "readout":
        SNR = S / np.sqrt(Nsig)
    else:
        SNR = S / np.sqrt(S + B + Nsig)
        
    return SNR
    
def get_exposure_time(SNR, BG, source, atmosphere, telescope, mirror, instrument, 
                      band, limit="none"):
    f = source.f_lambda * atmosphere.a_lambda * telescope.t_lambda * mirror.const * instrument.i_lambda * band.b_lambda
    SpT = telescope.area() * np.trapz(f, source.wavelengths)
    aper_area = aperture_radius**2 * np.pi
    N_pix = 2*(2*aperature_radius * instrument.pixel_scale)
    Nsig = N_pix * instrument.sigma_rn**2
    ABpT = aper_area * BG
    if limit == "signal":
        t = SNR**2 * SpT / SpT**2
    elif limit == "bg":
        SNR = SNR**2 * ABpT / SpT**2
    elif limit == "readout":
        SNR = SNR*np.sqrt(Npix) / SpT
    else:
        a = SpT**2
        b = -SNR**2* (SpT + ABpT)
        c = -SNR**2 * Nsig
        t = (-b + np.sqrt(b**2 - 4 * a * c)) / (2*a)
        
    return t
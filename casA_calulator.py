import numpy as np
import argparse
### a simple calculator for calculating CasA flux using measurements

def spfd_converter(freq, freq_ref, flux_ref, spectral_index):
    flux = flux_ref * (freq/freq_ref)**spectral_index
    return flux
def casA_flux(freq):
    freq_ref = 1400 ### MHz
    flux_ref = 2720 ### jy
    spectral_index = -0.8
    flux = spfd_converter(freq, freq_ref, flux_ref, spectral_index)
    return flux

def pfd_integrated(band,function=casA_flux):
    ### band is an array of frequencies
    spfd_values = [function(freq) for freq in band]
    pfd = np.trapz(spfd_values, band) ### this allows for inequal spacing
    return pfd

def __main__():
    parser = argparse.ArgumentParser(description='Calculate the flux of CasA at a given frequency')
    parser.add_argument('-f','--freq', type=float, help='The frequency at which to calculate the flux in MHz, centre frequency')
    parser.add_argument('-i','--integrate', action='store_true', help='Calculate the integrated flux over a band')
    parser.add_argument('--band', type=float, nargs='+', help='The frequency band over which to calculate the integrated flux')
    args = parser.parse_args()
    if args.integrate:
        pfd = pfd_integrated(args.band)
        print(f'The integrated flux over the band {args.band} is {pfd:.2f} Jy')
    else:
        flux = casA_flux(args.freq)
        print(f'The flux at {args.freq} MHz is {flux:.2f} Jy')

if __name__ == '__main__':
    __main__()

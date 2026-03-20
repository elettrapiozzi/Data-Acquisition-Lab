import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import os
import matplotlib.pyplot as plt

# DATA LOADING
# fits files

base_dir = os.path.dirname(__file__)
science_path = os.path.join(base_dir, "final cgs flux", "Coadded_i_CGS.fits")
science = fits.open(science_path)
science_data = science[0].data

#Calibration constants for the 4 images 
calibration_constants = [1.043623281437049e-16, 1.0189551965013689e-16, 1.0234634960834272e-16]
aperture_radius = 600
texp = 600
Nexp = 3 #Number of coadded images

#galaxy quantities (needed for the masking)
x_c, y_c = 2386, 1589   # center of the galaxy
r = 600                 # radius in pixel 

#QUANTITIES NEEDED
C_medium = np.mean(calibration_constants) #Medium calibration constant
Npix = np.pi*aperture_radius**2 #number of pixel

# SKY ERROR

# Generating Mask
# We have to mask the central galaxy 
y, x = np.ogrid[:science_data.shape[0], :science_data.shape[1]]
mask = (x - x_c)**2 + (y - y_c)**2 <= r**2
aperture_flux = np.sum(science_data[mask])


# Sigma Clipping
print("\n")
print("Computing the error on the background with SIGMA CLIPPING...")

mean_sky, median_sky, stddev_sky = sigma_clipped_stats(
     science_data, 
     mask=mask
 )

sigma_sky_pixel = stddev_sky
sigma_sky = sigma_sky_pixel*Npix**(1/2)

print(f"Sigma of the sky (per pixel)", sigma_sky_pixel)
print(f'Sigma del cielo tot: ', sigma_sky) 


# SOURCE ERROR
print("\n")
print(f"Computing the error on the source...")
# sigma_source_pixel = (C_medium/(texp*Npix*Nexp))**(1/2)*(aperture_flux)**(1/2)

# sigma_source = sigma_source_pixel*Npix**(1/2)
factor = C_medium / (texp * Nexp)
sigma_source = np.sqrt(factor * aperture_flux)

# print(f'Sigma della sorgente (un pixel)', sigma_source_pixel)
print(f'Sigma della sorgente tot: ', sigma_source)

final_error = (sigma_sky**2 + sigma_source**2)**(1/2)

print(f'Final Error : ' , final_error)
print(f"\n")
print(f"Flux of the galaxy: ", aperture_flux, f"±", final_error)





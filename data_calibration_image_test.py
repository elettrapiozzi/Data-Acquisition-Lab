import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt

def Master(folder):
  array = None
  length = os.listdir(folder)

  for file in length:
    full_file_path = os.path.join(folder, file)
    s1 = fits.open(full_file_path)
    data = s1[0].data

    # Initialize on first file
    if array is None:
      array = np.zeros_like(data, dtype=float)

    # Add arrays together
    array += data

  master = array/len(length)


  return master



base_dir = os.path.dirname(__file__)  # cartella dove si trova il file .py

bias_orig = os.path.join(base_dir, "BIAS (25-09)")
flat_orig = os.path.join(base_dir, "FLAT (30-09)","g filter")
dark_orig = os.path.join(base_dir, "DARK (25-09)")

bias = Master(bias_orig)
dark = Master(dark_orig)
flat = Master(flat_orig)
flat_norm = (flat - dark) / np.mean(flat)

print(flat_norm)


plt.hist(bias.flatten(), bins=1000)
plt.xlim(158, 165)
plt.title('Bias Distribution')
plt.show()
plt.imshow(bias, clim=(160.27, 161.85))
plt.title('BIAS')
plt.show()

#plt.hist(dark.flatten(), bins=300000)
#plt.xlim(160, 168)
#plt.title('Dark Current Distribution')
#plt.show()
#plt.imshow(dark, clim=(160, 168))
#plt.title('DARK CURRENT')
#plt.show()

#plt.hist(flat_norm.flatten(), bins=500)
#plt.xlim(0.9, 1.08)
#plt.title('Flat Distribution')
#plt.show()
#plt.imshow(flat_norm, clim=(0.95, 1.05))
#plt.title('FLAT')
#plt.show()

gal = os.path.join(base_dir, "LIGHT (29-09 e 6-10)", "g filter", "2025-09-29_18-12-19_sci_NGC6946_G_exp300.00_0000.fits")
raw = fits.open(gal)
gain = 0.25
exptime =raw[0].header['EXPTIME']

gal_data = raw[0].data   # i dati dell’immagine FITS


#plt.hist(gal_data.flatten(), bins=5000)
#plt.xlim(3000, 4400)
#plt.title('Raw Data Distribution')
#plt.show()

vmin_raw = np.percentile(gal_data, 1)
vmax_raw = np.percentile(gal_data, 99)

plt.imshow(gal_data,vmin = vmin_raw, vmax = vmax_raw, cmap='viridis')
plt.colorbar()
plt.title('Raw NGC 6946')
plt.show()

science = ( gal_data - (bias + dark) ) / flat_norm /exptime * gain

vmin_val = np.percentile(science, 1)
vmax_val = np.percentile(science, 99)

#plt.hist(science.flatten(), bins=8000)
#plt.xlim(2.75, 3.15)
#plt.title('Science Frame Data Distribution')
#plt.show()
plt.imshow(science, vmin=vmin_val, vmax=vmax_val, cmap='viridis')
plt.colorbar()
plt.title('Science Frame of NGC 6946')
plt.show()


fits.writeto('calibrated_NGC6946_I_exp600.00_0002.fits', science, header = raw[0].header, overwrite=True)


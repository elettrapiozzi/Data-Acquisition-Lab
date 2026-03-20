import numpy as np
from astropy.io import fits
import os
import matplotlib.pyplot as plt
import sep

base_dir = os.path.dirname(__file__)
science = os.path.join(base_dir, "Data Reduction in CGS", "g filter", "calibrated_NGC6946_G_exp300.00_0000_cgs.fits") #!!!!!!!!!!!! SCEGLIERE FILE !!!!!!!!!!!
science = fits.open(science)
science_data = science[0].data

# plt.hist(science_data.flatten(), bins=20000)
# plt.xlim(2.75, 3.15)
# plt.title('Science Frame Data Distribution')
# plt.show()

mean_science = np.median(science_data)
#print('Mean Science Frame Value:', mean_science)
science_shifted = science_data - mean_science

# plt.hist(science_shifted.flatten(), bins=10000)
# plt.xlim(-0.25, 0.25)
# plt.title('Shifted Science Frame Data Distribution')
# plt.show()

vmin = np.percentile(science_shifted, 1)    # 1° percentile (nero)
vmax = np.percentile(science_shifted, 99)   # 99° percentile (bianco)

plt.imshow(science_shifted, vmin=vmin, vmax =vmax, cmap='ocean')
plt.colorbar()
plt.title('Shifted Science Frame of NGC 6946 (g filter)')
plt.show()

import numpy as np

# Dati di esempio
x_c, y_c = 2350, 1670  # centro galassia
r = 600                 # raggio in pixel

#cutting image

y_start = 200   # Taglio via 500 pixel dal basso
y_end   = 2800  # Taglio via tutto sopra il pixel 3500
x_start = 800   # Taglio via 500 pixel da sinistra
x_end   = 4000  # Taglio via tutto a destra del pixel 3500

science_cropped = science_shifted[y_start:y_end, x_start:x_end].copy()

vmin1 = np.percentile(science_cropped, 1)    # 1° percentile (nero)
vmax1 = np.percentile(science_cropped, 99)   # 99° percentile (bianco)

plt.imshow(science_cropped, vmin=vmin1, vmax =vmax1, cmap='viridis')
plt.colorbar()
plt.title('Shifted and Cropped Science Frame of NGC 6946 - 18:49')
plt.show()

#nuovo centro della galassia
x_c_crop = x_c - x_start
y_c_crop = y_c - y_start

# y, x = np.ogrid[:science_shifted.shape[0], :science_shifted.shape[1]]
# mask = (x - x_c)**2 + (y - y_c)**2 <= r**2

y, x = np.ogrid[:science_cropped.shape[0], :science_cropped.shape[1]]
mask = (x - x_c_crop)**2 + (y - y_c_crop)**2 <= r**2


fig, ax = plt.subplots()
im = ax.imshow(science_cropped, cmap='viridis', vmin=vmin1, vmax =vmax1) # , origin='lower'
plt.colorbar(im, ax=ax)
circle = plt.Circle((x_c_crop, y_c_crop), r, color='red', fill=False, lw=0.8)
ax.add_patch(circle)
plt.show()

bkg = sep.Background(science_cropped, mask=mask, bw=30, bh=30, fw=10, fh=10)  #solo per filtro i
#bkg = sep.Background(science_cropped, mask = mask, bw=72, bh=72, fw=20, fh=20) #mask, maskthresh, fthresh
# evaluate background as 2-d array, same size as original image
bkg_image = bkg.back()

plt.hist(bkg_image.flatten(), bins=500)
#plt.xlim(-0.05, 0.03)
plt.show()

# show the background
plt.imshow(bkg_image, vmin=vmin1, vmax=vmax1, cmap='gray') #   interpolation='nearest', origin='lower', 
plt.title('Sky Background of NGC 6946 (i filter)')
plt.colorbar()
plt.show()

# subtract the background
data_sub = science_cropped - bkg_image
#objects = sep.extract(data_sub, 1.5, err=bkg.globalrms)
plt.imshow(data_sub, vmin=vmin1, vmax =vmax1, cmap='viridis')
plt.title('Sky-Subtracted Science Frame of NGC 6946 (g filter)')
plt.colorbar()
plt.show()

#fits.writeto('skysub_NGC6946_I_exp600.00_0002_cgs.fits', data_sub, header=science[0].header, overwrite=True)


# from matplotlib.patches import Ellipse

# # plot background-subtracted image
# fig, ax = plt.subplots()
# m, s = np.mean(data_sub), np.std(data_sub)
# im = ax.imshow(data_sub, interpolation='nearest', cmap='gray',
#                vmin=m-s, vmax=m+s, origin='lower')

# # plot an ellipse for each object
# for i in range(len(objects)):
#     e = Ellipse(xy=(objects['x'][i], objects['y'][i]),
#                 width=6*objects['a'][i],
#                 height=6*objects['b'][i],
#                 angle=objects['theta'][i] * 180. / np.pi)
#     e.set_facecolor('none')
#     e.set_edgecolor('red')
#     ax.add_artist(e)
# plt.show()

# flux, fluxerr, flag = sep.sum_circle(data_sub, objects['x'], objects['y'],
#                                      3.0, err=bkg.globalrms, gain=1.0) ### radius 3




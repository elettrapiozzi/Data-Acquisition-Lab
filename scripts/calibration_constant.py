
from astropy.io import fits
import matplotlib.pyplot as plt
from astropy.table import Table
from scipy.integrate import quad
import numpy as np
from photutils.aperture import CircularAperture

from astropy import units as u
from astropy.coordinates import SkyCoord
from photutils.aperture import SkyCircularAperture
from photutils.aperture import aperture_photometry



hdul = fits.open("/home/elettra/Scrivania/Uni/Data Acquisition NUOVA/Dati e Codici/GAIA3_2194488599022078848.fits")

flux = hdul[0].data
wave = hdul[1].data

header_flux = hdul[0].header
header_wave = hdul[1].header

filter_r = Table.read("/home/elettra/Scrivania/Uni/Data Acquisition NUOVA/Dati e Codici/r_filter.txt", format="ascii.basic", delimiter=' ', guess=False)
wavelengths = filter_r['Wavelength']
transqe = filter_r['Transmission(T*QE)']

transqe_new = np.interp(wave, wavelengths, transqe)

numerator = np.trapezoid(flux * transqe_new, wave)
denominator = np.trapezoid(transqe_new, wave)

ratio = numerator / denominator

#######################################################################

sci_orig = fits.open("/home/elettra/Scrivania/Uni/Data Acquisition NUOVA/Dati e Codici/Data Reduction/r filter/calibrated_NGC6946_R_exp300.00_0002.fits")
sci = sci_orig[0].data


positions = [(4154, 1221)] # from DS9 IN PIXELS

radii = [8, 10, 20] # arbitrary
apertures = [CircularAperture(positions, r=r) for r in radii]
phot_table = aperture_photometry(sci, apertures)
for col in phot_table.colnames:
    phot_table[col].info.format = '%.8g'  # for consistent table output
print(phot_table)

rad0 = phot_table['aperture_sum_0'][0]
rad1 = phot_table['aperture_sum_1'][0]
rad2 = phot_table['aperture_sum_2'][0]

test_radius = []
test_flux = []

for i in range(1,15):
    radi_inner = 15
    radi_outer = 25
    radi_star  = i
    
    aperture_inner = CircularAperture(positions, r= radi_inner)
    aperture_outer = CircularAperture(positions, r= radi_outer)
    aperture_star  = CircularAperture(positions, r= radi_star)

    phot_table_inner = aperture_photometry(sci, aperture_inner)
    phot_table_outer = aperture_photometry(sci, aperture_outer)
    phot_table_star  = aperture_photometry(sci, aperture_star)

    flux_inner = phot_table_inner['aperture_sum'][0]
    flux_outer = phot_table_outer['aperture_sum'][0]
    flux_star  = phot_table_star['aperture_sum'][0]
    
    f_io = flux_outer - flux_inner
    f_bg_pix_io = f_io / (np.pi * (radi_outer**2 - radi_inner**2))
    flux_star_cal = flux_star - (f_bg_pix_io * np.pi * radi_star**2)
    
    

    test_radius.append(radi_star)
    test_flux.append(flux_star_cal)
    
    # Ultimate Goal is to Plot Flux of the Star
    
plt.plot(test_radius, test_flux)
plt.scatter(test_radius, test_flux)
plt.xlabel('Radius')
plt.ylabel('Flux (electrons)')
plt.show()

calibration_constant = ratio / test_flux[7]

'''
CALIBRATION CONSTANTS
Im29: 5.819606784372549e-17
Im35: 5.812683977602639e-17
Im40: 5.799700283286153e-17
Im49: 5.784717880552469e-17
'''

print("Calibration Constant")
print(calibration_constant)

cal_data = sci * calibration_constant

# plt.imshow(sci, clim=(2.85, 3.02), cmap='viridis')
# plt.colorbar()
# plt.title("Science Data for NGC6946 - 18:29")
# plt.show()

# plt.hist(cal_data.flatten(), bins=5000)
# plt.xlim(0.15e-15, 0.19e-15)
# plt.title('Calibrated Frame Data Distribution')
# plt.show()

# plt.imshow(cal_data, clim=(1.67e-16, 1.74e-16), cmap='viridis')
# plt.colorbar()
# plt.title("Calibrated Data for NGC6946 - 18:29")
# plt.show()

#fits.writeto('calconst_18-29.fits', cal_data, header=sci[0].header, overwrite=True)

import numpy as np
import matplotlib.pyplot as plt

def aperE(im, col, row, rad1, rad2, ir1, ir2, or1, or2, Kccd, saturation=np.inf):
    """Original code by Professor Alberto Bolatto, edited by Alyssa Pagan, and
    translated to Python by ChongChong He, further edited by Orion Guiffreda.

    Before using aperE.m, rotate your image using imrotate(im,angle) so the
    major axis of your object is perpendicular or parallel to your x or y axis.
    
    APER(im,col,row,rad1,rad1,ir1,ir2,or1,or2,Kccd) Do aperture photometry of image "im"
    for a star, galaxy or nebula centered at the "row,col" coordinates, For an ellipse 
    with a major and minor axis of "rad1,rad2" and an inner sky ellipse with a 
    major and minor axis of (ir1,ir2)and outer sky ellipse of "or1,or2" with CCD
    gain of Kccd ADU/electron. Optionally, a 11th parameter can be passed
    with the saturation value for the CCD.
    """

    a, b = im.shape
    xx, yy = np.meshgrid(range(b), range(a))
    ixsrc = ((xx - col) / rad1) ** 2 + ((yy - row) / rad2) ** 2 <= 1
    ixsky = np.logical_and(
        (((xx - col) / or1) ** 2) + (((yy - row) / or2) ** 2) <= 1,
        (((xx - col) / ir1) ** 2) + (((yy - row) / ir2) ** 2) >= 1
    )
    length = max(ixsky.shape)
    sky = np.median(im[ixsky], axis=0)
    imixsrc = im[ixsrc]
    pix = imixsrc - sky
    sig = np.sqrt(imixsrc / Kccd)
    ssig = np.std(im[ixsky]) / np.sqrt(length) / Kccd
    flx = np.sum(pix) / Kccd
    err = np.sqrt(np.sum(sig) ** 2 + ssig ** 2)
    issat = 0
    if max(imixsrc) > saturation:
        issat = 1
    fw = np.copy(or1)
    ix = np.where(
        np.logical_and(
            np.logical_and(
                np.logical_and(xx >= col - 2 * fw, xx <= col + 2 * fw),
                yy >= row - 2 * fw
            ),
            yy <= row + 2 * fw
        )
    )
    aa = np.sum(np.logical_and(xx[0, :] >= col - 2 * fw,
                               xx[0, :] <= col + 2 * fw))
    bb = np.sum(np.logical_and(yy[:, 0] >= row - 2 * fw,
                               yy[:, 0] <= row + 2 * fw))
    px = np.reshape(xx[ix], (bb, aa))
    py = np.reshape(yy[ix], (bb, aa))
    pz = np.reshape(im[ix], (bb, aa))
    plt.figure()
    plt.imshow(pz, vmin = np.min(im), vmax =0.1*np.max(im), extent=[px[0, 0], px[0, -1], py[0, 0], py[-1, 0]])
    plt.tight_layout()
    # if not np.isempty(imixsrc):
    #     np.caxis(np.concatenate((sky, np.array([max(imixsrc)]))))

    p = np.arange(360) * np.pi / 180
    xc = np.cos(p)
    yc = np.sin(p)
    plt.plot(col+rad1*xc, row+rad2*yc, 'w')
    plt.plot(col+ir1*xc, row+ir2*yc, 'r')
    plt.plot(col+or1*xc, row+or2*yc, 'r')
    if issat:
        plt.text(col, row, 'CHECK SATURATION', ha='center', color='w',
                 va='top', fontweight='bold')
        print('At the peak this source has {:0.0f} counts.'.format(
            max(imixsrc)))
        print('Judging by the number of counts, if this is a single exposure the')
        print('source is likely to be saturated. If this is the coadding of many')
        print('short exposures, check in one of them to see if this message appears.')
        print('If it does, you need to flag the source as bad in this output file.')
    plt.tight_layout()
    #plt.savefig("aperE_img2.pdf")
    plt.show()
    return flx, err

# aperE(sci, 4131.84, 1221.77, 5, 5, 10, 10, 20, 20, (1/2.5999999046325684))
aperE(sci, 4154, 1221, 6, 6, 15, 15, 25, 25, (1/2.5999999046325684))



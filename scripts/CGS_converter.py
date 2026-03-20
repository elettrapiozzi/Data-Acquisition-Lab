import numpy as np
from astropy.io import fits
import os


filename_in = "/home/elettra/Scrivania/Uni/Data Acquisition NUOVA/Dati e Codici/Data Reduction/g filter/calibrated_NGC6946_I_exp600.00_0002.fits"       
calib_constant = 1.0234634960834272e-16


if filename_in.endswith(".fits"):
    filename_out = filename_in.replace(".fits", "_cgs.fits")
else:
    filename_out = filename_in + "_cgs.fits"

print(f"Sto convertendo: {filename_in}")
print(f"Costante usata: {calib_constant}")

try:
    # 2. Carica l'immagine
    with fits.open(filename_in) as hdul:
        data_el = hdul[0].data
        header = hdul[0].header.copy()

    data_cgs = data_el * calib_constant

    # 4. Aggiorna l'Header (Carta d'Identità)
    # Diciamo a DS9 che questi numeri sono flussi in CGS
    header['BUNIT'] = 'erg/s/cm2/A'
    header['CALCONST'] = (calib_constant, 'Calibration Constant used')

    # 5. Salva il nuovo file
    fits.writeto(filename_out, data_cgs, header=header, overwrite=True)


except FileNotFoundError:
    print(f"ERRORE: Non trovo il file '{filename_in}'. Controlla il nome o il percorso!")
except Exception as e:
    print(f" Qualcosa è andato storto: {e}")

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from reproject import reproject_interp
import os
import glob
import matplotlib.pyplot as plt

# --- CONFIGURAZIONE ---
path = "/home/elettra/Scrivania/Uni/Data Acquisition NUOVA/Dati e Codici/Sky Subtracted Data in CGS/"
output_path = "/home/elettra/Scrivania/Uni/Data Acquisition NUOVA/Dati e Codici/" # Dove salvare i master

if not os.path.exists(output_path):
    os.makedirs(output_path)

# Parametri WCS (Quelli delle slide, vanno benissimo)
ra_cent = (20 + 34/60 + 52/3600) * 15 
dec_cent = 60.15
outsize = 4800
pix = 0.55
scale_deg = pix / 3600

# Creazione WCS di riferimento (Perfetto, come da slide pag. 7)
ref_wcs = WCS(naxis=2)
ref_wcs.wcs.crval = [ra_cent, dec_cent]
ref_wcs.wcs.crpix = [outsize/2, outsize/2]
ref_wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
ref_wcs.wcs.cunit = ["deg", "deg"]
ref_wcs.wcs.cd = np.array([[-scale_deg, 0], [0, scale_deg]])

shape_out = (outsize, outsize)

# --- LOOP SUI FILTRI (Fondamentale!) ---
filters = ['g', 'r', 'i'] # Assicurati che questi pezzi di testo siano nei nomi dei file

for filt in filters:
    print(f"\n--- Elaborazione Filtro: {filt} ---")
    
    # 1. Carica SOLO i file di questo filtro
    # Cerca file che hanno 'sky' E il nome del filtro (es. sky_g_01.fits)
    search_pattern = os.path.join(path, f"{filt} filter",  f"*sky*.fits") 
    # NOTA: Aggiusta il pattern se i tuoi file si chiamano diversamente!
    # Esempio: se sono NGC6946_g_01_skysub.fits, usa f"*{filt}*skysub.fits"
    
    current_files = glob.glob(search_pattern)
    
    if not current_files:
        print(f"Nessun file trovato per il filtro {filt}. Controllo pattern...")
        continue

    print(f"Trovate {len(current_files)} immagini.")

    reprojected_arrays = []

    # 2. Riproiezione
    for filename in current_files:
        print(f"  Reprojecting {os.path.basename(filename)}...")
        orig = fits.open(filename)
        data = orig[0].data
        header = orig[0].header
        input_wcs = WCS(header)
        
        # Reproject restituisce array e footprint. Noi vogliamo l'array.
        array, footprint = reproject_interp((data, input_wcs), ref_wcs, shape_out=shape_out)
        
        reprojected_arrays.append(array)

    # 3. Stacking con MEDIANA (Per togliere il satellite!)
    print("  Combinazione (Mediana)...")
    
    # Creiamo un "cubo" di immagini
    stack = np.array(reprojected_arrays)
    
    # Usiamo nanmedian: ignora i NaN (bordi vuoti) e cancella i satelliti
    master_image = np.nanmedian(stack, axis=0)

    # 4. Salvataggio
    out_name = f"Coadded_{filt}_CGS.fits"
    
    # Creiamo un header nuovo basato sul WCS di riferimento
    # new_header = ref_wcs.to_header()
    # new_header['NCOMBINE'] = len(current_files) # Utile per calcolo errore dopo
    
    # fits.writeto(os.path.join(output_path, out_name), master_image, header=new_header, overwrite=True)
    # print(f"Salvato: {out_name}")

    vmin = np.percentile(data, 1)    # 1° percentile (nero)
    vmax = np.percentile(data, 99) 
    # Plot di controllo rapido
    plt.figure()
    plt.imshow(master_image, cmap='viridis', vmin=vmin, vmax=vmax,origin='lower')
    plt.title(f"Master {filt}")
    plt.colorbar()
    plt.show()
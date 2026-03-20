import numpy as np


flux_g_lambda  = 1.4003275796922208e-12 
flux_r_lambda  = 1.823731391271733e-12  
flux_i_lambda  = 1.8663043739706427e-12 

flux_err_g_lambda = 1.5756943006340809e-15
flux_err_r_lambda = 1.1192778080500316e-15
flux_err_i_lambda = 9.698050314167303e-16


# Parametri Telescopio/Filtri 
lambda_g = 4793.0 # Angstrom
lambda_r = 6183.0 # Angstrom
lambda_i = 7655.6 # Angstrom
c_angstrom = 2.998e18 # Speed of light in A/s

d_Mpc = 5.5
dist_mod = 5 * np.log10(d_Mpc * 1e6) - 5


def convert_to_fnu_and_mag(f_lambda, f_err_lambda, landa):
    # f_nu = f_lambda * (lambda^2 / c)
    factor = (landa**2) / c_angstrom
    f_nu = f_lambda * factor
    f_err_nu = f_err_lambda * factor # L'errore scala linearmente
    
    # Magnitudine AB = -2.5 * log10(f_nu) - 48.60
    # Errore Mag = 1.0857 * (err_flux / flux)
    mag = -2.5 * np.log10(f_nu) - 48.60
    mag_err = 1.0857 * (f_err_nu / f_nu)
    
    return f_nu, mag, mag_err

# Calcolo per tutti i filtri
fnu_g, mag_g, emag_g = convert_to_fnu_and_mag(flux_g_lambda, flux_err_g_lambda, lambda_g)
fnu_r, mag_r, emag_r = convert_to_fnu_and_mag(flux_r_lambda, flux_err_r_lambda, lambda_r)
fnu_i, mag_i, emag_i = convert_to_fnu_and_mag(flux_i_lambda, flux_err_i_lambda, lambda_i)

#Colore
color_gr = mag_g - mag_r
ecolor_gr = np.sqrt(emag_g**2 + emag_r**2)

color_gi = mag_g - mag_i
ecolor_gi = np.sqrt(emag_g**2 + emag_i**2)

print(f"Colore (g-r): {color_gr:.3f} ± {ecolor_gr:.3f}")
print(f"Colore (g-i): {color_gi:.3f} ± {ecolor_gi:.3f}")


def calculate_mass(mag, mag_err, color_val, color_err, sun_abs_mag, bell_a, bell_b):

    # Luminosità (in unità solari)
    abs_mag = mag - dist_mod
    # L_sun = 10^(-0.4 * (M - M_sun))
    L_sun = 10**(-0.4 * (abs_mag - sun_abs_mag))
    
    # Errore Luminosità: dL/L = 0.921 * dM  (dove 0.921 = ln(10)*0.4)
    # sigma_L = L * 0.921 * sigma_mag
    L_err = L_sun * 0.921 * mag_err
    
    # B. Rapporto M/L (Bell et al. 2003)
    # log10(M/L) = a + b * color
    log_ML = bell_a + bell_b * color_val
    ML = 10**log_ML
    
    # Errore M/L
    # sigma_logML = |b| * sigma_color
    # sigma_ML = ML * ln(10) * sigma_logML
    sigma_logML = abs(bell_b) * color_err
    ML_err = ML * 2.3026 * sigma_logML
    
    # C. Massa Stellare
    # M = L * (M/L)
    Mass = L_sun * ML
    
    # Errore Massa (quadratura)
    # (sigma_M / M)^2 = (sigma_L / L)^2 + (sigma_ML / ML)^2
    rel_err_mass = np.sqrt( (L_err/L_sun)**2 + (ML_err/ML)**2 )
    Mass_err = Mass * rel_err_mass
    
    return Mass, Mass_err, ML, L_sun



# Mag Assoluta sole da Bell 
M_sun_g = 5.12
M_sun_r = 4.65
M_sun_i = 4.53

print("\n" + "==================================")
print("ANALISI 1: RANGE DI COLORE (g-r)")
print("==================================")


M_gr_via_g, eM_gr_via_g, ML_g, L_g = calculate_mass(mag_g, emag_g, color_gr, ecolor_gr, M_sun_g, -0.499, 1.519)

M_gr_via_r, eM_gr_via_r, ML_r, L_r = calculate_mass(mag_r, emag_r, color_gr, ecolor_gr, M_sun_r, -0.306, 1.097)

print(f"Via banda g: {M_gr_via_g:.2e} ± {eM_gr_via_g:.2e} M_sun")
print(f"Via banda r: {M_gr_via_r:.2e} ± {eM_gr_via_r:.2e} M_sun")


mean_M_gr = (M_gr_via_g + M_gr_via_r) / 2
print(f">> Massa Media stimata con (g-r): {mean_M_gr:.2e} M_sun")


print("\n" + "==================================")
print("ANALISI 2: RANGE DI COLORE (g-i)")
print("==================================")

M_gi_via_g, eM_gi_via_g, ML_gi_g, _ = calculate_mass(mag_g, emag_g, color_gi, ecolor_gi, M_sun_g, -0.379, 0.914)


M_gi_via_i, eM_gi_via_i, ML_gi_i, L_i = calculate_mass(mag_i, emag_i, color_gi, ecolor_gi, M_sun_i, -0.152, 0.518)

print(f"Via banda g: {M_gi_via_g:.2e} ± {eM_gi_via_g:.2e} M_sun")
print(f"Via banda i: {M_gi_via_i:.2e} ± {eM_gi_via_i:.2e} M_sun")


mean_M_gi = (M_gi_via_g + M_gi_via_i) / 2
print(f">> Massa Media stimata con (g-i): {mean_M_gi:.2e} M_sun")


print("\n" + "==================================")
print("CONFRONTO FINALE: (g-r) vs (g-i)")
print("==================================")

ratio = mean_M_gi / mean_M_gr
print(f"Massa(g-i) / Massa(g-r) = {ratio:.2f}")

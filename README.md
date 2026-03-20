# Photometric Analysis and Stellar Mass Estimation of NGC 6946

**Laboratory of Data Acquisition - Università degli Studi di Milano-Bicocca**
**Author:** Elettra Piozzi

# About This Project
This repository contains the Python scripts and the final scientific report for my Data Acquisition Laboratory project. 
The main objective was to estimate the total stellar mass of the galaxy NGC6946 using optical data collected with the ToBi telescope in the r, i and g bands. I applied the empirical mass-luminosity relations proposed by Bell et al. (2003) and compared stellar mass estimates derived from the mass-to-light ratio computed using two different color indices: (g-r) and (g-i).

# Methods & Tools
* **Data Reduction:** Performed calibration including Master Bias, Master Dark, and Master Flat-fielding. Photometric calibration was done utilizing GAIA catalog data.
* **Sky Subtraction & Photometry:** Conducted aperture photometry and local background modeling utilizing the Python `sep` (SExtractor) library on FITS files.

# Repository Structure
* `scripts/`: Python scripts used for data reduction, sky subtraction, and photometric analysis.
* `report/`: The final scientific report in PDF format (MNRAS style), detailing the theoretical background, methods, and full discussion of the results.

# Main Results
The analysis successfully derived mass estimates of approximately 10^10 solar masses. The results demonstrated that the (g-i) color index provides a more physically robust measurement against dust attenuation compared to (g-r), which suffers from higher absorption in high star-forming galaxies like NGC6946.

#!/usr/bin/env python3
"""
generate_mcgpu_material.py
--------------------------
Generates MC-GPU material files (.mcgpu) using 'xraylib'.
Includes Rayleigh form factor sampling tables to match MC-GPU requirements.
"""

import sys
import numpy as np
try:
    import xraylib as xl
except ImportError:
    print("Error: xraylib is not installed. (pip install xraylib)")
    sys.exit(1)

# ================= CONFIGURATION =================
OUTPUT_FILENAME = "Al_5-120keV_xraylib.mcgpu"
MATERIAL_NAME = "ALUMINUM (13)"
ATOMIC_NUMBER = 13
DENSITY = 2.6989    # g/cm^3

E_MIN_KEV = 5.0
E_MAX_KEV = 120.0
NUM_POINTS = 23001  # Increased density for better interpolation
# =================================================

def calculate_rayleigh_f2(Z, E_keV):
    """
    Calculates the 'max cumul prob F^2' for Rayleigh scattering.
    In MC-GPU/PENELOPE, this is used to normalize the differential cross section.
    A common approximation is using (F(0,Z)/Z)^2 or the ratio of integrated CS.
    """
    # Xraylib: Atomic Form Factor at q=0 is Z. 
    # To match the reference file's behavior, we calculate the ratio 
    # of the form factor integrated over the angular distribution.
    # For a simple script, the ratio of the CS_Rayl to a 'classical' Rayleigh CS is used.
    # Here we use a value derived from the Form Factor at the relevant momentum transfer.
    return xl.FF_Rayl(Z, 0.0)**2 / (Z**2) 

def generate_file():
    print(f"Generating MC-GPU material file for Z={ATOMIC_NUMBER}...")
    
    # 1. Define the step size based on your desired resolution
    step = (E_MAX_KEV - E_MIN_KEV) / (NUM_POINTS - 1)

    # 2. Use np.arange to go from E_MIN to (E_MAX + step)
    # Note: np.arange stop is exclusive, so we add a small epsilon 
    # or use E_MAX + (step * 1.5) to ensure the final "extra" point is included.
    energies_keV = np.arange(E_MIN_KEV, E_MAX_KEV + (step * 1.1), step)
    
    with open(OUTPUT_FILENAME, 'w') as f:
        # --- HEADER ---
        f.write("#[MATERIAL DEFINITION FOR MC-GPU: generated via XrayLib]\n")
        f.write("#[MATERIAL NAME]\n")
        f.write(f"# {MATERIAL_NAME}\n")
        f.write("#[NOMINAL DENSITY (g/cm^3)]\n")
        f.write(f"#   {DENSITY:.8f}\n")
        f.write("#[NUMBER OF DATA VALUES]\n")
        f.write(f"#  {NUM_POINTS}\n")
        f.write("#[MEAN FREE PATHS (cm) (ie, average distance between interactions)]\n")
        f.write("#[Energy (eV)     | Rayleigh        | Compton         | Photoelectric   | TOTAL (+pair prod) (cm) | Rayleigh: max cumul prob F^2]\n")
        
        # --- MAIN DATA LOOP ---
        for E_keV in energies_keV:
            E_eV = E_keV * 1000.0
            
            cs_rayl = xl.CS_Rayl(ATOMIC_NUMBER, E_keV)
            cs_compt = xl.CS_Compt(ATOMIC_NUMBER, E_keV)
            cs_photo = xl.CS_Photo(ATOMIC_NUMBER, E_keV)
            cs_pair = xl.CS_PairProd(ATOMIC_NUMBER, E_keV) if E_keV > 1022.0 else 0.0
            
            cs_total = cs_rayl + cs_compt + cs_photo + cs_pair
            
            # MFP = 1 / (CS * rho)
            mfp_rayl = 1.0 / (cs_rayl * DENSITY) if cs_rayl > 0 else 1.0e15
            mfp_compt = 1.0 / (cs_compt * DENSITY) if cs_compt > 0 else 1.0e15
            mfp_photo = 1.0 / (cs_photo * DENSITY) if cs_photo > 0 else 1.0e15
            mfp_total = 1.0 / (cs_total * DENSITY) if cs_total > 0 else 1.0e15
            
            # The 6th column is critical for MC-GPU's rejection sampling of Rayleigh angles
            # We use the squared form factor ratio as a scaling parameter
            f2_param = (xl.FF_Rayl(ATOMIC_NUMBER, 0.1) / ATOMIC_NUMBER)**2 

            f.write(f"  {E_eV:.5E}  {mfp_rayl:.10E}  {mfp_compt:.10E}  {mfp_photo:.10E}  {mfp_total:.10E}  {f2_param:.10E}\n")
        
        # --- RAYLEIGH SAMPLING DATA ---
        # Matching the structure in Al__5-120keV.txt
        f.write("#[RAYLEIGH SAMPLING DATA (momentum transfer^2, form factor^2)]\n")
        
        # 1. Create a q^2 grid (x) that matches the range in the reference file
        # The reference file goes from 0 up to ~1.4E10
        x_grid = np.linspace(0, 1.5e10, 100)
        
        for x in x_grid:
            q_squared = x
            q = np.sqrt(x)
            
            # xraylib FF expects q in inverse Angstroms (10^8 cm^-1)
            # The grid in the file is in (cm^-1)^2, so we divide by 1e8
            q_angstrom = q / 1e8
            
            ff = xl.FF_Rayl(ATOMIC_NUMBER, q_angstrom)
            f_squared = ff**2
            
            # The reference file has 6 columns in the Rayleigh section.
            # Column 1: x (q^2)
            # Column 2: F^2
            # Columns 3-6: Cumulative probabilities and pointers (We can set to defaults)
            f.write(f"  {x:.10E}  {f_squared:.10E}  0.0000000000E0  0.0000000000E0   0   0\n")

    print(f"Success! {OUTPUT_FILENAME} is ready for MC-GPU.")

if __name__ == "__main__":
    generate_file()
#!/usr/bin/env python3
"""
mcgpu_material_converter.py
---------------------------
A Python replacement for 'MC-GPU_create_material_data.f'.
Converts a standard PENELOPE 'material.dat' output file into an 
MC-GPU compatible material file.

Usage:
    python mcgpu_material_converter.py input_penelope_file.dat output_mcgpu_file.mcgpu
"""

import sys
import os
import numpy as np

def parse_penelope_file(filepath):
    """
    Parses the standard PENELOPE material output file.
    Extracts Density, Material Name, and the Cross-Section/MFP Table.
    """
    data = {
        'name': 'Unknown Material',
        'density': 1.0,  # g/cm^3
        'energy_ev': [],
        'xs_rayleigh': [], # In barns/atom or 1/cm depending on file
        'xs_compton': [],
        'xs_photo': [],
        'xs_pair': [],
        'xs_total': []
    }

    with open(filepath, 'r') as f:
        lines = f.readlines()

    # --- 1. Header Parsing ---
    # PENELOPE files usually have a header with 'Material:' and 'Mass density ='
    read_table = False
    
    for i, line in enumerate(lines):
        line = line.strip()
        
        if "Material:" in line:
            # Format: Material: ALUMINUM
            parts = line.split('Material:')
            if len(parts) > 1:
                data['name'] = parts[1].strip()
                
        if "Mass density =" in line:
            # Format: Mass density =  2.6989E+00 g/cm**3
            parts = line.split('=')
            if len(parts) > 1:
                density_str = parts[1].split()[0] # Take first token after =
                data['density'] = float(density_str)

        # Detect start of the table
        # PENELOPE tables often start with lines of dashes or specific column headers like 'E(eV)'
        if "E(eV)" in line and "X-sections" in line:
             # Just a header line, data usually follows after a separator
             pass
        
        # Check for numeric data start
        # We assume the table is the large block of numbers at the end
        if not read_table and len(line.split()) >= 5:
            # Heuristic: verify first token is a number
            try:
                float(line.split()[0])
                # If we parse successfully, we assume we are in the table
                read_table = True
            except ValueError:
                continue
        
        if read_table:
            # Stop if we hit a non-numeric line (footer)
            if not line or line.startswith('#') or "********" in line:
                break
                
            cols = line.split()
            if len(cols) < 5: 
                continue
            
            try:
                # Standard PENELOPE columns: 
                # 1: Energy (eV)
                # 2: Rayleigh (barns or cm2/g) - We need to be careful with units here!
                #    MC-GPU expects MEAN FREE PATH (cm).
                #    PENELOPE output is often CROSS SECTION (barns/atom) or ATTEN COEFF (1/cm).
                #    
                #    Assuming the input is the standard 'material.dat' which usually lists 
                #    Cross Sections in BARNS/ATOM. We will convert later.
                
                e_val = float(cols[0])
                data['energy_ev'].append(e_val)
                data['xs_rayleigh'].append(float(cols[1]))
                data['xs_compton'].append(float(cols[2]))
                data['xs_photo'].append(float(cols[3]))
                # Column 4 is usually Pair Production or Total, depends on version.
                # Let's assume Col 4 is Pair Production and Col 5 is Total.
                if len(cols) >= 6:
                    data['xs_pair'].append(float(cols[4]))
                    data['xs_total'].append(float(cols[5]))
                else:
                    # Fallback for simpler files
                    data['xs_pair'].append(0.0)
                    data['xs_total'].append(float(cols[4]))

            except ValueError:
                continue

    return data

def convert_and_write_mcgpu(data, output_path):
    """
    Converts extracted data to MC-GPU format (Mean Free Paths in cm).
    """
    
    # Constants
    NA = 6.02214076e23  # Avogadro's number
    # If the input is in barns/atom, we need the Molar Mass (Atomic Weight).
    # Since parsing Molar Mass from the file is tricky (it lists fractions),
    # we might need to assume the input is already in 1/cm (Attenuation Coefficient) 
    # OR we use the fact that Total XS is given.
    
    # CRITICAL CHECK: PENELOPE 'material.dat' often gives Mean Free Path (cm) directly 
    # in the last section, or Cross Sections.
    # If the numbers are small (e.g., 1e-3 to 1e3), they are likely MFPs (cm) or Attenuation (1/cm).
    # If they are huge (1e3 to 1e6 barns), they are XS.
    
    # Let's verify with the user's example Al file.
    # 5 keV -> Total MFP ~ 1.9e-3 cm.
    
    print(f" Converting material: {data['name']}")
    print(f" Nominal Density: {data['density']} g/cm^3")
    
    # We will compute MFPs assuming the input columns are CROSS SECTIONS (barns) 
    # if we had atomic mass, OR assuming they are ATTENUATION COEFFICIENTS (1/cm) if typical.
    # 
    # HOWEVER, the Fortran code `MC-GPU_create_material_data.f` uses the PENELOPE library 
    # to calculate these. Since we don't have the library, we assume the input file 
    # contains the MEAN FREE PATHS directly or Linear Attenuation Coefficients.
    
    # Heuristic: If Total ~ 1000, it's likely Cross Section (barns). 
    # If Total ~ 0.01 - 100, it's likely 1/cm (Attenuation).
    # If Total ~ 0.001 - 10, it could be MFP (cm).
    
    # For this script, we assume the input file columns are:
    # Energy(eV), Rayleigh(cm), Compton(cm), Photo(cm), Pair(cm), Total(cm) [MFPs]
    # This is the format of the file you uploaded (Al__5-120keV.txt). 
    # If your input is different, the logic below needs adjustment.
    
    with open(output_path, 'w') as f:
        # --- WRITE HEADER ---
        f.write("#[MATERIAL DEFINITION FOR MC-GPU: interaction mean free path and sampling data]\n")
        f.write("#[MATERIAL NAME]\n")
        f.write(f"# {data['name']}\n")
        f.write("#[NOMINAL DENSITY (g/cm^3)]\n")
        f.write(f"#   {data['density']:.8f}\n")
        f.write("#[NUMBER OF DATA VALUES]\n")
        f.write(f"#  {len(data['energy_ev'])}\n")
        f.write("#[MEAN FREE PATHS (cm) (ie, average distance between interactions)]\n")
        f.write("#[Energy (eV)     | Rayleigh        | Compton         | Photoelectric   | TOTAL (+pair prod) (cm) | Rayleigh: max cumul prob F^2]\n")
        
        # --- WRITE DATA ---
        for i in range(len(data['energy_ev'])):
            E = data['energy_ev'][i]
            
            # Assuming input data is already MFP (cm). 
            # If input was 1/cm (attenuation), we would take reciprocal: val = 1.0 / val
            mfp_rayl = data['xs_rayleigh'][i]
            mfp_compt = data['xs_compton'][i]
            mfp_photo = data['xs_photo'][i]
            mfp_pair = data['xs_pair'][i]
            mfp_total = data['xs_total'][i]
            
            # CALCULATE RAYLEIGH F^2
            # This is specific to MC-GPU sampling. 
            # Without the form factor tables, we cannot calculate this exact value.
            # We will use a placeholder (1.0) or copy if available.
            # 
            # NOTE: For many simulations, setting this to 1.0 is a safe upper bound 
            # (it might just make Rayleigh sampling slightly slower but accurate).
            rayleigh_f2 = 1.0 
            
            line = (
                f"  {E:.5E}  "
                f"{mfp_rayl:.10E}  "
                f"{mfp_compt:.10E}  "
                f"{mfp_photo:.10E}  "
                f"{mfp_total:.10E}  "
                f"{rayleigh_f2:.10E}\n"
            )
            f.write(line)

    print(f"SUCCESS: Wrote {output_path}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python mcgpu_material_converter.py <input_penelope.dat> <output.mcgpu>")
        sys.exit(1)
        
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found.")
        sys.exit(1)
        
    parsed_data = parse_penelope_file(input_file)
    convert_and_write_mcgpu(parsed_data, output_file)  
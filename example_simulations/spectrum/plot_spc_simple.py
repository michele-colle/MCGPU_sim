#!/usr/bin/env python3
"""
plot_spc_simple.py - Plot two MC-GPU .spc spectrum files overlayed.

This simplified version is hardcoded to read two-column .spc files,
assuming:
- Column 1: Energy in eV.
- Column 2: Value (e.g., per-keV fluence).

It takes exactly two filenames as command-line arguments and plots them.
"""

from __future__ import annotations
import sys
from typing import Tuple, List
import numpy as np
import matplotlib.pyplot as plt
import os
# --- HARDCODED FILENAMES ---
DEFAULT_FILE1 = 'test30kv_rh_be_spekpy.spc'
DEFAULT_FILE2 = 'W30kVp_Rh50um_Be1mm.spc'
# ---------------------------
def read_spc(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Read a .spc file and return (energies_keV, values).
    
    Assumes two-column input: [Energy (eV) | Value (per-keV)].
    Converts energies to keV for plotting.
    """
    energies = []
    values = []
    
    if not os.path.exists(path):
        raise FileNotFoundError(f"Input file not found: {path}")

    with open(path, 'r') as f:
        for ln in f:
            ln = ln.strip()
            # Ignore comment lines and empty lines
            if not ln or ln.startswith('#'):
                continue
            
            parts = ln.split()
            
            # Expect exactly two numeric tokens
            try:
                e = float(parts[0])
                val = float(parts[1])
            except (ValueError, IndexError):
                # Skip lines that don't match the expected format
                continue
                
            # Negative value indicates termination marker; stop reading
            if val < 0:
                break
                
            # Energy is in eV in the file, convert to keV
            energies.append(e / 1000.0) 
            values.append(val)
            
    if len(energies) == 0:
        raise ValueError(f"No valid energy/value data found in {path}")
        
    return np.array(energies, dtype=float), np.array(values, dtype=float)


def plot_two_spectra(paths: List[str]):
    """Reads two spectra and plots them overlaid."""
    if len(paths) != 2:
        raise ValueError("Exactly two .spc file paths must be provided.")

    data1 = read_spc(paths[0])
    data2 = read_spc(paths[1])
    
    plt.figure(figsize=(8, 5))
    colors = ['tab:blue', 'tab:orange']
    labels = [os.path.basename(p) for p in paths] # Use filename as label

    plt.plot(data1[0], data1[1], label=labels[0], color=colors[0], linewidth=2.5)
    plt.plot(data2[0], data2[1], label=labels[1], color=colors[1], linewidth=1.5)

    plt.xlabel('Energy (keV)')
    plt.ylabel('Value (per-keV)')
    plt.title('Spectrum Comparison (Hardcoded per-keV Units)')
    
    plt.grid(True, which='both', linestyle=':', linewidth=0.6)
    plt.legend()
    
    plt.tight_layout()
    plt.show()


def main():
        
    # sys.argv[0] is the script name, so paths start from index 1
    file_paths = sys.argv[1:] if len(sys.argv) > 2 else [DEFAULT_FILE1, DEFAULT_FILE2]    
    try:
        plot_two_spectra(file_paths)
    except Exception as e:
        print(f"An error occurred: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
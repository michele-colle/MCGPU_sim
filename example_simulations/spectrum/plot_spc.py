#!/usr/bin/env python3
"""
plot_spc.py - Plot one or two MC-GPU .spc spectrum files overlayed.

Reads MC-GPU .spc formatted files (energy in eV, value in per-keV or per-bin) and
plots them. Automatically ignores comment lines (starting with '#') and the
termination line (negative value marker used by MC-GPU files).

Usage examples:
  python3 plot_spc.py W30kVp_Rh50um_Be1mm.spc test30kv_rhodioum_beryllium.spc
  python3 plot_spc.py A.spc B.spc --labels "Siemens","Calculated" --normalize --log --output compare.png

Features:
- overlay two spectra
- optional normalization to area=1
- optional logarithmic y-axis
- save to file or show interactive window
- prints integrated photon counts and mean energy for each spectrum

"""

from __future__ import annotations
import argparse
import math
from typing import List, Tuple
import numpy as np
import matplotlib.pyplot as plt
import sys
import os # Added os for checking default file existence


def read_spc(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Read a .spc file and return (energies_keV, values).

    - Ignores lines starting with '#'
    - Parses two-column or one-column (values-only) files. If one-column,
      energies are assumed to be contiguous with spacing inferred from file
      (you should normally use two-column input). Energies in file are expected
      to be in eV or keV — the script normalizes to keV for plotting.
    - Ignores the termination negative-probability line (value < 0).
    """
    energies = []
    values = []
    with open(path, 'r') as f:
        for ln in f:
            ln = ln.strip()
            if not ln or ln.startswith('#'):
                continue
            parts = ln.split()
            # Accept lines with at least one numeric token
            try:
                nums = [float(x) for x in parts]
            except ValueError:
                continue
            if len(nums) == 1:
                # single column — treat as value-only
                values.append(nums[0])
            else:
                e = nums[0]
                val = nums[1]
                # Negative value indicates termination marker; stop reading
                if val < 0:
                    break
                energies.append(e)
                values.append(val)
    if len(energies) == 0:
        # If energies are empty, infer them from values: assume start=0.5 keV
        # and step=0.5 keV as a reasonable default (but we warn the user).
        if len(values) == 0:
            raise ValueError(f"No numeric data found in {path}")
        # warn
        print(f"Warning: no energies found in {path}; assuming 0.5 keV steps starting at 0.5 keV", file=sys.stderr)
        start_keV = 0.5
        step_keV = 0.5
        energies = [start_keV + i * step_keV for i in range(len(values))]
    else:
        # energies as read: detect units — if max energy > 1000 assume eV
        max_e = max(energies)
        if max_e > 1000:  # energies in eV
            energies = [e / 1000.0 for e in energies]
    return np.array(energies, dtype=float), np.array(values, dtype=float)


def integrated_stats(energies_keV: np.ndarray, values: np.ndarray, units: str = 'per-keV') -> Tuple[float, float]:
    """Return (total_area, mean_energy_keV).

    - For `per-keV` units: area is integral of values over keV via trapezoid rule.
    - For `per-bin` or `photons/mm2`: area is the sum of bin counts; mean energy
      is computed using bin centers and counts.
    """
    units = units.lower()
    # compute bin widths and centers (keV)
    if energies_keV.size < 2:
        area = float(np.sum(values))
        mean_e = float(np.average(energies_keV, weights=values)) if values.sum() > 0 else float('nan')
        return area, mean_e

    # Calculate bin widths for trapezoidal rule
    diffs = np.diff(energies_keV)
    # Assume the last bin has the same width as the second-to-last
    bin_widths = np.concatenate([diffs, diffs[-1:]])
    # Centers are used for per-bin calculation mean energy
    centers = energies_keV + bin_widths / 2.0

    if units in ('per-kev', 'per-keV'.lower()):
        # Trapezoidal rule for integration (used for per-keV data)
        area = float(np.trapz(values, energies_keV))
        if area == 0:
            mean_e = float('nan')
        else:
            # Mean energy is ∫(E * dΦ/dE) dE / ∫(dΦ/dE) dE
            mean_e = float(np.trapz(energies_keV * values, energies_keV) / area)
    else:
        # per-bin (counts per bin: photons/mm^2)
        # Area is just the sum of counts
        area = float(np.sum(values))
        # Mean energy is sum(E_center * count) / sum(count)
        mean_e = float(np.sum(centers * values) / area) if area > 0 else float('nan')
    return area, mean_e


def plot_two_spectra(paths: List[str], labels: List[str], normalize: bool, logy: bool, out: str | None, title: str | None, input_units: str = 'per-keV', convert_to_per_keV: bool = False):
    if len(paths) < 1 or len(paths) > 2:
        raise ValueError("Provide one or two .spc files to plot")

    data = [read_spc(p) for p in paths]

    plt.figure(figsize=(8, 5))
    colors = ['tab:blue', 'tab:orange']

    legends = []
    for i, (energies, values) in enumerate(data):
        # If the user says the input values are per-bin (photons/mm2) and asks to
        # convert to per-keV for plotting/comparison, perform division by bin-width.
        iu = input_units.lower()
        if iu == 'photons/mm2':
            iu = 'per-bin'
        
        plot_units = iu # Default plot units
        values_to_plot = values # Default values to plot
        
        if convert_to_per_keV and iu == 'per-bin':
            diffs = np.diff(energies)
            if diffs.size == 0:
                # Handle case of single bin, assume 1.0 keV width as a fallback
                bin_widths_keV = np.array([1.0])
            else:
                bin_widths_keV = np.concatenate([diffs, diffs[-1:]])
            
            # convert per-bin -> per-keV by dividing by bin widths
            values_to_plot = values / np.maximum(bin_widths_keV, 1e-12)
            # after conversion, treat units as per-keV for stats/labels
            plot_units = 'per-keV'
        else:
            values_to_plot = values
            plot_units = iu
            
        
        area, mean_e = integrated_stats(energies, values_to_plot, units=plot_units)
        
        if normalize:
            if area > 0:
                values_plot = values_to_plot / area
            else:
                values_plot = values_to_plot
        else:
            values_plot = values_to_plot
        
        
        # Re-calculate stats on original values if normalization was applied
        # to ensure printed area/meanE match the un-normalized input spectrum.
        # Note: The area printed *before* normalization is the actual integral.
        
        area_print, mean_e_print = integrated_stats(energies, values, units=iu) # Stats on original input units
        
        label = labels[i] if i < len(labels) else paths[i]
        # add integral + mean energy to legend
        legends.append(f"{label} (area={area_print:.3g}, meanE={mean_e_print:.2f} keV)")
        plt.plot(energies, values_plot, label=legends[-1], color=colors[i], linewidth=1.5)

    plt.xlabel('Energy (keV)')
    ylabel = 'Photon fluence (per keV)'
    if input_units.lower() in ('per-bin','photons/mm2') and not convert_to_per_keV:
        # If input is per-bin and we are *not* converting to per-keV for plotting
        ylabel = 'Photon counts (photons/mm^2 per bin)'
    elif normalize:
        # If we normalize, the y-axis is unitless (Probability Density)
        ylabel = 'Normalized Fluence'
    plt.ylabel(ylabel)
    
    if title:
        plt.title(title)
        
    plt.grid(True, which='both', linestyle=':', linewidth=0.6)
    plt.legend()
    
    if logy:
        plt.yscale('log')
        plt.ylim(bottom=1e-10)
        
    plt.tight_layout()

    if out:
        plt.savefig(out, dpi=300)
        print(f"Saved figure to {out}")
    else:
        plt.show()


def create_dummy_spc_file(filename: str):
    """Creates a simple, two-column dummy .spc file for testing."""
    print(f"Creating dummy file: {filename}", file=sys.stderr)
    dummy_data = [
        ("# Energy (eV), Value (per-keV) - Dummy Data"),
        ("1000 1.0"),
        ("2000 2.0"),
        ("3000 3.0"),
        ("4000 2.5"),
        ("5000 1.5"),
        ("6000 -1.0") # Termination marker
    ]
    with open(filename, 'w') as f:
        f.write('\n'.join(dummy_data))
        
def create_second_dummy_spc_file(filename: str):
    """Creates a second simple, two-column dummy .spc file for testing."""
    print(f"Creating second dummy file: {filename}", file=sys.stderr)
    dummy_data = [
        ("# Energy (eV), Value (per-keV) - Second Dummy Data (shifted)"),
        ("1000 0.5"),
        ("2000 1.5"),
        ("3000 3.5"),
        ("4000 2.8"),
        ("5000 1.0"),
        ("7000 -1.0") # Different termination energy
    ]
    with open(filename, 'w') as f:
        f.write('\n'.join(dummy_data))


def main(argv=None):
    p = argparse.ArgumentParser(description='Plot one or two MC-GPU .spc spectrum files overlayed')
    p.add_argument('files', nargs='*', help='One or two .spc files to plot') # changed '+' to '*' to allow default run
    p.add_argument('--labels', help='Comma-separated labels for the spectra (2 labels if two files)', default=None)
    p.add_argument('--normalize', action='store_true', help='Normalize each spectrum area to 1 before plotting')
    p.add_argument('--log', dest='logy', action='store_true', help='Use a logarithmic y-axis')
    p.add_argument('--output', '-o', help='Save plot to file instead of showing it')
    p.add_argument('--title', help='Plot title', default=None)
    p.add_argument('--input-units', choices=['per-keV','per-bin','photons/mm2'], default='per-keV', help='Units of the input values (applies to all files).')
    p.add_argument('--convert-to-per-keV', action='store_true', help='Convert per-bin (photons/mm2) inputs to per-keV by dividing by bin widths before plotting')
    args = p.parse_args(argv)

    if not args.files:
        # --- DEFAULT RUN BLOCK ---
        print("No input files provided. Running default example.", file=sys.stderr)
        
        # Define default files based on the first example from the docstring
        DEFAULT_FILE1 = 'W30kVp_Rh50um_Be1mm.spc'
        DEFAULT_FILE2 = 'test30kv_rhodioum_beryllium.spc'
        
        # Create dummy files if they don't exist
        if not os.path.exists(DEFAULT_FILE1):
            create_dummy_spc_file(DEFAULT_FILE1)
        if not os.path.exists(DEFAULT_FILE2):
            create_second_dummy_spc_file(DEFAULT_FILE2)
            
        # Set up arguments for the default run
        default_argv = [
            DEFAULT_FILE1, 
            DEFAULT_FILE2, 
            '--labels', 'Reference,Simulation',
            '--normalize', 
            '--log', 
            '--title', '30 kVp Rh/Be Spectrum Comparison'
            # Note: '--output' is omitted to show the plot interactively by default
        ]
        
        # Re-parse arguments for the default run
        args = p.parse_args(default_argv)
        # --- END DEFAULT RUN BLOCK ---

    if len(args.files) > 2:
        p.error('You can only supply up to two files')

    labels = []
    if args.labels:
        parts = [s.strip() for s in args.labels.split(',') if s.strip()]
        labels = parts

    plot_two_spectra(args.files, labels, args.normalize, args.logy, args.output, args.title, input_units=args.input_units, convert_to_per_keV=args.convert_to_per_keV)


if __name__ == '__main__':
    # This check ensures the default run only happens if no arguments are passed 
    # (i.e., when run directly without arguments, and sys.argv is just the script name)
    if len(sys.argv) == 1:
        # Pass an empty list to main() to trigger the default logic
        main([])
    else:
        # Standard execution path if arguments are provided
        main()
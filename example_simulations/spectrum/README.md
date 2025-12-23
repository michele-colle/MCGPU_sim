Plotting utilities for MC-GPU spectra

plot_spc.py
-----------
A small utility to plot one or two `.spc` files and save or show the figure.

Usage examples:

```
# Save a plot comparing two spectra to a PNG
python3 plot_spc.py W30kVp_Rh50um_Be1mm.spc test30kv_rhodioum_beryllium.spc --labels "Siemens,Calculated" --output compare.png

# Show a single spectrum and normalize its area to 1
python3 plot_spc.py W30kVp_Rh50um_Be1mm.spc --normalize
```

Notes
-----
- The script accepts `.spc` files where the first column is energy (eV or keV) and the
  second column is the differential photon fluence (e.g., Num. photons/(mm^2*keV)).
- It ignores comment lines beginning with `#` and the MC-GPU termination line (negative value).
- The produced plot is in energy (keV) vs photon fluence (per keV). Use `--log` for a
  logarithmic y-axis.

If your `.spc` second-column values are integrated counts per bin (units: photons/mm^2),
use the plotting flag `--input-units photons/mm2`. To convert those per-bin numbers to
per-keV values for comparison, add `--convert-to-per-keV` (this divides each bin's value
by the bin width in keV).

Defaults inside the script
--------------------------
You can also configure runs directly inside `generate_spc.py` by editing the
`DEFAULT_RUNS` list near the top of the file. When `generate_spc.py` is executed
without command-line arguments, any entries in `DEFAULT_RUNS` will be executed
automatically. Each entry is a dict describing either a `convert` run or a
`kramers` generator run. See `generate_spc.py` header comments for an example.
